"""openFDA drug/event API 客户端。

设计要点
--------
* 使用 ``count=`` 聚合端点,一次请求即得到某药物的全部反应分布,
  避免逐反应请求(否则 N 药物 x M 反应 会产生 N*M 次调用)。
* 全库"药物-反应"报告对总数 N 通过一次独立的 count 查询(按药物聚合后求和)获得,
  而不是获取全部 2000 万条原始报告。
* 限流:调用之间固定 sleep(默认 1.0s),并对 429/5xx 做指数退避重试。
* openFDA 在查询无匹配时返回 HTTP 404(语义为"无结果",不是错误),
  本模块将其转换为空结果而非异常。

API 文档: https://open.fda.gov/apis/drug/event/
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

__all__ = ["OpenFDAError", "OpenFDAClient", "NoResultsError"]

BASE_URL = "https://api.fda.gov/drug/event.json"

# openFDA 建议的礼貌间隔;过快会被 429 限流
DEFAULT_MIN_INTERVAL = 1.0


class OpenFDAError(RuntimeError):
    """openFDA 返回了无法自动恢复的错误。"""


class NoResultsError(OpenFDAError):
    """查询无匹配结果(HTTP 404)。"""


@dataclass
class CountResult:
    """一次 count 聚合查询的结果。"""

    term: str
    count: int
    total: Optional[int] = None          # meta.results.total,该检索式的总匹配数
    raw: dict = field(default_factory=dict)


class OpenFDAClient:
    """带限流、重试与本地缓存的 openFDA 客户端。

    参数
    ----
    api_key : 可选,无 key 时限流为 240 req/min,本工具默认 1 req/s 远低于该上限。
    min_interval : 两次网络请求之间的最小间隔秒数。
    timeout : 单次请求超时秒数。
    max_retries : 429/5xx/网络错误的退避重试次数。
    use_cache : 是否在同一进程内缓存相同 URL 的响应(重复分析同一药物时省请求)。
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        min_interval: float = DEFAULT_MIN_INTERVAL,
        timeout: float = 45.0,
        max_retries: int = 4,
        use_cache: bool = True,
        verbose: bool = True,
        sleep_fn: Callable[[float], None] = time.sleep,
    ) -> None:
        self.api_key = api_key
        self.min_interval = max(0.0, float(min_interval))
        self.timeout = timeout
        self.max_retries = max_retries
        self.use_cache = use_cache
        self.verbose = verbose
        self._sleep = sleep_fn
        self._last_request_ts = 0.0
        self._cache: Dict[str, dict] = {}
        self.request_count = 0
        # grand_total() 实际使用的聚合条数,写入报告以说明数据口径
        self.grand_total_limit: Optional[int] = None

    # ------------------------------------------------------------------
    # 底层请求
    # ------------------------------------------------------------------
    def _throttle(self) -> None:
        """保证两次真实网络请求之间至少间隔 min_interval 秒。"""
        elapsed = time.monotonic() - self._last_request_ts
        wait = self.min_interval - elapsed
        if wait > 0:
            self._sleep(wait)

    def _log(self, msg: str) -> None:
        if self.verbose:
            print(f"[openfda] {msg}", flush=True)

    def fetch(self, params: Dict[str, object]) -> Optional[dict]:
        """执行一次 GET,返回解析后的 JSON;无结果(404)时返回 None。

        出错策略:404 -> None(空结果);429/5xx/URLError -> 指数退避重试;
        其它 4xx -> 抛 OpenFDAError(通常是查询语法问题,重试无意义)。
        """
        if self.api_key:
            params = {**params, "api_key": self.api_key}
        url = BASE_URL + "?" + urllib.parse.urlencode(params, safe=':"+()')

        if self.use_cache and url in self._cache:
            return self._cache[url]

        last_err: Optional[Exception] = None
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "FAERSSignal/1.0"})
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    payload = json.loads(resp.read().decode("utf-8"))
                self._last_request_ts = time.monotonic()
                self.request_count += 1
                if self.use_cache:
                    self._cache[url] = payload
                return payload
            except urllib.error.HTTPError as e:
                self._last_request_ts = time.monotonic()
                self.request_count += 1
                if e.code == 404:
                    # openFDA 以 404 表示"该检索式无匹配",属于正常空结果
                    if self.use_cache:
                        self._cache[url] = {}
                    return None
                if e.code == 429 or 500 <= e.code < 600:
                    last_err = e
                    backoff = min(30.0, 2.0 ** attempt) + self.min_interval
                    self._log(f"HTTP {e.code},第 {attempt + 1} 次退避重试,等待 {backoff:.1f}s")
                    self._sleep(backoff)
                    continue
                raise OpenFDAError(f"openFDA 返回 HTTP {e.code}: {e.reason}") from e
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
                self._last_request_ts = time.monotonic()
                self.request_count += 1
                last_err = e
                backoff = min(30.0, 2.0 ** attempt) + self.min_interval
                self._log(f"网络错误 {type(e).__name__},第 {attempt + 1} 次退避重试,等待 {backoff:.1f}s")
                self._sleep(backoff)

        raise OpenFDAError(f"重试 {self.max_retries} 次后仍失败: {last_err}")

    # ------------------------------------------------------------------
    # 高层查询
    # ------------------------------------------------------------------
    def count_reactions_for_drug(
        self,
        drug: str,
        limit: int = 200,
        field: str = "patient.reaction.reactionmeddrapt.exact",
    ) -> List[CountResult]:
        """某药物报告的全部反应(MedDRA 首选术语)及其报告数。

        检索式对药物名做精确短语匹配 ``patient.drug.medicinalproduct:"<DRUG>"``。
        openFDA 的药物名已规范为大写,大小写不敏感。
        """
        search = f'patient.drug.medicinalproduct:"{drug}"'
        payload = self.fetch({"search": search, "count": field, "limit": limit})
        if not payload or "results" not in payload:
            return []

        total = payload.get("meta", {}).get("results", {}).get("total")
        out: List[CountResult] = []
        for row in payload["results"]:
            term = str(row.get("term", "")).strip()
            if not term:
                continue
            out.append(CountResult(term=term, count=int(row.get("count", 0)), total=total))
        return out

    def count_drugs_for_reaction(
        self,
        reaction: str,
        limit: int = 800,
        field: str = "patient.drug.medicinalproduct.exact",
    ) -> List[CountResult]:
        """某反应报告涉及的全部药物及其报告数(用于估计反应的背景报告率)。

        注意:该 count 聚合受响应体大小限制(limit=1000 会被拒绝),
        因此加总得到的反应总数是**下限估计**。若要精确的反应总数,
        应改用 ``reaction_report_total()``,它直接读取 meta.results.total。
        """
        search = f'patient.reaction.reactionmeddrapt:"{reaction}"'
        payload = self.fetch({"search": search, "count": field, "limit": limit})
        if not payload or "results" not in payload:
            return []

        total = payload.get("meta", {}).get("results", {}).get("total")
        out: List[CountResult] = []
        for row in payload["results"]:
            term = str(row.get("term", "")).strip()
            if not term:
                continue
            out.append(CountResult(term=term, count=int(row.get("count", 0)), total=total))
        return out

    def drug_report_total(self, drug: str) -> int:
        """某药物的报告总数 n1 = a + b。"""
        search = f'patient.drug.medicinalproduct:"{drug}"'
        payload = self.fetch({"search": search, "limit": 1})
        if not payload:
            return 0
        return int(payload.get("meta", {}).get("results", {}).get("total", 0))

    def reaction_report_total(self, reaction: str) -> int:
        """某反应的全部报告数 m1 = a + c(直接由 meta.total 得到,不受 limit 截断)。"""
        search = f'patient.reaction.reactionmeddrapt:"{reaction}"'
        payload = self.fetch({"search": search, "limit": 1})
        if not payload:
            return 0
        return int(payload.get("meta", {}).get("results", {}).get("total", 0))

    def grand_total(self, limit: int = 1000) -> int:
        """全库药物-反应报告对总数 N。

        通过对全部药物名做 count 聚合后求和得到,即
        N = sum_over_drugs(该药物的报告数)。

        实际观测:openFDA 对 ``count`` 聚合的响应体大小有限制,
        ``limit=1000`` 在该端点上返回 HTTP 403(而 limit<=800 正常)。
        因此这里按递减序列尝试,取第一个成功的 limit,并把实际用了多少项
        记录在 ``self.grand_total_limit`` 中,供报告说明口径。

        由于只对报告数最高的前若干药物名求和,该值是**下限估计**。
        """
        for lim in (min(limit, 800), 500, 200, 100):
            try:
                payload = self.fetch(
                    {"count": "patient.drug.medicinalproduct.exact", "limit": lim}
                )
            except OpenFDAError as e:
                # 403 属于"该 limit 不被接受",继续尝试更小的 limit
                self._log(f"limit={lim} 失败({e}),尝试更小的 limit")
                continue
            if payload and "results" in payload:
                self.grand_total_limit = lim
                return int(sum(int(r.get("count", 0)) for r in payload["results"]))
        return 0
