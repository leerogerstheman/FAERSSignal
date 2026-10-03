"""SQLite 持久化层。

存放三类数据:
    analysis_run : 每次分析的元信息(时间、参数、数据源)
    signal_result: 每行一个药物-反应组合的完整统计结果
    api_cache    : 可选的 API 响应缓存,便于离线复算

使用 WAL 与索引以支持后续按药物/反应/信号强弱查询。
"""

from __future__ import annotations

import json
import sqlite3
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Optional

from .stats import SignalResult

__all__ = ["SignalStore"]

SCHEMA = """
CREATE TABLE IF NOT EXISTS analysis_run (
    run_id        INTEGER PRIMARY KEY AUTOINCREMENT,
    started_at    TEXT NOT NULL,
    finished_at   TEXT,
    drugs         TEXT NOT NULL,
    n_reactions   INTEGER,
    grand_total   INTEGER,
    data_source   TEXT,
    stats_version TEXT,
    notes         TEXT
);

CREATE TABLE IF NOT EXISTS signal_result (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id        INTEGER NOT NULL REFERENCES analysis_run(run_id),
    drug          TEXT NOT NULL,
    reaction      TEXT NOT NULL,
    a             REAL NOT NULL,
    b             REAL NOT NULL,
    c             REAL NOT NULL,
    d             REAL NOT NULL,
    n_total       REAL NOT NULL,
    expected      REAL,
    prr           REAL,
    prr_ci_low    REAL,
    prr_ci_high   REAL,
    chi2          REAL,
    chi2_p        REAL,
    ror           REAL,
    ror_ci_low    REAL,
    ror_ci_high   REAL,
    ic            REAL,
    ic025         REAL,
    ic975         REAL,
    ic_var        REAL,
    corrected     INTEGER NOT NULL DEFAULT 0,
    is_signal_prr INTEGER NOT NULL DEFAULT 0,
    is_signal_ror INTEGER NOT NULL DEFAULT 0,
    is_signal_ic  INTEGER NOT NULL DEFAULT 0,
    is_signal     INTEGER NOT NULL DEFAULT 0
);

CREATE INDEX IF NOT EXISTS idx_signal_drug     ON signal_result(drug);
CREATE INDEX IF NOT EXISTS idx_signal_reaction ON signal_result(reaction);
CREATE INDEX IF NOT EXISTS idx_signal_flag     ON signal_result(is_signal);
CREATE INDEX IF NOT EXISTS idx_signal_run      ON signal_result(run_id);

CREATE TABLE IF NOT EXISTS api_cache (
    url        TEXT PRIMARY KEY,
    payload    TEXT NOT NULL,
    fetched_at TEXT NOT NULL
);
"""

_RESULT_COLUMNS = [
    "drug", "reaction", "a", "b", "c", "d", "n_total", "expected",
    "prr", "prr_ci_low", "prr_ci_high", "chi2", "chi2_p",
    "ror", "ror_ci_low", "ror_ci_high",
    "ic", "ic025", "ic975", "ic_var",
    "corrected", "is_signal_prr", "is_signal_ror", "is_signal_ic", "is_signal",
]


class SignalStore:
    """SQLite 存储封装。既可作为上下文管理器,也可手动 open/close。"""

    def __init__(self, db_path: str | Path) -> None:
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row
        self.conn.executescript(SCHEMA)
        self.conn.commit()

    # ------------------------------------------------------------------
    def __enter__(self) -> "SignalStore":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def close(self) -> None:
        try:
            self.conn.close()
        except sqlite3.Error:
            pass

    # ------------------------------------------------------------------
    def start_run(
        self,
        drugs: Iterable[str],
        n_reactions: int,
        grand_total: int,
        data_source: str = "openFDA drug/event API",
        stats_version: str = "1.0",
        notes: str = "",
    ) -> int:
        """写入一次分析的元信息,返回 run_id。"""
        cur = self.conn.execute(
            """INSERT INTO analysis_run
               (started_at, drugs, n_reactions, grand_total, data_source, stats_version, notes)
               VALUES (?, ?, ?, ?, ?, ?, ?)""",
            (
                datetime.now(timezone.utc).isoformat(timespec="seconds"),
                json.dumps(list(drugs), ensure_ascii=False),
                int(n_reactions),
                int(grand_total),
                data_source,
                stats_version,
                notes,
            ),
        )
        self.conn.commit()
        return int(cur.lastrowid)

    def finish_run(self, run_id: int) -> None:
        self.conn.execute(
            "UPDATE analysis_run SET finished_at = ? WHERE run_id = ?",
            (datetime.now(timezone.utc).isoformat(timespec="seconds"), run_id),
        )
        self.conn.commit()

    def save_results(self, run_id: int, results: Iterable[SignalResult]) -> int:
        """批量写入统计结果,返回写入行数。"""
        rows = []
        for r in results:
            d = r.as_dict()
            rows.append(tuple(
                int(v) if isinstance(v, bool) else v for v in (d[c] for c in _RESULT_COLUMNS)
            ))
        if not rows:
            return 0
        placeholders = ", ".join("?" * (len(_RESULT_COLUMNS) + 1))
        sql = (
            f"INSERT INTO signal_result (run_id, {', '.join(_RESULT_COLUMNS)}) "
            f"VALUES ({placeholders})"
        )
        self.conn.executemany(sql, [(run_id, *row) for row in rows])
        self.conn.commit()
        return len(rows)

    def latest_run_id(self) -> Optional[int]:
        row = self.conn.execute("SELECT MAX(run_id) AS rid FROM analysis_run").fetchone()
        return int(row["rid"]) if row and row["rid"] is not None else None

    def top_signals(self, run_id: Optional[int] = None, limit: int = 10) -> List[sqlite3.Row]:
        """按 PRR 降序取信号行,便于快速核查。"""
        run_id = run_id if run_id is not None else self.latest_run_id()
        return list(self.conn.execute(
            """SELECT * FROM signal_result
               WHERE run_id = ? AND is_signal = 1
               ORDER BY prr DESC, ic025 DESC LIMIT ?""",
            (run_id, limit),
        ))

    # ------------------------------------------------------------------
    def cache_api_response(self, url: str, payload: dict) -> None:
        self.conn.execute(
            "INSERT OR REPLACE INTO api_cache (url, payload, fetched_at) VALUES (?, ?, ?)",
            (url, json.dumps(payload, ensure_ascii=False),
             datetime.now(timezone.utc).isoformat(timespec="seconds")),
        )
        self.conn.commit()
