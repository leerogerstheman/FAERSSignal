#!/usr/bin/env python
"""FAERSSignal 入口脚本。

用法示例::

    python faers_signal.py --drugs ASPIRIN,IBUPROFEN --top-reactions 30
    python faers_signal.py --drugs METFORMIN --csv output/metformin.csv

把项目根目录加入 sys.path 后从 src 包导入,使脚本既可在项目根目录直接运行,
也可通过 ``python -m src.cli`` 运行。
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.cli import main  # noqa: E402  (必须在 sys.path 调整之后导入)

if __name__ == "__main__":
    raise SystemExit(main())
