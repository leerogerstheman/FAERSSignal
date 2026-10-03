"""FAERSSignal —— 基于 openFDA 数据的药物警戒不相称性分析工具包。

模块划分:
    openfda : openFDA REST API 客户端(带限流与重试)
    stats   : PRR / ROR / BCPNN(IC) 等不相称性统计量,自实现,无 scipy 依赖
    store   : SQLite 持久化
    report  : Markdown 报告生成
    cli     : 命令行入口
"""

__version__ = "1.0.0"
__author__ = "leerogerstheman"

__all__ = ["__version__", "__author__"]
