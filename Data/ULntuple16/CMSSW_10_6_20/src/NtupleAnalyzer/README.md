# NtupleAnalyzer 使用入口

本页只帮助人快速定位核心内容；当前物理状态、正式 tag 和待确认项以 `Data/docs/handoff.md` 为准。

## 正式结果

正式结果只从 `fit_results/_CURRENT/` 进入，不要根据历史目录名猜测 nominal：

- `total_fit/total/total_fit_summary.csv`：total 4D fit 中央值和统计误差。
- `systematics/cross_sections.csv`：正式 total/differential 截面及 systematic，当前采用 H019 recombination。
- `differential_cross_sections/`：35-bin nominal/reference fit 和 base fit table。
- `chi2/chi2_binning_scan.csv`：projection chi2 binning scan。
- `toy_gof/gof_summary.txt`：正式 toy GOF 摘要。
- `fitter_stability/run.metadata.txt`：正式 fitter-stability 配方和状态。

`differential_cross_sections/cross_sections.csv` 中的旧 systematic 栏只作 fit-stage provenance；正式 systematic 必须读取 `systematics/cross_sections.csv`。

## 核心源码

- `Fit_4D_tot.cpp`：total corrected-error 4D fit。
- `Fit_4D_diff.cpp`：differential 4D fit。
- `Fit_Check.cpp`：projection chi2 和 toy GOF worker。
- `Plot_4D.hpp`：4D fit 共用绘图。
- `scripts/collect/recombine_systematics.py`：正式 systematic recombination。

## 新运行约定

- 新 campaign 写入 `fit_results/runs/<stage>/<stage>-YYYYMMDD-rNN/`。
- 临时 attempt、job 和日志写入 campaign 内部或 `fit_results/work/`，不得新增顶层 job/retry 目录。
- tag 只表达 stage、日期和 revision；物理参数与运行细节写入 metadata。
- 只有用户确认且 marker、metadata、checksum 完整的结果才能进入 `_CURRENT`。

现有顶层 ROOT、日志、ACLiC 文件和长名历史目录为 legacy compatibility/provenance，不代表当前正式结果。
当前正式结论和下一步请继续查看 `Data/docs/handoff.md`。
