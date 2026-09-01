# Data-driven SPS/DPS analysis

本目录是当前 data-driven nominal `f_DPS` 路线的独立代码区。它由历史
`tests/atlas_data_driven_sps/` 中 pre-updates 实际使用的功能选择性提取而来，
不复制旧大型 ROOT/PDF/PNG、ACLiC 产物或 bootstrap replica。旧
`tests/atlas_data_driven_sps/` 按用户指示暂时保留，但本目录的运行入口不依赖它；
待新链完整完成后由用户手动清理。

主要功能：

- 从统一 selection 构建 pair Data 和两个 single-J/psi candidate trees；
- 使用冻结的 4D model 做 single-J sPlot 和 12-cell prompt-prompt fits；
- event mixing、CR normalization 和 `f_DPS` 提取；
- fixed-seed full source-event bootstrap；
- single-J category、CR definition/purity 和 pure-DPS closure systematics；
- 四个一维 projection cross-check；
- `sigma_DPS`、`sigma_eff` 确定性传播及 pre-updates 图表。

当前状态是 `migration_complete_needs_full_root640_rerun`。旧 pre-updates 数值只作
regression reference；新的 nominal fitting model 改为 ROOT 6.40.02 corrected-error
workspace 后，central、bootstrap、systematics 和下游截面必须全部重跑并重新确认。

开始工作前阅读 `AGENTS.md`、`docs/handoff.md` 和
`inputs/input_manifest.txt`。所有新结果写入 `results/<unique-tag>/`，不得覆盖。

常用入口：

```bash
# 必须先确认 root-config --version 为 6.40.02
bash run_nominal.sh <unique-central-tag>
bash run_build_dps_reference.sh <unique-dps-reference-tag>
bash run_fdps_global_nominal_systematics.sh \
  <confirmed-central-tag> <dps-reference-tag> <unique-systematics-tag>
bash submit_full_bootstrap_fixed500_global.sh \
  <confirmed-central-tag> <unique-bootstrap-tag>
bash run_postprocess.sh <completed-systematics-tag> <unique-postprocess-tag>

# 只复现旧 pre-updates regression reference
bash presentation/scripts/run_all.sh
python3 postprocess/scripts/calculate_results.py
python3 postprocess/scripts/plot_figures.py
```

H015 F04 是唯一 runtime 例外：脚本自动进入 EL7/ROOT 6.14/09，仅从旧保存的
workspace/fit result 重绘投影，不执行 fit。其余新分析链固定为 ROOT 6.40.02。
