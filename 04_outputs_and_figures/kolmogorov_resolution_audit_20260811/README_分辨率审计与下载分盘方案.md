# 科莫格罗夫尺度审计与下载分盘方案 (2026-08-11)

## 1. 结论

当前不能把“之前一直处理的所有 case”统一判为满足 DNS 分辨率，也不应立即把所有远端 movie/HDF5 全量下载。

- 旧 `beta1` 初筛共 96 条：`resolved=23`、`marginal=10`、`under-resolved=27`、`not_checked=36`。
- 最近新程序的 14 个 `Ra=8e6` continuation 快照，在加入近壁点数和时间分辨率后：
  - `resolved=0`
  - `marginal=4`
  - `under-resolved=9`
  - `not_confirmed_nonrotating_BL=1`
- 小 Ek 的 `Ek=1.5e-4`、`2e-4` 和 `7e-4` 在当前单快照上通过耗散波数判据，但名义 Ekman 层内只有 1、1、3 个 z 网格点，不能宣称边界层已经解析。
- `Ek=1.5e-4` 和 `2e-4` 的 continuation 快照处于极低能量阶段，必须补取羽流爆发期三分量场，按最强耗散时段复查。
- 因严格条件尚未满足，本轮没有新增“全量 movie 下载”；已验证的 continuation、程序和诊断备份仍保存在迁移包中。

## 2. 论文依据

参考：Zhang et al. (2017), *Statistics of kinetic and thermal energy dissipation rates in two-dimensional turbulent Rayleigh-Benard convection*.

论文使用的动能耗散率为

\[
\epsilon_u=\frac{\nu}{2}\sum_{i,j}
\left(\frac{\partial u_j}{\partial x_i}+\frac{\partial u_i}{\partial x_j}\right)^2
=2\nu S_{ij}S_{ij}.
\]

科莫格罗夫尺度为

\[
\eta_K=\left(\frac{\nu^3}{\epsilon_u}\right)^{1/4}.
\]

对以 `Pr` 表示扩散率比的热标量，Batchelor 尺度为

\[
\eta_B=\frac{\eta_K}{\sqrt{Pr}}.
\]

该文报告的分辨率检查包括：

- `Pr=0.7` 时黏性边界层至少约 8 个点；
- 热边界层至少 10 个点；
- `Delta_g/eta_K <= 0.57`；
- `Delta_g/eta_B <= 0.48`；
- `Delta_t/tau_eta < 0.01`；
- 用全局精确耗散关系复核计算一致性。

注意：论文是二维、均匀网格、干 RB。当前三维旋转湿对流采用水平周期和非均匀 z 网格，因此同时保留两套判据：

1. 项目既有谱截断初筛：`kmax*eta_min >= 1`，近壁 `max(Delta z/eta) <= 1`。
2. Zhang 2017 的物理网格间距与边界层点数检查。

水汽和湿静能的最终 Batchelor 检查必须使用程序中的实际水汽扩散率/Schmidt 数；CSV 中的 `eta_B_min_Pr` 目前只是以 `Pr` 为扩散率比的热标量代理值。

## 3. 新程序 case 的当前状态

| case | 谱截断初筛 | Ekman层底部点数 | 时间分辨率 | 综合状态 |
| --- | --- | ---: | --- | --- |
| Ek1e-1 | marginal | 22 | pass | marginal |
| Ek1e-2 | marginal | 9 | fail | marginal |
| Ek1e-3 | under-resolved | 3 | fail | under-resolved |
| Ek1p5e-4 | resolved (低能快照) | 1 | pass | under-resolved |
| Ek2e-3 | under-resolved | 4 | fail | under-resolved |
| Ek2e-4 | resolved (低能快照) | 1 | pass | under-resolved |
| Ek3e-2 | marginal | 14 | fail | marginal |
| Ek3e-3 | under-resolved | 5 | fail | under-resolved |
| Ek5e-2 | marginal | 17 | marginal | marginal |
| Ek5e-3 | under-resolved | 6 | fail | under-resolved |
| Ek5e-4 | under-resolved | 2 | fail | under-resolved |
| Ek7e-3 | under-resolved | 7 | fail | under-resolved |
| Ek7e-4 | resolved (低能快照) | 3 | pass | under-resolved |
| norotating | marginal | 未定义 | marginal | 未确认非旋转黏性BL |

详细数值见 `kolmogorov_resolution_current_cases_restart.csv`。

## 4. 数据来源与限制

- 新程序快照来源：
  `C:\Users\jiasenzhang\Desktop\修改程序\rotating_moist_migration_bundle_20260808\03_current_cases\<case>\run`
- 使用文件：`continua_q1.h5`、`continua_q2.h5`、`continua_q3.h5`、`field_gridc.h5`、`bou.in`。
- 对重复的水平周期端点先去重，再把交错速度插值到单元中心。
- 水平导数用周期 Fourier 导数，垂直导数使用实际拉伸 z 网格。
- 同时计算 `2 nu Sij Sij` 和 `nu sum_j,i (partial_j u_i)^2`；全局相对差异约为千分量级，可作为插值和不可压缩一致性检查。
- 当前是 continuation 单快照筛查，不是稳态多帧最坏耗散统计。
- 旧表 `kolmogorov_resolution_all_cases_latest_field_clean.csv` 同样主要是最新单帧筛查。

## 5. 下一轮必须补的审计数据

优先只下载用于分辨率审计的少量完整三分量场，而不是先下载全部 movie：

1. 间歇 case：每个 case 至少取 3 个羽流爆发峰值时刻和 3 个低谷时刻。
2. 稳态 case：最新稳定区间等间隔取 5-10 帧。
3. 对每帧计算 `eta_K(z,t)`，最终使用时间窗口内最小 `eta_K` 和最大 `Delta/eta_K`。
4. 从速度/标量剖面定义实际边界层厚度，统计层内网格点数，不只使用 `sqrt(Ek)`。
5. 用实际水汽扩散率检查 q 和 m 的 Batchelor/Obukhov-Corrsin 小尺度。

## 6. 本机磁盘与 FileZilla 状态

| 盘符 | 可用空间 | 用途建议 |
| --- | ---: | --- |
| C: | 56.95 GB | 系统盘，只放程序、CSV、continuation 小备份 |
| D: | 26.36 GB | 不放原始 HDF5 |
| E: | 180.62 GB | 已有分析结果与 `Ra8e6` 本地数据，继续放 reduced data |
| F: | 109.78 GB | 放一个中等完整 case，例如 `Ek1e-3` |
| G: | 0 GB | 当前已满，禁止继续写入 |
| H: | 292.04 GB | 放一个大完整 case，例如 `Ek1e-1` |

现有主要数据：

- `G:\moist convection`: 511.43 GB
- `E:\moist RB\moist\result\data\ns\transition_study\rotating_self_aggregation\Ra8e6`: 110.31 GB
- migration bundle: 4.63 GB

FileZilla 当前队列约 `3356` 文件、`102.6 GiB`，目标仍包含 `G:`。队列未被本次审计清空或改写。

建议按完整 case 搬移，不拆散同一个 case：

- `G:\moist convection\Ek1e-1` (246.72 GB) -> `H:\rotating_moist_archive\Pr0p7\Ra8e6\Ek1e-1`
- `G:\moist convection\Ek1e-3` (67.17 GB) -> `F:\rotating_moist_archive\Pr0p7\Ra8e6\Ek1e-3`

完成复制、文件数/字节数校验和抽样 HDF5 打开验证后，才能删除 G: 原副本；随后 G: 可释放约 314 GB，足以容纳当前 102.6 GiB 队列。不要在 FileZilla 正在写这些目录时搬移。

## 7. 目录命名标准

新下载保持远端参数层级：

```text
<drive>:\rotating_case\Pr0p7\Ra8e6\<Ek-or-norotating>\AR16\Beta1p02\qbot0p5_qtop0p004978\N385x385x65\<segment>\run
```

同一 case 不跨盘。只有单个 case 大于任一盘可用空间时，才按完整子目录 `restart/`、`movie/`、`diagnostics/` 分盘，并在根目录保存映射 CSV。

## 8. 可重复运行

审计脚本：

`E:\moist RB\rotating_case_inventory\02_inventory_and_plot_scripts\check_kolmogorov_resolution_local.py`

输出：

- `kolmogorov_resolution_current_cases_restart.csv`
- `kolmogorov_resolution_current_cases_restart.md`
- `beta1_resolved_candidates_20260711.csv`

