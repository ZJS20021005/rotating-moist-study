# 14 个新程序 case：本地存储与下载状态（2026-08-12）

## 已完成

- 14 个目标 case 的主存储位置已经固定，任何一个 case 都不跨盘：
  - `Ek1e-1` 完整放在 `H:\rotating_case\Pr0p7\Ra8e6\Ek1e-1`。
  - 其余 13 个 case 完整放在 `G:\moist convection` 下各自目录。
  - `F:` 保持为空，作为备用盘。
- `Ek1e-1` 搬运完成：3021 个文件，136865412267 bytes（127.465848 GiB）。
- G 上原 `Ek1e-1` 源目录已经不存在；robocopy 退出码为 3（成功，且小于 8）。
- 抽检 5 个代表性 HDF5 文件均能由 `h5py` 打开并读取数据集。
- 两组已确认重复数据已删除，共释放约 119.329 GiB；`dry` 没有删除。
- E 盘散落的 `Ek1p5e-4/AR4` 四个非重复场文件已搬到 G 盘同一 case 下：
  `supplemental_fields_from_E_20260801`。

## FileZilla 队列

- 原队列：3356 项，102.584662 GiB。
- 保留目标队列：2697 项，74.690359 GiB，只包含：
  - `Ek1e-3`, `qbot=0.5`
  - `Ek1p5e-4`, `qbot=0.5`
- 已移除的非目标待下载项：659 项，27.894303 GiB，包括：
  - `Ek1e-3`, `qbot=1`
  - `Ek1e-4`, `qbot=0.5/1`
  - `Ek1p5e-3`, `qbot=1`
- 仅从待下载队列移除；没有删除这些参数在本地已经存在的文件。
- 原始队列数据库已备份为：
  `filezilla_queue_before_target_filter_20260811.sqlite3`。
- FileZilla 已启动“处理队列”，目前停在 `phfile.hnaicc.cn:2100` 密码输入框；输入密码并确认后即可继续。

## 明确保留

- `G:\moist convection\dry`：6.887135 GiB，按用户要求保留。
- E 盘 AR10 旧研究数据：约 2.033 GiB，未删除。
- E 盘 `nor\botsat`：约 0.509 GiB，未删除。

## 尚待决定的已下载非目标内容

- `G:\moist convection\Ek1e-2\AR16\Beta1p02\qbot1_qtop0p004978`：27.886185 GiB。
- 少量 qbot=1、`Ek1e-4`、`Ek1p5e-3` 配置文件与空目录。
- 这些内容不是已确认重复项，因此本轮未删除。
