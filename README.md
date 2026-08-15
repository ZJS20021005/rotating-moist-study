# rotating_case_inventory

这是旋转湿对流算例的本地管理入口。新建算例默认放在：

```text
/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case
```

## 新建单个算例

在 VS Code 中打开：

```text
E:\moist RB\rotating_case_inventory
```

运行：

```powershell
.\create_case.ps1
```

程序依次询问 `Ra、Pr、Ek、beta、AR、gamma、alpha、tau、qbot、
qtop、vortex_lc、n1、n2、n3、TMAX`。直接回车沿用上一次输入值。

也可以直接传参：

```powershell
.\create_case.ps1 --ra 8e6 --pr 0.7 --ek 3e-3 `
  --beta 1.02 --gamma 1.1 --alpha 3 --tau 1e-3 `
  --qbot 0.5 --qtop 0.004978 --aspect-ratio 16 `
  --n1 257 --n2 257 --n3 65 --time 2000
```

## Drizzle 初始条件：必须分两步

新建或更新后的每个 case，其 `run` 目录包含：

```text
prepare_drizzle_initial_condition.sh
generate_drizzle_initial_condition.py
stability_solver.py
check_drizzle_before_submit.py
simexec
```

每次开始一个新的 `nread=0` 计算前，先进入该 case 的 `run` 目录：

```bash
./prepare_drizzle_initial_condition.sh
```

该程序会：

1. 重新读取当前 `bou.in`；
2. 使用用户提供的 `stability_solver.moist_base_state` 求解该 case
   自己的一维 drizzle 基态；
3. 禁止在线性解 fallback 后继续；
4. 使用 `saturation_width=1e-8`，与 DNS 的
   `tanh(1e8*(q-qs))` 一致；
5. 生成并覆盖 `drizzle_init.dat` 与
   `drizzle_init_meta.json`；
6. 对参数、边界值和扰动幅值进行检查。

然后再运行或提交 `simexec`。在 `nread=0` 时，DNS 读取这个文件并设置：

```text
q(x,y,z) = q_drizzle(z)
b(x,y,z) = b_drizzle(z) + 1e-4*sin(pi*z/H)*N(0,1)
u = v = w = 0
```

因此，浮力扰动叠加在 drizzle 浮力基态上，而不是线性浮力剖面上。
水汽基态不加扰动。

若 `drizzle_init.dat` 与当前 `bou.in` 不一致，提交脚本中的检查和
DNS 启动检查都会停止计算。建 case 程序只放置工具和最新版
`simexec`，不会替你静默生成 drizzle，也不会提交作业。

## 提交

确认刚执行过 drizzle 生成程序后，可以在 case 根目录运行：

```bash
./submit_after_check.sh
```

或者在 `run` 目录中运行：

```bash
csub < subjob.sh
```

`subjob.sh` 会在启动 DNS 前再次检查 drizzle 文件，但不会重新求解。

## 批量创建

编辑：

```text
01_case_builder_remote_upload\batch_cases.json
```

执行：

```powershell
.\create_cases_batch.ps1
```

只预览、不创建：

```powershell
.\create_cases_batch.ps1 --dry-run
```

批量提交前先更新清单并生成命令：

```powershell
.\update_inventory.ps1
.\prepare_batch_submit.ps1
```

检查命令文件后再执行：

```powershell
.\prepare_batch_submit.ps1 -Execute
```

批量提交不会替代逐 case 的 drizzle 生成步骤。

## 参数映射

```text
invRo = sqrt(Pr/Ra)/Ek
norotating: invRo = 0

beta -> betaqs
dsalbot = 0
dsaltop = beta - 1
qvapbot = qbot
qvaptop = qtop

UBCBOT = 1  # 下壁 no-slip
UBCTOP = 0  # 上壁 free-slip
```

`vortex_lc` 由用户在 `bou.in` 中给定，不由 `k_peak` 自动换算。

## 主要文件

```text
00_latest_program\
  source\fluid_solver\inqpr.f90
  source\simexec

01_case_builder_remote_upload\
  create_rotating_case.py
  prepare_drizzle_initial_condition.sh
  generate_drizzle_initial_condition.py
  check_drizzle_before_submit.py
  linear_stability_reference\stability_solver.py

03_inventory_tables\
  rotating_case_inventory.xlsx
  rotating_case_inventory_latest.json
```
