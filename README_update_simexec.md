# 更新已有 case 的 simexec 与 drizzle 工具

远端最新版可执行文件位于：

```text
/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_program/source/simexec
```

当前初始化采用“手动生成、DNS 读取”的两步流程。每个 `run` 目录应有：

```text
prepare_drizzle_initial_condition.sh
generate_drizzle_initial_condition.py
stability_solver.py
check_drizzle_before_submit.py
simexec
```

更新工具或 `simexec` 后，不要自动提交。对于每个要从 `nread=0`
开始的 case，先运行：

```bash
cd "该 case 的 run 目录"
./prepare_drizzle_initial_condition.sh
```

再执行：

```bash
csub < subjob.sh
```

`prepare_drizzle_initial_condition.sh` 每次都会重新读取当前
`bou.in`，求解并覆盖该 case 自己的 drizzle 基态。`simexec`
不会自行求解 drizzle，只读取并验证 `drizzle_init.dat`。

2026-07-28 的全量更新报告：

```text
/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/manual_drizzle_workflow_update_report_20260728.json
```
