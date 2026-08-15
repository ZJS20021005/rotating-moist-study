param(
    [switch]$DryRun,
    [string]$Root = "",
    [string]$Simexec = ""
)

$remoteScript = "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/postprocess/strict_force_balance_20260805/prepare_strict_force_continuations.py"
$remoteArgs = @()
if ($DryRun) { $remoteArgs += "--dry-run" }
if ($Root -ne "") { $remoteArgs += "--root `"$Root`"" }
if ($Simexec -ne "") { $remoteArgs += "--simexec `"$Simexec`"" }

$command = "python3 `"$remoteScript`" " + ($remoteArgs -join " ")
ssh c01n0011 $command
