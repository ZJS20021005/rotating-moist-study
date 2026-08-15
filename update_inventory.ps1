Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Script = Join-Path $Root "02_inventory_and_plot_scripts\update_rotating_case_inventory.py"
$BundledPython = Join-Path $env:USERPROFILE ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
$Python = if (Test-Path -LiteralPath $BundledPython) { $BundledPython } else { "python" }

& $Python $Script @args
