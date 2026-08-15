param(
    [switch]$Execute,
    [string]$Status = "CREATED_ONLY"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$JsonPath = Join-Path $Root "03_inventory_tables\rotating_case_inventory_latest.json"
if (-not (Test-Path -LiteralPath $JsonPath)) {
    throw "Inventory JSON not found. Run .\update_inventory.ps1 first."
}

$Rows = Get-Content -LiteralPath $JsonPath -Raw | ConvertFrom-Json
$Targets = @($Rows | Where-Object { $_.run_status -eq $Status })

$Stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$OutFile = Join-Path $Root ("submit_commands_{0}_{1}.txt" -f $Status, $Stamp)
$Commands = @()
foreach ($row in $Targets) {
    $run = [string]$row.run_path
    if ($run.Length -eq 0) { continue }
    $Commands += (
        "cd '$run' && " +
        "python3 ./check_drizzle_before_submit.py " +
        "--bou-in ./bou.in --profile ./drizzle_init.dat " +
        "--expected-perturb 1e-4 && " +
        "csub < subjob.sh"
    )
}

$Commands | Set-Content -LiteralPath $OutFile -Encoding UTF8
Write-Host ("Prepared {0} command(s): {1}" -f $Commands.Count, $OutFile)

if (-not $Execute) {
    Write-Host "Preview only. No job submitted."
    Write-Host "After checking the command file, run:"
    Write-Host ".\prepare_batch_submit.ps1 -Execute"
    exit 0
}

foreach ($cmd in $Commands) {
    Write-Host "Submitting: $cmd"
    ssh c01n0006 $cmd
}
