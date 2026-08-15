Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Builder = Join-Path $Root "01_case_builder_remote_upload"

Push-Location $Builder
try {
    & ".\update_bouin_vortex_lc.ps1" @args
}
finally {
    Pop-Location
}
