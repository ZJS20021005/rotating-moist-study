Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Builder = Join-Path $Root "01_case_builder_remote_upload"

Push-Location $Builder
try {
    & ".\create_cases_batch.ps1" @args
}
finally {
    Pop-Location
}
