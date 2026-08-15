$ErrorActionPreference = "Stop"

$repo = Split-Path -Parent $MyInvocation.MyCommand.Path
$codexHome = if ($env:CODEX_HOME) { $env:CODEX_HOME } else { Join-Path $HOME "CodexHome" }
$skillRoot = Join-Path $codexHome "skills"

New-Item -ItemType Directory -Force -Path $skillRoot | Out-Null

foreach ($name in @("rotating-moist-study", "jiasen-scientific-plot-style")) {
    $source = Join-Path $repo "05_source_reference\codex_skills\$name"
    $target = Join-Path $skillRoot $name
    if (-not (Test-Path -LiteralPath $source)) {
        throw "Missing skill mirror: $source"
    }
    if (Test-Path -LiteralPath $target) {
        Remove-Item -LiteralPath $target -Recurse -Force
    }
    Copy-Item -LiteralPath $source -Destination $target -Recurse -Force
    Write-Output "Installed $name -> $target"
}

