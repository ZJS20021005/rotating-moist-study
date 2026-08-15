param(
    [string]$BundleRoot = (Split-Path -Parent $PSScriptRoot)
)

$templatePath = Join-Path $PSScriptRoot "subjob_xh5.template.sh"
$configSource = Join-Path $PSScriptRoot "platform_config.xh5.sh"
$configTarget = Join-Path $PSScriptRoot "platform_config.sh"
$caseRoot = Join-Path $BundleRoot "03_current_cases"
$utf8NoBom = [System.Text.UTF8Encoding]::new($false)

if (-not (Test-Path -LiteralPath $templatePath)) {
    throw "Missing PBS template: $templatePath"
}

$template = [System.IO.File]::ReadAllText($templatePath)
$cases = Get-ChildItem -LiteralPath $caseRoot -Directory | Sort-Object Name
if ($cases.Count -ne 14) {
    throw "Expected 14 cases, found $($cases.Count)"
}

foreach ($case in $cases) {
    $runDir = Join-Path $case.FullName "run"
    $bouPath = Join-Path $runDir "bou.in"
    if (-not (Test-Path -LiteralPath $bouPath)) {
        throw "Missing bou.in: $bouPath"
    }

    $bou = [System.IO.File]::ReadAllText($bouPath)
    if ($bou -notmatch "(?ms)Restart flag.*?\r?\n.*?\r?\n\s*1\s+") {
        throw "Continuation NREAD is not 1: $bouPath"
    }

    $lines = $bou -split "\r?\n"
    $tmaxUpdated = $false
    for ($i = 0; $i -lt $lines.Count - 1; $i++) {
        if ($lines[$i] -match "NTST" -and $lines[$i] -match "TMAX") {
            $values = $lines[$i + 1] -split "\s+" | Where-Object { $_ -ne "" }
            if ($values.Count -lt 6) {
                throw "Invalid timing row in $bouPath"
            }
            $values[4] = "500d0"
            $lines[$i + 1] = $values -join "     "
            $tmaxUpdated = $true
            break
        }
    }
    if (-not $tmaxUpdated) {
        throw "Could not locate TMAX row in $bouPath"
    }
    $bou = ($lines -join "`n").TrimEnd() + "`n"
    [System.IO.File]::WriteAllText($bouPath, $bou, $utf8NoBom)

    $jobCase = $case.Name -replace "[^A-Za-z0-9]", ""
    $jobName = "Ra8e6Pr07${jobCase}C500"
    $subjob = $template.Replace("__JOB_NAME__", $jobName)
    [System.IO.File]::WriteAllText((Join-Path $runDir "subjob.sh"), $subjob, $utf8NoBom)
}

Copy-Item -LiteralPath $configSource -Destination $configTarget -Force
Write-Host "Prepared xh5 PBS/qsub subjob.sh for all 14 continuation cases."
Write-Host "Set TMAX=500d0 and retained NREAD=1 for all cases."
Write-Host "Installed xh5 platform_config.sh. No job was submitted."
