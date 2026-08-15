$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
python "$ScriptDir\create_rotating_case.py" --config "$ScriptDir\case_config.json" --batch "$ScriptDir\batch_cases.json" @args
