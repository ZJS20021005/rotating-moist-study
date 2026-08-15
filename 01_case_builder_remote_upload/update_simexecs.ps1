$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
python "$ScriptDir\update_existing_case_simexecs.py" --config "$ScriptDir\case_config.json" @args
