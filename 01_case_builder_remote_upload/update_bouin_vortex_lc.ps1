$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
python "$ScriptDir\update_existing_bouin_vortex_lc.py" --config "$ScriptDir\case_config.json" @args
