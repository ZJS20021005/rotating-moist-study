$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
if ($args.Count -eq 0) {
    python "$ScriptDir\create_rotating_case.py" --config "$ScriptDir\case_config.json" --interactive
} else {
    python "$ScriptDir\create_rotating_case.py" --config "$ScriptDir\case_config.json" @args
}
