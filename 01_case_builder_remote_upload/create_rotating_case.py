from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from itertools import product
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "case_config.json"
DEFAULT_BATCH_CONFIG = HERE / "batch_cases.json"
DEFAULT_INTERACTIVE_CONFIG = HERE / "last_interactive_config.json"
DRIZZLE_GENERATOR = HERE / "generate_drizzle_initial_condition.py"
DRIZZLE_CHECKER = HERE / "check_drizzle_before_submit.py"
DRIZZLE_PREPARE_SCRIPT = HERE / "prepare_drizzle_initial_condition.sh"
DRIZZLE_REFERENCE_SOLVER = (
    HERE / "linear_stability_reference" / "stability_solver.py"
)
CLI_OVERRIDE_KEYS = [
    "platform",
    "ssh_alias",
    "remote_root",
    "template_run",
    "program_bin",
    "ra",
    "pr",
    "ek",
    "gamma",
    "beta",
    "alpha",
    "tau",
    "qvapbot",
    "qvaptop",
    "n1",
    "n2",
    "n3",
    "aspect_ratio",
    "rext1",
    "rext2",
    "dt",
    "tmax",
    "nread",
    "vortex_lc",
    "continue_case",
    "continue_time",
]


REMOTE_SCRIPT = r"""
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


def parse_payload() -> dict:
    if "PAYLOAD_JSON" in globals():
        return json.loads(PAYLOAD_JSON)
    return json.loads(__import__("sys").stdin.read())


def fortran_float(value) -> str:
    v = float(value)
    if v == 0.0:
        return "0"
    av = abs(v)
    if av >= 1.0e4 or av < 1.0e-3:
        exp = int(math.floor(math.log10(av)))
        mant = v / (10.0 ** exp)
        s = f"{mant:.12g}".rstrip("0").rstrip(".")
        return f"{s}d{exp}"
    s = f"{v:.12g}"
    return f"{s}d0" if "e" not in s and "E" not in s else s.replace("e", "d").replace("E", "d")


def parse_fortran_float(value) -> float:
    return float(str(value).replace("D", "E").replace("d", "e"))


def inv_ro(ra: float, pr: float, ek) -> float:
    if ek is None or float(ek) == 0.0:
        return 0.0
    return math.sqrt(float(pr) / float(ra)) / float(ek)


def value_row(lines: list[str], token: str) -> tuple[int, list[str]]:
    for i, line in enumerate(lines[:-1]):
        if token in line:
            return i + 1, lines[i + 1].split()
    raise RuntimeError(f"Could not find bou.in header containing {token!r}")


def ensure_row_length(row: list[str], expected: int, defaults: list[str], token: str) -> list[str]:
    if len(defaults) != expected:
        raise RuntimeError(f"Internal defaults for {token!r} have wrong length")
    if len(row) == expected:
        return row
    if len(row) > expected:
        raise RuntimeError(
            f"Header {token!r}: expected {expected} values, got {len(row)}"
        )
    return row + defaults[len(row):]


def set_row(lines: list[str], token: str, values: list[str]) -> None:
    idx, old = value_row(lines, token)
    if len(old) != len(values):
        raise RuntimeError(
            f"Header {token!r}: expected {len(old)} values, got {len(values)}"
        )
    lines[idx] = "     ".join(values)


def set_job_name(text: str, job_name: str, comment: str, submit_style: str) -> str:
    if submit_style == "sbatch":
        line = f"#SBATCH -J {job_name}           # {comment}"
        pattern = r"(?m)^#SBATCH\s+-J\s+.*$"
    else:
        line = f"#CSUB -J {job_name}           # {comment}"
        pattern = r"(?m)^#CSUB\s+-J\s+.*$"
    updated, count = re.subn(pattern, line, text, count=1)
    if count != 1:
        raise RuntimeError("Expected exactly one #CSUB -J line in subjob.sh")
    return updated


def install_drizzle_job_check(text: str, perturb_amp: float) -> str:
    marker = "# DRIZZLE_INITIAL_CONDITION_CHECK"
    block = (
        marker
        + "\n"
        + "python3 ./check_drizzle_before_submit.py "
        + "--bou-in ./bou.in --profile ./drizzle_init.dat "
        + f"--expected-perturb {perturb_amp:.17e}\n"
    )
    if marker in text:
        start = text.index(marker)
        end = text.find("\n", start)
        if end < 0:
            end = len(text)
        second_end = text.find("\n", end + 1)
        if second_end < 0:
            second_end = len(text)
        return text[:start] + block + text[second_end + 1 :]

    anchors = [
        "rm -r fact flowmov data contstr movie",
        "mpirun -launcher fork",
        "mpirun",
        "impi-mpirun ./simexec",
    ]
    for anchor in anchors:
        position = text.find(anchor)
        if position >= 0:
            return text[:position] + block + text[position:]
    return text.rstrip() + "\n\n" + block


def chmod_exec(path: Path) -> None:
    if path.exists():
        path.chmod(path.stat().st_mode | 0o755)


def copy_file_if_present(src: Path, dst: Path, required: bool = False) -> None:
    if src.exists():
        shutil.copy2(src, dst)
    elif required:
        raise FileNotFoundError(str(src))


def prepare_case(cfg: dict) -> dict:
    remote_root = Path(cfg["remote_root"])
    template_run = Path(cfg["template_run"])
    program_bin = Path(cfg["program_bin"]) if cfg.get("program_bin") else None
    case_dir = Path(cfg["case_dir"])
    continue_case = bool(cfg.get("continue_case", False))
    continuation_number = None
    source_run = None
    if continue_case:
        if not case_dir.is_dir():
            raise RuntimeError(
                f"Cannot continue missing case directory: {case_dir}"
            )
        candidates = []
        existing_numbers = []
        for child in case_dir.iterdir():
            match = re.fullmatch(r"conti(\d+)", child.name, flags=re.IGNORECASE)
            if match and child.is_dir():
                number = int(match.group(1))
                existing_numbers.append(number)
                if (child / "run").is_dir():
                    candidates.append((number, child))
        continuation_number = max(existing_numbers, default=0) + 1
        source_case = max(candidates, default=None, key=lambda item: item[0])
        source_run = (source_case[1] / "run") if source_case else (case_dir / "run")
        if not source_run.is_dir():
            raise RuntimeError(f"Cannot find restart source run: {source_run}")
        run_dir = case_dir / f"conti{continuation_number}" / "run"
    else:
        run_dir = case_dir / "run"

    ra = float(cfg["ra"])
    pr = float(cfg["pr"])
    ek_value = cfg.get("ek")
    ek = None if ek_value in [None, "", "norotating", "NR", "no"] else float(ek_value)
    beta = float(cfg["beta"])
    qvaptop = float(cfg["qvaptop"])
    qvapbot = float(cfg["qvapbot"])
    drizzle_perturb_amp = 1.0e-4
    expected_dsaltop = beta - 1.0
    expected_dsalbot = 0.0
    if continue_case:
        cfg["nread"] = 1
        cfg["tmax"] = float(cfg.get("continue_time", cfg.get("tmax", 500.0)))

    dry_run = bool(cfg.get("dry_run", False))
    allow_existing = bool(cfg.get("allow_existing", False))
    if dry_run:
        return {
            "status": "dry_run",
            "case_dir": str(case_dir),
            "run_dir": str(run_dir),
            "continue_case": continue_case,
            "continuation_number": continuation_number,
            "job_name": cfg["job_name"],
            "computed_invRo": inv_ro(ra, pr, ek),
            "dsaltop": expected_dsaltop,
            "dsalbot": expected_dsalbot,
            "qvaptop": qvaptop,
            "qvapbot": qvapbot,
            "drizzle_perturb_amp": drizzle_perturb_amp,
            "vortex_lc": float(cfg.get("vortex_lc", 0.0)),
            "boundary_rule": "dsaltop=beta-1, dsalbot=0",
            "submit_command": (
                f"cd '{run_dir}' && sbatch subjob.sh"
                if cfg.get("submit_style") == "sbatch"
                else f"cd '{run_dir}' && csub < subjob.sh"
            ),
        }

    if not continue_case and not (template_run / "bou.in").is_file():
        raise RuntimeError(f"template_run is invalid: {template_run}")
    if program_bin is not None and not program_bin.is_file():
        raise RuntimeError(f"program_bin is missing: {program_bin}")

    if not continue_case and case_dir.exists() and not allow_existing:
        raise RuntimeError(
            f"Destination already exists: {case_dir}. "
            "Use --allow-existing only if you want to refresh editable files."
        )

    remote_root.mkdir(parents=True, exist_ok=True)
    case_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)

    input_run = source_run if continue_case else template_run
    for name in [
        "bou.in", "subjob.sh", "runsim.sh", "field_gridc.h5",
        "field_gridd.h5", "field_ms.xmf", "simexec",
        "drizzle_init.dat", "drizzle_init_meta.json",
        "generate_drizzle_initial_condition.py",
        "check_drizzle_before_submit.py",
        "prepare_drizzle_initial_condition.sh", "stability_solver.py",
    ]:
        copy_file_if_present(
            input_run / name,
            run_dir / name,
            required=name in {"bou.in", "subjob.sh"},
        )
    if continue_case:
        for restart_file in source_run.glob("continua_*"):
            if restart_file.is_file():
                shutil.copy2(restart_file, run_dir / restart_file.name)

    if program_bin is not None:
        shutil.copy2(program_bin, run_dir / "simexec")
        shutil.copy2(program_bin, run_dir / program_bin.name)
    else:
        copy_file_if_present(input_run / "simexec", run_dir / "simexec", required=True)

    if (template_run / "postprocess").is_dir() and not (run_dir / "postprocess").exists():
        shutil.copytree(template_run / "postprocess", run_dir / "postprocess")

    for folder in ["fact", "flowmov", "data", "contstr", "movie"]:
        (run_dir / folder).mkdir(exist_ok=True)

    drizzle_tools = {
        "generate_drizzle_initial_condition.py": cfg.get(
            "_drizzle_generator_source"
        ),
        "check_drizzle_before_submit.py": cfg.get("_drizzle_checker_source"),
        "prepare_drizzle_initial_condition.sh": cfg.get(
            "_drizzle_prepare_source"
        ),
        "stability_solver.py": cfg.get("_drizzle_reference_solver_source"),
    }
    missing_tools = [
        name for name, content in drizzle_tools.items() if not content
    ]
    if missing_tools:
        raise RuntimeError(
            "Missing required drizzle preparation tool(s): "
            + ", ".join(missing_tools)
        )
    for name, content in drizzle_tools.items():
        (run_dir / name).write_text(str(content), encoding="utf-8")

    bou_path = run_dir / "bou.in"
    lines = bou_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    grid_i, grid = value_row(lines, "N1      N2")
    grid = ensure_row_length(grid, 4, ["129", "129", "65", "3"], "N1      N2")
    grid[0:3] = [str(int(cfg["n1"])), str(int(cfg["n2"])), str(int(cfg["n3"]))]
    lines[grid_i] = "     ".join(grid)

    time_i, time_row = value_row(lines, "NTST      TRESTART")
    time_row = ensure_row_length(time_row, 6, ["5d5", "10d0", "0.1d0", "10d0", "800d0", "0"], "NTST      TRESTART")
    time_row[0] = cfg.get("ntst", time_row[0])
    time_row[1] = fortran_float(cfg.get("trestart", 10.0))
    time_row[2] = fortran_float(cfg.get("tpin", 0.1))
    time_row[3] = fortran_float(cfg.get("tframe", 10.0))
    time_row[4] = fortran_float(cfg.get("tmax", 800.0))
    time_row[5] = str(int(cfg.get("idtv", time_row[5])))
    lines[time_i] = "     ".join(time_row)

    restart_i, restart = value_row(lines, "NREAD     IRESET")
    restart = ensure_row_length(restart, 3, ["0", "1", "0"], "NREAD     IRESET")
    restart[0] = str(int(cfg.get("nread", 0)))
    restart[1] = str(int(cfg.get("ireset", restart[1])))
    restart[2] = str(int(cfg.get("pread", restart[2])))
    lines[restart_i] = "     ".join(restart)

    geo_i, geo = value_row(lines, "ALX3D     REXT1")
    geo = ensure_row_length(geo, 6, ["1d0", "10d0", "10d0", "1", "16d0", "1"], "ALX3D     REXT1")
    geo[0] = fortran_float(cfg.get("alx3d", 1.0))
    geo[1] = fortran_float(cfg.get("rext1", 10.0))
    geo[2] = fortran_float(cfg.get("rext2", 10.0))
    geo[3] = str(int(cfg.get("istr3", geo[3])))
    geo[4] = fortran_float(cfg.get("str3", 16.0))
    geo[5] = str(int(cfg.get("lmax", geo[5])))
    lines[geo_i] = "     ".join(geo)

    ubc_i, ubc = value_row(lines, "UBCBOT    UBCTOP")
    ubc = ensure_row_length(ubc, 2, ["1", "0"], "UBCBOT    UBCTOP")
    ubc[0] = str(int(cfg.get("ubcbot", ubc[0])))
    ubc[1] = str(int(cfg.get("ubctop", ubc[1])))
    lines[ubc_i] = "     ".join(ubc)

    ctl_i, ctl = value_row(lines, "Ra       Pr     invRo")
    ctl = ensure_row_length(ctl, 9, ["1d8", "0.7d0", "0", "0", "1.1", "1.0", "3.0", "1.0", "1d-3"], "Ra       Pr     invRo")
    ctl[0] = fortran_float(ra)
    ctl[1] = fortran_float(pr)
    ctl[2] = fortran_float(inv_ro(ra, pr, ek))
    ctl[3] = fortran_float(cfg.get("control_alpha", 0.0))
    ctl[4] = fortran_float(cfg["gamma"])
    ctl[5] = fortran_float(cfg.get("sm", 1.0))
    ctl[6] = fortran_float(cfg["alpha"])
    ctl[7] = fortran_float(beta)
    ctl[8] = fortran_float(cfg["tau"])
    lines[ctl_i] = "     ".join(ctl)

    mod_i, mod = value_row(lines, "A_stopmod")
    mod = ensure_row_length(mod, 8, ["0", "2d0", "0", "0", "0", "0", "0.004978d0", "1d0"], "A_stopmod")
    mod[0] = fortran_float(cfg.get("a_stopmod", 0.0))
    mod[1] = fortran_float(cfg.get("k_stopmod", 2.0))
    mod[2] = fortran_float(cfg.get("a_sbotmod", 0.0))
    mod[3] = fortran_float(cfg.get("k_sbotmod", 0.0))
    mod[4] = fortran_float(expected_dsaltop)
    mod[5] = fortran_float(expected_dsalbot)
    mod[6] = fortran_float(qvaptop)
    mod[7] = fortran_float(qvapbot)
    lines[mod_i] = "     ".join(mod)

    cfl_i, cfl = value_row(lines, "DTMAX(dt var.)")
    cfl = ensure_row_length(cfl, 5, ["7d-3", "5d0", "1.2d0", "1d-2", "5d0"], "DTMAX(dt var.)")
    cfl[0] = fortran_float(cfg.get("dtmax", cfg.get("dt", 0.007)))
    cfl[1] = fortran_float(cfg.get("resid", 5.0))
    cfl[2] = fortran_float(cfg.get("cflmax", 1.2))
    cfl[3] = fortran_float(cfg.get("dt", 0.007))
    cfl[4] = fortran_float(cfg.get("cfllim", 5.0))
    lines[cfl_i] = "     ".join(cfl)

    movie_i, movie = value_row(lines, "mov_zcut_k")
    movie = ensure_row_length(movie, 4, ["5", "2d2", "600", "0"], "mov_zcut_k")
    movie[0] = str(int(cfg.get("mov_zcut_k", movie[0])))
    movie[1] = fortran_float(cfg.get("tframe_me", movie[1].replace("d", "e")))
    movie[2] = str(int(cfg.get("stat_me", movie[2])))
    movie[3] = fortran_float(cfg.get("vortex_lc", 0.0))
    lines[movie_i] = "     ".join(movie)

    bou_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    written_lines = bou_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    for token, expected in [
        ("N1      N2", 4),
        ("NTST      TRESTART", 6),
        ("NREAD     IRESET", 3),
        ("ALX3D     REXT1", 6),
        ("UBCBOT    UBCTOP", 2),
        ("Ra       Pr     invRo", 9),
        ("A_stopmod", 8),
        ("DTMAX(dt var.)", 5),
        ("mov_zcut_k", 4),
    ]:
        _, written = value_row(written_lines, token)
        if len(written) != expected:
            raise RuntimeError(
                f"Written bou.in row {token!r} has {len(written)} values, expected {expected}"
            )
    _, written_ubc = value_row(written_lines, "UBCBOT    UBCTOP")
    if written_ubc != [str(int(cfg.get("ubcbot", 1))), str(int(cfg.get("ubctop", 0)))]:
        raise RuntimeError(
            "Velocity boundary check failed: wrote "
            + " ".join(written_ubc)
            + f", expected {int(cfg.get('ubcbot', 1))} {int(cfg.get('ubctop', 0))}"
        )
    _, written_mod = value_row(written_lines, "A_stopmod")
    dsaltop = parse_fortran_float(written_mod[4])
    dsalbot = parse_fortran_float(written_mod[5])
    qvaptop_written = parse_fortran_float(written_mod[6])
    qvapbot_written = parse_fortran_float(written_mod[7])
    checks = [
        ("dsaltop", dsaltop, expected_dsaltop),
        ("dsalbot", dsalbot, expected_dsalbot),
        ("qvaptop", qvaptop_written, qvaptop),
        ("qvapbot", qvapbot_written, qvapbot),
    ]
    for name, actual, expected in checks:
        if abs(actual - expected) > 1.0e-12 * max(1.0, abs(expected)):
            raise RuntimeError(
                f"Boundary check failed for {name}: wrote {actual}, expected {expected}"
            )

    job_name = cfg["job_name"]
    comment = cfg["job_comment"]
    subjob_path = run_dir / "subjob.sh"
    subjob_text = set_job_name(
        subjob_path.read_text(encoding="utf-8", errors="ignore"),
        job_name,
        comment,
        cfg.get("submit_style", "csub"),
    )
    subjob_text = install_drizzle_job_check(
        subjob_text, drizzle_perturb_amp
    )
    subjob_path.write_text(subjob_text, encoding="utf-8")

    chmod_exec(run_dir / "simexec")
    chmod_exec(run_dir / "subjob.sh")
    chmod_exec(run_dir / "runsim.sh")
    chmod_exec(run_dir / "generate_drizzle_initial_condition.py")
    chmod_exec(run_dir / "check_drizzle_before_submit.py")
    chmod_exec(run_dir / "prepare_drizzle_initial_condition.sh")

    info = {
        key: value
        for key, value in cfg.items()
        if not key.startswith("_drizzle_")
    }
    info["computed_invRo"] = inv_ro(ra, pr, ek)
    info["dsaltop"] = dsaltop
    info["dsalbot"] = dsalbot
    info["qvaptop"] = qvaptop_written
    info["qvapbot"] = qvapbot_written
    info["boundary_row"] = {
        "A_stopmod": written_mod[0],
        "k_stopmod": written_mod[1],
        "A_sbotmod": written_mod[2],
        "k_sbotmod": written_mod[3],
        "dsaltop": written_mod[4],
        "dsalbot": written_mod[5],
        "qvaptop": written_mod[6],
        "qvapbot": written_mod[7],
    }
    submit_command = (
        f"cd '{run_dir}' && sbatch subjob.sh"
        if cfg.get("submit_style") == "sbatch"
        else f"cd '{run_dir}' && csub < subjob.sh"
    )
    info["submit_command"] = submit_command
    info["continue_case"] = continue_case
    info["continuation_number"] = continuation_number
    info["restart_source_run"] = str(source_run) if source_run else None
    info_path = run_dir.parent / "case_info.json" if continue_case else case_dir / "case_info.json"
    info_path.write_text(json.dumps(info, indent=2), encoding="utf-8")
    helper_path = run_dir.parent / "submit_after_check.sh" if continue_case else case_dir / "submit_after_check.sh"
    submit_line = "sbatch subjob.sh\n" if cfg.get("submit_style") == "sbatch" else "csub < subjob.sh\n"
    helper_path.write_text(
        "#!/bin/bash\n"
        "set -euo pipefail\n"
        f"cd {repr(str(run_dir))}\n"
        "python3 ./check_drizzle_before_submit.py "
        "--bou-in ./bou.in --profile ./drizzle_init.dat "
        f"--expected-perturb {drizzle_perturb_amp:.17e}\n"
        + submit_line,
        encoding="utf-8",
    )
    chmod_exec(helper_path)

    return {
        "status": "prepared",
        "case_dir": str(case_dir),
        "run_dir": str(run_dir),
        "continue_case": continue_case,
        "continuation_number": continuation_number,
        "job_name": job_name,
        "computed_invRo": inv_ro(ra, pr, ek),
        "dsaltop": dsaltop,
        "dsalbot": dsalbot,
        "qvaptop": qvaptop_written,
        "qvapbot": qvapbot_written,
        "drizzle_profile": str(run_dir / "drizzle_init.dat"),
        "drizzle_metadata": str(run_dir / "drizzle_init_meta.json"),
        "drizzle_perturb_amp": drizzle_perturb_amp,
        "prepare_drizzle_command": (
            f"cd '{run_dir}' && ./prepare_drizzle_initial_condition.sh"
        ),
        "boundary_row": info["boundary_row"],
        "submit_command": submit_command,
    }


cfg = parse_payload()
result = prepare_case(cfg)
print(json.dumps(result, indent=2))
"""


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def load_interactive_defaults(config_path: Path) -> dict[str, Any]:
    cfg = load_json(config_path)
    if DEFAULT_INTERACTIVE_CONFIG.exists():
        saved = load_json(DEFAULT_INTERACTIVE_CONFIG)
        cfg.update(saved)
    return cfg


PLATFORM_DEFAULTS: dict[str, dict[str, Any]] = {
    "hainan": {
        "ssh_alias": "c01n0006",
        "remote_root": "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case",
        "template_run": "/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study/Ra1e8/AR10/norotating/norotating/run",
        "program_bin": "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_program/source/simexec",
        "submit_style": "csub",
    },
    "shuguang": {
        "ssh_alias": "xh5",
        "remote_root": "/work/home/jiasenzhang/rotating_case",
        "template_run": "/work/home/jiasenzhang/rotating_case/template/run",
        "program_bin": "/work/home/jiasenzhang/rotating_case/latest_program/source/simexec",
        "submit_style": "sbatch",
    },
}


def apply_platform_defaults(cfg: dict[str, Any]) -> dict[str, Any]:
    platform = str(cfg.get("platform", "hainan")).strip().lower()
    aliases = {"hna": "hainan", "海南": "hainan", "xh5": "shuguang", "曙光": "shuguang"}
    platform = aliases.get(platform, platform)
    if platform not in PLATFORM_DEFAULTS:
        raise SystemExit(
            f"Unknown platform {platform!r}; choose hainan or shuguang."
        )
    cfg["platform"] = platform
    for key, value in PLATFORM_DEFAULTS[platform].items():
        cfg[key] = value
    return cfg


def save_interactive_defaults(cfg: dict[str, Any]) -> None:
    skip_keys = {
        "case_dir",
        "run_dir",
        "job_name",
        "job_comment",
        "dry_run",
        "allow_existing",
        "computed_invRo",
        "boundary_row",
        "submit_command",
        "batch_index",
        "case_dir",
        "run_dir",
    }
    clean = {key: value for key, value in cfg.items() if key not in skip_keys}
    clean = {
        key: value for key, value in clean.items() if not key.startswith("_")
    }
    DEFAULT_INTERACTIVE_CONFIG.write_text(
        json.dumps(clean, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def prompt_value(label: str, current: Any, *, cast=str, aliases: dict[str, Any] | None = None) -> Any:
    aliases = aliases or {}
    if current is None:
        shown = ""
    else:
        shown = str(current)
    raw = input(f"{label} [{shown}]: ").strip()
    if raw == "":
        return current
    lowered = raw.lower()
    if lowered in aliases:
        return aliases[lowered]
    try:
        return cast(raw)
    except Exception as exc:
        raise SystemExit(f"无法解析 {label}: {raw!r} ({exc})")


def prompt_interactive_cfg(cfg: dict[str, Any]) -> dict[str, Any]:
    print("")
    print("Input case parameters. Press Enter to use the default in brackets.")
    print("For platform, use hainan or shuguang/xh5.")
    print("For Ek, use norotating / NR / 0 for a nonrotating case.")
    print("For continuation, enter y to create contiN without overwriting older runs.")
    print("")
    selected_platform = prompt_value(
        "platform",
        cfg.get("platform", "hainan"),
        cast=str,
        aliases={"hna": "hainan", "海南": "hainan", "xh5": "shuguang", "曙光": "shuguang"},
    )
    cfg["platform"] = selected_platform
    cfg = apply_platform_defaults(cfg)
    cfg["ra"] = prompt_value("Ra", cfg.get("ra"), cast=float)
    cfg["pr"] = prompt_value("Pr", cfg.get("pr"), cast=float)
    cfg["ek"] = prompt_value(
        "Ek",
        cfg.get("ek", "norotating"),
        cast=str,
        aliases={"0": "norotating", "nr": "norotating", "nonrotating": "norotating", "no": "norotating"},
    )
    cfg["beta"] = prompt_value("beta", cfg.get("beta"), cast=float)
    cfg["aspect_ratio"] = prompt_value(
        "AR / aspect_ratio",
        cfg.get("aspect_ratio", cfg.get("rext1", 10.0)),
        cast=float,
    )
    cfg["gamma"] = prompt_value("gamma", cfg.get("gamma"), cast=float)
    cfg["alpha"] = prompt_value("alpha / alphaqs", cfg.get("alpha"), cast=float)
    cfg["tau"] = prompt_value("tau / tau_cond", cfg.get("tau"), cast=float)
    cfg["qvapbot"] = prompt_value("qbot", cfg.get("qvapbot"), cast=float)
    cfg["qvaptop"] = prompt_value("qtop", cfg.get("qvaptop"), cast=float)
    cfg["vortex_lc"] = prompt_value("vortex_lc", cfg.get("vortex_lc", 0.0), cast=float)
    cfg["n1"] = prompt_value("n1", cfg.get("n1"), cast=int)
    cfg["n2"] = prompt_value("n2", cfg.get("n2"), cast=int)
    cfg["n3"] = prompt_value("n3", cfg.get("n3"), cast=int)
    cfg["tmax"] = prompt_value("time / TMAX", cfg.get("tmax"), cast=float)
    cfg["continue_case"] = prompt_value(
        "continue previous case? (y/n)",
        cfg.get("continue_case", False),
        cast=lambda value: value.lower() in {"y", "yes", "1", "true"},
    )
    if cfg["continue_case"]:
        cfg["continue_time"] = prompt_value(
            "continuation time",
            cfg.get("continue_time", 500.0),
            cast=float,
        )
    print("")
    return cfg


def label_number(value: float, *, beta: bool = False, sci: bool = False) -> str:
    if beta:
        return f"{float(value):.2f}".replace(".", "p")
    v = float(value)
    if sci and v != 0.0:
        s = f"{v:.6e}"
        mant, exp = s.split("e")
        mant = mant.rstrip("0").rstrip(".")
        exp_int = int(exp)
        return f"{mant}e{exp_int}".replace(".", "p")
    s = f"{v:.8g}"
    return s.replace("-", "m").replace(".", "p").replace("+", "")


def normalize_ek(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"", "none", "no", "nr", "norotating", "nonrotating", "0"}:
            return None
        return float(text)
    if float(value) == 0.0:
        return None
    return float(value)


def build_case_layout(cfg: dict[str, Any]) -> dict[str, str]:
    ra = float(cfg["ra"])
    pr = float(cfg["pr"])
    ek = normalize_ek(cfg.get("ek"))
    beta_value = float(cfg["beta"])
    qbot = float(cfg["qvapbot"])
    qtop = float(cfg["qvaptop"])
    n1, n2, n3 = int(cfg["n1"]), int(cfg["n2"]), int(cfg["n3"])
    ar = float(cfg.get("aspect_ratio", cfg.get("rext1", 10.0)))

    pr_layer = "Pr" + label_number(pr)
    ra_layer = "Ra" + label_number(ra, sci=True)
    ek_layer = "norotating" if ek is None else "Ek" + label_number(ek, sci=True)
    ar_layer = "AR" + label_number(ar)
    beta_layer = "Beta" + label_number(beta_value, beta=True)
    moisture_layer = "qbot" + label_number(qbot) + "_qtop" + label_number(qtop)
    grid_layer = f"N{n1}x{n2}x{n3}"

    short_ek = "NR" if ek is None else "Ek" + label_number(ek, sci=True)
    job_name = (
        "Ra"
        + label_number(ra, sci=True)
        + "Pr"
        + label_number(pr)
        + short_ek
        + "B"
        + label_number(beta_value, beta=True)
        + f"N{n1}"
    )
    job_name = re_safe_job_name(job_name)[:80]

    remote_root = cfg["remote_root"].rstrip("/")
    case_dir = "/".join(
        [
            remote_root,
            pr_layer,
            ra_layer,
            ek_layer,
            ar_layer,
            beta_layer,
            moisture_layer,
            grid_layer,
        ]
    )
    return {
        "pr_layer": pr_layer,
        "ra_layer": ra_layer,
        "ek_layer": ek_layer,
        "ar_layer": ar_layer,
        "beta_layer": beta_layer,
        "moisture_layer": moisture_layer,
        "grid_layer": grid_layer,
        "case_dir": case_dir,
        "job_name": job_name,
    }


def re_safe_job_name(text: str) -> str:
    keep = []
    for ch in text:
        keep.append(ch if ch.isalnum() or ch in {"_", "-"} else "_")
    return "".join(keep)


def apply_cli_overrides(cfg: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    for key in CLI_OVERRIDE_KEYS:
        value = getattr(args, key, None)
        if value is not None:
            cfg[key] = value
    cfg["dry_run"] = bool(args.dry_run)
    cfg["allow_existing"] = bool(args.allow_existing)
    return cfg


def validate_cfg(cfg: dict[str, Any]) -> None:
    required = [
        "ssh_alias",
        "remote_root",
        "template_run",
        "ra",
        "pr",
        "gamma",
        "beta",
        "alpha",
        "tau",
        "qvapbot",
        "qvaptop",
        "n1",
        "n2",
        "n3",
    ]
    missing = [name for name in required if cfg.get(name) in [None, ""]]
    if missing:
        raise SystemExit("Missing required config keys: " + ", ".join(missing))


def expanded_case_specs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    list_keys = [key for key, value in spec.items() if isinstance(value, list)]
    if not list_keys:
        return [dict(spec)]
    fixed = {key: value for key, value in spec.items() if key not in list_keys}
    expanded = []
    for values in product(*(spec[key] for key in list_keys)):
        row = dict(fixed)
        row.update(dict(zip(list_keys, values)))
        expanded.append(row)
    return expanded


def finalize_cfg(cfg: dict[str, Any]) -> dict[str, Any]:
    cfg = apply_platform_defaults(cfg)
    validate_cfg(cfg)
    ek = normalize_ek(cfg.get("ek"))
    cfg["ek"] = None if ek is None else ek
    if cfg.get("aspect_ratio") is not None:
        cfg["rext1"] = float(cfg["aspect_ratio"])
        cfg["rext2"] = float(cfg["aspect_ratio"])
    layout = build_case_layout(cfg)
    cfg.update(layout)
    cfg["job_comment"] = (
        f"{layout['ra_layer']}, {layout['pr_layer']}, {layout['ek_layer']}, "
        f"{layout['ar_layer']}, {layout['beta_layer']}, "
        f"{layout['moisture_layer']}, {layout['grid_layer']}"
    )
    return cfg


def attach_drizzle_tools(cfg: dict[str, Any]) -> dict[str, Any]:
    for required_path in [
        DRIZZLE_GENERATOR,
        DRIZZLE_CHECKER,
        DRIZZLE_PREPARE_SCRIPT,
        DRIZZLE_REFERENCE_SOLVER,
    ]:
        if not required_path.is_file():
            raise SystemExit(f"Missing drizzle workflow file: {required_path}")

    cfg["_drizzle_generator_source"] = DRIZZLE_GENERATOR.read_text(
        encoding="utf-8"
    )
    cfg["_drizzle_checker_source"] = DRIZZLE_CHECKER.read_text(
        encoding="utf-8"
    )
    cfg["_drizzle_prepare_source"] = DRIZZLE_PREPARE_SCRIPT.read_text(
        encoding="utf-8"
    )
    cfg["_drizzle_reference_solver_source"] = (
        DRIZZLE_REFERENCE_SOLVER.read_text(encoding="utf-8")
    )
    return cfg


def run_remote(cfg: dict[str, Any]) -> dict[str, Any]:
    payload = json.dumps(cfg)
    bootstrap = f"PAYLOAD_JSON = {payload!r}\n" + REMOTE_SCRIPT
    process = subprocess.run(
        ["ssh", str(cfg["ssh_alias"]), "python3 -"],
        input=bootstrap,
        text=True,
        capture_output=True,
        check=False,
        encoding="utf-8",
        errors="replace",
    )
    if process.returncode != 0:
        sys.stderr.write(process.stdout)
        sys.stderr.write(process.stderr)
        raise SystemExit(process.returncode)
    print(process.stdout.strip())
    return json.loads(process.stdout)


def run_single(args: argparse.Namespace) -> None:
    cfg = load_json(args.config)
    cfg = apply_platform_defaults(cfg)
    cfg = apply_cli_overrides(cfg, args)
    cfg = finalize_cfg(cfg)
    if not cfg.get("dry_run", False):
        cfg = attach_drizzle_tools(cfg)

    print("Local plan:")
    print(json.dumps({k: cfg[k] for k in ["ssh_alias", "case_dir", "job_name"]}, indent=2))
    run_remote(cfg)


def run_interactive(args: argparse.Namespace) -> None:
    cfg = load_interactive_defaults(args.config)
    cfg = prompt_interactive_cfg(cfg)
    cfg = apply_cli_overrides(cfg, args)
    cfg = finalize_cfg(cfg)
    if not cfg.get("dry_run", False):
        cfg = attach_drizzle_tools(cfg)

    print("Local plan:")
    print(json.dumps({
        "ssh_alias": cfg["ssh_alias"],
        "platform": cfg["platform"],
        "case_dir": cfg["case_dir"],
        "job_name": cfg["job_name"],
        "tmax": cfg["tmax"],
        "qvapbot": cfg["qvapbot"],
        "qvaptop": cfg["qvaptop"],
        "continue_case": cfg.get("continue_case", False),
        "continue_time": cfg.get("continue_time"),
    }, indent=2))
    run_remote(cfg)
    if not cfg.get("dry_run", False):
        save_interactive_defaults(cfg)
        print(f"Saved defaults: {DEFAULT_INTERACTIVE_CONFIG}")


def run_batch(args: argparse.Namespace) -> None:
    base_cfg = load_json(args.config)
    base_cfg = apply_platform_defaults(base_cfg)
    batch = load_json(args.batch)
    defaults = batch.get("defaults", {})
    case_specs = batch.get("cases", [])
    if not isinstance(case_specs, list) or not case_specs:
        raise SystemExit("Batch file must contain a non-empty 'cases' list.")

    prepared_cfgs: list[dict[str, Any]] = []
    for index, raw_spec in enumerate(case_specs, start=1):
        for expanded in expanded_case_specs(raw_spec):
            cfg = dict(base_cfg)
            cfg.update(defaults)
            cfg.update(expanded)
            cfg = apply_platform_defaults(cfg)
            cfg = apply_cli_overrides(cfg, args)
            cfg = finalize_cfg(cfg)
            if not cfg.get("dry_run", False):
                cfg = attach_drizzle_tools(cfg)
            cfg["batch_index"] = index
            prepared_cfgs.append(cfg)

    print(f"Batch plan: {len(prepared_cfgs)} case(s)")
    results = []
    for index, cfg in enumerate(prepared_cfgs, start=1):
        print(
            f"[{index}/{len(prepared_cfgs)}] {cfg['case_dir']} "
            f"job={cfg['job_name']}"
        )
        results.append(run_remote(cfg))

    summary_path = HERE / "last_batch_results.json"
    summary_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"Batch result summary: {summary_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare a rotating moist RB case on the remote cluster."
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--batch",
        type=Path,
        help="Create a batch from a JSON file. Defaults are still read from --config.",
    )
    parser.add_argument("--ssh-alias")
    parser.add_argument("--platform", choices=["hainan", "shuguang", "xh5"])
    parser.add_argument("--remote-root")
    parser.add_argument("--template-run")
    parser.add_argument("--program-bin")
    parser.add_argument("--ra", type=float)
    parser.add_argument("--pr", type=float)
    parser.add_argument("--ek", help="Use 'norotating' or 0 for invRo=0.")
    parser.add_argument("--gamma", type=float)
    parser.add_argument("--beta", type=float)
    parser.add_argument("--alpha", type=float, help="Mapped to alphaqs in bou.in.")
    parser.add_argument("--tau", type=float)
    parser.add_argument("--aspect-ratio", type=float)
    parser.add_argument("--qvapbot", "--qbot", dest="qvapbot", type=float)
    parser.add_argument("--qvaptop", "--qtop", dest="qvaptop", type=float)
    parser.add_argument("--n1", type=int)
    parser.add_argument("--n2", type=int)
    parser.add_argument("--n3", type=int)
    parser.add_argument("--rext1", type=float)
    parser.add_argument("--rext2", type=float)
    parser.add_argument("--dt", type=float)
    parser.add_argument("--tmax", "--time", dest="tmax", type=float)
    parser.add_argument("--nread", type=int)
    parser.add_argument("--vortex-lc", dest="vortex_lc", type=float)
    parser.add_argument("--continue-case", action="store_true", default=None, dest="continue_case")
    parser.add_argument("--continue-time", type=float, dest="continue_time")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--allow-existing", action="store_true")
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        help="Prompt for Ra, Pr, Ek, beta, AR, gamma, alpha, tau, qbot, qtop, n1, n2, n3, and TMAX. Empty input reuses the last value.",
    )
    args = parser.parse_args()

    if args.batch:
        if args.interactive:
            raise SystemExit("--interactive cannot be combined with --batch.")
        run_batch(args)
    elif args.interactive:
        run_interactive(args)
    else:
        run_single(args)


if __name__ == "__main__":
    main()
