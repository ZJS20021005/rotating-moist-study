#!/usr/bin/env python3
"""Stage unsubmitted 500-time continuation runs for the Ra=8e6 study.

This script is intended to run on c01n0011.  It finds the latest existing
qbot=0.5, beta=1.02 run for each requested Ek, copies restart files and
metadata, installs the newest pressure-enabled simexec, changes bou.in to
NREAD=1 and TMAX=500d0, and writes a continuation subjob.sh.

It never calls csub.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
from pathlib import Path


DEFAULT_ROOT = Path(
    "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/Pr0p7/Ra8e6"
)
DEFAULT_SIMEXEC = Path(
    "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_program/source/simexec"
)

# Keep the AR4 weak-rotation runs separate; the rest of this requested set is
# the established AR16 qbot=0.5 family.
TARGETS = {
    "1p5e-4": "AR4",
    "2e-4": "AR4",
    "5e-4": "AR16",
    "7e-4": "AR16",
    "1e-3": "AR16",
    "2e-3": "AR16",
    "3e-3": "AR16",
    "5e-3": "AR16",
    "7e-3": "AR16",
    "1e-2": "AR16",
    "3e-2": "AR16",
    "5e-2": "AR16",
    "1e-1": "AR16",
}

SKIP_COPY_DIRS = {
    "movie",
    "flowmov",
    "data",
    "contstr",
    "diagnostics",
    "__pycache__",
}


def parse_float(text: str) -> float:
    return float(text.replace("D", "e").replace("d", "e"))


def read_bou(path: Path) -> tuple[int, str, float]:
    lines = path.read_text(errors="ignore").splitlines()
    n1 = 0
    tmax = 0.0
    for i, line in enumerate(lines):
        if "N1" in line and "N2" in line and "N3" in line and i + 1 < len(lines):
            parts = lines[i + 1].split("!")[0].split()
            if len(parts) >= 3:
                n1 = int(parse_float(parts[0]))
        if "TMAX" in line and i + 1 < len(lines):
            parts = lines[i + 1].split("!")[0].split()
            if len(parts) >= 5:
                tmax = parse_float(parts[4])
    return n1, lines[0] if lines else "", tmax


def latest_time(run: Path) -> float:
    candidates = list(run.glob("data/avgvar.out"))
    if not candidates:
        return -1.0
    last = None
    for line in candidates[0].read_text(errors="ignore").splitlines():
        fields = line.split()
        if not fields:
            continue
        try:
            last = parse_float(fields[0])
        except ValueError:
            continue
    return -1.0 if last is None else last


def candidate_runs(ek_dir: Path, ar: str) -> list[tuple[Path, int, float]]:
    out: list[tuple[Path, int, float]] = []
    for bou in ek_dir.glob(f"{ar}/Beta1p02/qbot0p5_qtop0p004978/N*/**/run/bou.in"):
        run = bou.parent
        # Never use a staged strict-force continuation as the source for
        # another staging pass.  This prevents nested
        # conti_strict_force_500/conti_strict_force_500 directories.
        if "conti_strict_force_500" in run.parts:
            continue
        if not (run / "continua_q1.h5").exists():
            continue
        n1, _, _ = read_bou(bou)
        out.append((run, n1, latest_time(run)))
    return out


def patch_bou(text: str) -> str:
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if "TMAX" in line and i + 1 < len(lines):
            parts = lines[i + 1].split("!")
            values = parts[0].split()
            if len(values) >= 5:
                values[4] = "500d0"
                lines[i + 1] = "     ".join(values) + (
                    (" !" + parts[1]) if len(parts) > 1 else ""
                )
        if "NREAD" in line and "IRESET" in line and i + 1 < len(lines):
            parts = lines[i + 1].split("!")
            values = parts[0].split()
            if len(values) >= 3:
                values[0] = "1"
                lines[i + 1] = "     ".join(values) + (
                    (" !" + parts[1]) if len(parts) > 1 else ""
                )
    return "\n".join(lines) + "\n"


def make_subjob(source: str, job_name: str) -> str:
    lines = source.splitlines()
    out: list[str] = []
    saw_mkdir = False
    for line in lines:
        if line.startswith("#CSUB -J"):
            out.append(f"#CSUB -J {job_name}")
            continue
        if "check_drizzle_before_submit.py" in line:
            continue
        if re.search(r"rm\s+-r\s+fact\s+flowmov\s+data\s+contstr\s+movie", line):
            out.append("mkdir -p fact flowmov data contstr movie")
            saw_mkdir = True
            continue
        if re.search(r"mkdir\s+fact\s+flowmov\s+data\s+contstr\s+movie", line):
            if not saw_mkdir:
                out.append("mkdir -p fact flowmov data contstr movie")
                saw_mkdir = True
            continue
        if "impi-mpirun ./simexec" in line:
            out.append("mkdir -p fact flowmov data contstr movie")
            out.append("")
            out.append("impi-mpirun ./simexec")
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def stage_one(source_run: Path, simexec: Path, label: str) -> Path:
    parent = source_run.parent
    cont_name = "conti_strict_force_500"
    target = parent / cont_name / "run"
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing {target}")
    target.mkdir(parents=True)

    for item in source_run.iterdir():
        if item.name in SKIP_COPY_DIRS:
            continue
        if item.name in {"subjob.sh", "bou.in", "simexec"}:
            continue
        if item.name.endswith((".out", ".err")):
            continue
        destination = target / item.name
        if item.is_dir():
            shutil.copytree(item, destination)
        else:
            shutil.copy2(item, destination)

    shutil.copy2(simexec, target / "simexec")
    (target / "simexec").chmod(0o755)

    source_bou = source_run / "bou.in"
    (target / "bou.in").write_text(patch_bou(source_bou.read_text(errors="ignore")))

    source_subjob = source_run / "subjob.sh"
    job_name = f"Ra8e6Pr0p7{label}strict500"
    (target / "subjob.sh").write_text(make_subjob(source_subjob.read_text(errors="ignore"), job_name))
    (target / "subjob.sh").chmod(0o755)

    for dirname in ("fact", "flowmov", "data", "contstr", "movie"):
        (target / dirname).mkdir(exist_ok=True)
    digest = hashlib.sha256(simexec.read_bytes()).hexdigest()
    (target / "CONTINUATION_INFO.txt").write_text(
        "Prepared continuation for strict force-balance diagnostics\n"
        f"Source run: {source_run}\n"
        "Continuation duration: 500 physical time units (TMAX=500d0)\n"
        "Restart: NREAD=1\n"
        f"simexec SHA-256: {digest}\n"
        "Submission status: NOT SUBMITTED\n"
        "Submit manually from this run directory with: csub < subjob.sh\n"
    )
    return target


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--simexec", type=Path, default=DEFAULT_SIMEXEC)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not args.simexec.exists():
        raise SystemExit(f"missing simexec: {args.simexec}")

    report: list[str] = []
    for ek, ar in TARGETS.items():
        ek_dir = args.root / f"Ek{ek}"
        candidates = candidate_runs(ek_dir, ar)
        if not candidates:
            report.append(f"SKIP Ek{ek} ({ar}): no restartable qbot0.5 run")
            continue
        source, n1, tlast = max(candidates, key=lambda item: (item[1], item[2], str(item[0])))
        target = source.parent / "conti_strict_force_500" / "run"
        if args.dry_run:
            report.append(
                f"PLAN Ek{ek}: source={source} N={n1} tlast={tlast:g} target={target}"
            )
            continue
        try:
            staged = stage_one(source, args.simexec, f"Ek{ek}")
            report.append(
                f"READY Ek{ek}: source={source} N={n1} tlast={tlast:g} target={staged}"
            )
        except FileExistsError as exc:
            report.append(f"EXISTS Ek{ek}: {exc}")

    print("\n".join(report))
    print("No csub command was executed.")


if __name__ == "__main__":
    main()
