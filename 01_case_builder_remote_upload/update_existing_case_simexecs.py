from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "case_config.json"


REMOTE_UPDATE_SCRIPT = r"""
import hashlib
import json
import os
import shutil
import time
from pathlib import Path


def parse_payload() -> dict:
    if "PAYLOAD_JSON" in globals():
        return json.loads(PAYLOAD_JSON)
    return json.loads(__import__("sys").stdin.read())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def is_case_run(path: Path) -> bool:
    if path.name != "run":
        return False
    if not (path / "bou.in").is_file():
        return False
    parts = set(path.parts)
    skip_names = {
        "latest_program",
        "Rotating_Moist_RB3-ns-deltab727",
        "source",
        "postprocess_scripts_index",
    }
    return not bool(parts.intersection(skip_names))


def main(cfg: dict) -> dict:
    remote_root = Path(cfg["remote_root"])
    latest_simexec = Path(cfg["latest_simexec"])
    dry_run = bool(cfg.get("dry_run", False))
    include_existing_same_hash = bool(cfg.get("include_existing_same_hash", False))
    only_missing = bool(cfg.get("only_missing", False))

    if not latest_simexec.is_file():
        raise FileNotFoundError(f"latest simexec is missing: {latest_simexec}")
    if latest_simexec.stat().st_size <= 0:
        raise RuntimeError(f"latest simexec is empty: {latest_simexec}")

    latest_hash = sha256(latest_simexec)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    updated = []
    skipped = []
    backed_up = []

    run_dirs = sorted(p for p in remote_root.rglob("run") if is_case_run(p))
    for run_dir in run_dirs:
        target = run_dir / "simexec"
        item = {
            "run_dir": str(run_dir),
            "target": str(target),
        }
        if only_missing and target.exists():
            item["reason"] = "exists"
            skipped.append(item)
            continue
        if target.exists() and target.is_file():
            old_hash = sha256(target)
            item["old_sha256"] = old_hash
            if old_hash == latest_hash and not include_existing_same_hash:
                item["reason"] = "same_sha256"
                skipped.append(item)
                continue
            backup = run_dir / f"simexec_before_latest_{stamp}"
            item["backup"] = str(backup)
            if not dry_run:
                shutil.copy2(target, backup)
                backed_up.append(str(backup))
        else:
            item["old_sha256"] = None

        item["new_sha256"] = latest_hash
        if not dry_run:
            shutil.copy2(latest_simexec, target)
            os.chmod(target, target.stat().st_mode | 0o755)
        updated.append(item)

    report = {
        "status": "dry_run" if dry_run else "updated",
        "remote_root": str(remote_root),
        "latest_simexec": str(latest_simexec),
        "latest_sha256": latest_hash,
        "num_run_dirs_found": len(run_dirs),
        "num_updated": len(updated),
        "num_skipped": len(skipped),
        "num_backed_up": len(backed_up),
        "updated": updated,
        "skipped": skipped,
    }
    report_path = remote_root / f"latest_simexec_update_report_{stamp}.json"
    if not dry_run:
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    report["report_path"] = str(report_path)
    return report


cfg = parse_payload()
print(json.dumps(main(cfg), indent=2))
"""


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def run_remote(cfg: dict) -> dict:
    payload = json.dumps(cfg)
    bootstrap = f"PAYLOAD_JSON = {payload!r}\n" + REMOTE_UPDATE_SCRIPT
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Update existing remote case run/simexec files from the latest compiled simexec."
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--ssh-alias")
    parser.add_argument("--remote-root")
    parser.add_argument("--latest-simexec")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--only-missing", action="store_true")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Copy even when an existing run/simexec already has the same SHA256.",
    )
    args = parser.parse_args()

    base = load_json(args.config)
    cfg = {
        "ssh_alias": args.ssh_alias or base["ssh_alias"],
        "remote_root": args.remote_root or base["remote_root"],
        "latest_simexec": args.latest_simexec
        or base.get("program_bin")
        or (base["remote_root"].rstrip("/") + "/latest_program/source/simexec"),
        "dry_run": args.dry_run,
        "only_missing": args.only_missing,
        "include_existing_same_hash": args.force,
    }
    run_remote(cfg)


if __name__ == "__main__":
    main()
