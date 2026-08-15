from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "case_config.json"


REMOTE_BOUIN_SCRIPT = r"""
import json
import time
from pathlib import Path


def parse_payload() -> dict:
    if "PAYLOAD_JSON" in globals():
        return json.loads(PAYLOAD_JSON)
    return json.loads(__import__("sys").stdin.read())


def is_case_run(path: Path) -> bool:
    if path.name != "run":
        return False
    if not (path / "bou.in").is_file():
        return False
    skip_names = {
        "latest_program",
        "Rotating_Moist_RB3-ns-deltab727",
        "source",
        "postprocess_scripts_index",
    }
    return not bool(set(path.parts).intersection(skip_names))


def find_movie_row(lines: list[str]) -> tuple[int, list[str]]:
    for i, line in enumerate(lines[:-1]):
        if "mov_zcut_k" in line:
            return i + 1, lines[i + 1].split()
    raise RuntimeError("could not find header containing mov_zcut_k")


def main(cfg: dict) -> dict:
    remote_root = Path(cfg["remote_root"])
    vortex_lc = str(cfg.get("vortex_lc", "0"))
    dry_run = bool(cfg.get("dry_run", False))
    set_existing = bool(cfg.get("set_existing", False))
    stamp = time.strftime("%Y%m%d_%H%M%S")

    changed = []
    skipped = []
    errors = []

    run_dirs = sorted(p for p in remote_root.rglob("run") if is_case_run(p))
    for run_dir in run_dirs:
        bou = run_dir / "bou.in"
        item = {"run_dir": str(run_dir), "bou_in": str(bou)}
        try:
            lines = bou.read_text(encoding="utf-8", errors="ignore").splitlines()
            row_i, vals = find_movie_row(lines)
            item["old_values"] = vals
            new_vals = list(vals)
            if len(new_vals) == 4:
                if set_existing:
                    new_vals[3] = vortex_lc
                else:
                    item["reason"] = "already_has_vortex_lc"
                    item["new_values"] = new_vals
                    skipped.append(item)
                    continue
            elif len(new_vals) == 3:
                new_vals.append(vortex_lc)
            elif len(new_vals) == 2:
                new_vals.extend(["600", vortex_lc])
            elif len(new_vals) == 1:
                new_vals.extend(["200d0", "600", vortex_lc])
            else:
                item["reason"] = f"unexpected_value_count_{len(new_vals)}"
                errors.append(item)
                continue

            item["new_values"] = new_vals
            backup = run_dir / f"bou.in_before_vortex_lc_{stamp}"
            item["backup"] = str(backup)
            if not dry_run:
                backup.write_text("\n".join(lines) + "\n", encoding="utf-8")
                lines[row_i] = "     ".join(new_vals)
                bou.write_text("\n".join(lines) + "\n", encoding="utf-8")
            changed.append(item)
        except Exception as exc:
            item["error"] = str(exc)
            errors.append(item)

    report = {
        "status": "dry_run" if dry_run else "updated",
        "remote_root": str(remote_root),
        "vortex_lc_written": vortex_lc,
        "set_existing": set_existing,
        "num_run_dirs_found": len(run_dirs),
        "num_changed": len(changed),
        "num_skipped": len(skipped),
        "num_errors": len(errors),
        "changed": changed,
        "skipped": skipped,
        "errors": errors,
    }
    report_path = remote_root / f"bouin_vortex_lc_update_report_{stamp}.json"
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
    bootstrap = f"PAYLOAD_JSON = {payload!r}\n" + REMOTE_BOUIN_SCRIPT
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
        description="Ensure existing remote bou.in files expose vortex_lc in the mov_zcut_k row."
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--ssh-alias")
    parser.add_argument("--remote-root")
    parser.add_argument("--vortex-lc", default="0")
    parser.add_argument("--set-existing", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    base = load_json(args.config)
    cfg = {
        "ssh_alias": args.ssh_alias or base["ssh_alias"],
        "remote_root": args.remote_root or base["remote_root"],
        "vortex_lc": args.vortex_lc,
        "set_existing": args.set_existing,
        "dry_run": args.dry_run,
    }
    run_remote(cfg)


if __name__ == "__main__":
    main()
