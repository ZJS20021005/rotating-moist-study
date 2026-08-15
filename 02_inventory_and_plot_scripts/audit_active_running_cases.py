from __future__ import annotations

import json
import subprocess
from datetime import datetime
from pathlib import Path


REMOTE_HOST = "c01n0006"
REMOTE_ROOT = "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case"
OUTPUT = Path(
    r"E:\moist RB\rotating_case_inventory\03_inventory_tables"
) / "active_running_cases_latest.json"
OBSERVATION_SECONDS = 60


REMOTE_CODE = r'''
import json, math, re, time
from pathlib import Path

ROOT = Path(r"""__REMOTE_ROOT__""")
WAIT = __WAIT__


def as_float(value):
    try:
        return float(str(value).replace("D", "e").replace("d", "e"))
    except Exception:
        return None


def row_after(lines, marker):
    for index, line in enumerate(lines):
        if marker in line:
            for following in lines[index + 1:]:
                clean = following.split("!")[0].strip()
                if clean:
                    return clean.split()
    return []


def read_bou(run):
    try:
        lines = (run / "bou.in").read_text(errors="ignore").splitlines()
    except Exception:
        return None
    grid = row_after(lines, "N1      N2")
    geometry = row_after(lines, "ALX3D")
    control = row_after(lines, "Ra       Pr")
    boundary = row_after(lines, "A_stopmod")
    if len(grid) < 3 or len(geometry) < 3 or len(control) < 9 or len(boundary) < 8:
        return None
    Ra, Pr, invRo = as_float(control[0]), as_float(control[1]), as_float(control[2])
    Ek = None
    if Ra and Pr is not None and invRo is not None and abs(invRo) > 1.0e-15:
        Ek = math.sqrt(Pr / Ra) / invRo
    return {
        "Ra": Ra,
        "Pr": Pr,
        "Ek": Ek,
        "Nx": int(as_float(grid[0])),
        "Ny": int(as_float(grid[1])),
        "Nz": int(as_float(grid[2])),
        "AR": as_float(geometry[1]),
        "beta": as_float(control[7]),
        "qbot": as_float(boundary[7]),
        "qtop": as_float(boundary[6]),
    }


def last_numeric_time(path):
    try:
        with path.open("rb") as handle:
            handle.seek(0, 2)
            end = handle.tell()
            handle.seek(max(0, end - 131072))
            text = handle.read().decode("utf-8", "ignore")
    except Exception:
        return None
    answer = None
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        value = as_float(parts[0])
        if value is not None and math.isfinite(value):
            answer = value
    return answer


def monitored_files(run):
    paths = [
        run / "data" / "avgvar.out",
        run / "data" / "mse_aggregation.out",
        run / "data" / "mse_aggregation_scales.out",
        run / "diagnostics" / "thermo" / "mprime_square.dat",
        run / "diagnostics" / "scale" / "convective_scale.dat",
        run / "diagnostics" / "scale" / "moist_integral_scale.dat",
        run / "diagnostics" / "scale" / "peak_scale.dat",
    ]
    return [path for path in paths if path.exists()]


def scan():
    records = {}
    for bou in ROOT.rglob("bou.in"):
        run = bou.parent
        if run.name != "run":
            continue
        params = read_bou(run)
        if params is None:
            continue
        files = {}
        for path in monitored_files(run):
            try:
                stat = path.stat()
            except Exception:
                continue
            files[str(path.relative_to(run))] = {
                "size": stat.st_size,
                "mtime": stat.st_mtime,
                "last_time": last_numeric_time(path),
            }
        latest_required = [
            "diagnostics/thermo/mprime_square.dat",
            "diagnostics/scale/convective_scale.dat",
            "diagnostics/scale/moist_integral_scale.dat",
            "diagnostics/scale/peak_scale.dat",
        ]
        params.update({
            "run_path": str(run),
            "latest_program": all(name in files and files[name]["size"] > 0 for name in latest_required),
            "files": files,
        })
        records[str(run)] = params
    return records


first = scan()
time.sleep(WAIT)
second = scan()
now = time.time()
results = []
for run_path in sorted(set(first) | set(second)):
    before = first.get(run_path, {})
    after = second.get(run_path, {})
    row = {key: value for key, value in after.items() if key != "files"}
    row["run_path"] = run_path
    changes = []
    newest_mtime = None
    latest_times = []
    for name, current in after.get("files", {}).items():
        previous = before.get("files", {}).get(name)
        newest_mtime = current["mtime"] if newest_mtime is None else max(newest_mtime, current["mtime"])
        if current.get("last_time") is not None:
            latest_times.append(current["last_time"])
        if previous is None:
            changes.append({"file": name, "kind": "created", "delta_time": None})
            continue
        advanced = (
            current["size"] != previous["size"]
            or current["mtime"] > previous["mtime"] + 1.0e-6
            or (
                current.get("last_time") is not None
                and previous.get("last_time") is not None
                and current["last_time"] > previous["last_time"] + 1.0e-8
            )
        )
        if advanced:
            delta_time = None
            if current.get("last_time") is not None and previous.get("last_time") is not None:
                delta_time = current["last_time"] - previous["last_time"]
            changes.append({"file": name, "kind": "advanced", "delta_time": delta_time})
    row["active_during_observation"] = bool(changes)
    row["changes"] = changes
    row["latest_diagnostic_time"] = max(latest_times) if latest_times else None
    row["seconds_since_latest_write"] = None if newest_mtime is None else now - newest_mtime
    results.append(row)

print(json.dumps({
    "remote_host": "__REMOTE_HOST__",
    "remote_root": str(ROOT),
    "observation_seconds": WAIT,
    "remote_scan_epoch": now,
    "records": results,
}, ensure_ascii=False))
'''


def main() -> None:
    code = (
        REMOTE_CODE.replace("__REMOTE_ROOT__", REMOTE_ROOT)
        .replace("__REMOTE_HOST__", REMOTE_HOST)
        .replace("__WAIT__", str(OBSERVATION_SECONDS))
    )
    process = subprocess.run(
        ["ssh", REMOTE_HOST, "python3", "-"],
        input=code,
        text=True,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=OBSERVATION_SECONDS + 90,
    )
    if process.returncode != 0:
        raise RuntimeError(
            "Remote active-run audit failed.\n"
            f"STDOUT:\n{process.stdout}\nSTDERR:\n{process.stderr}"
        )
    payload = json.loads(process.stdout)
    payload["local_scanned_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    active = [row for row in payload["records"] if row["active_during_observation"]]
    print(f"saved {OUTPUT}")
    print(f"active during {OBSERVATION_SECONDS}s: {len(active)}")
    for row in active:
        ek = "NR" if row.get("Ek") is None else f"{row['Ek']:.6g}"
        print(
            f"Ra={row.get('Ra'):.6g} Ek={ek} AR={row.get('AR'):g} "
            f"beta={row.get('beta'):g} qbot={row.get('qbot'):g} "
            f"N={row.get('Nx')}x{row.get('Ny')}x{row.get('Nz')} "
            f"t={row.get('latest_diagnostic_time')}"
        )


if __name__ == "__main__":
    main()
