#!/usr/bin/env python3
"""Reduce remote movie VZ_me fields to a full-volume vertical-speed history."""

from __future__ import print_function

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


TIME_RE = re.compile(r'<Time\s+Value="([+\-0-9.eEdD]+)"')
FRAME_RE = re.compile(r"field(\d+)\.xmf$")


def read_time(xmf_path):
    text = xmf_path.read_text(errors="ignore")
    match = TIME_RE.search(text)
    if match is None:
        return None
    return float(match.group(1).replace("D", "e").replace("d", "e"))


def snapshots(run_paths):
    candidates = []
    for priority, run_text in enumerate(run_paths):
        movie = Path(run_text) / "movie"
        if not movie.is_dir():
            continue
        for xmf in movie.glob("field[0-9]*.xmf"):
            match = FRAME_RE.search(xmf.name)
            if match is None:
                continue
            frame = int(match.group(1))
            h5 = movie / ("field{:05d}.h5".format(frame))
            if not h5.exists():
                continue
            time = read_time(xmf)
            if time is None:
                time = 10.0 * frame
            candidates.append((float(time), priority, frame, h5))

    selected = {}
    for item in candidates:
        key = (item[0], item[2])
        if key not in selected or item[1] > selected[key][1]:
            selected[key] = item
    return sorted(selected.values(), key=lambda item: (item[0], item[2]))


def dump_vz(h5_path, output_path):
    subprocess.check_call(
        ["h5dump", "-d", "/VZ_me", "-b", "LE", "-o", str(output_path), str(h5_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return np.fromfile(str(output_path), dtype="<f8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases_json")
    parser.add_argument("output_csv")
    parser.add_argument("--stride", type=int, default=5)
    args = parser.parse_args()

    cases = json.loads(Path(args.cases_json).read_text())
    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="w_rms_"))
    rows = []
    errors = []
    try:
        for case in cases:
            label = case.get("Ek_label", "unknown")
            available = snapshots(case["run_paths"])
            stride = max(1, args.stride)
            selected = available[::stride]
            if selected and selected[-1] != available[-1]:
                selected.append(available[-1])
            print("[{}] snapshots={}".format(label, len(selected)), flush=True)
            for time, priority, frame, h5_path in selected:
                try:
                    raw = dump_vz(h5_path, work / "vz.bin")
                    if raw.size == 0 or not np.all(np.isfinite(raw)):
                        raise RuntimeError("empty or non-finite VZ_me data")
                    wrms = float(np.sqrt(np.mean(raw * raw)))
                    rows.append({
                        "Ra": case["Ra"],
                        "Pr": case["Pr"],
                        "Ek": case.get("Ek"),
                        "AR": case["AR"],
                        "Nx": case.get("Nx", ""),
                        "Ny": case.get("Ny", ""),
                        "Nz": case.get("Nz", ""),
                        "time": time,
                        "frame": frame,
                        "w_rms": wrms,
                        "field_path": str(h5_path),
                    })
                except Exception as exc:
                    errors.append({
                        "Ek": case.get("Ek"),
                        "time": time,
                        "frame": frame,
                        "field_path": str(h5_path),
                        "error": str(exc),
                    })
    finally:
        shutil.rmtree(str(work), ignore_errors=True)

    fields = [
        "Ra", "Pr", "Ek", "AR", "Nx", "Ny", "Nz",
        "time", "frame", "w_rms", "field_path",
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    error_path = output_csv.with_name(output_csv.stem + "_errors.json")
    error_path.write_text(json.dumps(errors, indent=2))
    print("completed={} errors={}".format(len(rows), len(errors)))


if __name__ == "__main__":
    main()
