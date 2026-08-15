#!/usr/bin/env python3
"""Inventory and preserve critical remote Rainy-Benard files before cleanup.

This tool is intentionally non-destructive.  It never removes or modifies a
remote file.  The default action only writes a remote file manifest locally.
Use --download-critical to additionally stream two tar.gz archives to the
local machine: programs/scripts and lightweight case metadata/diagnostics.
Large HDF5 movie and checkpoint files remain a separate reviewed download.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shlex
import subprocess
from pathlib import Path


DEFAULT_REMOTE_ROOT = "/share/org/SHUTUANL/shu_zhangjs/rainy model"
DEFAULT_LOCAL_ROOT = Path(r"H:\moist_RB_remote_archive_20260807")

PROGRAM_PATHS = (
    "rotating_case/latest_program",
    "rotating_case/postprocess",
    "ns/transition_study/beta1/postprocess_scripts_index",
)

STUDY_PATHS = (
    "rotating_case",
    "ns/transition_study",
)

LIGHTWEIGHT_EXCLUDES = (
    "*.h5",
    "*.bin",
    "*.o",
    "*.mod",
    "*.pyc",
    "__pycache__",
    "strict_force_*",
    "core.*",
    "simexec*",
)


def ssh(host: str, command: str, *, text: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["ssh", "-o", "BatchMode=yes", host, command],
        check=True,
        text=text,
        capture_output=True,
    )


def verify_connection(host: str, remote_root: str) -> None:
    command = (
        "set -eu; "
        "hostname; "
        f"test -d {shlex.quote(remote_root)}; "
        f"du -sh {shlex.quote(remote_root)}"
    )
    result = ssh(host, command)
    print(result.stdout.strip())


def write_manifest(host: str, remote_root: str, output: Path) -> None:
    command = (
        f"find {shlex.quote(remote_root)} -type f "
        "-printf '%p\\t%s\\t%T@\\n'"
    )
    result = ssh(host, command)
    output.write_text(
        "remote_path\tbytes\tmtime_epoch\n" + result.stdout,
        encoding="utf-8",
        newline="",
    )
    print(f"wrote {output}")


def existing_relative_paths(host: str, remote_root: str, paths: tuple[str, ...]) -> list[str]:
    existing = []
    for relative in paths:
        full = f"{remote_root.rstrip('/')}/{relative}"
        result = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", host, f"test -e {shlex.quote(full)}"],
            check=False,
        )
        if result.returncode == 0:
            existing.append(relative)
    return existing


def stream_tar(
    host: str,
    remote_root: str,
    relative_paths: list[str],
    output: Path,
    excludes: tuple[str, ...] = (),
) -> None:
    if not relative_paths:
        print(f"skip {output.name}: no requested remote path exists")
        return
    exclude_args = " ".join(
        f"--exclude={shlex.quote(pattern)}" for pattern in excludes
    )
    paths = " ".join(shlex.quote(path) for path in relative_paths)
    command = (
        "set -eu; "
        f"tar -czf - -C {shlex.quote(remote_root)} {exclude_args} {paths}"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as handle:
        process = subprocess.Popen(
            ["ssh", "-o", "BatchMode=yes", host, command],
            stdout=handle,
            stderr=subprocess.PIPE,
        )
        _, stderr = process.communicate()
    if process.returncode != 0:
        output.unlink(missing_ok=True)
        raise RuntimeError(stderr.decode(errors="replace"))
    print(f"wrote {output} ({output.stat().st_size / 2**20:.1f} MiB)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_local_checksums(root: Path) -> None:
    rows = []
    for path in sorted(root.glob("*.tar.gz")):
        rows.append(
            {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    with (root / "archive_checksums.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", required=True, help="Configured SSH alias")
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--local-root", type=Path, default=DEFAULT_LOCAL_ROOT)
    parser.add_argument("--download-critical", action="store_true")
    args = parser.parse_args()

    local_root = args.local_root.resolve()
    local_root.mkdir(parents=True, exist_ok=True)
    verify_connection(args.host, args.remote_root)
    write_manifest(
        args.host,
        args.remote_root,
        local_root / "remote_all_files_manifest.tsv",
    )

    if args.download_critical:
        programs = existing_relative_paths(
            args.host, args.remote_root, PROGRAM_PATHS
        )
        studies = existing_relative_paths(
            args.host, args.remote_root, STUDY_PATHS
        )
        stream_tar(
            args.host,
            args.remote_root,
            programs,
            local_root / "remote_programs_and_scripts.tar.gz",
        )
        stream_tar(
            args.host,
            args.remote_root,
            studies,
            local_root / "remote_case_metadata_and_timeseries_no_hdf5.tar.gz",
            LIGHTWEIGHT_EXCLUDES,
        )
        write_local_checksums(local_root)

    print("No remote file was modified or deleted.")


if __name__ == "__main__":
    main()
