#!/usr/bin/env python3
"""Validate that drizzle_init.dat matches the current case bou.in."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


PARAMETER_NAMES = [
    "Ra",
    "Pr",
    "invRo",
    "gamma",
    "Sm",
    "alphaqs",
    "betaqs",
    "tau_cond",
    "b_bot",
    "b_top",
    "q_bot",
    "q_top",
    "saturation_width",
]


def number(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def value_row(lines: list[str], token: str) -> list[str]:
    for index, line in enumerate(lines[:-1]):
        if token in line:
            return lines[index + 1].split()
    raise RuntimeError(f"Missing bou.in header containing {token!r}")


def close_enough(actual: float, expected: float) -> bool:
    return abs(actual - expected) <= 2.0e-10 * max(1.0, abs(expected))


def read_bou_parameters(path: Path) -> tuple[list[float], list[float]]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    control = [number(value) for value in value_row(lines, "Ra       Pr     invRo")]
    boundary = [number(value) for value in value_row(lines, "A_stopmod")]
    if len(control) != 9:
        raise RuntimeError(f"Expected 9 control values, found {len(control)}")
    if len(boundary) != 8:
        raise RuntimeError(f"Expected 8 boundary values, found {len(boundary)}")
    if abs(boundary[0]) > 1.0e-14 or abs(boundary[2]) > 1.0e-14:
        raise RuntimeError(
            "One-dimensional drizzle initialization requires "
            "A_stopmod=A_sbotmod=0"
        )
    expected = [
        control[0],
        control[1],
        control[2],
        control[4],
        control[5],
        control[6],
        control[7],
        control[8],
        boundary[5],
        boundary[4],
        boundary[7],
        boundary[6],
    ]
    return expected, boundary


def validate(
    bou_path: Path,
    profile_path: Path,
    expected_perturb: float,
) -> None:
    expected, boundary = read_bou_parameters(bou_path)
    raw_lines = [
        line.strip()
        for line in profile_path.read_text(encoding="ascii").splitlines()
        if line.strip()
    ]
    if len(raw_lines) < 4:
        raise RuntimeError("drizzle_init.dat is incomplete")
    first = raw_lines[0].split()
    if len(first) != 2:
        raise RuntimeError("First profile row must contain nprofile and perturb_amp")
    nprofile = int(first[0])
    perturb_amp = number(first[1])
    if nprofile < 33:
        raise RuntimeError(f"nprofile={nprofile} is too small")
    if not close_enough(perturb_amp, expected_perturb):
        raise RuntimeError(
            f"Perturbation mismatch: profile={perturb_amp:.16e}, "
            f"expected={expected_perturb:.16e}"
        )

    profile_parameters = [number(value) for value in raw_lines[1].split()]
    if len(profile_parameters) != len(PARAMETER_NAMES):
        raise RuntimeError(
            f"Expected {len(PARAMETER_NAMES)} profile parameters, "
            f"found {len(profile_parameters)}"
        )
    expected_with_width = expected + [profile_parameters[-1]]
    for name, actual, wanted in zip(
        PARAMETER_NAMES, profile_parameters, expected_with_width
    ):
        if not close_enough(actual, wanted):
            raise RuntimeError(
                f"{name} mismatch: profile={actual:.16e}, "
                f"bou.in={wanted:.16e}"
            )

    data_lines = raw_lines[2:]
    if len(data_lines) != nprofile:
        raise RuntimeError(
            f"Profile declares {nprofile} rows but contains {len(data_lines)}"
        )
    rows = [[number(value) for value in line.split()] for line in data_lines]
    if any(len(row) != 3 for row in rows):
        raise RuntimeError("Every profile data row must contain z, b and q")
    if any(not all(math.isfinite(value) for value in row) for row in rows):
        raise RuntimeError("Profile contains a non-finite value")
    if any(rows[index + 1][0] <= rows[index][0] for index in range(nprofile - 1)):
        raise RuntimeError("Profile z values are not strictly increasing")
    if not close_enough(rows[0][0], 0.0) or not close_enough(rows[-1][0], 1.0):
        raise RuntimeError("Profile must cover 0 <= z/H <= 1")

    endpoints = [
        ("b_bot", rows[0][1], boundary[5]),
        ("b_top", rows[-1][1], boundary[4]),
        ("q_bot", rows[0][2], boundary[7]),
        ("q_top", rows[-1][2], boundary[6]),
    ]
    for name, actual, wanted in endpoints:
        if not close_enough(actual, wanted):
            raise RuntimeError(
                f"{name} endpoint mismatch: profile={actual:.16e}, "
                f"bou.in={wanted:.16e}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bou-in", type=Path, default=Path("bou.in"))
    parser.add_argument(
        "--profile", type=Path, default=Path("drizzle_init.dat")
    )
    parser.add_argument("--expected-perturb", type=float, default=1.0e-4)
    args = parser.parse_args()
    validate(args.bou_in, args.profile, args.expected_perturb)
    print(
        "Drizzle initialization check passed: "
        f"{args.profile} matches {args.bou_in}"
    )


if __name__ == "__main__":
    main()
