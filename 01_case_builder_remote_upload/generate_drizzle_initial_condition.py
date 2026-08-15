#!/usr/bin/env python3
"""Generate a case-specific Rainy-Benard drizzle initial condition.

This wrapper deliberately reuses ``moist_base_state`` from the linear
stability package supplied by the user.  It converts one case configuration
into:

* ``drizzle_init.dat``: compact profile read by the DNS at ``nread=0``;
* ``drizzle_init_meta.json``: auditable parameters, residuals and provenance.

The generator refuses the linear-profile fallback in ``moist_base_state``.
If the nonlinear drizzle solve does not converge, case preparation fails.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
SOLVER_CANDIDATES = [
    HERE / "linear_stability_reference" / "stability_solver.py",
    HERE / "stability_solver.py",
]
SOLVER_PATH = next(
    (candidate for candidate in SOLVER_CANDIDATES if candidate.is_file()),
    SOLVER_CANDIDATES[0],
)
REFERENCE_DIR = SOLVER_PATH.parent
sys.path.insert(0, str(REFERENCE_DIR))

from stability_solver import (  # noqa: E402
    Params,
    fd_matrices,
    moist_base_state,
    smooth_heaviside,
)


PROFILE_FORMAT_VERSION = 1
DEFAULT_PROFILE_N = 401
DEFAULT_PERTURB_AMP = 1.0e-4
# Match compute_qsat_from_b.f90:
# H=0.5*(1+tanh(1e8*(q-qs))) = H_smooth(q-qs; width=1e-8).
DEFAULT_SATURATION_WIDTH = 1.0e-8


def parse_fortran_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def bou_value_row(lines: list[str], token: str) -> list[str]:
    for index, line in enumerate(lines[:-1]):
        if token in line:
            return lines[index + 1].split()
    raise ValueError(f"Could not find bou.in row after header {token!r}")


def config_from_bou(path: Path) -> dict[str, Any]:
    """Read the current case parameters directly from ``bou.in``."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    control = bou_value_row(lines, "Ra       Pr     invRo")
    boundary = bou_value_row(lines, "A_stopmod")
    velocity_bc = bou_value_row(lines, "UBCBOT    UBCTOP")
    if len(control) < 9:
        raise ValueError("The Ra/Pr/invRo row in bou.in needs 9 values")
    if len(boundary) < 8:
        raise ValueError("The scalar-boundary row in bou.in needs 8 values")
    if len(velocity_bc) < 2:
        raise ValueError("The velocity-boundary row in bou.in needs 2 values")

    values = [parse_fortran_float(value) for value in control[:9]]
    bounds = [parse_fortran_float(value) for value in boundary[:8]]
    ra, pr, inv_ro = values[:3]
    beta = values[7]
    ek = None if abs(inv_ro) < 1.0e-30 else math.sqrt(pr / ra) / inv_ro
    dsal_top, dsal_bot = bounds[4], bounds[5]
    if abs(dsal_bot) > 1.0e-12:
        raise ValueError(
            f"Drizzle initialization requires dsalbot=0, got {dsal_bot}"
        )
    expected_top = beta - 1.0
    if abs(dsal_top - expected_top) > 1.0e-12 * max(1.0, abs(expected_top)):
        raise ValueError(
            "Drizzle initialization requires dsaltop=beta-1; "
            f"got dsaltop={dsal_top}, beta-1={expected_top}"
        )

    return {
        "ra": ra,
        "pr": pr,
        "ek": ek,
        "gamma": values[4],
        "sm": values[5],
        "alpha": values[6],
        "beta": beta,
        "tau": values[8],
        "a_stopmod": bounds[0],
        "a_sbotmod": bounds[2],
        "qvaptop": bounds[6],
        "qvapbot": bounds[7],
        "ubcbot": int(velocity_bc[0]),
        "ubctop": int(velocity_bc[1]),
    }


def normalize_ek(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"", "none", "no", "nr", "norotating", "nonrotating", "0"}:
            return None
        return float(text)
    numeric = float(value)
    return None if numeric == 0.0 else numeric


def computed_inv_ro(ra: float, pr: float, ek_value: Any) -> float:
    ek = normalize_ek(ek_value)
    return 0.0 if ek is None else math.sqrt(pr / ra) / ek


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def drizzle_residual(
    z: np.ndarray,
    b: np.ndarray,
    q: np.ndarray,
    params: Params,
) -> tuple[float, float, float]:
    _, _, d2 = fd_matrices(z.size)
    qs = np.exp(params.alphaqs * (b - params.betaqs * z))
    supersat = q - qs
    condensation = (
        smooth_heaviside(supersat, params.saturation_width)
        * supersat
        / params.tau_cond
    )
    residual_b = d2 @ b + params.gamma * condensation
    residual_q = params.sm * (d2 @ q) - condensation
    residual_b[0] = b[0] - params.b_bot
    residual_b[-1] = b[-1] - params.b_top
    residual_q[0] = q[0] - params.q_bot
    residual_q[-1] = q[-1] - params.q_top
    return (
        float(np.max(np.abs(np.r_[residual_b, residual_q]))),
        float(np.max(np.abs(residual_b))),
        float(np.max(np.abs(residual_q))),
    )


def build_profile(cfg: dict[str, Any]) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    ra = float(cfg["ra"])
    pr = float(cfg["pr"])
    gamma = float(cfg["gamma"])
    sm = float(cfg.get("sm", 1.0))
    alphaqs = float(cfg["alpha"])
    betaqs = float(cfg["beta"])
    tau_cond = float(cfg["tau"])
    b_bot = 0.0
    b_top = betaqs - 1.0
    q_bot = float(cfg["qvapbot"])
    q_top = float(cfg["qvaptop"])
    profile_n = int(cfg.get("drizzle_profile_n", DEFAULT_PROFILE_N))
    perturb_amp = float(
        cfg.get("drizzle_perturb_amp", DEFAULT_PERTURB_AMP)
    )
    saturation_width = float(
        cfg.get("drizzle_saturation_width", DEFAULT_SATURATION_WIDTH)
    )

    if profile_n < 33:
        raise ValueError("drizzle_profile_n must be at least 33")
    if profile_n % 2 == 0:
        profile_n += 1
    if perturb_amp <= 0.0 or perturb_amp > 1.0e-3:
        raise ValueError(
            "drizzle_perturb_amp must be positive and no larger than 1e-3"
        )
    if saturation_width <= 0.0:
        raise ValueError("drizzle_saturation_width must be positive")

    for name in ("a_stopmod", "a_sbotmod"):
        if abs(float(cfg.get(name, 0.0))) > 1.0e-14:
            raise ValueError(
                "The one-dimensional drizzle base state requires "
                f"{name}=0, got {cfg.get(name)!r}"
            )

    params = Params(
        pr=pr,
        gamma=gamma,
        sm=sm,
        alphaqs=alphaqs,
        betaqs=betaqs,
        tau_cond=tau_cond,
        b_bot=b_bot,
        b_top=b_top,
        q_bot=q_bot,
        q_top=q_top,
        ubc_bot="noslip" if int(cfg.get("ubcbot", 1)) == 1 else "freeslip",
        ubc_top="noslip" if int(cfg.get("ubctop", 0)) == 1 else "freeslip",
        condensation_mode="smooth",
        saturation_width=saturation_width,
        linearize_switch=False,
    )

    state = moist_base_state(profile_n, params)
    converged = bool(round(float(state["converged"][0])))
    message = str(state["message"][0])
    if not converged:
        raise RuntimeError(
            "The nonlinear drizzle solve did not converge. "
            "The linear fallback was rejected. Solver message: " + message
        )

    z = np.asarray(state["z"], dtype=float)
    b = np.asarray(state["b"], dtype=float)
    q = np.asarray(state["q"], dtype=float)
    qs = np.asarray(state["qs"], dtype=float)
    supersat = np.asarray(state["supersat"], dtype=float)
    if not np.all(np.isfinite(np.r_[z, b, q, qs, supersat])):
        raise RuntimeError("The converged drizzle profile contains non-finite values")
    if not np.all(np.diff(z) > 0.0):
        raise RuntimeError("The drizzle profile z grid is not strictly increasing")

    residual_max, residual_b_max, residual_q_max = drizzle_residual(
        z, b, q, params
    )
    if residual_max > 1.0e-5:
        raise RuntimeError(
            f"Drizzle residual is too large ({residual_max:.6e}); "
            "case preparation was stopped"
        )

    inv_ro = computed_inv_ro(ra, pr, cfg.get("ek"))
    header_parameters = np.asarray(
        [
            ra,
            pr,
            inv_ro,
            gamma,
            sm,
            alphaqs,
            betaqs,
            tau_cond,
            b_bot,
            b_top,
            q_bot,
            q_top,
            saturation_width,
        ],
        dtype=float,
    )
    metadata: dict[str, Any] = {
        "profile_format_version": PROFILE_FORMAT_VERSION,
        "source": "user-supplied stability_solver.moist_base_state",
        "source_file": str(SOLVER_PATH),
        "source_sha256": sha256(SOLVER_PATH),
        "solver_converged": True,
        "solver_message": message,
        "profile_n": int(z.size),
        "drizzle_perturb_amp": perturb_amp,
        "perturbed_dns_variable": "b only",
        "perturbation_envelope": "sin(pi*z/H)",
        "perturbation_distribution": "zero-mean unit-variance Gaussian",
        "ra": ra,
        "pr": pr,
        "ek": cfg.get("ek"),
        "invRo": inv_ro,
        "gamma": gamma,
        "sm": sm,
        "alphaqs": alphaqs,
        "betaqs": betaqs,
        "tau_cond": tau_cond,
        "b_bot": b_bot,
        "b_top": b_top,
        "q_bot": q_bot,
        "q_top": q_top,
        "saturation_width": saturation_width,
        "residual_max": residual_max,
        "residual_b_max": residual_b_max,
        "residual_q_max": residual_q_max,
        "min_supersaturation": float(np.min(supersat)),
        "max_supersaturation": float(np.max(supersat)),
        "equations": {
            "buoyancy": "d2b/dz2 + gamma*C = 0",
            "moisture": "Sm*d2q/dz2 - C = 0",
            "condensation": "C=H_smooth(q-qs)*(q-qs)/tau",
            "saturation": "qs=exp(alphaqs*(b-betaqs*z))",
        },
        "header_parameter_order": [
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
        ],
    }
    arrays = {
        "z": z,
        "b": b,
        "q": q,
        "header_parameters": header_parameters,
    }
    return metadata, arrays


def write_outputs(
    cfg: dict[str, Any],
    profile_path: Path,
    metadata_path: Path,
) -> dict[str, Any]:
    metadata, arrays = build_profile(cfg)
    profile_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    with profile_path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(
            f"{arrays['z'].size:d} "
            f"{metadata['drizzle_perturb_amp']:.17e}\n"
        )
        stream.write(
            " ".join(f"{value:.17e}" for value in arrays["header_parameters"])
            + "\n"
        )
        for z_value, b_value, q_value in zip(
            arrays["z"], arrays["b"], arrays["q"]
        ):
            stream.write(
                f"{z_value:.17e} {b_value:.17e} {q_value:.17e}\n"
            )

    metadata["profile_file"] = profile_path.name
    metadata["profile_sha256"] = sha256(profile_path)
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a converged case-specific drizzle base state"
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--params-json", type=Path)
    source.add_argument("--bou-in", type=Path)
    parser.add_argument(
        "--profile", type=Path, default=Path("drizzle_init.dat")
    )
    parser.add_argument(
        "--metadata", type=Path, default=Path("drizzle_init_meta.json")
    )
    parser.add_argument("--profile-n", type=int, default=DEFAULT_PROFILE_N)
    parser.add_argument(
        "--perturb-amp", type=float, default=DEFAULT_PERTURB_AMP
    )
    parser.add_argument(
        "--saturation-width",
        type=float,
        default=DEFAULT_SATURATION_WIDTH,
    )
    args = parser.parse_args()

    if args.params_json is not None:
        cfg = json.loads(args.params_json.read_text(encoding="utf-8"))
    else:
        cfg = config_from_bou(args.bou_in)
    cfg["drizzle_profile_n"] = args.profile_n
    cfg["drizzle_perturb_amp"] = args.perturb_amp
    cfg["drizzle_saturation_width"] = args.saturation_width
    metadata = write_outputs(cfg, args.profile, args.metadata)
    print(json.dumps(metadata, ensure_ascii=False))


if __name__ == "__main__":
    main()
