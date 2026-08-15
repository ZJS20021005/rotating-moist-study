#!/usr/bin/env python3
"""Local linear-stability scans for dry and moist rotating RB models.

The solver uses a primitive-variable generalized eigenvalue problem for one
horizontal Fourier mode.  It is designed as a transparent local baseline:

  * dry RB / rotating dry RB sanity checks;
  * current rotating moist Rainy-Benard parameters;
  * critical wavelength scaling with Ek.

The moist model follows the local DNS variables:

    m = b + gamma q,
    qs = exp(alphaqs * (b - betaqs z)),
    C' = chi/tau_cond * (q' - alphaqs qs b').

Here chi is evaluated on a one-dimensional conductive/drizzle base state.  The
default base state solves steady 1-D diffusion plus condensation for b and q.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import numpy as np
from numpy.typing import NDArray
from scipy import linalg, optimize


@dataclass(frozen=True)
class Params:
    pr: float = 0.7
    gamma: float = 1.1
    sm: float = 1.0
    alphaqs: float = 3.0
    betaqs: float = 1.0
    tau_cond: float = 1.0e-3
    b_bot: float = 0.0
    b_top: float = 0.0
    q_bot: float = 1.0
    q_top: float = 0.004978
    ubc_bot: str = "noslip"
    ubc_top: str = "freeslip"
    condensation_mode: str = "smooth"
    saturation_width: float = 1.0e-4
    linearize_switch: bool = False


def fd_matrices(n: int) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Return z, first derivative, second derivative on a uniform grid."""
    z = np.linspace(0.0, 1.0, n)
    dz = z[1] - z[0]
    d1 = np.zeros((n, n), dtype=float)
    d2 = np.zeros((n, n), dtype=float)

    for i in range(1, n - 1):
        d1[i, i - 1] = -0.5 / dz
        d1[i, i + 1] = 0.5 / dz
        d2[i, i - 1] = 1.0 / dz**2
        d2[i, i] = -2.0 / dz**2
        d2[i, i + 1] = 1.0 / dz**2

    # Second-order one-sided derivatives at the walls.
    d1[0, 0] = -3.0 / (2.0 * dz)
    d1[0, 1] = 4.0 / (2.0 * dz)
    d1[0, 2] = -1.0 / (2.0 * dz)
    d1[-1, -1] = 3.0 / (2.0 * dz)
    d1[-1, -2] = -4.0 / (2.0 * dz)
    d1[-1, -3] = 1.0 / (2.0 * dz)

    # One-sided second derivatives, mostly used by derivative BC checks.
    d2[0, 0] = 2.0 / dz**2
    d2[0, 1] = -5.0 / dz**2
    d2[0, 2] = 4.0 / dz**2
    d2[0, 3] = -1.0 / dz**2
    d2[-1, -1] = 2.0 / dz**2
    d2[-1, -2] = -5.0 / dz**2
    d2[-1, -3] = 4.0 / dz**2
    d2[-1, -4] = -1.0 / dz**2
    return z, d1, d2


def smooth_heaviside(x: NDArray[np.float64], width: float) -> NDArray[np.float64]:
    if width <= 0:
        return (x > 0).astype(float)
    arg = np.clip(x / width, -80.0, 80.0)
    return 0.5 * (1.0 + np.tanh(arg))


def smooth_heaviside_derivative(x: NDArray[np.float64], width: float) -> NDArray[np.float64]:
    if width <= 0:
        return np.zeros_like(x)
    arg = np.clip(x / width, -80.0, 80.0)
    sech2 = 1.0 / np.cosh(arg) ** 2
    return 0.5 * sech2 / width


def moist_base_state(n: int, params: Params) -> dict[str, NDArray[np.float64]]:
    """Solve the 1-D conductive/drizzle base state for b and q."""
    z, d1, d2 = fd_matrices(n)
    kappa = 1.0  # steady equation can be divided by the common diffusivity
    sm = params.sm

    b_linear = params.b_bot + (params.b_top - params.b_bot) * z
    q_linear = params.q_bot + (params.q_top - params.q_bot) * z
    y0 = np.r_[b_linear, q_linear]

    def residual(y: NDArray[np.float64]) -> NDArray[np.float64]:
        b = y[:n]
        q = y[n:]
        qs = np.exp(params.alphaqs * (b - params.betaqs * z))
        x = q - qs
        h = smooth_heaviside(x, params.saturation_width)
        c = h * x / params.tau_cond
        r = np.empty(2 * n)
        r[:n] = kappa * (d2 @ b) + params.gamma * c
        r[n:] = sm * kappa * (d2 @ q) - c
        r[0] = b[0] - params.b_bot
        r[n - 1] = b[-1] - params.b_top
        r[n] = q[0] - params.q_bot
        r[-1] = q[-1] - params.q_top
        return r

    sol = optimize.root(residual, y0, method="hybr", options={"xtol": 1.0e-10, "maxfev": 20000})
    if not sol.success:
        # A sharp condensation layer can make the root solve finicky.  The
        # linear conductive profile is still a useful first approximation, but
        # record the failure so downstream interpretation stays honest.
        b = b_linear
        q = q_linear
        converged = False
        message = sol.message
    else:
        b = sol.x[:n]
        q = sol.x[n:]
        converged = True
        message = sol.message
    qs = np.exp(params.alphaqs * (b - params.betaqs * z))
    supersat = q - qs
    chi = smooth_heaviside(supersat, params.saturation_width)
    return {
        "z": z,
        "b": b,
        "q": q,
        "qs": qs,
        "supersat": supersat,
        "chi": chi,
        "dbdz": d1 @ b,
        "dqdz": d1 @ q,
        "converged": np.array([float(converged)]),
        "message": np.array([message], dtype=object),
    }


@lru_cache(maxsize=32)
def cached_moist_base_state(n: int, params: Params) -> dict[str, NDArray[np.float64]]:
    """Cache the 1-D moist base state across neutral-curve root solves."""
    return moist_base_state(n, params)


def idx(var: int, n: int, i: int) -> int:
    return var * n + i


def impose_value(row: int, col: int, a: NDArray[np.complex128], b: NDArray[np.complex128]) -> None:
    a[row, :] = 0.0
    b[row, :] = 0.0
    a[row, col] = 1.0


def impose_derivative(
    row: int,
    var: int,
    n: int,
    d1: NDArray[np.float64],
    wall: int,
    a: NDArray[np.complex128],
    b: NDArray[np.complex128],
) -> None:
    a[row, :] = 0.0
    b[row, :] = 0.0
    for j in range(n):
        a[row, idx(var, n, j)] = d1[wall, j]


def max_growth_dry(
    ra: float,
    k: float,
    ek: float | None,
    n: int,
    pr: float,
    ubc_bot: str,
    ubc_top: str,
) -> tuple[float, complex]:
    z, d1, d2 = fd_matrices(n)
    nu = math.sqrt(pr / ra)
    kap = 1.0 / math.sqrt(ra * pr)
    inv_ro = 0.0 if ek is None or math.isinf(ek) else math.sqrt(pr / ra) / ek
    lap = d2 - k * k * np.eye(n)
    nvar = 5
    size = nvar * n
    a = np.zeros((size, size), dtype=np.complex128)
    bmat = np.zeros((size, size), dtype=np.complex128)

    for i in range(n):
        # u
        row = idx(0, n, i)
        a[row, idx(0, n, i)] += nu * lap[i, i]
        for j in range(n):
            if j != i:
                a[row, idx(0, n, j)] += nu * lap[i, j]
        a[row, idx(1, n, i)] += inv_ro
        a[row, idx(3, n, i)] += -1j * k
        bmat[row, idx(0, n, i)] = 1.0

        # v
        row = idx(1, n, i)
        for j in range(n):
            a[row, idx(1, n, j)] += nu * lap[i, j]
        a[row, idx(0, n, i)] += -inv_ro
        bmat[row, idx(1, n, i)] = 1.0

        # w
        row = idx(2, n, i)
        for j in range(n):
            a[row, idx(2, n, j)] += nu * lap[i, j]
            a[row, idx(3, n, j)] += -d1[i, j]
        a[row, idx(4, n, i)] += 1.0
        bmat[row, idx(2, n, i)] = 1.0

        # continuity
        row = idx(3, n, i)
        a[row, idx(0, n, i)] += 1j * k
        for j in range(n):
            a[row, idx(2, n, j)] += d1[i, j]

        # thermal anomaly, background T_z=-1
        row = idx(4, n, i)
        a[row, idx(2, n, i)] += 1.0
        for j in range(n):
            a[row, idx(4, n, j)] += kap * lap[i, j]
        bmat[row, idx(4, n, i)] = 1.0

    # Boundary conditions.
    for wall in (0, n - 1):
        impose_value(idx(2, n, wall), idx(2, n, wall), a, bmat)  # w=0
        impose_value(idx(4, n, wall), idx(4, n, wall), a, bmat)  # theta=0

    # Bottom horizontal velocity.
    if ubc_bot == "noslip":
        impose_value(idx(0, n, 0), idx(0, n, 0), a, bmat)
        impose_value(idx(1, n, 0), idx(1, n, 0), a, bmat)
    else:
        impose_derivative(idx(0, n, 0), 0, n, d1, 0, a, bmat)
        impose_derivative(idx(1, n, 0), 1, n, d1, 0, a, bmat)

    # Top horizontal velocity.
    if ubc_top == "noslip":
        impose_value(idx(0, n, n - 1), idx(0, n, n - 1), a, bmat)
        impose_value(idx(1, n, n - 1), idx(1, n, n - 1), a, bmat)
    else:
        impose_derivative(idx(0, n, n - 1), 0, n, d1, n - 1, a, bmat)
        impose_derivative(idx(1, n, n - 1), 1, n, d1, n - 1, a, bmat)

    # Pressure gauge.
    impose_value(idx(3, n, 0), idx(3, n, 0), a, bmat)

    vals = linalg.eigvals(a, bmat)
    vals = vals[np.isfinite(vals)]
    vals = vals[np.abs(vals) < 1.0e6]
    if vals.size == 0:
        return -np.inf, np.nan + 0j
    imax = int(np.argmax(vals.real))
    return float(vals[imax].real), complex(vals[imax])


def max_growth_moist(ra: float, k: float, ek: float, n: int, params: Params) -> tuple[float, complex, dict]:
    base = cached_moist_base_state(n, params)
    z, d1, d2 = fd_matrices(n)
    nu = math.sqrt(params.pr / ra)
    kap = 1.0 / math.sqrt(ra * params.pr)
    inv_ro = 0.0 if math.isinf(ek) else math.sqrt(params.pr / ra) / ek
    lap = d2 - k * k * np.eye(n)
    nvar = 6
    size = nvar * n
    a = np.zeros((size, size), dtype=np.complex128)
    bmat = np.zeros((size, size), dtype=np.complex128)

    dbdz = np.asarray(base["dbdz"], dtype=float)
    dqdz = np.asarray(base["dqdz"], dtype=float)
    qs = np.asarray(base["qs"], dtype=float)
    chi = np.asarray(base["chi"], dtype=float)
    supersat = np.asarray(base["supersat"], dtype=float)
    if params.linearize_switch:
        chi_eff = chi + supersat * smooth_heaviside_derivative(supersat, params.saturation_width)
    else:
        chi_eff = chi
    cb = -chi_eff * params.alphaqs * qs / params.tau_cond
    cq = chi_eff / params.tau_cond

    for i in range(n):
        # u
        row = idx(0, n, i)
        for j in range(n):
            a[row, idx(0, n, j)] += nu * lap[i, j]
        a[row, idx(1, n, i)] += inv_ro
        a[row, idx(3, n, i)] += -1j * k
        bmat[row, idx(0, n, i)] = 1.0

        # v
        row = idx(1, n, i)
        for j in range(n):
            a[row, idx(1, n, j)] += nu * lap[i, j]
        a[row, idx(0, n, i)] += -inv_ro
        bmat[row, idx(1, n, i)] = 1.0

        # w
        row = idx(2, n, i)
        for j in range(n):
            a[row, idx(2, n, j)] += nu * lap[i, j]
            a[row, idx(3, n, j)] += -d1[i, j]
        a[row, idx(4, n, i)] += 1.0
        bmat[row, idx(2, n, i)] = 1.0

        # continuity
        row = idx(3, n, i)
        a[row, idx(0, n, i)] += 1j * k
        for j in range(n):
            a[row, idx(2, n, j)] += d1[i, j]

        # b
        row = idx(4, n, i)
        a[row, idx(2, n, i)] += -dbdz[i]
        for j in range(n):
            a[row, idx(4, n, j)] += kap * lap[i, j]
        a[row, idx(4, n, i)] += params.gamma * cb[i]
        a[row, idx(5, n, i)] += params.gamma * cq[i]
        bmat[row, idx(4, n, i)] = 1.0

        # q
        row = idx(5, n, i)
        a[row, idx(2, n, i)] += -dqdz[i]
        for j in range(n):
            a[row, idx(5, n, j)] += params.sm * kap * lap[i, j]
        a[row, idx(4, n, i)] += -cb[i]
        a[row, idx(5, n, i)] += -cq[i]
        bmat[row, idx(5, n, i)] = 1.0

    for wall in (0, n - 1):
        impose_value(idx(2, n, wall), idx(2, n, wall), a, bmat)  # w=0
        impose_value(idx(4, n, wall), idx(4, n, wall), a, bmat)  # b'=0
        impose_value(idx(5, n, wall), idx(5, n, wall), a, bmat)  # q'=0

    if params.ubc_bot == "noslip":
        impose_value(idx(0, n, 0), idx(0, n, 0), a, bmat)
        impose_value(idx(1, n, 0), idx(1, n, 0), a, bmat)
    else:
        impose_derivative(idx(0, n, 0), 0, n, d1, 0, a, bmat)
        impose_derivative(idx(1, n, 0), 1, n, d1, 0, a, bmat)

    if params.ubc_top == "noslip":
        impose_value(idx(0, n, n - 1), idx(0, n, n - 1), a, bmat)
        impose_value(idx(1, n, n - 1), idx(1, n, n - 1), a, bmat)
    else:
        impose_derivative(idx(0, n, n - 1), 0, n, d1, n - 1, a, bmat)
        impose_derivative(idx(1, n, n - 1), 1, n, d1, n - 1, a, bmat)

    impose_value(idx(3, n, 0), idx(3, n, 0), a, bmat)
    vals = linalg.eigvals(a, bmat)
    vals = vals[np.isfinite(vals)]
    vals = vals[np.abs(vals) < 1.0e7]
    if vals.size == 0:
        return -np.inf, np.nan + 0j, base
    imax = int(np.argmax(vals.real))
    return float(vals[imax].real), complex(vals[imax]), base


def neutral_ra_for_k(
    model: str,
    k: float,
    ek: float,
    n: int,
    params: Params,
    ra_min: float,
    ra_max: float,
) -> tuple[float, float]:
    def growth(log_ra: float) -> float:
        ra = 10.0**log_ra
        if model == "dry":
            g, _ = max_growth_dry(ra, k, None if math.isinf(ek) else ek, n, params.pr, params.ubc_bot, params.ubc_top)
        else:
            g, _, _ = max_growth_moist(ra, k, ek, n, params)
        return g

    lo = math.log10(ra_min)
    hi = math.log10(ra_max)
    glo = growth(lo)
    ghi = growth(hi)
    # Expand upper bound if needed.
    while ghi < 0.0 and hi < 12.5:
        hi += 0.5
        ghi = growth(hi)
    if glo > 0.0:
        return 10.0**lo, glo
    if ghi < 0.0:
        return math.nan, ghi
    root = optimize.brentq(growth, lo, hi, xtol=1.0e-4, rtol=1.0e-4, maxiter=50)
    return 10.0**root, growth(root)


def scan_critical(
    model: str,
    eks: Iterable[float],
    k_values: NDArray[np.float64],
    n: int,
    params: Params,
    ra_min: float,
    ra_max: float,
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for ek in eks:
        best = None
        for k in k_values:
            ra_c, g = neutral_ra_for_k(model, float(k), float(ek), n, params, ra_min, ra_max)
            row = {
                "model": model,
                "Ek": float(ek),
                "k": float(k),
                "Ra_neutral": float(ra_c),
                "growth_at_root": float(g),
            }
            rows.append(row)
            if np.isfinite(ra_c) and (best is None or ra_c < best["Ra_neutral"]):
                best = row
        if best is not None:
            print(
                f"{model} Ek={ek:g}: best k={best['k']:.6g}, "
                f"Ra_c={best['Ra_neutral']:.6g}, lambda={2*math.pi/best['k']:.6g}",
                flush=True,
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize_scan(rows: list[dict[str, float]]) -> list[dict[str, float]]:
    by_key: dict[tuple[str, float], list[dict[str, float]]] = {}
    for row in rows:
        if not np.isfinite(row["Ra_neutral"]):
            continue
        by_key.setdefault((row["model"], row["Ek"]), []).append(row)
    summary = []
    for (model, ek), group in sorted(by_key.items(), key=lambda x: (x[0][0], x[0][1])):
        best = min(group, key=lambda r: r["Ra_neutral"])
        summary.append(
            {
                "model": model,
                "Ek": ek,
                "k_c": best["k"],
                "lambda_c": 2.0 * math.pi / best["k"],
                "Ra_c": best["Ra_neutral"],
                "n_k": len(group),
            }
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", choices=("dry", "moist"), default="dry")
    parser.add_argument("--n", type=int, default=45)
    parser.add_argument("--eks", type=str, default="inf,1e-2,3e-3,1e-3,3e-4,1e-4")
    parser.add_argument("--kmin", type=float, default=1.0)
    parser.add_argument("--kmax", type=float, default=80.0)
    parser.add_argument("--nk", type=int, default=36)
    parser.add_argument("--ra-min", type=float, default=1.0e2)
    parser.add_argument("--ra-max", type=float, default=1.0e10)
    parser.add_argument("--ubc-bot", choices=("noslip", "freeslip"), default="noslip")
    parser.add_argument("--ubc-top", choices=("noslip", "freeslip"), default="freeslip")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    eks = [math.inf if token == "inf" else float(token) for token in args.eks.split(",")]
    # Geometric k spacing is better for Ek scaling across decades.
    k_values = np.geomspace(args.kmin, args.kmax, args.nk)
    params = Params(ubc_bot=args.ubc_bot, ubc_top=args.ubc_top)
    rows = scan_critical(args.model, eks, k_values, args.n, params, args.ra_min, args.ra_max)
    write_csv(args.output / f"{args.model}_neutral_scan.csv", rows)
    summary = summarize_scan(rows)
    write_csv(args.output / f"{args.model}_critical_summary.csv", summary)
    (args.output / f"{args.model}_run_config.json").write_text(
        json.dumps(
            {
                "model": args.model,
                "n": args.n,
                "eks": args.eks,
                "kmin": args.kmin,
                "kmax": args.kmax,
                "nk": args.nk,
                "ra_min": args.ra_min,
                "ra_max": args.ra_max,
                "ubc_bot": args.ubc_bot,
                "ubc_top": args.ubc_top,
                "params": asdict(params),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
