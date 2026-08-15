#!/usr/bin/env bash
set -u

PY=/share/org/SHUTUANL/shu_zhangjs/rainy\ model/rotating_case/postprocess/strict_force_balance_20260805/compute_fp3d_from_pressure_movies.py
OUT=/share/org/SHUTUANL/shu_zhangjs/rainy\ model/rotating_case/postprocess/strict_force_balance_20260805/fp3d
ROOT=/share/org/SHUTUANL/shu_zhangjs/rainy\ model/rotating_case/Pr0p7/Ra8e6

run_one() {
    ek="$1"
    run="$2"
    mkdir -p "$OUT/$ek"
    python3 "$PY" "$ROOT/$run" --output-dir "$OUT/$ek" --zmin 0.1 --zmax 0.9
}

run_one Ek1p5e-4 'Ek1p5e-4/AR4/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti1/conti_strict_force_500/run'
run_one Ek1e-3 'Ek1e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti1/conti_strict_force_500/run'
run_one Ek2e-3 'Ek2e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65/conti_strict_force_500/run'
run_one Ek3e-3 'Ek3e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti3/conti_strict_force_500/run'
run_one Ek5e-3 'Ek5e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti2/conti1/conti_strict_force_500/run'
run_one Ek7e-3 'Ek7e-3/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti2/conti1/conti_strict_force_500/run'
run_one Ek1e-2 'Ek1e-2/AR16/Beta1p02/qbot0p5_qtop0p004978/N257x257x65/conti2/conti1/conti_strict_force_500/run'
run_one Ek3e-2 'Ek3e-2/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65/conti2/conti1/conti_strict_force_500/run'
run_one Ek5e-2 'Ek5e-2/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65/conti_strict_force_500/run'
run_one Ek1e-1 'Ek1e-1/AR16/Beta1p02/qbot0p5_qtop0p004978/N385x385x65/conti1/conti_strict_force_500/run'
