#!/usr/bin/env python3
"""
primitive_xi_resolvent_audit.py

Computational audit for the primitive Xi-resolvent S-fraction program.

This script DOES NOT prove the Riemann Hypothesis. It provides a reproducible
finite-precision audit of the following theory candidate:

    Phi(u) = Xi(1/2 + i sqrt(u)) / Xi(1/2)
    R(u)   = - Phi'(u) / Phi(u)

Conjectural positive-growth Markov/S-fraction representation:

    R(u) = int_0^infty d pi(x) / (1 - u x),   d pi >= 0.

The script verifies finite computational gates:

    1. Taylor coefficient positivity:       r_n = [u^n] R(u) >= 0
    2. Hankel moment positivity:            det[r_{i+j}] >= 0
    3. shifted Hankel moment positivity:    det[r_{i+j+1}] >= 0
    4. positive-growth S-fraction stripping coefficients a_n > 0
    5. optional positive Pick/Herglotz sign: Im z > 0 => Im R(z) >= 0
    6. optional finite-zero comparison against first Riemann zeros

If all tested gates pass, this supports the primitive-resolvent S-fraction
program. It remains a finite numerical verification, not an unconditional proof.

Typical usage:

    python primitive_xi_resolvent_audit.py --dps 120 --max-gate 12 --depth 25 --tail-order 80 --radius 4 --pick --pick-trials 300

Stronger usage:

    python primitive_xi_resolvent_audit.py --dps 160 --max-gate 16 --depth 35 --tail-order 120 --radius 4 --pick --pick-trials 500 --zeros 10

Dependencies:

    pip install mpmath tqdm

Author: Buddhika Weerasooriya / generated with ChatGPT assistance
License suggestion: MIT
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
import time
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

import mpmath as mp

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    tqdm = None


# -----------------------------------------------------------------------------
# Utility
# -----------------------------------------------------------------------------


def progress(iterable: Iterable, total: Optional[int] = None, desc: str = ""):
    if tqdm is None:
        return iterable
    return tqdm(iterable, total=total, desc=desc, ncols=100)


def fmt(x, digits: int = 18) -> str:
    try:
        return mp.nstr(x, digits)
    except Exception:
        return str(x)


def complex_to_str(z, digits: int = 50) -> str:
    return fmt(z, digits)


def is_small_imag(z, tol: mp.mpf) -> bool:
    return abs(mp.im(z)) <= 1000 * tol


# -----------------------------------------------------------------------------
# Completed Riemann Xi
# -----------------------------------------------------------------------------


def xi(s: mp.mpc) -> mp.mpc:
    """Completed Riemann Xi function in the standard normalization."""
    return mp.mpf("0.5") * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


# -----------------------------------------------------------------------------
# Cauchy coefficient extraction
# -----------------------------------------------------------------------------


def xi_taylor_cauchy(max_order: int, radius: mp.mpf, samples: int) -> List[mp.mpc]:
    """
    Compute coefficients a_n of

        Xi(1/2 + q) = sum a_n q^n

    using the Cauchy integral formula on a circle |q|=radius.
    """
    s0 = mp.mpf("0.5")
    vals: List[mp.mpc] = []
    omegas: List[mp.mpc] = []

    for j in progress(range(samples), total=samples, desc="Cauchy samples"):
        theta = 2j * mp.pi * j / samples
        omega = mp.e ** theta
        omegas.append(omega)
        vals.append(xi(s0 + radius * omega))

    coeffs: List[mp.mpc] = []
    for n in progress(range(max_order + 1), total=max_order + 1, desc="Xi coeffs"):
        total = mp.mpc(0)
        for j in range(samples):
            total += vals[j] * (omegas[j] ** (-n))
        coeffs.append(total / samples / (radius ** n))

    return coeffs


# -----------------------------------------------------------------------------
# Phi and positive primitive resolvent R
# -----------------------------------------------------------------------------


def phi_series_from_xi_coeffs(a: Sequence[mp.mpc], max_order_u: int) -> List[mp.mpc]:
    """
    Build coefficients p_j of

        Phi(u) = Xi(1/2+i sqrt(u)) / Xi(1/2) = sum p_j u^j.

    Since q=i sqrt(u), q^(2j)=(-1)^j u^j, so

        p_j = a_{2j} (-1)^j / a_0.
    """
    a0 = a[0]
    return [a[2 * j] * ((-1) ** j) / a0 for j in range(max_order_u + 1)]


# -----------------------------------------------------------------------------
# Formal power series helpers
# -----------------------------------------------------------------------------


def trim_series(p: Sequence[mp.mpc], n: int) -> List[mp.mpc]:
    if len(p) >= n:
        return list(p[:n])
    return list(p) + [mp.mpc(0)] * (n - len(p))


def series_derivative(p: Sequence[mp.mpc]) -> List[mp.mpc]:
    if len(p) <= 1:
        return [mp.mpc(0)]
    return [(k + 1) * p[k + 1] for k in range(len(p) - 1)]


def series_inverse(p: Sequence[mp.mpc], n: int) -> List[mp.mpc]:
    p = trim_series(p, n)
    if abs(p[0]) == 0:
        raise ZeroDivisionError("series_inverse requires nonzero constant term")

    q = [mp.mpc(0)] * n
    q[0] = 1 / p[0]

    for k in range(1, n):
        s = mp.mpc(0)
        for j in range(1, k + 1):
            s += p[j] * q[k - j]
        q[k] = -s / p[0]

    return q


def series_mul(a: Sequence[mp.mpc], b: Sequence[mp.mpc], n: int) -> List[mp.mpc]:
    a = trim_series(a, n)
    b = trim_series(b, n)
    c = [mp.mpc(0)] * n

    for i in range(n):
        s = mp.mpc(0)
        for j in range(i + 1):
            s += a[j] * b[i - j]
        c[i] = s

    return c


def series_scale(p: Sequence[mp.mpc], c: mp.mpc) -> List[mp.mpc]:
    return [c * x for x in p]


def eval_series(coeffs: Sequence[mp.mpc], z: mp.mpc) -> mp.mpc:
    s = mp.mpc(0)
    for c in reversed(coeffs):
        s = s * z + c
    return s


def positive_resolvent_series(phi_coeffs: Sequence[mp.mpc], n: int) -> List[mp.mpc]:
    """
    Compute coefficients of

        R(u) = - Phi'(u)/Phi(u).

    Under the RH-formal product

        Phi(u) = prod_n (1 - u/gamma_n^2),

    this becomes

        R(u) = sum_n 1/(gamma_n^2-u)
             = sum_n x_n/(1-u x_n), x_n=1/gamma_n^2>0.
    """
    dphi = series_derivative(phi_coeffs)
    inv_phi = series_inverse(phi_coeffs, n)
    log_deriv = series_mul(dphi, inv_phi, n)
    return [-c for c in log_deriv]


# -----------------------------------------------------------------------------
# Moment and Hankel gates
# -----------------------------------------------------------------------------


@dataclass
class MomentRow:
    n: int
    value: mp.mpc
    status: str


@dataclass
class HankelRow:
    gate: int
    det_h0: mp.mpc
    det_h1: mp.mpc
    status_h0: str
    status_h1: str



def hankel_matrix(seq: Sequence[mp.mpc], shift: int, gate: int) -> mp.matrix:
    size = gate + 1
    M = mp.matrix(size)
    for i in range(size):
        for j in range(size):
            M[i, j] = seq[i + j + shift]
    return M


def test_positive_moments(
    R_coeffs: Sequence[mp.mpc], max_gate: int, tol: mp.mpf
) -> Tuple[List[mp.mpc], List[MomentRow], List[Tuple[int, mp.mpc]], List[HankelRow], List[Tuple[str, int, mp.mpc]]]:
    needed = 2 * max_gate + 2
    moments = list(R_coeffs[: needed + 1])

    moment_rows: List[MomentRow] = []
    moment_violations: List[Tuple[int, mp.mpc]] = []

    for n, val in enumerate(moments):
        re_val = mp.re(val)
        status = "PASS"
        if re_val < -tol or not is_small_imag(val, tol):
            status = "FAIL"
            moment_violations.append((n, val))
        moment_rows.append(MomentRow(n=n, value=val, status=status))

    hankel_rows: List[HankelRow] = []
    hankel_violations: List[Tuple[str, int, mp.mpc]] = []

    for gate in range(max_gate + 1):
        H0 = hankel_matrix(moments, 0, gate)
        H1 = hankel_matrix(moments, 1, gate)
        det0 = mp.det(H0)
        det1 = mp.det(H1)

        status0 = "PASS"
        status1 = "PASS"

        if mp.re(det0) < -tol or not is_small_imag(det0, tol):
            status0 = "FAIL"
            hankel_violations.append(("H0", gate, det0))
        if mp.re(det1) < -tol or not is_small_imag(det1, tol):
            status1 = "FAIL"
            hankel_violations.append(("H1", gate, det1))

        hankel_rows.append(
            HankelRow(
                gate=gate,
                det_h0=det0,
                det_h1=det1,
                status_h0=status0,
                status_h1=status1,
            )
        )

    return moments, moment_rows, moment_violations, hankel_rows, hankel_violations


# -----------------------------------------------------------------------------
# Positive-growth S-fraction stripping
# -----------------------------------------------------------------------------


@dataclass
class TailRow:
    depth: int
    a: mp.mpc
    status: str
    min_coeff: mp.mpc
    min_coeff_index: int
    coeff_violations: int
    error: str = ""



def strip_tail_positive_growth(F: Sequence[mp.mpc], n_out: int) -> Tuple[mp.mpc, List[mp.mpc]]:
    """
    Positive-growth S-fraction tail stripping.

    Given F_n(u)=a_n+b_n u+..., use

        F_{n+1}(u) = (1 - a_n/F_n(u))/u.

    Then F_{n+1}(0)=b_n/a_n if a_n != 0.
    """
    if abs(F[0]) == 0:
        raise ZeroDivisionError("Cannot strip tail: F_n(0)=0")

    needed = n_out + 1
    invF = series_inverse(F, needed)
    a = F[0]
    A_over_F = series_scale(invF, a)
    tail = [-A_over_F[k] for k in range(1, needed)]
    return a, tail


def strip_sfraction_positive_growth(
    F0: Sequence[mp.mpc], max_depth: int, tail_order: int, tol: mp.mpf
) -> List[TailRow]:
    tails: List[TailRow] = []
    coeffs = trim_series(F0, tail_order + max_depth + 2)

    for depth in progress(range(max_depth + 1), total=max_depth + 1, desc="S-strip"):
        a = coeffs[0]
        status = "PASS"
        if mp.re(a) <= tol or not is_small_imag(a, tol):
            status = "FAIL"

        coeff_violations = 0
        min_coeff = mp.inf
        min_coeff_index = -1
        check_len = min(len(coeffs), tail_order)

        for k in range(check_len):
            val = mp.re(coeffs[k])
            if val < min_coeff:
                min_coeff = val
                min_coeff_index = k
            if val < -tol:
                coeff_violations += 1

        row = TailRow(
            depth=depth,
            a=a,
            status=status,
            min_coeff=min_coeff,
            min_coeff_index=min_coeff_index,
            coeff_violations=coeff_violations,
        )
        tails.append(row)

        if depth == max_depth:
            break

        coeffs_needed_next = tail_order + (max_depth - depth) + 2
        try:
            _, coeffs = strip_tail_positive_growth(coeffs, coeffs_needed_next)
        except Exception as e:
            row.status = "FAIL"
            row.error = str(e)
            break

        if len(coeffs) < 2:
            break

    return tails


# -----------------------------------------------------------------------------
# Pick/Herglotz gate
# -----------------------------------------------------------------------------


@dataclass
class PickResult:
    violations: int
    min_im: mp.mpf
    worst_z: mp.mpc
    worst_value: mp.mpc



def pick_test_positive_growth(
    coeffs: Sequence[mp.mpc],
    trials: int,
    xmin: mp.mpf,
    xmax: mp.mpf,
    ymin: mp.mpf,
    ymax: mp.mpf,
    tol: mp.mpf,
    seed: int,
) -> PickResult:
    """
    Positive-growth Markov orientation:

        R(z)=int d pi(x)/(1-zx), d pi >= 0

    should satisfy

        Im z > 0 => Im R(z) >= 0.
    """
    random.seed(seed)
    violations = 0
    min_im = mp.inf
    worst_z = mp.mpc(0)
    worst_value = mp.mpc(0)

    for _ in progress(range(trials), total=trials, desc="Pick R"):
        xr = mp.e ** (mp.log(xmin) + random.random() * (mp.log(xmax) - mp.log(xmin)))
        yr = mp.e ** (mp.log(ymin) + random.random() * (mp.log(ymax) - mp.log(ymin)))
        if random.random() < 0.5:
            xr = -xr
        z = mp.mpc(xr, yr)
        val = eval_series(coeffs, z)
        im_val = mp.im(val)
        if im_val < min_im:
            min_im = im_val
            worst_z = z
            worst_value = val
        if im_val < -tol:
            violations += 1

    return PickResult(violations=violations, min_im=min_im, worst_z=worst_z, worst_value=worst_value)


# -----------------------------------------------------------------------------
# Optional zero comparison
# -----------------------------------------------------------------------------


@dataclass
class ZeroComparisonRow:
    n: int
    gamma: mp.mpf
    x: mp.mpf
    r_pred: List[mp.mpf]



def finite_zero_moments(num_zeros: int, max_n: int) -> List[mp.mpf]:
    """
    Approximate moments from the first num_zeros critical zeros:

        R(u) approx sum_k 1/(gamma_k^2-u)
             = sum_n (sum_k gamma_k^{-2n-2}) u^n.

    This is only a finite partial comparison, not a proof.
    """
    moments = [mp.mpf("0") for _ in range(max_n + 1)]
    for k in progress(range(1, num_zeros + 1), total=num_zeros, desc="Zeros"):
        gamma = mp.im(mp.zetazero(k))
        for n in range(max_n + 1):
            moments[n] += gamma ** (-(2 * n + 2))
    return moments


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------


def write_csv(
    path: str,
    phi_coeffs: Sequence[mp.mpc],
    R_coeffs: Sequence[mp.mpc],
    moment_rows: Sequence[MomentRow],
    hankel_rows: Sequence[HankelRow],
    tails: Sequence[TailRow],
    pick_res: Optional[PickResult],
    zero_moments: Optional[Sequence[mp.mpf]],
):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["section", "index", "value", "status", "extra"])

        for n, c in enumerate(phi_coeffs):
            w.writerow(["Phi_coeff", n, complex_to_str(c, 80), "", ""])
        for n, c in enumerate(R_coeffs):
            w.writerow(["R_coeff", n, complex_to_str(c, 80), "", ""])
        for row in moment_rows:
            w.writerow(["positive_moment_r_n", row.n, complex_to_str(row.value, 80), row.status, ""])
        for row in hankel_rows:
            w.writerow(["H0_det_r_i+j", row.gate, complex_to_str(row.det_h0, 80), row.status_h0, ""])
            w.writerow(["H1_det_r_i+j+1", row.gate, complex_to_str(row.det_h1, 80), row.status_h1, ""])
        for row in tails:
            w.writerow(["s_fraction_a_n", row.depth, complex_to_str(row.a, 80), row.status, row.error])
            w.writerow([
                "tail_coeff_min",
                row.depth,
                complex_to_str(row.min_coeff, 80),
                "",
                f"index={row.min_coeff_index}; coeff_viol={row.coeff_violations}",
            ])
        if pick_res is not None:
            w.writerow(["pick_positive_growth", "violations", pick_res.violations, "", ""])
            w.writerow(["pick_positive_growth", "min_im", fmt(pick_res.min_im, 80), "", ""])
            w.writerow(["pick_positive_growth", "worst_z", complex_to_str(pick_res.worst_z, 80), "", ""])
            w.writerow(["pick_positive_growth", "worst_value", complex_to_str(pick_res.worst_value, 80), "", ""])
        if zero_moments is not None:
            for n, val in enumerate(zero_moments):
                w.writerow(["finite_zero_moment", n, fmt(val, 80), "", "partial zero sum"])


def write_json_summary(path: str, summary: dict):
    def convert(obj):
        if isinstance(obj, mp.mpf):
            return fmt(obj, 80)
        if isinstance(obj, mp.mpc):
            return complex_to_str(obj, 80)
        return obj

    with open(path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=convert)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Finite computational audit for the primitive Xi positive-resolvent S-fraction program."
    )

    parser.add_argument("--dps", type=int, default=120, help="mpmath decimal precision")
    parser.add_argument("--max-gate", type=int, default=12, help="largest Hankel gate N")
    parser.add_argument("--depth", type=int, default=25, help="S-fraction stripping depth")
    parser.add_argument("--tail-order", type=int, default=80, help="number of tail coefficients to retain/check")
    parser.add_argument("--radius", type=str, default="4.0", help="Cauchy extraction radius in q-plane")
    parser.add_argument("--samples", type=int, default=None, help="Cauchy sample count; auto if omitted")
    parser.add_argument("--tol", type=str, default="1e-40", help="failure tolerance")

    parser.add_argument("--pick", action="store_true", help="run optional Pick/Herglotz sign test")
    parser.add_argument("--pick-trials", type=int, default=300)
    parser.add_argument("--pick-xmin", type=str, default="1e-6")
    parser.add_argument("--pick-xmax", type=str, default="0.05")
    parser.add_argument("--pick-ymin", type=str, default="1e-6")
    parser.add_argument("--pick-ymax", type=str, default="0.05")
    parser.add_argument("--seed", type=int, default=1234)

    parser.add_argument("--zeros", type=int, default=0, help="optional compare with first N zeta zeros")
    parser.add_argument("--csv", type=str, default="primitive_xi_resolvent_audit.csv")
    parser.add_argument("--json", type=str, default="primitive_xi_resolvent_audit_summary.json")

    args = parser.parse_args()
    mp.mp.dps = args.dps
    radius = mp.mpf(args.radius)
    tol = mp.mpf(args.tol)

    max_moment_needed = 2 * args.max_gate + 2
    max_coeff = max(args.tail_order + args.depth + 8, max_moment_needed + 8)
    max_xi_order = 2 * max_coeff

    if args.samples is None:
        samples = 4
        while samples < 8 * (max_xi_order + 1):
            samples *= 2
    else:
        samples = args.samples

    print("=" * 120)
    print("PRIMITIVE XI POSITIVE-RESOLVENT S-FRACTION AUDIT")
    print("=" * 120)
    print("This is a finite numerical audit, not an unconditional proof of RH.")
    print()
    print("Primitive kernel:")
    print("  Phi(u) = Xi(1/2+i sqrt(u)) / Xi(1/2)")
    print("Positive primitive resolvent:")
    print("  R(u) = - Phi'(u)/Phi(u)")
    print()
    print("Conjectural positive-growth representation:")
    print("  R(u) = int_0^infty d pi(x)/(1-u x), d pi >= 0")
    print("-" * 120)
    print(f"dps          = {args.dps}")
    print(f"max_gate     = {args.max_gate}")
    print(f"depth        = {args.depth}")
    print(f"tail_order   = {args.tail_order}")
    print(f"max_coeff    = {max_coeff}")
    print(f"max_xi_order = {max_xi_order}")
    print(f"radius       = {radius}")
    print(f"samples      = {samples}")
    print(f"tol          = {tol}")
    print(f"pick         = {args.pick}")
    print(f"zeros        = {args.zeros}")
    print(f"csv          = {args.csv}")
    print(f"json         = {args.json}")
    print("=" * 120)

    t0 = time.time()

    print("\nComputing Xi Taylor coefficients...")
    xi_coeffs = xi_taylor_cauchy(max_xi_order, radius, samples)

    print("Building Phi(u) coefficients...")
    phi_coeffs = phi_series_from_xi_coeffs(xi_coeffs, max_coeff)

    print("Computing R(u)=-Phi'/Phi series...")
    R_coeffs = positive_resolvent_series(phi_coeffs, max_coeff)

    print("\n" + "=" * 120)
    print("PHI / R CHECK")
    print("=" * 120)
    print(f"Phi(0) coeff p_0 = {fmt(phi_coeffs[0], 60)}")
    print(f"Phi'(0) p_1      = {fmt(phi_coeffs[1], 60)}")
    print(f"R(0)=-Phi'(0)    = {fmt(R_coeffs[0], 60)}")

    print("\n" + "=" * 120)
    print("GATE 1/2: POSITIVE MOMENTS AND HANKEL POSITIVITY")
    print("=" * 120)
    moments, moment_rows, moment_violations, hankel_rows, hankel_violations = test_positive_moments(
        R_coeffs, args.max_gate, tol
    )

    print("\nMoments r_n = [u^n] R(u):")
    for row in moment_rows:
        print(
            f"r_{row.n:<3d} = {fmt(mp.re(row.value), 36):>46s} "
            f"imag={fmt(mp.im(row.value), 8):>14s} {row.status}"
        )
    print("-" * 120)
    print(f"moment violations: {len(moment_violations)}")

    print("\nHankel gates:")
    for row in hankel_rows:
        print("-" * 120)
        print(f"gate N={row.gate}")
        print(
            f"det[r_i+j]     = {fmt(mp.re(row.det_h0), 36):>46s} "
            f"imag={fmt(mp.im(row.det_h0), 8):>14s} {row.status_h0}"
        )
        print(
            f"det[r_i+j+1]   = {fmt(mp.re(row.det_h1), 36):>46s} "
            f"imag={fmt(mp.im(row.det_h1), 8):>14s} {row.status_h1}"
        )
    print("-" * 120)
    print(f"hankel violations: {len(hankel_violations)}")

    print("\n" + "=" * 120)
    print("GATE 3: POSITIVE-GROWTH S-FRACTION TAIL STRIPPING")
    print("=" * 120)
    tails = strip_sfraction_positive_growth(R_coeffs, args.depth, args.tail_order, tol)

    sfailures = 0
    coeff_failures = 0
    for row in tails:
        if row.status != "PASS":
            sfailures += 1
        coeff_failures += row.coeff_violations
        print(
            f"n={row.depth:<3d} a_n={fmt(mp.re(row.a), 42):>50s} "
            f"imag={fmt(mp.im(row.a), 8):>14s} "
            f"min_coeff={fmt(row.min_coeff, 18):>24s} "
            f"coeff_viol={row.coeff_violations:<3d} {row.status}"
        )
        if row.error:
            print(f"    error: {row.error}")
    print("-" * 120)
    print(f"S-fraction coefficient failures: {sfailures}")
    print(f"tail coefficient sign violations: {coeff_failures}")

    pick_res: Optional[PickResult] = None
    if args.pick:
        print("\n" + "=" * 120)
        print("GATE 4: OPTIONAL POSITIVE PICK TEST")
        print("=" * 120)
        pick_res = pick_test_positive_growth(
            coeffs=R_coeffs[: args.tail_order],
            trials=args.pick_trials,
            xmin=mp.mpf(args.pick_xmin),
            xmax=mp.mpf(args.pick_xmax),
            ymin=mp.mpf(args.pick_ymin),
            ymax=mp.mpf(args.pick_ymax),
            tol=tol,
            seed=args.seed,
        )
        print(f"Pick violations : {pick_res.violations}")
        print(f"worst z         : {fmt(pick_res.worst_z, 30)}")
        print(f"worst R(z)      : {fmt(pick_res.worst_value, 40)}")
        print(f"min Im R        : {fmt(pick_res.min_im, 40)}")

    zero_moments: Optional[List[mp.mpf]] = None
    if args.zeros > 0:
        print("\n" + "=" * 120)
        print("OPTIONAL FINITE ZERO COMPARISON")
        print("=" * 120)
        compare_n = min(10, len(R_coeffs) - 1)
        zero_moments = finite_zero_moments(args.zeros, compare_n)
        for n in range(compare_n + 1):
            print(
                f"n={n:<3d} R_coeff={fmt(mp.re(R_coeffs[n]), 30):>38s} "
                f"partial_zero_sum={fmt(zero_moments[n], 30):>38s}"
            )
        print("Note: finite zero sums are partial sums and are expected to be below full coefficients.")

    write_csv(args.csv, phi_coeffs, R_coeffs, moment_rows, hankel_rows, tails, pick_res, zero_moments)

    total_failures = (
        len(moment_violations)
        + len(hankel_violations)
        + sfailures
        + coeff_failures
        + (pick_res.violations if pick_res is not None else 0)
    )

    elapsed = time.time() - t0

    summary = {
        "status": "PASS" if total_failures == 0 else "FAIL_OR_WARNING",
        "total_failures": total_failures,
        "moment_violations": len(moment_violations),
        "hankel_violations": len(hankel_violations),
        "sfraction_coefficient_failures": sfailures,
        "tail_coefficient_sign_violations": coeff_failures,
        "pick_violations": pick_res.violations if pick_res is not None else None,
        "dps": args.dps,
        "max_gate": args.max_gate,
        "depth": args.depth,
        "tail_order": args.tail_order,
        "radius": fmt(radius, 30),
        "samples": samples,
        "elapsed_seconds": elapsed,
        "important_note": "Finite numerical audit only; this does not prove RH or the global Markov/S-fraction representation.",
        "theory_target": "Prove globally and unconditionally that R(u)=-Phi'(u)/Phi(u) admits R(u)=int d pi(x)/(1-u x), d pi>=0.",
    }
    write_json_summary(args.json, summary)

    print("\n" + "=" * 120)
    print("SUMMARY")
    print("=" * 120)
    if total_failures == 0:
        print("PASS: no numerical violations detected.")
        print()
        print("Interpretation:")
        print("  This supports that the corrected primitive resolvent")
        print("    R(u)=-Phi'(u)/Phi(u)")
        print("  behaves as a positive-growth S-fraction / Markov-type object")
        print("  over the tested finite gates and sampled region.")
    else:
        print(f"FAIL/WARNING: {total_failures} possible violations detected.")
        print()
        print("Interpretation:")
        print("  The tested finite gates did not all pass. Inspect the CSV/JSON outputs.")

    print()
    print("Mathematical caution:")
    print("  A finite numerical pass is not a proof of RH. The analytic target is to prove")
    print("  the global positive Markov/S-fraction representation for R without assuming RH.")
    print()
    print(f"Elapsed seconds: {elapsed:.2f}")
    print(f"CSV written to : {args.csv}")
    print(f"JSON written to: {args.json}")
    print("=" * 120)

    return 0 if total_failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
