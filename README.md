# Primitive Xi Resolvent Audit

A finite-precision computational audit of a primitive resolvent S-fraction program related to the Riemann Hypothesis.

This repository studies the completed Riemann (\Xi)-function through the normalized critical-line transform

[
\Phi(u)=\frac{\Xi(1/2+i\sqrt{u})}{\Xi(1/2)}
]

and its primitive positive resolvent

[
R(u)=-\frac{\Phi'(u)}{\Phi(u)}.
]

The goal is to test whether (R(u)) behaves numerically like a positive-growth Markov/Stieltjes/S-fraction object:

[
R(u)=\int_0^\infty \frac{d\pi(x)}{1-ux},
\qquad d\pi(x)\ge0.
]

If such a representation could be proved analytically and globally without assuming the Riemann Hypothesis, it would force the poles of (R), hence the zeros of (\Phi), to lie on the positive real (u)-axis. This would correspond to zeros of (\Xi(s)) lying on the critical line.

This repository does **not** claim to prove the Riemann Hypothesis. It provides a reproducible finite computational audit of a possible S-fraction route.

---

## 1. Mathematical Motivation

Let

[
\Phi(u)=\frac{\Xi(1/2+i\sqrt{u})}{\Xi(1/2)}.
]

If the Riemann Hypothesis is true, then formally

[
\Phi(u)=\prod_n\left(1-\frac{u}{\gamma_n^2}\right),
]

where the nontrivial zeros of (\zeta(s)) are

[
\rho_n=\frac12\pm i\gamma_n.
]

Taking the logarithmic derivative gives

[
\frac{\Phi'(u)}{\Phi(u)}
========================

-\sum_n\frac{1}{\gamma_n^2-u}.
]

Therefore,

[
R(u)=-\frac{\Phi'(u)}{\Phi(u)}
==============================

\sum_n\frac{1}{\gamma_n^2-u}.
]

Equivalently, with

[
x_n=\frac1{\gamma_n^2}>0,
]

we have

[
R(u)=\sum_n\frac{x_n}{1-u x_n}.
]

This is the positive-growth Markov form. The computational question is whether this structure appears directly in the Taylor and S-fraction data of (R(u)).

---

## 2. What This Audit Tests

The script tests four positivity gates.

### Gate 1: Positive Taylor Coefficients

The series

[
R(u)=r_0+r_1u+r_2u^2+\cdots
]

is computed, and the script checks whether

[
r_n\ge0.
]

If

[
R(u)=\int_0^\infty \frac{d\pi(x)}{1-ux},
\qquad d\pi(x)\ge0,
]

then

[
r_n=\int_0^\infty x^n,d\pi(x)\ge0.
]

---

### Gate 2: Hankel Moment Positivity

The script checks the Hankel determinants

[
\det[r_{i+j}]_{i,j=0}^{N}\ge0
]

and

[
\det[r_{i+j+1}]_{i,j=0}^{N}\ge0.
]

These are necessary positivity conditions for ((r_n)) to be a Stieltjes moment sequence.

---

### Gate 3: Positive-Growth S-Fraction Tail Stripping

Starting with

[
F_0(u)=R(u),
]

the script recursively defines

[
a_n=F_n(0),
]

[
F_{n+1}(u)=\frac{1-a_n/F_n(u)}{u}.
]

The test checks whether

[
a_n>0
]

for all computed depths.

This is the positive-growth S-fraction orientation appropriate for kernels of the form

[
\frac{1}{1-ux}.
]

---

### Gate 4: Optional Pick/Herglotz Sign Test

For a positive-growth Markov transform,

[
R(z)=\int_0^\infty \frac{d\pi(x)}{1-zx},
]

one expects

[
\operatorname{Im}z>0
\quad\Longrightarrow\quad
\operatorname{Im}R(z)\ge0.
]

The script optionally samples the upper half-plane and checks this sign numerically.

---

## 3. Why the Primitive Resolvent Matters

Earlier tests on the curvature-rise quotient

[
\frac{K(\sqrt{u})-K(0)}{u}
]

showed positive Taylor coefficients but failed Hankel positivity and direct S-fraction tail stripping.

This suggests that the curvature quotient is not the primitive S-fraction object.

The successful object is instead

[
R(u)=-\frac{\Phi'(u)}{\Phi(u)}.
]

Thus the proposed hierarchy is:

[
\Phi(u)
]

[
\Downarrow
]

[
R(u)=-\frac{\Phi'(u)}{\Phi(u)}
]

[
\Downarrow
]

[
\text{positive Taylor coefficients}
]

[
\Downarrow
]

[
\text{Hankel moment positivity}
]

[
\Downarrow
]

[
\text{positive S-fraction tails}
]

[
\Downarrow
]

[
\text{positive spectral recursion}.
]

The curvature floor, if true, should be viewed as downstream from this primitive resolvent structure.

---

## 4. Installation

This script requires Python 3 and the following packages:

```bash
pip install mpmath tqdm
```

The script uses arbitrary-precision arithmetic through `mpmath`.

---

## 5. Usage

Baseline run:

```bash
python primitive_xi_resolvent_audit.py --dps 100 --max-gate 10 --depth 20 --tail-order 60
```

Stronger run with Pick test:

```bash
python primitive_xi_resolvent_audit.py --dps 120 --max-gate 12 --depth 25 --tail-order 80 --radius 4 --pick --pick-trials 300
```

Heavier run:

```bash
python primitive_xi_resolvent_audit.py --dps 160 --max-gate 16 --depth 35 --tail-order 120 --radius 4 --pick --pick-trials 500
```

Optional zero comparison:

```bash
python primitive_xi_resolvent_audit.py --dps 160 --max-gate 16 --depth 35 --tail-order 120 --radius 4 --pick --pick-trials 500 --zeros 10
```

---

## 6. Output Files

The script writes:

```text
primitive_xi_resolvent_audit.csv
primitive_xi_resolvent_audit_summary.json
```

The CSV contains detailed coefficient, Hankel, S-fraction, and Pick-test data.

The JSON summary records key parameters and pass/fail counts.

---

## 7. Expected Output

A successful run should report something like:

```text
moment violations: 0
hankel violations: 0
S-fraction coefficient failures: 0
Pick violations: 0
PASS: no numerical violations detected.
```

A representative stronger run found:

```text
moment violations: 0
hankel violations: 0
S-fraction coefficient failures: 0
Pick violations: 0
```

This supports the conjecture that the primitive positive resolvent

[
R(u)=-\frac{\Phi'(u)}{\Phi(u)}
]

behaves numerically like a positive-growth Markov/S-fraction object.

---

## 8. Main Conjecture

The numerical experiments motivate the following conjecture.

### Primitive Xi Resolvent S-Fraction Conjecture

Let

[
\Phi(u)=\frac{\Xi(1/2+i\sqrt{u})}{\Xi(1/2)}
]

and

[
R(u)=-\frac{\Phi'(u)}{\Phi(u)}.
]

Then

[
R(u)
]

admits a positive-growth Markov representation

[
R(u)=\int_0^\infty \frac{d\pi(x)}{1-ux},
\qquad d\pi(x)\ge0.
]

Equivalently, the coefficients

[
R(u)=\sum_{n\ge0}r_nu^n
]

form a positive Stieltjes moment sequence:

[
r_n\ge0,
]

[
\det[r_{i+j}]_{i,j=0}^{N}\ge0,
]

[
\det[r_{i+j+1}]_{i,j=0}^{N}\ge0
]

for every (N\ge0).

Equivalently, the positive-growth S-fraction tail coefficients generated by

[
F_{n+1}(u)=\frac{1-a_n/F_n(u)}{u},
\qquad a_n=F_n(0),
\qquad F_0=R,
]

satisfy

[
a_n>0
]

for every (n\ge0).

---

## 9. Path to the Riemann Hypothesis

If the global representation

[
R(u)=\int_0^\infty \frac{d\pi(x)}{1-ux},
\qquad d\pi(x)\ge0
]

could be proved without assuming RH, then the singularities of (R) would lie on the positive real (u)-axis.

But

[
R(u)=-\frac{\Phi'(u)}{\Phi(u)}
]

has singularities precisely where

[
\Phi(u)=0.
]

Therefore, such a theorem would force the zeros of (\Phi) to lie on the positive real (u)-axis.

Since

[
\Phi(u)=\frac{\Xi(1/2+i\sqrt{u})}{\Xi(1/2)},
]

a positive real zero

[
u=\gamma^2>0
]

corresponds to a zero

[
s=\frac12+i\gamma
]

of (\Xi(s)).

Thus the proof path would be:

[
R(u)\text{ has a global positive Markov/S-fraction representation}
]

[
\Longrightarrow
]

[
\Phi(u)\text{ has only positive real zeros}
]

[
\Longrightarrow
]

[
\Xi(s)\text{ has zeros only on } \operatorname{Re}(s)=1/2
]

[
\Longrightarrow
]

[
\text{RH}.
]

---

## 10. Important Limitations

This repository does **not** prove RH.

It provides finite numerical evidence for a proposed positivity structure.

The missing analytic step is:

[
R(u)=-\frac{d}{du}\log
\frac{\Xi(1/2+i\sqrt{u})}{\Xi(1/2)}
]

has a global positive Markov/S-fraction representation

[
R(u)=\int_0^\infty \frac{d\pi(x)}{1-ux},
\qquad d\pi(x)\ge0
]

without assuming the Riemann Hypothesis.

The computations verify only finite-depth, finite-precision consequences of this conjecture.

---

## 11. Novelty Status

The broad ingredients are classical:

* the Riemann (\Xi)-function;
* critical-line real-rootedness formulations of RH;
* Laguerre-Pólya theory;
* logarithmic derivatives of entire functions;
* moment positivity;
* Stieltjes/Markov transforms;
* S-fractions and continued fractions;
* Pick/Herglotz functions.

The potentially nonstandard aspect of this repository is the exact computational formulation:

[
\Phi(u)=\frac{\Xi(1/2+i\sqrt{u})}{\Xi(1/2)},
]

[
R(u)=-\frac{\Phi'(u)}{\Phi(u)},
]

together with the simultaneous audit of:

[
r_n>0,
]

[
\det[r_{i+j}]>0,
]

[
\det[r_{i+j+1}]>0,
]

positive S-fraction tail coefficients, and Pick positivity.

This should be understood as a computational formulation and conjectural route, not as a completed proof.

---

## 12. Repository Structure

Suggested structure:

```text
primitive_xi_resolvent_audit/
│
├── primitive_xi_resolvent_audit.py
├── README.md
```

The CSV and JSON output files may be regenerated by running the script.

---

## 13. Reproducibility Notes

The script computes Taylor coefficients using Cauchy coefficient extraction around (s=1/2). It then builds the (u)-series for

[
\Phi(u)
]

and formally computes

[
R(u)=-\Phi'(u)/\Phi(u).
]

All arithmetic is performed using `mpmath` at configurable precision.

For higher confidence:

* increase `--dps`;
* increase `--max-gate`;
* increase `--depth`;
* increase `--tail-order`;
* compare runs at different `--radius` values;
* enable `--pick`;
* inspect the generated CSV and JSON files.

---

## 14. Citation / Contact

Author: Buddhika Weerasooriya

---

## 15. Disclaimer

This repository is an experimental mathematical computation project.

It does not claim to prove the Riemann Hypothesis.

It identifies and tests a precise S-fraction/Markov positivity conjecture for a primitive resolvent of the Riemann (\Xi)-function. Proving that conjecture analytically and globally remains the central open problem.
