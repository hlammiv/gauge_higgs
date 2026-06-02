#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Ferrenberg-Swendsen / WHAM / MBAR histogram reweighting (pure numpy).

PURE post-processing of Monte Carlo time series -- no lattice, no group, no
dimension. Given measurements made at one or several reference couplings, it
estimates expectation values <O> at arbitrary TARGET couplings by reweighting,
so a small set of simulated (beta, kappa, ...) points can trace a transition
line / ridge and localize a triple point by interpolation.

This is the Python mirror of ``src/measure/reweight.hpp``. It implements the
SAME math (same log-sum-exp stabilization, same WHAM/MBAR fixed point, same
delete-1 jackknife) and is intended to be numerically equivalent.

MODEL / CONVENTION (fixed; the self-tests pin the exact analytic identities)
----------------------------------------------------------------------------
The Monte Carlo weight of a configuration is exp(-S) with the action LINEAR in
the couplings,

    S(config; lambda) = sum_i lambda_i E_i(config),

where ``lambda`` is the vector of couplings and ``E_i`` are the conjugate
"energy" observables. The CALLER fixes the signs so that this holds.

For the U(1) gauge+Higgs campaign consumed here (see ``src/u1_scan.cpp`` and
the time-series header), the columns store ``A`` and ``B`` per trajectory with

    weight = exp(-S),  S = beta*A + kappa*(-B) + on-site const,
    A = plaq_energy_sum = sum_plaq (1 - cos theta_pl)            [conj to beta]
    B = hop_energy_sum  = sum_{x,mu} 2 Re[conj(phi) e^{i q theta} phi]
        so that the energy conjugate to kappa is  E_kappa = -B.

Therefore, in the generic ``lambda / E_i`` framework,

    lambda = (beta, kappa),   E_1 = A,   E_2 = -B.

Reweighting a sample from reference (beta0, kappa0) to target (beta, kappa)
multiplies its RELATIVE weight by

    exp( -(S(lambda) - S(lambda0)) )
      = exp( -[(beta - beta0) * A + (kappa - kappa0) * (-B)] )
      = exp( -(beta - beta0) * A + (kappa - kappa0) * B ).

The on-site terms are fixed across the scan and cancel in every ratio below.

NUMERICAL STABILITY
-------------------
The reweighting exponents can be large and positive, so exp() of them would
overflow. EVERY sum-of-exponentials below uses the log-sum-exp trick: shift by
the maximum exponent, exponentiate the (now <= 0) shifted exponents, sum, and
add the max back in log space. We never form exp() of a large positive number
directly.

PUBLIC API
----------
    log_sum_exp(a)
    reweighted_mean(O, A, B, beta0, kappa0, beta, kappa)
    reweighted_susceptibility(O, A, B, V, beta0, kappa0, beta, kappa)
    reweighted_mean_jackknife(O, A, B, beta0, kappa0, beta, kappa, nblocks=16)
    Ensemble(beta0, kappa0, A, B, O)
    multi_histogram(ensembles, beta, kappa, observable=...)
    chi_surface(ensembles, V, beta_mesh, kappa_mesh, observable=...)
    load_ts(path)

Run ``python3 scripts/reweight.py --selftest`` for the analytic checks.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence, Union

import numpy as np

ArrayLike = Union[Sequence[float], np.ndarray]


# ---------------------------------------------------------------------------
#  Numerically stable log-sum-exp:  log( sum_t exp(a_t) ).
# ---------------------------------------------------------------------------
def log_sum_exp(a: ArrayLike) -> float:
    """Return ``log(sum_t exp(a_t))`` computed stably; ``-inf`` for empty input.

    Mirrors ``gh::log_sum_exp`` in ``reweight.hpp``: shift by ``max(a)`` so every
    exponentiated term lies in ``(0, 1]`` (no overflow), sum, add the max back in
    log space. An empty input returns ``-inf`` (the additive identity of a
    log-sum). If the max is non-finite (all ``-inf``, or a ``+inf``), it is passed
    through unchanged.
    """
    a = np.asarray(a, dtype=float).ravel()
    if a.size == 0:
        return -np.inf
    m = np.max(a)
    if not np.isfinite(m):
        return float(m)  # all -inf -> -inf;  any +inf -> +inf
    s = np.sum(np.exp(a - m))
    return float(m + np.log(s))


def _weighted_mean_exp(o: np.ndarray, a: np.ndarray) -> float:
    """Stable ``( sum_t o_t exp(a_t) ) / ( sum_t exp(a_t) )``.

    Shifts BOTH sums by ``max(a)`` so the common factor ``exp(max)`` cancels.
    ``o_t`` may have any sign. Returns NaN for empty input. Mirrors
    ``reweight_detail::weighted_mean_exp``.
    """
    o = np.asarray(o, dtype=float).ravel()
    a = np.asarray(a, dtype=float).ravel()
    if a.size == 0:
        return float("nan")
    m = np.max(a)
    if not np.isfinite(m):
        # Degenerate weights; fall back to plain mean to avoid 0/0 nonsense.
        return float(np.mean(o)) if o.size else float("nan")
    w = np.exp(a - m)  # in (0, 1]
    den = np.sum(w)
    num = np.sum(o * w)
    return float(num / den)


# ---------------------------------------------------------------------------
#  Convention helper: per-sample log relative weight a_t for a single ensemble.
#
#  In the (beta, kappa) -> (A, -B) framework,
#      a_t = -(beta - beta0) * A_t - (kappa - kappa0) * (-B_t)
#          = -(beta - beta0) * A_t + (kappa - kappa0) * B_t .
#  This is exactly  a_t = -sum_i (lambda_i - lambda0_i) E_i(t)  with
#  E_1 = A (conj to beta) and E_2 = -B (conj to kappa).
# ---------------------------------------------------------------------------
def _log_weights(A: np.ndarray, B: np.ndarray,
                 beta0: float, kappa0: float,
                 beta: float, kappa: float) -> np.ndarray:
    A = np.asarray(A, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    if A.shape != B.shape:
        raise ValueError("A and B must have the same shape")
    return -(beta - beta0) * A + (kappa - kappa0) * B


# ---------------------------------------------------------------------------
#  SINGLE-HISTOGRAM (Ferrenberg-Swendsen) reweighting of a mean.
# ---------------------------------------------------------------------------
def reweighted_mean(O: ArrayLike, A: ArrayLike, B: ArrayLike,
                    beta0: float, kappa0: float,
                    beta: float, kappa: float) -> float:
    """Single-histogram FS estimate of ``<O>`` at target ``(beta, kappa)``.

    ``O``, ``A``, ``B`` are 1-D arrays over samples taken at the reference
    coupling ``(beta0, kappa0)``. The estimator is

        a_t     = -(beta - beta0) * A_t + (kappa - kappa0) * B_t
        <O>     = ( sum_t O_t e^{a_t} ) / ( sum_t e^{a_t} ),

    evaluated with log-sum-exp (shift by ``max_t a_t``). Returns NaN on empty /
    mismatched input.
    """
    O = np.asarray(O, dtype=float).ravel()
    A = np.asarray(A, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    if O.size == 0 or O.shape != A.shape or O.shape != B.shape:
        return float("nan")
    a = _log_weights(A, B, beta0, kappa0, beta, kappa)
    return _weighted_mean_exp(O, a)


# ---------------------------------------------------------------------------
#  SINGLE-HISTOGRAM reweighted susceptibility (variance * volume).
# ---------------------------------------------------------------------------
def reweighted_susceptibility(O: ArrayLike, A: ArrayLike, B: ArrayLike, V: float,
                              beta0: float, kappa0: float,
                              beta: float, kappa: float) -> float:
    """Reweighted susceptibility ``V * (<O^2> - <O>^2)`` at ``(beta, kappa)``.

    Both moments are reweighted with the SAME per-sample weights, then combined.
    This is the connected variance of ``O`` times the volume ``V`` -- the
    standard observable used to locate transitions (peaks / ridges). Returns NaN
    on empty / mismatched input.
    """
    O = np.asarray(O, dtype=float).ravel()
    A = np.asarray(A, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    if O.size == 0 or O.shape != A.shape or O.shape != B.shape:
        return float("nan")
    a = _log_weights(A, B, beta0, kappa0, beta, kappa)
    mean = _weighted_mean_exp(O, a)
    mean2 = _weighted_mean_exp(O * O, a)
    var = mean2 - mean * mean
    return float(V * var)


# ---------------------------------------------------------------------------
#  SINGLE-HISTOGRAM mean + delete-1 BLOCK jackknife error.
# ---------------------------------------------------------------------------
def reweighted_mean_jackknife(O: ArrayLike, A: ArrayLike, B: ArrayLike,
                              beta0: float, kappa0: float,
                              beta: float, kappa: float,
                              nblocks: int = 16) -> tuple[float, float]:
    """Return ``(mean, stderr)`` for the reweighted mean via a block jackknife.

    The samples are partitioned into ``nblocks`` contiguous blocks (autocorrelated
    MC data is best decimated in blocks rather than single samples). Each
    leave-one-block-out subset gives a replica ``theta_k``; the delete-1
    jackknife error is

        err^2 = (m - 1) / m * sum_k (theta_k - theta_bar)^2 ,

    with ``m`` the number of (non-empty) blocks. The reported ``mean`` is the
    full-sample reweighted mean. With fewer than 2 usable blocks the error is 0.
    Returns ``(nan, 0)`` on empty / mismatched input.
    """
    O = np.asarray(O, dtype=float).ravel()
    A = np.asarray(A, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    if O.size == 0 or O.shape != A.shape or O.shape != B.shape:
        return (float("nan"), 0.0)

    a = _log_weights(A, B, beta0, kappa0, beta, kappa)
    mean = _weighted_mean_exp(O, a)

    n = O.size
    m = int(max(1, min(nblocks, n)))
    if m < 2:
        return (mean, 0.0)

    # Contiguous block index for every sample (numpy array_split is balanced).
    block_of = np.empty(n, dtype=int)
    for bidx, idx in enumerate(np.array_split(np.arange(n), m)):
        block_of[idx] = bidx

    replicas = []
    for k in range(m):
        keep = block_of != k
        if not np.any(keep):
            continue
        replicas.append(_weighted_mean_exp(O[keep], a[keep]))
    replicas = np.asarray(replicas, dtype=float)
    mm = replicas.size
    if mm < 2:
        return (mean, 0.0)
    tbar = float(np.mean(replicas))
    sw = float(np.sum((replicas - tbar) ** 2))
    err = float(np.sqrt((mm - 1) / mm * sw))
    return (mean, err)


# ---------------------------------------------------------------------------
#  Ensemble container: one reference ensemble for the multi-histogram solver.
# ---------------------------------------------------------------------------
@dataclass
class Ensemble:
    """One reference ensemble sampled at ``(beta0, kappa0)``.

    ``A``, ``B``, ``O`` are 1-D time series of equal length. ``A`` is the energy
    conjugate to beta and ``-B`` the energy conjugate to kappa (see module
    docstring). ``O`` is the observable to reweight; if not supplied it defaults
    to zeros (still usable to solve the WHAM free energies, but ``<O>`` would be
    0 unless an ``observable=`` override is passed to the evaluator).
    """
    beta0: float
    kappa0: float
    A: np.ndarray
    B: np.ndarray
    O: np.ndarray = field(default=None)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        self.A = np.asarray(self.A, dtype=float).ravel()
        self.B = np.asarray(self.B, dtype=float).ravel()
        if self.A.shape != self.B.shape:
            raise ValueError("Ensemble A and B must have the same length")
        if self.O is None:
            self.O = np.zeros_like(self.A)
        else:
            self.O = np.asarray(self.O, dtype=float).ravel()
            if self.O.shape != self.A.shape:
                raise ValueError("Ensemble O must match the length of A, B")

    @property
    def n_samples(self) -> int:
        return self.A.size


# Type for an `observable` selector: either a name (resolved against the
# Ensemble fields A/B/O) or a callable mapping (A, B, O) flattened arrays -> O.
ObservableSpec = Union[str, Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray], None]


@dataclass
class WhamSolution:
    """Result of the WHAM/MBAR self-consistency solve.

    f         : dimensionless free energies f_k, with f_0 anchored to 0.
    g         : per-sample log mixture denominator g(n), length M_total.
    A, B, O   : flattened energies / default observable over all M_total samples.
    residual  : achieved max_k |df_k| at the last sweep.
    iters     : number of sweeps performed.
    converged : residual < tol within max_iter.
    """
    f: np.ndarray
    g: np.ndarray
    A: np.ndarray
    B: np.ndarray
    O: np.ndarray
    residual: float
    iters: int
    converged: bool

    @property
    def n_total(self) -> int:
        return self.O.size


def _resolve_observable(spec: ObservableSpec,
                        A: np.ndarray, B: np.ndarray, O: np.ndarray) -> np.ndarray:
    """Map an `observable` spec to a flattened observable array over all samples.

    ``None``        -> the default (flattened ensemble O).
    a field name    -> "A", "B", or "O".
    a callable      -> spec(A, B, O), returning a length-M array.
    """
    if spec is None or spec == "O":
        return O
    if isinstance(spec, str):
        table = {"A": A, "B": B, "O": O}
        if spec not in table:
            raise ValueError(f"unknown observable name {spec!r}; use 'A','B','O' or a callable")
        return table[spec]
    if callable(spec):
        out = np.asarray(spec(A, B, O), dtype=float).ravel()
        if out.shape != O.shape:
            raise ValueError("observable callable must return one value per sample")
        return out
    raise TypeError("observable must be None, a field name, or a callable")


# ---------------------------------------------------------------------------
#  MULTI-HISTOGRAM (WHAM / MBAR).
#
#  Combine K reference ensembles sampled at couplings {(beta0_k, kappa0_k)} with
#  N_k samples each into one optimal estimator. Index every sample globally by
#  n = 0..M-1. Let
#      u_k(n) = beta0_k * A(n) + kappa0_k * (-B(n))                            (S of n at k)
#  be the reduced action of sample n evaluated at ensemble k's couplings.
#
#  The dimensionless free energies f_k (defined up to an additive constant)
#  satisfy the MBAR/WHAM fixed point
#      exp(-f_k) = sum_n  exp(-u_k(n)) / ( sum_l N_l exp(f_l - u_l(n)) ).
#  With the per-sample log mixture denominator
#      g(n) = log_sum_exp_l [ log N_l + f_l - u_l(n) ],
#  the update is
#      f_k <- -log_sum_exp_n [ -u_k(n) - g(n) ],
#  iterated to convergence, re-anchoring f_0 := 0 each sweep. Stop when
#  max_k |df_k| < tol or after max_iter sweeps.
#
#  A target (beta, kappa) has reduced action u_t(n) = beta*A(n) + kappa*(-B(n))
#  and MBAR weights log W(n) = -u_t(n) - g(n); <O> = weighted_mean_exp(O, logW).
# ---------------------------------------------------------------------------
def _reduced_action(beta: float, kappa: float,
                    A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """u(n) = beta*A(n) + kappa*(-B(n)) = beta*A - kappa*B  (the action S)."""
    return beta * A - kappa * B


def solve_wham(ensembles: Sequence[Ensemble],
               tol: float = 1e-10, max_iter: int = 100000) -> WhamSolution:
    """Solve the MBAR/WHAM free energies f_k for K reference ensembles.

    Returns a :class:`WhamSolution` carrying the converged ``f``, the per-sample
    mixture denominator ``g``, the flattened ``A/B/O`` time series, and
    convergence diagnostics. The expensive solve is done once; any number of
    targets can then be evaluated cheaply via :func:`multi_histogram_eval`.
    """
    K = len(ensembles)
    if K == 0:
        raise ValueError("need at least one ensemble")

    # Flatten all ensembles into global sample arrays.
    A = np.concatenate([e.A for e in ensembles]) if K else np.array([])
    B = np.concatenate([e.B for e in ensembles])
    O = np.concatenate([e.O for e in ensembles])
    Nk = np.array([float(e.n_samples) for e in ensembles], dtype=float)
    M = A.size
    if M == 0:
        raise ValueError("ensembles contain no samples")

    # U[k][n] = reduced action of sample n at ensemble k's couplings.
    U = np.empty((K, M), dtype=float)
    for k, e in enumerate(ensembles):
        U[k] = _reduced_action(e.beta0, e.kappa0, A, B)

    logNk = np.where(Nk > 0.0, np.log(np.where(Nk > 0.0, Nk, 1.0)), -np.inf)

    f = np.zeros(K, dtype=float)
    residual = np.inf
    it = 0
    g = None
    for it in range(1, max_iter + 1):
        # g(n) = log_sum_exp_l [ logN_l + f_l - U[l][n] ]   (column-wise LSE).
        terms = (logNk + f)[:, None] - U            # shape (K, M)
        g = _logsumexp_axis0(terms)                 # length M
        # f_k(new) = -log_sum_exp_n [ -U[k][n] - g(n) ].
        fnew = -_logsumexp_axis1(-U - g[None, :])   # length K
        fnew -= fnew[0]                             # re-anchor gauge: f_0 := 0
        residual = float(np.max(np.abs(fnew - f)))
        f = fnew
        if residual < tol:
            break

    # Final g consistent with the converged f.
    terms = (logNk + f)[:, None] - U
    g = _logsumexp_axis0(terms)

    return WhamSolution(f=f, g=g, A=A, B=B, O=O,
                        residual=residual, iters=it, converged=residual < tol)


def _logsumexp_axis0(x: np.ndarray) -> np.ndarray:
    """Stable log-sum-exp over axis 0 (rows), returning a length-ncols vector."""
    m = np.max(x, axis=0)
    finite = np.isfinite(m)
    out = np.array(m, dtype=float)  # non-finite columns pass max through
    if np.any(finite):
        shifted = x[:, finite] - m[finite]
        out[finite] = m[finite] + np.log(np.sum(np.exp(shifted), axis=0))
    return out


def _logsumexp_axis1(x: np.ndarray) -> np.ndarray:
    """Stable log-sum-exp over axis 1 (cols), returning a length-nrows vector."""
    m = np.max(x, axis=1)
    finite = np.isfinite(m)
    out = np.array(m, dtype=float)
    if np.any(finite):
        shifted = x[finite, :] - m[finite][:, None]
        out[finite] = m[finite] + np.log(np.sum(np.exp(shifted), axis=1))
    return out


def multi_histogram_eval(sol: WhamSolution, beta: float, kappa: float,
                         observable: ObservableSpec = None) -> float:
    """Reweight a target ``(beta, kappa)`` using a precomputed WHAM solution.

        a(n) = -u_t(n) - g(n),   u_t(n) = beta*A(n) - kappa*B(n)
        <O>  = ( sum_n O(n) e^{a(n)} ) / ( sum_n e^{a(n)} ),

    evaluated with :func:`_weighted_mean_exp`. ``observable`` selects which
    quantity to average (default: the flattened ensemble O).
    """
    O = _resolve_observable(observable, sol.A, sol.B, sol.O)
    a = -_reduced_action(beta, kappa, sol.A, sol.B) - sol.g
    return _weighted_mean_exp(O, a)


def multi_histogram(ensembles: Sequence[Ensemble], beta: float, kappa: float,
                    observable: ObservableSpec = None,
                    tol: float = 1e-10, max_iter: int = 100000,
                    return_solution: bool = False):
    """One-shot WHAM/MBAR ``<O>`` at target ``(beta, kappa)``.

    Combines the K reference ensembles (so a target between simulated points uses
    all overlapping data), solving the free energies once and reweighting the
    target. With ``return_solution=True`` also returns the :class:`WhamSolution`
    (for convergence assertions / reuse).
    """
    sol = solve_wham(ensembles, tol=tol, max_iter=max_iter)
    val = multi_histogram_eval(sol, beta, kappa, observable=observable)
    if return_solution:
        return val, sol
    return val


def multi_histogram_susceptibility(sol: WhamSolution, V: float,
                                   beta: float, kappa: float,
                                   observable: ObservableSpec = None) -> float:
    """``V * (<O^2> - <O>^2)`` at ``(beta, kappa)`` from a precomputed solution."""
    O = _resolve_observable(observable, sol.A, sol.B, sol.O)
    a = -_reduced_action(beta, kappa, sol.A, sol.B) - sol.g
    mean = _weighted_mean_exp(O, a)
    mean2 = _weighted_mean_exp(O * O, a)
    return float(V * (mean2 - mean * mean))


# ---------------------------------------------------------------------------
#  Susceptibility surface over a (beta, kappa) mesh, via the multi-histogram.
# ---------------------------------------------------------------------------
def chi_surface(ensembles: Sequence[Ensemble], V: float,
                beta_mesh: ArrayLike, kappa_mesh: ArrayLike,
                observable: ObservableSpec = None,
                tol: float = 1e-10, max_iter: int = 100000):
    """Reweighted susceptibility surface over a ``(beta, kappa)`` mesh.

    Solves the WHAM free energies ONCE from ``ensembles`` and evaluates
    ``V * (<O^2> - <O>^2)`` at every mesh point. ``beta_mesh`` and ``kappa_mesh``
    are 1-D coordinate vectors.

    Returns ``(BETA, KAPPA, CHI)`` where ``BETA, KAPPA = np.meshgrid(beta_mesh,
    kappa_mesh)`` and ``CHI`` has shape ``(len(kappa_mesh), len(beta_mesh))`` --
    directly usable by ``matplotlib`` ``contour`` / ``pcolormesh`` to trace ridge
    lines off the simulated grid.
    """
    beta_mesh = np.asarray(beta_mesh, dtype=float).ravel()
    kappa_mesh = np.asarray(kappa_mesh, dtype=float).ravel()
    sol = solve_wham(ensembles, tol=tol, max_iter=max_iter)

    BETA, KAPPA = np.meshgrid(beta_mesh, kappa_mesh)
    CHI = np.empty_like(BETA, dtype=float)
    for j in range(KAPPA.shape[0]):
        for i in range(BETA.shape[1]):
            CHI[j, i] = multi_histogram_susceptibility(
                sol, V, BETA[j, i], KAPPA[j, i], observable=observable)
    return BETA, KAPPA, CHI


# ---------------------------------------------------------------------------
#  Time-series loader.
# ---------------------------------------------------------------------------
def load_ts(path: str) -> dict:
    """Load a campaign time-series ``.dat`` file into ``{column_name: array}``.

    File format (see ``src/u1_scan.cpp``):
        # ... arbitrary comment lines ...
        # columns: traj  A  B  avg_plaquette  higgs_length  link_energy  monopole_density
        0 322.1 2359.2 0.79 1.60 1.15 0.156
        ...

    Comment lines start with ``#``. The ``# columns:`` line names the columns; if
    absent, generic ``c0, c1, ...`` names are used. Robust to extra whitespace and
    to a column header with more/fewer names than data columns (extra names are
    dropped, missing names are filled with ``cN``).
    """
    names: Optional[list[str]] = None
    rows: list[list[float]] = []
    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                body = s.lstrip("#").strip()
                low = body.lower()
                if low.startswith("columns:"):
                    names = body[len("columns:"):].split()
                continue
            parts = s.split()
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                # Not a numeric data row (e.g. a stray text line); skip it.
                continue

    if not rows:
        ncol = len(names) if names else 0
        return {(names[i] if names else f"c{i}"): np.array([], dtype=float)
                for i in range(ncol)}

    ncol = max(len(r) for r in rows)
    # Pad ragged rows defensively with NaN so the array is rectangular.
    data = np.full((len(rows), ncol), np.nan, dtype=float)
    for r_i, r in enumerate(rows):
        data[r_i, :len(r)] = r

    if names is None:
        names = [f"c{i}" for i in range(ncol)]
    else:
        names = list(names[:ncol]) + [f"c{i}" for i in range(len(names), ncol)]

    return {names[i]: data[:, i] for i in range(ncol)}


# ===========================================================================
#  SELF-TEST
# ===========================================================================
def _selftest() -> int:
    """Run analytic self-tests; print PASS/FAIL; return 0 on all-pass else 1."""
    rng = np.random.default_rng(20240530)
    fails = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal fails
        status = "PASS" if ok else "FAIL"
        if not ok:
            fails += 1
        msg = f"[{status}] {name}"
        if detail:
            msg += f"   ({detail})"
        print(msg)

    # ---- log_sum_exp sanity -------------------------------------------------
    check("log_sum_exp empty == -inf", log_sum_exp([]) == -np.inf)
    a = np.array([1000.0, 1000.5, 999.5])  # would overflow naively
    lse = log_sum_exp(a)
    lse_ref = 1000.5 + np.log(np.exp(-0.5) + np.exp(0.0) + np.exp(-1.0))
    check("log_sum_exp stable & correct", abs(lse - lse_ref) < 1e-9,
          f"|d|={abs(lse - lse_ref):.2e}")

    # ---- (a) Zero-shift identity -------------------------------------------
    n = 4000
    A = rng.normal(300.0, 30.0, n)
    B = rng.normal(2400.0, 200.0, n)
    O = rng.normal(0.77, 0.02, n)
    beta0, kappa0, V = 1.0, 0.3, 256.0
    m_rw = reweighted_mean(O, A, B, beta0, kappa0, beta0, kappa0)
    m_plain = float(np.mean(O))
    check("(a) zero-shift mean == plain mean",
          abs(m_rw - m_plain) < 1e-10, f"|d|={abs(m_rw - m_plain):.2e}")
    chi_rw = reweighted_susceptibility(O, A, B, V, beta0, kappa0, beta0, kappa0)
    chi_plain = float(V * np.var(O))  # population variance (1/N) matches estimator
    check("(a) zero-shift susceptibility == V*var",
          abs(chi_rw - chi_plain) < 1e-9, f"|d|={abs(chi_rw - chi_plain):.2e}")

    # ---- (b) Two-state analytic model --------------------------------------
    # A system with two macrostates s in {0, 1}. State s has FIXED conjugate
    # energies A_s = alpha_s, B_s = gamma_s and observable value O_s. Sampled at
    # the reference (beta0, kappa0), the ensemble visits state s exactly N_s
    # times. The FS single-histogram estimator at a target (beta, kappa) is then
    #   <O> = sum_s O_s N_s e^{-(dbeta A_s - dkappa B_s)}
    #         / sum_s N_s e^{-(dbeta A_s - dkappa B_s)} ,   d* = (* - *0),
    # which is the EXACT two-state Boltzmann mean whose reference weights are the
    # N_s. We build the analytic curve from those same N_s, so agreement is exact
    # (machine precision) and does NOT depend on Monte-Carlo sampling noise -- the
    # estimator is a deterministic function of the (A, B, O) arrays.
    alpha = np.array([10.0, 16.0])     # A_s
    gamma = np.array([4.0, 9.0])       # B_s
    Ovals = np.array([0.2, 0.9])       # O_s
    counts = np.array([7, 3])          # N_s, exact integer reference multiplicities
    br, kr = 0.5, 0.7                  # reference couplings

    def analytic_mean(b, k):
        # Effective reference weight of state s is exactly N_s; reweighting to
        # (b, k) multiplies it by exp(-((b-br)A_s - (k-kr)B_s)).
        logw = (np.log(counts.astype(float))
                - ((b - br) * alpha - (k - kr) * gamma))
        w = np.exp(logw - logw.max())
        return float(np.sum(Ovals * w) / np.sum(w))

    s_idx = np.concatenate([np.full(c, s) for s, c in enumerate(counts)])
    A2 = alpha[s_idx]
    B2 = gamma[s_idx]
    O2 = Ovals[s_idx]

    worst = 0.0
    for b in (0.3, 0.5, 0.8, 1.2):
        for k in (0.4, 0.7, 1.0):
            est = reweighted_mean(O2, A2, B2, br, kr, b, k)
            ref = analytic_mean(b, k)
            worst = max(worst, abs(est - ref))
    check("(b) two-state reweighted mean matches analytic curve",
          worst < 1e-6, f"worst|d|={worst:.2e}")

    # ---- (c) single- vs multi-histogram with ONE ensemble ------------------
    ens1 = Ensemble(br, kr, A2, B2, O2)
    worst_cm = 0.0
    for b in (0.4, 0.6, 0.9):
        for k in (0.5, 0.8):
            single = reweighted_mean(O2, A2, B2, br, kr, b, k)
            multi, sol1 = multi_histogram([ens1], b, k, return_solution=True)
            worst_cm = max(worst_cm, abs(single - multi))
    check("(c) single == multi (one ensemble)",
          worst_cm < 1e-9, f"worst|d|={worst_cm:.2e}")

    # ---- (d) cross-check vs brute-force direct (no LSE) for small dlambda ---
    # Small shift so exp() of the raw exponents cannot overflow.
    db, dk = 0.001, 0.001
    bt, kt = beta0 + db, kappa0 + dk
    # Direct, non-log-sum-exp computation (centered to avoid trivial overflow,
    # but using a plain ratio of exp-sums rather than the LSE path).
    a_direct = -(bt - beta0) * A + (kt - kappa0) * B
    a_direct = a_direct - a_direct.mean()  # harmless common shift; small range
    w_direct = np.exp(a_direct)
    mean_direct = float(np.sum(O * w_direct) / np.sum(w_direct))
    mean_lse = reweighted_mean(O, A, B, beta0, kappa0, bt, kt)
    check("(d) LSE path == brute-force direct (small dlambda)",
          abs(mean_lse - mean_direct) < 1e-10,
          f"|d|={abs(mean_lse - mean_direct):.2e}")

    # ---- (e) jackknife error positive and shrinks ~1/sqrt(n) ---------------
    # IID synthetic data, zero shift so the reweighted mean is the plain mean and
    # the jackknife error should match the naive standard error sigma/sqrt(n).
    def jk_err(nn):
        Oi = rng.normal(0.0, 1.0, nn)
        Ai = rng.normal(300.0, 30.0, nn)
        Bi = rng.normal(2400.0, 200.0, nn)
        _, err = reweighted_mean_jackknife(Oi, Ai, Bi, beta0, kappa0,
                                           beta0, kappa0, nblocks=nn)
        return err

    e_small = jk_err(1000)
    e_large = jk_err(16000)  # 16x the samples -> ~4x smaller error
    ratio = e_small / e_large if e_large > 0 else np.inf
    check("(e) jackknife error positive", e_small > 0 and e_large > 0,
          f"err(1k)={e_small:.4f}")
    check("(e) jackknife error shrinks ~1/sqrt(n)",
          3.0 < ratio < 5.0, f"ratio(16x more data)={ratio:.2f} (expect ~4)")

    print()
    if fails == 0:
        print("ALL SELF-TESTS PASSED")
    else:
        print(f"{fails} SELF-TEST(S) FAILED")
    return 1 if fails else 0


def _main(argv: list[str]) -> int:
    if "--selftest" in argv:
        return _selftest()
    print(__doc__)
    print("Run with --selftest to execute the analytic checks.")
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
