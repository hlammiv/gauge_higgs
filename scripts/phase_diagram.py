#!/usr/bin/env python3
"""Group-agnostic (beta,kappa) phase-diagram analyzer + plotter for a u1_scan campaign.

Turns a u1_scan_campaign output directory into labeled (beta,kappa) phase diagrams
with susceptibility-ridge phase-boundary lines, one figure per (L,q), plus a
multi-q junction-trajectory summary. The CSV/time-series format is the only thing
this depends on, so the same tool serves the eventual SU(N) discrete-group scans.

WHAT IT COMPUTES (numpy)
  Points are enumerated from <campaign>/manifest.txt (fallback: walk the tree for
  ts_*.dat) and grouped by (L,q). For each point we read its per-trajectory time
  series ts_*.dat and, for every intensive observable O present among
      {avg_plaquette, link_energy, monopole_density},
  compute the mean and the SUSCEPTIBILITY
      chi_O = V * Var(O),   V = L**ndim   (ndim from --ndim, default 4).
  NOTE: the phase BOUNDARY is the LOCATION of the chi PEAK along a scan line. The
  volume prefactor V only RESCALES chi by a (beta,kappa)-independent constant; it
  does NOT move the argmax. Boundary-finding is therefore prefactor-robust -- the
  detected ridge is identical whether or not V (or any constant) multiplies Var.
  The error on chi is a simple BLOCKING error: split the series into blocks, form
  chi on each block, and take V * std(block variances)  (autocorrelation-aware).

BOUNDARY DETECTION (susceptibility ridges) -- per (L,q):
  We lay the points on the (beta,kappa) grid.
   - Coulomb<->confinement line: use chi[monopole_density] if that column exists,
     else chi[avg_plaquette]. For each fixed-kappa ROW, the beta of maximum chi is a
     boundary point (kappa, beta*).
   - Higgs<->confinement line:   chi[link_energy]. For each fixed-beta COLUMN, the
     kappa of maximum chi is a boundary point (beta, kappa*).
  On a coarse grid "peak" is just the argmax over the available points (noted as
  coarse). Each ridge is collected as an ordered point set; where the two ridges
  meet is the triple-point region.

OUTPUTS
  <figdir>/phase_q<q>_L<L>.png   per (L,q): chi heatmap + ridge points/lines + region
                                 labels + triple-point marker.
  <figdir>/trajectory.png        multi-q: all ridge lines overlaid, and the junction
                                 (beta_t,kappa_t) vs q trajectory.
  <campaign>/boundaries.csv      ALWAYS written, columns:
        L,q,line_type,fixed_coord,fixed_value,peak_coord,peak_value,chi_max,chi_err
  (<figdir> defaults to <campaign>/figs.)

ROBUSTNESS
  Header-driven columns (an added ts column just works); ridges whose observable is
  absent are skipped; coarse grids and single-row / single-column groups do not
  crash. If matplotlib import fails we still write boundaries.csv + a text summary
  and exit 0 with a notice.

CLI
  python3 scripts/phase_diagram.py <campaign_dir> [--ndim 4] [--figdir figs]
  python3 scripts/phase_diagram.py --selftest
"""
import os
import sys
import csv
import math
import argparse
import tempfile
import shutil
from collections import defaultdict, OrderedDict

import numpy as np

# Observables we know how to treat as intensive order-parameter proxies. The boundary
# rules below reference these names; any are skipped silently if the column is absent.
OBSERVABLES = ("avg_plaquette", "link_energy", "monopole_density")
# Coulomb<->confinement ridge prefers the monopole density (the U(1) order parameter)
# and falls back to the plaquette; Higgs<->confinement ridge uses the link energy.
COULOMB_PREF = ("monopole_density", "avg_plaquette")
HIGGS_OBS = "link_energy"


# --------------------------------------------------------------------------- I/O
def _parse_columns_line(path):
    """Return the header-declared column NAME->index map from a ts_*.dat file.

    Reads the comment block and parses the single line of the exact form
        '# columns: traj  A  B  avg_plaquette ...'
    Header-driven by design: appending a new column (e.g. monopole_density) to the
    time series is picked up automatically. Returns {} if no such line is found.
    """
    try:
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    break  # comment block ended; no columns line
                low = line.lower()
                if "columns:" in low:
                    after = line.split("columns:", 1)[1]
                    names = after.split()
                    return {name: i for i, name in enumerate(names)}
    except OSError:
        return {}
    return {}


def load_timeseries(path):
    """Load a ts_*.dat into a dict NAME->1D float array using its '# columns:' header.

    Falls back to the documented default column order if no header line is present.
    Rows that do not parse as all-float are dropped. Returns {} if nothing loaded.
    """
    colmap = _parse_columns_line(path)
    if not colmap:
        colmap = {n: i for i, n in enumerate(
            ["traj", "A", "B", "avg_plaquette", "higgs_length", "link_energy"])}
    try:
        # Empty / still-being-written ts files (live campaigns) trigger a loadtxt
        # "no data" UserWarning; that case is expected and handled by the size check.
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = np.loadtxt(path, comments="#", ndmin=2)
    except (OSError, ValueError):
        return {}
    if data.size == 0 or data.shape[0] == 0:
        return {}
    out = {}
    ncol = data.shape[1]
    for name, idx in colmap.items():
        if 0 <= idx < ncol:
            out[name] = data[:, idx]
    return out


def _strip_comments(line):
    return line.split("#", 1)[0].strip()


def enumerate_points(campaign):
    """Enumerate campaign points as a list of dicts.

    Each point: {L, q, beta, kappa, tsfile}. Prefer manifest.txt (rows
    'L q beta kappa base_seed outdir'); fall back to walking the tree for ts_*.dat
    and inferring (L,q,beta,kappa) from the directory layout / filename.
    """
    points = []
    manifest = os.path.join(campaign, "manifest.txt")
    if os.path.isfile(manifest):
        with open(manifest) as f:
            for raw in f:
                row = _strip_comments(raw)
                if not row:
                    continue
                p = row.split()
                if len(p) < 6:
                    continue
                try:
                    L = int(float(p[0])); q = int(float(p[1]))
                    beta = float(p[2]); kappa = float(p[3])
                except ValueError:
                    continue
                outdir = p[5]
                ts = _find_ts(campaign, outdir, beta, kappa)
                if ts:
                    points.append(dict(L=L, q=q, beta=beta, kappa=kappa, tsfile=ts))
    if points:
        return points
    # ---- fallback: walk the tree ------------------------------------------------
    for root, _dirs, files in os.walk(campaign):
        for fn in files:
            if not (fn.startswith("ts_") and fn.endswith(".dat")):
                continue
            ts = os.path.join(root, fn)
            beta, kappa = _parse_bk_from_name(fn)
            L, q = _parse_lq_from_path(root)
            if beta is None or kappa is None:
                continue
            points.append(dict(L=L if L is not None else 0,
                               q=q if q is not None else 0,
                               beta=beta, kappa=kappa, tsfile=ts))
    return points


def _find_ts(campaign, outdir, beta, kappa):
    """Locate the ts_*.dat for a manifest point. outdir may be relative to repo or
    to the campaign dir; the ts file is named ts_b<beta>_k<kappa>.dat inside it."""
    fname = "ts_b%.6f_k%.6f.dat" % (beta, kappa)
    candidates = [
        os.path.join(outdir, fname),                       # repo-root relative
        os.path.join(campaign, outdir, fname),             # campaign relative
        os.path.join(os.path.dirname(campaign), outdir, fname),
    ]
    # outdir often starts with the campaign basename; also try stripping that prefix.
    base = os.path.basename(os.path.normpath(campaign))
    parts = outdir.replace("\\", "/").split("/")
    if parts and parts[0] == base:
        candidates.append(os.path.join(campaign, *parts[1:], fname))
    for c in candidates:
        if os.path.isfile(c):
            return c
    # last resort: any ts file in the resolved dir
    for d in (outdir, os.path.join(campaign, outdir)):
        if os.path.isdir(d):
            for f in os.listdir(d):
                if f.startswith("ts_") and f.endswith(".dat"):
                    return os.path.join(d, f)
    return None


def _parse_bk_from_name(fn):
    """ts_b0.800000_k0.450000.dat -> (0.8, 0.45)."""
    try:
        stem = fn[len("ts_"):-len(".dat")]
        bpart, kpart = stem.split("_k")
        beta = float(bpart[len("b"):])
        kappa = float(kpart)
        return beta, kappa
    except (ValueError, IndexError):
        return None, None


def _parse_lq_from_path(path):
    """Infer (L,q) from a path containing .../L<L>/q<q>/..."""
    L = q = None
    for seg in path.replace("\\", "/").split("/"):
        if seg.startswith("L") and seg[1:].isdigit():
            L = int(seg[1:])
        elif seg.startswith("q") and seg[1:].isdigit():
            q = int(seg[1:])
    return L, q


# ------------------------------------------------------- susceptibility & blocking
def susceptibility(series, volume):
    """chi = V * Var(O). Population variance; the V prefactor only rescales (does not
    move the argmax), so ridge LOCATIONS are prefactor-robust."""
    if series is None or len(series) < 2:
        return float("nan")
    return float(volume) * float(np.var(series))


def chi_blocking_error(series, volume, nblocks=16):
    """Simple blocking error on chi: split the series into blocks, compute chi on each
    block, return V * std(block variances) (autocorrelation-aware). 0 if too short."""
    n = 0 if series is None else len(series)
    if n < 4:
        return 0.0
    nb = min(nblocks, max(2, n // 2))
    bs = n // nb
    if bs < 2:
        return 0.0
    block_var = [float(np.var(series[b * bs:(b + 1) * bs])) for b in range(nb)]
    return float(volume) * float(np.std(block_var, ddof=1)) if len(block_var) > 1 else 0.0


# --------------------------------------------------------------- group / grid model
class Group:
    """All points for one (L,q): per-point means + chi for every observable present,
    laid out on the (beta,kappa) grid."""

    def __init__(self, L, q, ndim):
        self.L = L
        self.q = q
        self.volume = float(L) ** ndim if L else 1.0
        # keyed by (beta,kappa)
        self.mean = defaultdict(dict)   # (b,k) -> {obs: mean}
        self.chi = defaultdict(dict)    # (b,k) -> {obs: chi}
        self.chi_err = defaultdict(dict)
        self.betas = set()
        self.kappas = set()
        self.present = set()            # observable names actually seen

    def add_point(self, beta, kappa, ts):
        bk = (round(beta, 9), round(kappa, 9))
        self.betas.add(bk[0]); self.kappas.add(bk[1])
        for obs in OBSERVABLES:
            if obs in ts and len(ts[obs]) >= 2:
                s = ts[obs]
                self.mean[bk][obs] = float(np.mean(s))
                self.chi[bk][obs] = susceptibility(s, self.volume)
                self.chi_err[bk][obs] = chi_blocking_error(s, self.volume)
                self.present.add(obs)

    def beta_grid(self):
        return sorted(self.betas)

    def kappa_grid(self):
        return sorted(self.kappas)

    def chi_at(self, beta, kappa, obs):
        return self.chi.get((round(beta, 9), round(kappa, 9)), {}).get(obs, float("nan"))

    def chi_err_at(self, beta, kappa, obs):
        return self.chi_err.get((round(beta, 9), round(kappa, 9)), {}).get(obs, 0.0)


def _coulomb_obs(group):
    """Pick the Coulomb<->confinement observable that is present (monopole preferred)."""
    for o in COULOMB_PREF:
        if o in group.present:
            return o
    return None


def detect_boundaries(group):
    """Return ridge point sets for this group.

    ridges = {
      'coulomb': {'obs': name, 'points': [ (kappa_fixed, beta_peak, chi, chi_err), ...]},
      'higgs':   {'obs': name, 'points': [ (beta_fixed, kappa_peak, chi, chi_err), ...]},
    }
    Coulomb ridge: for each fixed kappa ROW, beta of max chi.
    Higgs ridge:   for each fixed beta COLUMN, kappa of max chi.
    A ridge is omitted entirely if its observable is absent.
    """
    ridges = {}
    betas = group.beta_grid()
    kappas = group.kappa_grid()

    cobs = _coulomb_obs(group)
    if cobs is not None and len(betas) >= 1:
        pts = []
        for k in kappas:
            best = None  # (chi, beta, err)
            for b in betas:
                c = group.chi_at(b, k, cobs)
                if not math.isnan(c) and (best is None or c > best[0]):
                    best = (c, b, group.chi_err_at(b, k, cobs))
            if best is not None:
                pts.append((k, best[1], best[0], best[2]))
        if pts:
            ridges["coulomb"] = dict(obs=cobs, points=pts)

    if HIGGS_OBS in group.present and len(kappas) >= 1:
        pts = []
        for b in betas:
            best = None  # (chi, kappa, err)
            for k in kappas:
                c = group.chi_at(b, k, HIGGS_OBS)
                if not math.isnan(c) and (best is None or c > best[0]):
                    best = (c, k, group.chi_err_at(b, k, HIGGS_OBS))
            if best is not None:
                pts.append((b, best[1], best[0], best[2]))
        if pts:
            ridges["higgs"] = dict(obs=HIGGS_OBS, points=pts)

    return ridges


def junction(ridges):
    """Estimate the triple-point / ridge-intersection (beta_t, kappa_t).

    The coulomb ridge maps kappa->beta*(kappa); the higgs ridge maps beta->kappa*(beta).
    Their crossing is the self-consistent point closest to both. We pick the (beta,kappa)
    from the discrete ridge points that minimizes the combined mismatch. Returns
    (beta_t, kappa_t) or None if either ridge is missing.
    """
    if "coulomb" not in ridges or "higgs" not in ridges:
        return None
    c = ridges["coulomb"]["points"]  # (kappa, beta*, ...)
    h = ridges["higgs"]["points"]    # (beta, kappa*, ...)
    if not c or not h:
        return None
    # interpolators on the (sparse) ridge samples
    ck = np.array([p[0] for p in c]); cb = np.array([p[1] for p in c])
    hb = np.array([p[0] for p in h]); hk = np.array([p[1] for p in h])
    cs = np.argsort(ck); ck, cb = ck[cs], cb[cs]
    hs = np.argsort(hb); hb, hk = hb[hs], hk[hs]
    # candidate kappas: union of both ridges' kappa coordinates
    cand_k = np.unique(np.concatenate([ck, hk]))
    best = None
    for kq in cand_k:
        beta_from_coul = float(np.interp(kq, ck, cb)) if len(ck) > 1 else float(cb[0])
        # kappa predicted by higgs ridge at that beta:
        kappa_from_higgs = float(np.interp(beta_from_coul, hb, hk)) if len(hb) > 1 else float(hk[0])
        mismatch = abs(kappa_from_higgs - kq)
        if best is None or mismatch < best[0]:
            best = (mismatch, beta_from_coul, 0.5 * (kq + kappa_from_higgs))
    if best is None:
        return None
    return (best[1], best[2])


# ------------------------------------------------------------------------ outputs
def write_boundaries_csv(path, groups):
    """Always write the boundary point set. Columns:
    L,q,line_type,fixed_coord,fixed_value,peak_coord,peak_value,chi_max,chi_err."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["L", "q", "line_type", "fixed_coord", "fixed_value",
                    "peak_coord", "peak_value", "chi_max", "chi_err"])
        for g in groups:
            ridges = detect_boundaries(g)
            if "coulomb" in ridges:
                for (kfix, bpk, chi, err) in ridges["coulomb"]["points"]:
                    w.writerow([g.L, g.q, "coulomb_confinement", "kappa",
                                "%.6f" % kfix, "beta", "%.6f" % bpk,
                                "%.8g" % chi, "%.8g" % err])
            if "higgs" in ridges:
                for (bfix, kpk, chi, err) in ridges["higgs"]["points"]:
                    w.writerow([g.L, g.q, "higgs_confinement", "beta",
                                "%.6f" % bfix, "kappa", "%.6f" % kpk,
                                "%.8g" % chi, "%.8g" % err])


def write_text_summary(path, groups):
    """Text fallback when matplotlib is unavailable."""
    with open(path, "w") as f:
        f.write("# phase_diagram text summary (matplotlib unavailable)\n")
        for g in groups:
            ridges = detect_boundaries(g)
            f.write("L=%d q=%d  observables=%s  grid=%dx%d (beta x kappa)\n" % (
                g.L, g.q, ",".join(sorted(g.present)),
                len(g.beta_grid()), len(g.kappa_grid())))
            for name in ("coulomb", "higgs"):
                if name in ridges:
                    r = ridges[name]
                    f.write("  %-20s ridge (obs=%s): %d points\n" %
                            (name, r["obs"], len(r["points"])))
            jt = junction(ridges)
            if jt:
                f.write("  triple-point region ~ (beta=%.4f, kappa=%.4f)\n" % jt)


# ------------------------------------------------------------------------ plotting
def _heatmap_arrays(group, obs):
    """Build (beta-grid, kappa-grid, chi 2D array [nk,nb]) for pcolormesh. NaN where
    a point is missing (partial campaigns)."""
    betas = group.beta_grid()
    kappas = group.kappa_grid()
    Z = np.full((len(kappas), len(betas)), np.nan)
    for i, k in enumerate(kappas):
        for j, b in enumerate(betas):
            Z[i, j] = group.chi_at(b, k, obs)
    return np.array(betas), np.array(kappas), Z


def plot_group(group, figdir):
    """Per-(L,q) phase diagram. Returns the saved path or None on skip."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ridges = detect_boundaries(group)
    betas = group.beta_grid()
    kappas = group.kappa_grid()
    if not betas or not kappas:
        return None

    fig, ax = plt.subplots(figsize=(7.5, 6.0))

    # Background heatmap: prefer the link-energy chi (the Higgs-line driver), else the
    # coulomb-line observable. Use the chi that is most informative across the grid.
    bg_obs = HIGGS_OBS if HIGGS_OBS in group.present else _coulomb_obs(group)
    mappable = None
    if bg_obs is not None:
        B, K, Z = _heatmap_arrays(group, bg_obs)
        if np.isfinite(Z).any():
            if len(B) >= 2 and len(K) >= 2:
                # pcolormesh wants cell edges; build them from midpoints.
                be = _edges(B); ke = _edges(K)
                mappable = ax.pcolormesh(be, ke, np.ma.masked_invalid(Z),
                                         shading="auto", cmap="viridis")
            else:
                # degenerate (single row/col) -> scatter colored by chi
                bb, kk, cc = [], [], []
                for i, k in enumerate(K):
                    for j, b in enumerate(B):
                        if np.isfinite(Z[i, j]):
                            bb.append(b); kk.append(k); cc.append(Z[i, j])
                if cc:
                    mappable = ax.scatter(bb, kk, c=cc, s=200, marker="s",
                                          cmap="viridis")
            if mappable is not None:
                cb = fig.colorbar(mappable, ax=ax)
                cb.set_label(r"$\chi_{\mathrm{%s}}=V\,\mathrm{Var}$" %
                             bg_obs.replace("_", r"\_"))

    # Overlay ridges.
    if "coulomb" in ridges:
        pts = sorted(ridges["coulomb"]["points"], key=lambda p: p[0])  # by kappa
        rb = [p[1] for p in pts]; rk = [p[0] for p in pts]
        ax.plot(rb, rk, "o-", color="white", mec="black", mew=1.2, ms=8,
                lw=2, label="Coulomb-confine (%s)" % ridges["coulomb"]["obs"])
    if "higgs" in ridges:
        pts = sorted(ridges["higgs"]["points"], key=lambda p: p[0])    # by beta
        hb = [p[0] for p in pts]; hk = [p[1] for p in pts]
        ax.plot(hb, hk, "s--", color="red", mec="black", mew=1.0, ms=8,
                lw=2, label="Higgs-confine (%s)" % ridges["higgs"]["obs"])

    # Triple-point marker.
    jt = junction(ridges)
    if jt is not None:
        ax.plot([jt[0]], [jt[1]], marker="*", ms=22, color="gold", mec="black",
                mew=1.4, lw=0, label="triple-point region")

    # Phase-region text labels (standard U(1)+Higgs layout).
    bmin, bmax = min(betas), max(betas)
    kmin, kmax = min(kappas), max(kappas)

    def _frac(lo, hi, f):
        return lo + (hi - lo) * f if hi > lo else lo

    ax.text(_frac(bmin, bmax, 0.10), _frac(kmin, kmax, 0.12), "Confined",
            ha="left", va="bottom", fontsize=12, fontweight="bold",
            bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    ax.text(_frac(bmin, bmax, 0.92), _frac(kmin, kmax, 0.10), "Coulomb",
            ha="right", va="bottom", fontsize=12, fontweight="bold",
            bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    ax.text(_frac(bmin, bmax, 0.50), _frac(kmin, kmax, 0.92), "Higgs",
            ha="center", va="top", fontsize=12, fontweight="bold",
            bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    coarse = (len(betas) <= 8 or len(kappas) <= 8)
    ax.set_xlabel(r"$\beta$ (gauge coupling)")
    ax.set_ylabel(r"$\kappa$ (hopping coupling)")
    ttl = "U(1) charge q=%d, L=%d" % (group.q, group.L)
    if coarse:
        ttl += "  (coarse grid: peak = argmax)"
    ax.set_title(ttl)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.85)
    fig.tight_layout()
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, "phase_q%d_L%d.png" % (group.q, group.L))
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def _edges(centers):
    """Cell edges from sorted centers for pcolormesh."""
    c = np.asarray(centers, dtype=float)
    if len(c) == 1:
        d = 1.0
        return np.array([c[0] - d / 2, c[0] + d / 2])
    mids = 0.5 * (c[:-1] + c[1:])
    first = c[0] - (mids[0] - c[0])
    last = c[-1] + (c[-1] - mids[-1])
    return np.concatenate([[first], mids, [last]])


def plot_trajectory(groups, figdir):
    """Multi-q summary: overlay every group's ridges on one (beta,kappa) plot, and
    plot the junction (beta_t,kappa_t) vs q. Returns saved path or None."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 5.5))
    cmap = plt.get_cmap("tab10")
    junctions = []  # (q, L, beta_t, kappa_t)

    for idx, g in enumerate(sorted(groups, key=lambda gg: (gg.q, gg.L))):
        ridges = detect_boundaries(g)
        col = cmap(idx % 10)
        lab = "q=%d L=%d" % (g.q, g.L)
        if "coulomb" in ridges:
            pts = sorted(ridges["coulomb"]["points"], key=lambda p: p[0])
            ax0.plot([p[1] for p in pts], [p[0] for p in pts], "o-", color=col,
                     ms=5, lw=1.6, label="%s Coul" % lab)
        if "higgs" in ridges:
            pts = sorted(ridges["higgs"]["points"], key=lambda p: p[0])
            ax0.plot([p[0] for p in pts], [p[1] for p in pts], "s--", color=col,
                     ms=5, lw=1.6, label="%s Higgs" % lab)
        jt = junction(ridges)
        if jt is not None:
            ax0.plot([jt[0]], [jt[1]], marker="*", ms=16, color=col, mec="black")
            junctions.append((g.q, g.L, jt[0], jt[1]))

    ax0.set_xlabel(r"$\beta$"); ax0.set_ylabel(r"$\kappa$")
    ax0.set_title("phase boundaries, all (q,L)")
    ax0.legend(fontsize=7, ncol=2, framealpha=0.85)

    if junctions:
        js = sorted(junctions, key=lambda t: t[0])
        qs = [t[0] for t in js]
        bt = [t[2] for t in js]
        kt = [t[3] for t in js]
        ax1.plot(qs, bt, "o-", color="tab:blue", label=r"$\beta_t$")
        ax1.plot(qs, kt, "s-", color="tab:red", label=r"$\kappa_t$")
        ax1.set_xlabel("charge q"); ax1.set_ylabel("triple-point coord")
        ax1.set_title(r"junction $(\beta_t,\kappa_t)$ vs $q$")
        ax1.legend()
    else:
        ax1.text(0.5, 0.5, "no junctions\n(need both ridges)", ha="center",
                 va="center", transform=ax1.transAxes)
        ax1.set_title(r"junction vs $q$ (none)")

    fig.tight_layout()
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, "trajectory.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- driver
def build_groups(campaign, ndim):
    """Enumerate + load all points, grouped by (L,q)."""
    points = enumerate_points(campaign)
    groups = OrderedDict()  # (L,q) -> Group
    n_loaded = 0
    for p in points:
        ts = load_timeseries(p["tsfile"])
        if not ts:
            continue
        key = (p["L"], p["q"])
        if key not in groups:
            groups[key] = Group(p["L"], p["q"], ndim)
        groups[key].add_point(p["beta"], p["kappa"], ts)
        n_loaded += 1
    return list(groups.values()), len(points), n_loaded


def run(campaign, ndim, figdir):
    """Analyze a campaign dir: write boundaries.csv (+ figures if possible)."""
    if not os.path.isdir(campaign):
        print("ERROR: campaign dir not found: %s" % campaign, file=sys.stderr)
        return 2
    groups, n_points, n_loaded = build_groups(campaign, ndim)
    print("# enumerated %d points, loaded %d time series, %d (L,q) groups" %
          (n_points, n_loaded, len(groups)))
    for g in groups:
        print("#   L=%d q=%d : %d beta x %d kappa, observables={%s}" %
              (g.L, g.q, len(g.beta_grid()), len(g.kappa_grid()),
               ",".join(sorted(g.present))))

    bpath = os.path.join(campaign, "boundaries.csv")
    write_boundaries_csv(bpath, groups)
    print("# wrote %s" % bpath)

    if not os.path.isabs(figdir):
        figdir = os.path.join(campaign, figdir)

    try:
        import matplotlib  # noqa: F401
    except Exception as e:  # pragma: no cover - environment dependent
        tsum = os.path.join(campaign, "phase_summary.txt")
        write_text_summary(tsum, groups)
        print("# NOTICE: matplotlib unavailable (%s); wrote text summary %s" % (e, tsum))
        return 0

    figs = []
    for g in groups:
        try:
            p = plot_group(g, figdir)
            if p:
                figs.append(p); print("# wrote %s" % p)
        except Exception as e:
            print("# WARN: figure for L=%d q=%d failed: %s" % (g.L, g.q, e),
                  file=sys.stderr)
    if groups:
        try:
            t = plot_trajectory(groups, figdir)
            if t:
                figs.append(t); print("# wrote %s" % t)
        except Exception as e:
            print("# WARN: trajectory figure failed: %s" % e, file=sys.stderr)
    print("# done: %d figures, boundaries.csv" % len(figs))
    return 0


# ---------------------------------------------------------------------- self-test
def _synth_series(rng, mean, var, n):
    """Gaussian series with the requested mean/variance (so chi=V*Var is controlled)."""
    if var <= 0:
        return np.full(n, mean)
    s = rng.normal(0.0, 1.0, size=n)
    s = (s - s.mean()) / (s.std() + 1e-300)   # zero-mean unit-var exactly
    return mean + math.sqrt(var) * s


def build_synthetic_campaign(root, rng):
    """Construct a campaign whose ts files have KNOWN ridges:
      chi[avg_plaquette] peaks at beta_peak(kappa)  (a ROW argmax over beta)
      chi[link_energy]  peaks at kappa_peak(beta)   (a COLUMN argmax over kappa)
    We bake the target variances into per-point synthetic series. Returns
    (betas, kappas, beta_peak_of_kappa, kappa_peak_of_beta) for assertions.
    """
    L = 4
    q = 2
    betas = [0.80, 0.90, 1.00, 1.10, 1.20]
    kappas = [0.10, 0.20, 0.30, 0.40]
    npts = 400

    # Known ridges (must land ON grid nodes so argmax recovers them exactly):
    #  Coulomb line: beta_peak rises with kappa -> index into betas.
    def beta_peak_idx(ik):  # one beta index per kappa row
        return min(ik + 1, len(betas) - 1)          # 1,2,3,3
    #  Higgs line: kappa_peak rises with beta -> index into kappas.
    def kappa_peak_idx(ib):  # one kappa index per beta column
        return min(ib, len(kappas) - 1)             # 0,1,2,3,3

    beta_peak_of_kappa = {kappas[ik]: betas[beta_peak_idx(ik)] for ik in range(len(kappas))}
    kappa_peak_of_beta = {betas[ib]: kappas[kappa_peak_idx(ib)] for ib in range(len(betas))}

    base_var = 1.0e-4
    peak_boost = 6.0   # variance multiplier exactly on the ridge node

    pdir_root = os.path.join(root, "L%d" % L, "q%d" % q)
    manifest_rows = []
    for ib, b in enumerate(betas):
        for ik, k in enumerate(kappas):
            d = os.path.join(pdir_root, "b%.6f_k%.6f" % (b, k))
            os.makedirs(d, exist_ok=True)
            fn = os.path.join(d, "ts_b%.6f_k%.6f.dat" % (b, k))

            # avg_plaquette variance: maximal when ib == beta_peak_idx(ik)
            var_plaq = base_var * (peak_boost if ib == beta_peak_idx(ik) else 1.0)
            # link_energy variance: maximal when ik == kappa_peak_idx(ib)
            var_link = base_var * (peak_boost if ik == kappa_peak_idx(ib) else 1.0)

            plaq = _synth_series(rng, 0.4 + 0.05 * b, var_plaq, npts)
            link = _synth_series(rng, 0.01 + 0.02 * k, var_link, npts)
            higgs = _synth_series(rng, 0.8, base_var, npts)
            A = _synth_series(rng, 1.0e4, 1.0, npts)
            B = _synth_series(rng, 10.0 * k, 1.0, npts)

            with open(fn, "w") as f:
                f.write("# synthetic U(1)+Higgs ts. L=%d q=%d beta=%.6f kappa=%.6f\n" % (L, q, b, k))
                f.write("# columns: traj  A  B  avg_plaquette  higgs_length  link_energy\n")
                for t in range(npts):
                    f.write("%d %.10g %.10g %.10g %.10g %.10g\n" %
                            (t, A[t], B[t], plaq[t], higgs[t], link[t]))
            manifest_rows.append((L, q, b, k, os.path.relpath(d, os.path.dirname(root))
                                  if False else d))

    # manifest with outdir paths relative to the campaign's parent (exercise resolution)
    with open(os.path.join(root, "manifest.txt"), "w") as f:
        f.write("# synthetic manifest (NDIM=4)\n")
        f.write("# columns: L q beta kappa base_seed outdir\n")
        seed = 1
        for (L_, q_, b_, k_, d_) in manifest_rows:
            f.write("%d %d %.6f %.6f %d %s\n" % (L_, q_, b_, k_, seed, d_))
            seed += 1

    return betas, kappas, beta_peak_of_kappa, kappa_peak_of_beta


def selftest(ndim=4):
    """Build a synthetic campaign with known ridges, run detection, ASSERT recovery
    within one grid spacing, assert boundaries.csv + >=1 figure written."""
    tmp = tempfile.mkdtemp(prefix="phasediag_selftest_")
    try:
        rng = np.random.default_rng(20260529)
        campaign = os.path.join(tmp, "campaign")
        os.makedirs(campaign, exist_ok=True)
        betas, kappas, beta_peak_of_kappa, kappa_peak_of_beta = \
            build_synthetic_campaign(campaign, rng)

        groups, n_points, n_loaded = build_groups(campaign, ndim)
        if len(groups) != 1:
            print("SELFTEST FAIL: expected 1 (L,q) group, got %d" % len(groups))
            return 1
        g = groups[0]
        ridges = detect_boundaries(g)

        bgrid = sorted(betas)
        kgrid = sorted(kappas)
        bspace = min(bgrid[i + 1] - bgrid[i] for i in range(len(bgrid) - 1))
        kspace = min(kgrid[i + 1] - kgrid[i] for i in range(len(kgrid) - 1))

        # --- Coulomb ridge: per fixed kappa, beta_peak must match known within one space.
        if "coulomb" not in ridges:
            print("SELFTEST FAIL: coulomb ridge not detected")
            return 1
        for (kfix, bpk, chi, err) in ridges["coulomb"]["points"]:
            want = beta_peak_of_kappa[round(kfix, 9)]
            if abs(bpk - want) > bspace + 1e-9:
                print("SELFTEST FAIL: coulomb ridge at kappa=%.3f recovered beta=%.3f, "
                      "expected %.3f (>1 grid spacing %.3f)" % (kfix, bpk, want, bspace))
                return 1

        # --- Higgs ridge: per fixed beta, kappa_peak must match known within one space.
        if "higgs" not in ridges:
            print("SELFTEST FAIL: higgs ridge not detected")
            return 1
        for (bfix, kpk, chi, err) in ridges["higgs"]["points"]:
            want = kappa_peak_of_beta[round(bfix, 9)]
            if abs(kpk - want) > kspace + 1e-9:
                print("SELFTEST FAIL: higgs ridge at beta=%.3f recovered kappa=%.3f, "
                      "expected %.3f (>1 grid spacing %.3f)" % (bfix, kpk, want, kspace))
                return 1

        # --- prefactor robustness: argmax must not move when V changes. -----------
        g2 = Group(g.L, g.q, ndim)
        g2.volume = g.volume * 1000.0  # rescale prefactor
        # recompute chi with the rescaled volume by reloading
        for p in enumerate_points(campaign):
            ts = load_timeseries(p["tsfile"])
            if ts:
                g2.add_point(p["beta"], p["kappa"], ts)
        r2 = detect_boundaries(g2)
        for name in ("coulomb", "higgs"):
            a = [(round(x[0], 9), round(x[1], 9)) for x in ridges[name]["points"]]
            b = [(round(x[0], 9), round(x[1], 9)) for x in r2[name]["points"]]
            if a != b:
                print("SELFTEST FAIL: ridge '%s' moved under V rescale (not prefactor-robust)" % name)
                return 1

        # --- outputs ----------------------------------------------------------------
        rc = run(campaign, ndim, "figs")
        if rc != 0:
            print("SELFTEST FAIL: run() returned %d" % rc)
            return 1
        bpath = os.path.join(campaign, "boundaries.csv")
        if not os.path.isfile(bpath):
            print("SELFTEST FAIL: boundaries.csv not written")
            return 1
        with open(bpath) as f:
            nrows = sum(1 for ln in f if ln and not ln.startswith("L,"))
        if nrows < 1:
            print("SELFTEST FAIL: boundaries.csv has no data rows")
            return 1

        figdir = os.path.join(campaign, "figs")
        figs = [x for x in (os.listdir(figdir) if os.path.isdir(figdir) else [])
                if x.endswith(".png")]
        # matplotlib may legitimately be absent; only require a figure if it imported.
        have_mpl = True
        try:
            import matplotlib  # noqa: F401
        except Exception:
            have_mpl = False
        if have_mpl and len(figs) < 1:
            print("SELFTEST FAIL: no figure written despite matplotlib present")
            return 1

        # junction sanity: should fall inside the grid box.
        jt = junction(ridges)
        if jt is not None:
            if not (min(bgrid) - bspace <= jt[0] <= max(bgrid) + bspace and
                    min(kgrid) - kspace <= jt[1] <= max(kgrid) + kspace):
                print("SELFTEST FAIL: junction (%.3f,%.3f) outside grid" % jt)
                return 1

        print("SELFTEST PASS (groups=%d, coulomb_pts=%d, higgs_pts=%d, figs=%d, "
              "junction=%s)" % (
                  len(groups), len(ridges["coulomb"]["points"]),
                  len(ridges["higgs"]["points"]),
                  len(figs), ("(%.3f,%.3f)" % jt) if jt else "none"))
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="U(1)/SU(N) (beta,kappa) phase-diagram analyzer + plotter")
    ap.add_argument("campaign", nargs="?", help="campaign output dir (from u1_scan_campaign)")
    ap.add_argument("--ndim", type=int, default=4, help="spacetime dim for V=L**ndim (default 4)")
    ap.add_argument("--figdir", default="figs", help="figure subdir (default: <campaign>/figs)")
    ap.add_argument("--selftest", action="store_true", help="run the synthetic self-test and exit")
    args = ap.parse_args(argv)

    if args.selftest:
        sys.exit(selftest(args.ndim))
    if not args.campaign:
        ap.error("campaign dir required (or use --selftest)")
    sys.exit(run(args.campaign, args.ndim, args.figdir))


if __name__ == "__main__":
    main()
