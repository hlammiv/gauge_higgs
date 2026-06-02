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

PHASE CLASSIFICATION -- per (L,q):
  PRIMARY PATH (when the per-point summary.csv carries the transverse photon mass
  m_gamma, the new Coulomb discriminator):
    m_gamma = log(C0/|C1|) of the transverse photon correlator is SMALL with a TIGHT
    error ONLY in the Coulomb phase (massless photon); it is LARGE (and often noisy)
    in BOTH the confined and Higgs phases. m_gamma ALONE marks the Coulomb 'valley'.
    rho_M (monopole density) is a CONFINEMENT order parameter: HIGH in confined, LOW
    in both Coulomb and Higgs -- so it separates Confined from Higgs (both heavy-m_gamma).
    Each (beta,kappa) cell is labeled:
        Coulomb  if m_gamma is small AND well-resolved (small relative error);
        else Confined if rho_M is high;
        else Higgs.
    The Coulomb-confinement LINE is the boundary of the small-m_gamma region against
    the high-rho_M region; the Coulomb-Higgs LINE is the boundary of the small-m_gamma
    region against the low-rho_M large-m_gamma region (the m_gamma rise with kappa at
    fixed beta). If the small-m_gamma region reaches kappa_max without closing off, we
    report status 'no_junction' ("Coulomb wedge persists") -- the q>=5 signature.
    The (old) susceptibility ridges are still computed for the figures' overlays, but
    the PHASE LABELS are driven by m_gamma + rho_M, not by chi_plaq (which is BROKEN at
    high q: the strong Higgs transition's plaquette susceptibility truncates/mislabels
    the genuine pure-gauge Coulomb line).

  FALLBACK PATH (older runs WITHOUT m_gamma, e.g. stage-C data): we classify with
  rho_M (confinement) + link_energy (Higgs proxy) ONLY, and CLEARLY flag that
  "Coulomb vs Higgs is unresolved without m_gamma" -- the q>=5 intermediate-Coulomb
  phase has an ordered scalar (high link energy) and CANNOT be told apart from Higgs by
  link energy. The fallback never claims a Coulomb phase or a wedge: it shows only
  confined-vs-deconfined and a link-energy crossover.

  THRESHOLDS are calibrated from the per-group DATA (max-gap split + relative-error
  cut), not hard-coded magic numbers, and are printed / written to the CSV / shown on
  the figure (see _mgamma_threshold / _rho_threshold).

BOUNDARY DETECTION (legacy susceptibility ridges, used for figure overlays) -- per (L,q):
  We lay the points on the (beta,kappa) grid.
   - Coulomb<->confinement line: chi[avg_plaquette] = V*Var(plaq), the plaquette
     specific heat (PEAKS at beta_c). For each fixed-kappa ROW, the beta of maximum
     chi is a boundary point (kappa, beta*). (V*Var(monopole_density) is monotone in
     beta -- an order parameter, not a peaking susceptibility -- so it is NOT used.)
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
# Coulomb<->confinement ridge uses the PLAQUETTE SPECIFIC HEAT chi = V*Var(plaq),
# which PEAKS at beta_c. (V*Var(monopole_density) is NOT a peaking susceptibility:
# the monopole density is an ORDER PARAMETER -- large and noisy deep in the confined
# phase, dropping sharply across beta_c -- so its variance is monotone in beta and
# argmax would wrongly land on the scan EDGE. The monopole density's sharp drop
# corroborates beta_c as an order parameter, but is not used for the ridge argmax.)
# Higgs<->confinement ridge uses the link-energy susceptibility.
COULOMB_PREF = ("avg_plaquette",)
HIGGS_OBS = "link_energy"
# A transition line exists only where its susceptibility peak is SIGNIFICANT. Drop
# ridge points whose chi falls below this fraction of the ridge's max chi -- this
# truncates a line where the transition ceases (e.g. the Coulomb-confinement line in
# the deep-Higgs region, where chi_plaq -> 0 and the per-row argmax is just noise).
RIDGE_SIGNIF_FRAC = 0.15

# --- m_gamma-based phase classification (PRIMARY path) ----------------------------
# The Coulomb-cell test on m_gamma = log(C0/|C1|) (transverse photon mass proxy):
#   (1) m_gamma must be SMALL: below a threshold DERIVED FROM THE DATA -- the midpoint
#       of the largest relative gap in the sorted m_gamma values, which separates the
#       small-tight (Coulomb) cluster from the large (confined/Higgs) cluster. This is
#       a 1-D max-gap split (Jenks-style), NOT a magic constant; see _mgamma_threshold.
#   (2) m_gamma must be WELL-RESOLVED: its error must be a small fraction of its value
#       (a massless photon is measured precisely). A cell with a huge relative error is
#       a heavy/noisy photon masquerading as small -- it is NOT Coulomb. We require
#       m_gamma_err <= MGAMMA_REL_ERR_MAX * |m_gamma| (with a tiny abs-error floor so a
#       genuinely near-zero m_gamma with near-zero error still qualifies).
# MGAMMA_REL_ERR_MAX: error at most this fraction of the value. 0.5 means the value is
#   at least ~2 sigma away from zero AND, more importantly, a heavy-phase cell (large
#   m_gamma with large absolute error) fails because its error swamps the small-cluster
#   scale. Documented, conservative, and not tuned to a target answer.
MGAMMA_REL_ERR_MAX = 0.5
# Absolute error floor (lattice units): a cell with m_gamma_err <= this is "resolved"
# regardless of relative error, so a true Coulomb cell with m_gamma ~ 0 isn't excluded
# by a 0/0 relative test. Small compared with the heavy-phase m_gamma scale (~O(1)).
MGAMMA_ABS_ERR_FLOOR = 0.05
# Sanity FLOOR for the data-derived m_gamma threshold: in a degenerate all-small group
# (every cell light, e.g. a pure Coulomb/Higgs slab with no confined corner) the max-gap
# can land at a tiny value; we never let the Coulomb cut drop below this floor (a
# massless-photon log(C0/|C1|) proxy is O(0.1-0.5) on these short-Lt volumes). The cut
# is otherwise pinned to the DATA: it is always kept strictly inside the observed
# [min,max] m_gamma range (see _mgamma_threshold), so a CLEAN data-driven gap is never
# overridden by an arbitrary ceiling -- only the floor guards the degenerate case.
MGAMMA_THR_FLOOR = 0.05


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


def load_summary_row(tsfile, beta, kappa):
    """Read the per-point summary.csv that sits beside a ts_*.dat and return the row
    for (beta,kappa) as a NAME->float dict, or {} if absent/unparsable.

    The summary.csv header is a CSV column line
        beta,kappa,L,q,plaq,...,rho_M,rho_M_err[,m_gamma,m_gamma_err,m_gamma_cosh,...]
    The NEW (photon-enabled) runs append m_gamma,m_gamma_err,m_gamma_cosh,m_gamma_cosh_err;
    OLDER runs stop at rho_M_err. Both are handled: we just map whatever columns exist.
    Comment lines ('#...') above the header are skipped. We match the data row whose
    (beta,kappa) equals this point (a per-point summary may also be a 1-row file)."""
    d = os.path.join(os.path.dirname(tsfile), "summary.csv")
    if not os.path.isfile(d):
        return {}
    header = None
    best = {}
    try:
        with open(d, newline="") as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                cells = [c.strip() for c in line.split(",")]
                if header is None:
                    # the CSV header line begins with the literal 'beta'
                    if cells and cells[0].lower() == "beta":
                        header = cells
                    continue
                if len(cells) < 2:
                    continue
                try:
                    rb = float(cells[0]); rk = float(cells[1])
                except ValueError:
                    continue
                row = {}
                for name, val in zip(header, cells):
                    try:
                        row[name] = float(val)
                    except ValueError:
                        pass  # non-numeric (shouldn't happen past beta,kappa) -> skip
                if abs(rb - beta) < 1e-6 and abs(rk - kappa) < 1e-6:
                    return row  # exact match
                if not best:
                    best = row  # fall back to the sole/first row if no exact match
    except OSError:
        return {}
    return best


def _strip_comments(line):
    return line.split("#", 1)[0].strip()


def enumerate_points(campaign):
    """Enumerate campaign points as a list of dicts.

    Each point: {L, q, beta, kappa, tsfile}. PRIMARY source = walking the tree for
    ts_*.dat (the on-disk data is the ground truth). manifest.txt is only a FALLBACK:
    a RESUMED campaign overwrites manifest.txt with just its newly-run points, so the
    manifest is NOT a reliable complete index (that truncation silently dropped the
    kappa<=0.6 points of a resumed scan). Enumerating from the ts tree is robust to it.
    """
    points = []
    # ---- primary: walk the tree for every ts_*.dat that actually has data ----------
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
    if points:
        return points
    # ---- fallback: manifest.txt (rows 'L q beta kappa base_seed outdir') -----------
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
        # summary.csv-derived per-point scalars (the NEW Coulomb discriminator + the
        # confinement order parameter at the summary's full statistics).
        self.mgamma = {}                # (b,k) -> m_gamma  (log(C0/|C1|): small=Coulomb)
        self.mgamma_err = {}            # (b,k) -> m_gamma_err
        self.rho_sum = {}               # (b,k) -> rho_M (monopole density mean, summary)
        self.has_mgamma = False         # True once any cell carries a finite m_gamma
        self.has_rho_sum = False        # True once any cell carries a summary rho_M

    def add_point(self, beta, kappa, ts, summary=None):
        bk = (round(beta, 9), round(kappa, 9))
        self.betas.add(bk[0]); self.kappas.add(bk[1])
        for obs in OBSERVABLES:
            if obs in ts and len(ts[obs]) >= 2:
                s = ts[obs]
                self.mean[bk][obs] = float(np.mean(s))
                self.chi[bk][obs] = susceptibility(s, self.volume)
                self.chi_err[bk][obs] = chi_blocking_error(s, self.volume)
                self.present.add(obs)
        # absorb the per-point summary.csv scalars (m_gamma is the Coulomb locator).
        if summary:
            mg = summary.get("m_gamma")
            mge = summary.get("m_gamma_err")
            if mg is not None and not math.isnan(mg):
                self.mgamma[bk] = float(mg)
                self.mgamma_err[bk] = float(mge) if (mge is not None and
                                                     not math.isnan(mge)) else 0.0
                self.has_mgamma = True
            rho = summary.get("rho_M")
            if rho is not None and not math.isnan(rho):
                self.rho_sum[bk] = float(rho)
                self.has_rho_sum = True

    def mgamma_at(self, beta, kappa):
        return self.mgamma.get((round(beta, 9), round(kappa, 9)), float("nan"))

    def mgamma_err_at(self, beta, kappa):
        return self.mgamma_err.get((round(beta, 9), round(kappa, 9)), float("nan"))

    def rho_at(self, beta, kappa):
        """Confinement order parameter at (beta,kappa): prefer the summary.csv rho_M
        (full-statistics mean) and fall back to the ts monopole_density mean."""
        bk = (round(beta, 9), round(kappa, 9))
        if bk in self.rho_sum:
            return self.rho_sum[bk]
        return self.mean.get(bk, {}).get("monopole_density", float("nan"))

    def has_rho(self):
        return self.has_rho_sum or ("monopole_density" in self.present)

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


# ---------------------------------------------- m_gamma / rho_M phase classification
def _max_gap_split(values):
    """1-D two-cluster split by the LARGEST RELATIVE GAP in sorted unique values.

    Sort the finite values; the cut is the midpoint of the consecutive pair whose gap
    (v[i+1]-v[i]) is the largest *relative to the lower endpoint* -- this is a scale-
    aware max-gap (Jenks-style) split that cleanly separates a tight low cluster from
    a spread-out high cluster, which is exactly the m_gamma small(Coulomb)/large(heavy)
    structure. Returns (cut, lo_cluster_max, hi_cluster_min) or (None,None,None) if
    there is no gap (one cluster). Tie-breaks toward the largest ABSOLUTE gap.
    """
    v = sorted(float(x) for x in values if x is not None and not math.isnan(x))
    if len(v) < 2:
        return None, None, None
    best_i, best_score, best_abs = None, -1.0, -1.0
    for i in range(len(v) - 1):
        gap = v[i + 1] - v[i]
        if gap <= 1e-12:
            continue
        denom = abs(v[i]) + 1e-9
        score = gap / denom
        if score > best_score + 1e-12 or (abs(score - best_score) <= 1e-12 and gap > best_abs):
            best_score, best_abs, best_i = score, gap, i
    if best_i is None:
        return None, None, None
    lo_max = v[best_i]
    hi_min = v[best_i + 1]
    return 0.5 * (lo_max + hi_min), lo_max, hi_min


def _mgamma_threshold(group):
    """Data-derived m_gamma cut separating the small-tight (Coulomb) cluster from the
    large (confined/Higgs) cluster, with a sanity clamp. PRINCIPLED, not hand-tuned:
    it is the midpoint of the largest relative gap in this group's own m_gamma values
    (max-gap 1-D split). Returns None if m_gamma is unavailable for the group.

    If the split degenerates (all values one cluster) we fall back to a fraction of the
    m_gamma range, then clamp into [MGAMMA_THR_MIN, MGAMMA_THR_MAX]."""
    if not group.has_mgamma:
        return None
    vals = [group.mgamma[bk] for bk in group.mgamma]
    vals = [v for v in vals if not math.isnan(v)]
    if not vals:
        return None
    lo, hi = min(vals), max(vals)
    cut, _lo, _hi = _max_gap_split(vals)
    if cut is None:
        cut = lo + 0.5 * (hi - lo)  # no gap: just split the range
    # Pin to the DATA: the cut must lie strictly inside the observed range (a clean gap
    # is respected as-is) and never below the physical floor for a massless photon.
    eps = 1e-6 * (hi - lo + 1.0)
    cut = min(hi - eps, max(lo + eps, cut))
    return max(MGAMMA_THR_FLOOR, cut)


def _rho_threshold(group):
    """Data-derived rho_M cut separating HIGH-rho_M (confined) from LOW-rho_M
    (deconfined: Coulomb or Higgs). Uses the same max-gap split on this group's rho_M
    values -- the monopole density drops sharply across the confinement transition, so
    the steepest-drop midpoint is the natural high/low divider. Falls back to half the
    range when there is no clear gap. Returns None if rho_M is unavailable."""
    if not group.has_rho():
        return None
    betas = group.beta_grid(); kappas = group.kappa_grid()
    vals = []
    for b in betas:
        for k in kappas:
            r = group.rho_at(b, k)
            if not math.isnan(r):
                vals.append(r)
    if not vals:
        return None
    cut, _lo, _hi = _max_gap_split(vals)
    if cut is None:
        lo, hi = min(vals), max(vals)
        cut = lo + 0.5 * (hi - lo)
    return cut


def _is_coulomb_cell(mg, mge, mthr):
    """True iff cell (m_gamma=mg, err=mge) is in the Coulomb valley: m_gamma SMALL
    (below the data-derived threshold) AND WELL-RESOLVED (error a small fraction of the
    value, or below the absolute floor). A heavy/noisy photon (large mg or large mge)
    is NOT Coulomb."""
    if mthr is None or math.isnan(mg):
        return False
    if not (mg < mthr):
        return False
    if math.isnan(mge):
        return True  # no error reported: rely on the small-value test alone
    return (mge <= MGAMMA_ABS_ERR_FLOOR) or (mge <= MGAMMA_REL_ERR_MAX * abs(mg))


def classify_cells(group):
    """PRIMARY classification grid using m_gamma (Coulomb locator) + rho_M (confinement).

    Returns (P, info) where P[i,j] in {0 Confined, 1 Coulomb, 2 Higgs, nan missing}
    indexed [kappa_i, beta_j], and info carries the thresholds + a per-row Coulomb-kappa
    summary used to build the boundary lines and the junction. Each cell:
        Coulomb  if m_gamma is small & well-resolved;
        else Confined if rho_M is high;
        else Higgs.
    Caller must check group.has_mgamma first (fallback path handles its absence)."""
    betas = group.beta_grid(); kappas = group.kappa_grid()
    mthr = _mgamma_threshold(group)
    rthr = _rho_threshold(group)
    P = np.full((len(kappas), len(betas)), np.nan)
    for i, k in enumerate(kappas):
        for j, b in enumerate(betas):
            mg = group.mgamma_at(b, k)
            mge = group.mgamma_err_at(b, k)
            rho = group.rho_at(b, k)
            if math.isnan(mg) and math.isnan(rho):
                continue
            if _is_coulomb_cell(mg, mge, mthr):
                P[i, j] = 1.0                                   # Coulomb
            elif rthr is not None and not math.isnan(rho) and rho >= rthr:
                P[i, j] = 0.0                                   # Confined
            else:
                P[i, j] = 2.0                                   # Higgs
    return P, {"mgamma_threshold": mthr, "rho_threshold": rthr,
               "betas": betas, "kappas": kappas}


def mgamma_lines(group):
    """Build the two m_gamma-based boundary LINES + a junction status from the cell
    classification. Returns a dict:
        {'coulomb_confine': [(beta, kappa_lo_of_coulomb), ...],   # Coulomb vs Confined
         'coulomb_higgs':   [(beta, kappa_hi_of_coulomb), ...],   # Coulomb vs Higgs
         'junction': <junction-status dict (see junction())>,
         'mgamma_threshold':.., 'rho_threshold':..,
         'P':.., 'betas':.., 'kappas':..}

    For each fixed-beta COLUMN we find the contiguous Coulomb band in kappa (cells with
    label 1). The LOWER edge of that band borders the Confined region (rising into
    Coulomb as kappa increases off the confined floor); the UPPER edge borders the
    Higgs region (m_gamma rising again as the scalar orders). The Coulomb-confinement
    line is the set of lower edges; the Coulomb-Higgs line is the set of upper edges.

    JUNCTION / WEDGE: if for EVERY beta column that has a Coulomb band the band runs up
    to kappa_max (no Higgs cap), the Coulomb wedge persists -> status 'no_junction'
    ("Coulomb wedge persists: m_gamma stays small to kappa_max"). Otherwise the triple
    point is the highest-beta column where the band is CAPPED below kappa_max meeting
    the lowest-beta column where Coulomb first appears (the corner where both lines
    meet); we report it as a genuine 'junction'."""
    P, info = classify_cells(group)
    betas = info["betas"]; kappas = info["kappas"]
    kmax = max(kappas) if kappas else float("nan")
    ktol = (_min_step(kappas) or 0.0) * 0.5

    lower_edges = []   # (beta, kappa)  Coulomb<->Confined
    upper_edges = []   # (beta, kappa)  Coulomb<->Higgs (only where the band is capped)
    capped_betas = []  # betas whose Coulomb band is capped BELOW kappa_max (Higgs above)
    coulomb_betas = []  # betas that have any Coulomb cell at all
    for j, b in enumerate(betas):
        col = [(i, k) for i, k in enumerate(kappas) if P[i, j] == 1.0]
        if not col:
            continue
        coulomb_betas.append(b)
        k_lo = min(k for _i, k in col)
        k_hi = max(k for _i, k in col)
        lower_edges.append((b, k_lo))
        # Is there a Higgs cell ABOVE the band in this column? (label 2 at kappa>k_hi)
        higgs_above = any(P[i, j] == 2.0 and kappas[i] > k_hi + 1e-12
                          for i in range(len(kappas)))
        band_reaches_top = (k_hi >= kmax - ktol - 1e-12)
        if higgs_above and not band_reaches_top:
            upper_edges.append((b, k_hi))
            capped_betas.append(b)

    lower_edges.sort(); upper_edges.sort()

    def _no(reason, beta_t=None, kappa_t=None):
        return {"status": "no_junction", "reason": reason,
                "beta_t": beta_t, "kappa_t": kappa_t, "mismatch": None}

    if not coulomb_betas:
        jstat = _no("no Coulomb (small, well-resolved m_gamma) cells in scanned range")
    elif not capped_betas:
        # Coulomb band reaches kappa_max in every column that has it: the wedge persists.
        jstat = _no("Coulomb wedge persists: m_gamma stays small to kappa_max "
                    "(no Higgs cap on the Coulomb band; q>=5-like)")
    else:
        # Triple point: the corner where the Coulomb-Higgs cap (upper line) meets the
        # Coulomb-confinement onset. Take the lowest-beta capped column's cap kappa and
        # that column's beta -- the lower-left corner of the capped Coulomb wedge, where
        # the rising Coulomb-confinement line and the Coulomb-Higgs line converge.
        b_t = min(capped_betas)
        # kappa at the meeting: the cap of the lowest capped column.
        cap_at_bt = [k for (bb, k) in upper_edges if abs(bb - b_t) < 1e-9]
        k_t = min(cap_at_bt) if cap_at_bt else max(k for (_b, k) in lower_edges)
        jstat = {"status": "junction", "beta_t": float(b_t), "kappa_t": float(k_t),
                 "mismatch": 0.0}

    return {"coulomb_confine": lower_edges, "coulomb_higgs": upper_edges,
            "junction": jstat, "mgamma_threshold": info["mgamma_threshold"],
            "rho_threshold": info["rho_threshold"],
            "P": P, "betas": betas, "kappas": kappas}


def classify_fallback(group):
    """FALLBACK classification when m_gamma is ABSENT (older/stage-C data). We can ONLY
    distinguish confined (HIGH rho_M) from deconfined (LOW rho_M), and within deconfined
    flag a link-energy crossover -- but we CANNOT separate Coulomb from Higgs, because
    the q>=5 intermediate-Coulomb phase has an ordered scalar (HIGH link energy) yet a
    still-massless photon, indistinguishable from Higgs by link energy alone.

    Returns (P, info): P[i,j] in {0 Confined, 3 Deconfined(Coulomb/Higgs UNRESOLVED),
    nan missing}. We deliberately do NOT emit a Coulomb (1) or Higgs (2) label and do
    NOT report a wedge -- doing so would fabricate a Coulomb phase the data can't show.
    info carries the rho_M threshold, a link-energy crossover threshold (for the
    crossover line only), and the honest caveat string."""
    betas = group.beta_grid(); kappas = group.kappa_grid()
    rthr = _rho_threshold(group)
    L_thr = _range_mid(group, HIGGS_OBS)  # link-energy crossover midpoint (line only)
    P = np.full((len(kappas), len(betas)), np.nan)
    for i, k in enumerate(kappas):
        for j, b in enumerate(betas):
            rho = group.rho_at(b, k)
            link = group.mean.get((round(b, 9), round(k, 9)), {}).get(HIGGS_OBS, float("nan"))
            if math.isnan(rho) and math.isnan(link):
                continue
            if rthr is not None and not math.isnan(rho) and rho >= rthr:
                P[i, j] = 0.0                       # Confined
            else:
                P[i, j] = 3.0                       # Deconfined: Coulomb/Higgs unresolved
    caveat = ("Coulomb vs Higgs is UNRESOLVED without m_gamma (the q>=5 intermediate-"
              "Coulomb phase has an ordered scalar and cannot be distinguished from "
              "Higgs by link energy). No Coulomb phase or wedge is claimed.")
    return P, {"rho_threshold": rthr, "link_threshold": L_thr,
               "caveat": caveat, "betas": betas, "kappas": kappas}


def _significant(pts, frac=RIDGE_SIGNIF_FRAC):
    """Keep only ridge points whose chi (tuple index 2) is >= frac * the ridge's max
    chi. Truncates a ridge where the transition ceases (e.g. the Coulomb-confinement
    line in the deep-Higgs region, where chi_plaq -> 0 and the argmax is noise)."""
    if not pts:
        return pts
    cmax = max(p[2] for p in pts)
    if not (cmax > 0):
        return pts
    return [p for p in pts if p[2] >= frac * cmax]


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
        pts = _significant(pts)
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
        pts = _significant(pts)
        if pts:
            ridges["higgs"] = dict(obs=HIGGS_OBS, points=pts)

    return ridges


def _min_step(values):
    """Smallest positive consecutive gap among sorted unique `values`, or None."""
    u = sorted(set(round(float(v), 9) for v in values))
    if len(u) < 2:
        return None
    gaps = [u[i + 1] - u[i] for i in range(len(u) - 1)]
    gaps = [g for g in gaps if g > 1e-12]
    return min(gaps) if gaps else None


def junction(ridges):
    """Decide whether the Coulomb and Higgs ridges genuinely CROSS, and where.

    PHYSICS: for compact U(1) + charge-q Higgs there is a triple point at finite
    (beta_t,kappa_t) only for low charge (q<=4); for high charge (q>=5) the Coulomb
    phase persists as a WEDGE up to kappa->infinity and the two transition lines never
    meet -- there is NO triple point in range. So this MUST be able to report "no
    junction", not always manufacture one.

    The coulomb ridge maps kappa->beta*(kappa) (the gauge transition, defined only on
    the kappa rows where chi_plaq actually peaks); the higgs ridge maps beta->kappa*(beta).
    We find the self-consistent crossing that minimizes the combined mismatch, then
    ACCEPT it only if it is a REAL intersection inside the scanned box:
      (i)   the minimal mismatch is small -- at most ~one kappa grid step (otherwise the
            "crossing" is just the least-bad extrapolation of two lines that miss);
      (ii)  the crossing kappa lies within [min,max] of the Higgs ridge kappa range AND
            the Coulomb ridge has SUPPORT near that kappa, i.e. the crossing kappa is
            within ~one kappa step of the Coulomb ridge's own kappa support. This is the
            key wedge test: for q>=5 the gauge (Coulomb) line only exists at HIGH kappa
            (its low-kappa support has receded), so a crossing extrapolated DOWN to where
            the Higgs line sits has no gauge transition under it -> reject;
      (iii) the crossing beta lies within the Coulomb ridge's beta range (the beta grid
            it was sampled on), not extrapolated past the scanned betas.

    RETURNS a dict (callers branch on 'status'):
      {'status':'junction','beta_t':..,'kappa_t':..,'mismatch':..}    -- genuine crossing
      {'status':'no_junction','reason':..,'beta_t':..,'kappa_t':..,'mismatch':..}
            -- no triple point in range; beta_t/kappa_t carry the best (rejected)
               crossing estimate (or None) purely for diagnostics/annotation.
    """
    def _no(reason, beta_t=None, kappa_t=None, mismatch=None):
        return {"status": "no_junction", "reason": reason,
                "beta_t": beta_t, "kappa_t": kappa_t, "mismatch": mismatch}

    if "coulomb" not in ridges or "higgs" not in ridges:
        return _no("missing ridge (need both coulomb and higgs)")
    c = ridges["coulomb"]["points"]  # (kappa, beta*, ...)
    h = ridges["higgs"]["points"]    # (beta, kappa*, ...)
    if not c or not h:
        return _no("empty ridge")
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
        return _no("no crossing candidate")
    mismatch, beta_t, kappa_t = best

    # ---- grid-step estimates from the ridges' OWN coordinates --------------------
    # Coulomb fixed kappas ARE grid rows; Higgs fixed betas ARE grid columns. Their
    # spacing recovers the scan steps without needing the Group handed in here.
    kstep = _min_step(ck)
    if kstep is None:
        kstep = _min_step(hk)
    if kstep is None:
        kstep = 0.1  # last-resort default; a single-kappa scan can't constrain a wedge
    # tolerance for "support near": half a kappa step is enough to keep an on-grid
    # crossing while rejecting one extrapolated a full step beyond the Coulomb support.
    ktol = 0.5 * kstep

    # (i) mismatch must be small -- at most ~one kappa grid step.
    if not (mismatch <= kstep + 1e-9):
        return _no("ridges do not cross (mismatch %.3f > one kappa step %.3f)"
                   % (mismatch, kstep), beta_t, kappa_t, mismatch)

    # (ii) crossing kappa must sit within the Higgs ridge kappa range AND within ~one
    #      step of the Coulomb ridge's kappa support (the wedge test).
    hk_lo, hk_hi = float(hk.min()), float(hk.max())
    if not (hk_lo - ktol <= kappa_t <= hk_hi + ktol):
        return _no("crossing kappa %.3f outside Higgs ridge kappa range [%.3f,%.3f]"
                   % (kappa_t, hk_lo, hk_hi), beta_t, kappa_t, mismatch)
    ck_lo, ck_hi = float(ck.min()), float(ck.max())
    if not (ck_lo - ktol <= kappa_t <= ck_hi + ktol):
        # Coulomb (gauge) transition has no support at the crossing kappa: the wedge
        # persists -- the gauge line has receded to higher kappa than where the Higgs
        # line sits, so they never meet. This is the q>=5 "no triple point" signature.
        return _no("Coulomb wedge persists: no gauge transition near kappa=%.3f "
                   "(Coulomb support kappa in [%.3f,%.3f])"
                   % (kappa_t, ck_lo, ck_hi), beta_t, kappa_t, mismatch)

    # (iii) crossing beta must be within the Coulomb ridge's sampled beta range.
    cb_lo, cb_hi = float(cb.min()), float(cb.max())
    bstep = _min_step(hb) or _min_step(cb) or 0.1
    if not (cb_lo - 0.5 * bstep <= beta_t <= cb_hi + 0.5 * bstep):
        return _no("crossing beta %.3f outside Coulomb ridge beta range [%.3f,%.3f]"
                   % (beta_t, cb_lo, cb_hi), beta_t, kappa_t, mismatch)

    return {"status": "junction", "beta_t": beta_t, "kappa_t": kappa_t,
            "mismatch": mismatch}


# Convenience accessors so callers can stay terse.
def _is_junction(j):
    return isinstance(j, dict) and j.get("status") == "junction"


def _junction_bk(j):
    """(beta_t, kappa_t) tuple for a genuine junction, else None."""
    if _is_junction(j):
        return (j["beta_t"], j["kappa_t"])
    return None


def group_junction(group):
    """The ACTIVE junction status for a group and the single source of truth the CSV,
    figures, and trajectory all drive off of.

    PRIMARY (m_gamma present): the m_gamma-based triple-point / wedge status.
    FALLBACK (no m_gamma): we CANNOT honestly locate a Coulomb-bounded triple point --
    a chi_plaq/chi_link ridge crossing is NOT a Coulomb signature without the photon
    mass (the q>=5 ordered-scalar phase masquerades as Higgs). So we ALWAYS return
    'no_junction' with the honest 'unresolved' reason -- never a fabricated triple
    point. (The legacy junction() is still used for the figure overlays, not for any
    Coulomb/triple-point CLAIM.)"""
    if group.has_mgamma:
        return mgamma_lines(group)["junction"]
    return {"status": "no_junction",
            "reason": "Coulomb vs Higgs unresolved without m_gamma (no Coulomb-bounded "
                      "triple point can be claimed from chi ridges alone)",
            "beta_t": None, "kappa_t": None, "mismatch": None}


def transition_strength(ridge):
    """Heuristic transition-strength summary for one ridge (coulomb or higgs).

    Returns a dict {peak_chi, median_chi, peak_over_median, label} or None if the
    ridge is missing/empty. The label is:
      'transition-like'  -- the chi peak along the scan is TALL and LOCALIZED
                            (peak/median large), suggestive of a sharp transition;
      'crossover-like'   -- the chi profile is broad / weak (peak/median near 1),
                            suggestive of an analytic crossover.

    HEAVY CAVEAT (documented on purpose): a SINGLE lattice size CANNOT distinguish a
    true (thermodynamic) phase transition from a smooth crossover. That distinction
    requires FINITE-SIZE SCALING -- e.g. the susceptibility peak height growing ~ the
    volume (and width shrinking) as L increases. With one L the peak/median ratio only
    tells us how pronounced the feature is on THIS volume; it is a flag, not a verdict.
    """
    if not ridge or not ridge.get("points"):
        return None
    chis = [float(p[2]) for p in ridge["points"] if not math.isnan(float(p[2]))]
    if not chis:
        return None
    peak = max(chis)
    med = float(np.median(chis))
    ratio = (peak / med) if med > 0 else float("inf")
    # Heuristic cut: a peak >~2x the line's median chi is "transition-like" on this
    # volume; flatter profiles are "crossover-like". Tune is intentionally loose --
    # see the finite-size-scaling caveat above; this is only a flag.
    label = "transition-like" if ratio >= 2.0 else "crossover-like"
    return {"peak_chi": peak, "median_chi": med,
            "peak_over_median": ratio, "label": label}


# ------------------------------------------------------------------------ outputs
def write_boundaries_csv(path, groups):
    """Always write the boundary point set + per-line transition-strength flags + the
    ACTIVE (m_gamma-based when available, else legacy) junction-status summary + the
    chosen classification thresholds, one block per group. Columns:
        L,q,line_type,fixed_coord,fixed_value,peak_coord,peak_value,chi_max,chi_err,
        line_peak_chi,line_peak_over_median,line_strength_label,
        classifier,mgamma_threshold,rho_threshold,
        junction_status,beta_t,kappa_t,junction_note

    line_type values:
      coulomb_confinement / higgs_confinement -- LEGACY chi-ridge points (kept for the
        figure overlays + back-compat; on the m_gamma path these are diagnostic only).
      mgamma_coulomb_confine / mgamma_coulomb_higgs -- the NEW m_gamma-based phase-
        boundary lines (peak_value = the boundary kappa at fixed beta).
    classifier is 'm_gamma+rho_M' (primary) or 'rho_M-only (Coulomb/Higgs unresolved)'
    (fallback); the thresholds are the data-derived cuts (see _mgamma_threshold /
    _rho_threshold). The strength label is a single-volume HEURISTIC flag only --
    distinguishing a true transition from a crossover needs finite-size scaling."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["L", "q", "line_type", "fixed_coord", "fixed_value",
                    "peak_coord", "peak_value", "chi_max", "chi_err",
                    "line_peak_chi", "line_peak_over_median", "line_strength_label",
                    "classifier", "mgamma_threshold", "rho_threshold",
                    "junction_status", "beta_t", "kappa_t", "junction_note"])
        for g in groups:
            ridges = detect_boundaries(g)
            if g.has_mgamma:
                ml = mgamma_lines(g)
                j = ml["junction"]
                classifier = "m_gamma+rho_M"
                mthr = ml["mgamma_threshold"]; rthr = ml["rho_threshold"]
            else:
                ml = None
                j = group_junction(g)  # honest 'unresolved' no_junction in fallback
                _, finfo = classify_fallback(g)
                classifier = "rho_M-only (Coulomb/Higgs unresolved)"
                mthr = None; rthr = finfo["rho_threshold"]
            jstatus = j.get("status", "no_junction")
            if _is_junction(j):
                jbt = "%.6f" % j["beta_t"]; jkt = "%.6f" % j["kappa_t"]
                jnote = "mismatch=%.4g" % (j.get("mismatch") or 0.0)
            else:
                jbt = jkt = ""
                jnote = j.get("reason", "")
            mthr_s = ("%.6g" % mthr) if mthr is not None else ""
            rthr_s = ("%.6g" % rthr) if rthr is not None else ""
            strengths = {nm: transition_strength(ridges.get(nm))
                         for nm in ("coulomb", "higgs")}

            def _emit(line_type, nm, fixed_coord, peak_coord, fixed_v, peak_v, chi, err):
                s = strengths.get(nm) or {}
                w.writerow([
                    g.L, g.q, line_type, fixed_coord, "%.6f" % fixed_v,
                    peak_coord, "%.6f" % peak_v, "%.8g" % chi, "%.8g" % err,
                    ("%.8g" % s["peak_chi"]) if s else "",
                    ("%.4g" % s["peak_over_median"]) if s else "",
                    s.get("label", "") if s else "",
                    classifier, mthr_s, rthr_s,
                    jstatus, jbt, jkt, jnote])

            if "coulomb" in ridges:
                for (kfix, bpk, chi, err) in ridges["coulomb"]["points"]:
                    _emit("coulomb_confinement", "coulomb", "kappa", "beta",
                          kfix, bpk, chi, err)
            if "higgs" in ridges:
                for (bfix, kpk, chi, err) in ridges["higgs"]["points"]:
                    _emit("higgs_confinement", "higgs", "beta", "kappa",
                          bfix, kpk, chi, err)
            # m_gamma-based boundary lines (primary path only): chi columns left 0.
            if ml is not None:
                for (bfix, kbnd) in ml["coulomb_confine"]:
                    _emit("mgamma_coulomb_confine", None, "beta", "kappa",
                          bfix, kbnd, 0.0, 0.0)
                for (bfix, kbnd) in ml["coulomb_higgs"]:
                    _emit("mgamma_coulomb_higgs", None, "beta", "kappa",
                          bfix, kbnd, 0.0, 0.0)


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
                    s = transition_strength(r)
                    slab = ("  [%s peak_chi=%.4g peak/med=%.3g]"
                            % (s["label"], s["peak_chi"], s["peak_over_median"])) if s else ""
                    f.write("  %-20s ridge (obs=%s): %d points%s\n" %
                            (name, r["obs"], len(r["points"]), slab))
            if g.has_mgamma:
                ml = mgamma_lines(g)
                j = ml["junction"]
                f.write("  classifier=m_gamma+rho_M  m_gamma_threshold=%.4g  "
                        "rho_threshold=%s\n" % (
                            ml["mgamma_threshold"],
                            ("%.4g" % ml["rho_threshold"]) if ml["rho_threshold"] is not None else "n/a"))
            else:
                j = group_junction(g)
                _, finfo = classify_fallback(g)
                f.write("  classifier=rho_M-only (NO m_gamma)  rho_threshold=%s\n" % (
                    ("%.4g" % finfo["rho_threshold"]) if finfo["rho_threshold"] is not None else "n/a"))
                f.write("  CAVEAT: %s\n" % finfo["caveat"])
            if _is_junction(j):
                f.write("  triple-point region ~ (beta=%.4f, kappa=%.4f)  [mismatch=%.4g]\n"
                        % (j["beta_t"], j["kappa_t"], j.get("mismatch") or 0.0))
            else:
                f.write("  NO triple point in scanned range: %s\n" % j.get("reason", ""))


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

    # Triple-point marker (only for a GENUINE junction; otherwise the Coulomb wedge
    # persists and there is no triple point to mark).
    jt = _junction_bk(junction(ridges))
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


def _mean_grid(group, obs):
    """(betas, kappas, Z[nk,nb]) of the MEAN of obs over the grid (NaN where absent)."""
    betas = group.beta_grid(); kappas = group.kappa_grid()
    Z = np.full((len(kappas), len(betas)), np.nan)
    for i, k in enumerate(kappas):
        for j, b in enumerate(betas):
            Z[i, j] = group.mean.get((round(b, 9), round(k, 9)), {}).get(obs, np.nan)
    return np.array(betas), np.array(kappas), Z


def _overlay_coulomb(ax, ridges):
    if "coulomb" in ridges:
        pts = sorted(ridges["coulomb"]["points"], key=lambda p: p[0])
        ax.plot([p[1] for p in pts], [p[0] for p in pts], "o-", color="white",
                mec="black", mew=1.2, ms=7, lw=2, label="Coulomb-confine")


def _overlay_higgs(ax, ridges):
    if "higgs" in ridges:
        pts = sorted(ridges["higgs"]["points"], key=lambda p: p[0])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "s--", color="red",
                mec="black", mew=1.0, ms=7, lw=2, label="Higgs-confine")


def _range_mid(group, obs):
    """Midpoint of obs's mean-value range across the grid (None if absent)."""
    vals = [m[obs] for m in group.mean.values() if obs in m and not math.isnan(m[obs])]
    return 0.5 * (min(vals) + max(vals)) if vals else None


def _rho_grid(group):
    """(betas, kappas, Z[nk,nb]) of the confinement order parameter rho_M (summary.csv
    rho_M preferred, ts monopole_density fallback). NaN where absent."""
    betas = group.beta_grid(); kappas = group.kappa_grid()
    Z = np.full((len(kappas), len(betas)), np.nan)
    for i, k in enumerate(kappas):
        for j, b in enumerate(betas):
            Z[i, j] = group.rho_at(b, k)
    return np.array(betas), np.array(kappas), Z


def _mgamma_grid(group):
    """(betas, kappas, Zval[nk,nb], Zerr[nk,nb]) of m_gamma and its error. NaN where
    absent. Heavy-phase cells (large m_gamma, large error) are left as ordinary values
    here; the plot masks/annotates them using the resolved-ness test."""
    betas = group.beta_grid(); kappas = group.kappa_grid()
    Zv = np.full((len(kappas), len(betas)), np.nan)
    Ze = np.full((len(kappas), len(betas)), np.nan)
    for i, k in enumerate(kappas):
        for j, b in enumerate(betas):
            Zv[i, j] = group.mgamma_at(b, k)
            Ze[i, j] = group.mgamma_err_at(b, k)
    return np.array(betas), np.array(kappas), Zv, Ze


def _overlay_mgamma_lines(ax, ml):
    """Overlay the two m_gamma-based phase-boundary lines (each as (beta,kappa) points)."""
    cc = ml.get("coulomb_confine") or []
    ch = ml.get("coulomb_higgs") or []
    if cc:
        cc = sorted(cc)
        ax.plot([p[0] for p in cc], [p[1] for p in cc], "o-", color="white",
                mec="black", mew=1.2, ms=7, lw=2, label="Coulomb-confine (m_gamma)")
    if ch:
        ch = sorted(ch)
        ax.plot([p[0] for p in ch], [p[1] for p in ch], "s--", color="red",
                mec="black", mew=1.0, ms=7, lw=2, label="Coulomb-Higgs (m_gamma)")


def _ridge_interp(ridge):
    """f(x)->y interpolating a ridge's (p[0],p[1]) points (clamped to the endpoints),
    or None if the ridge is missing/empty. Used to classify cells by their position
    relative to the drawn boundary lines, so the panel-C colors match the lines."""
    if not ridge or not ridge.get("points"):
        return None
    pts = sorted(ridge["points"], key=lambda p: p[0])
    xs = np.array([p[0] for p in pts], dtype=float)
    ys = np.array([p[1] for p in pts], dtype=float)
    if len(xs) == 1:
        y0 = float(ys[0])
        return lambda x: y0
    return lambda x: float(np.interp(x, xs, ys))


def plot_group_panels(group, figdir):
    """3-panel per-(L,q) figure:
        (a) rho_M  -- confinement order parameter (HIGH confined, LOW Coulomb/Higgs);
        (b) m_gamma -- transverse photon mass (heatmap; SMALL+tight = Coulomb valley;
            large-error heavy-phase cells are HATCHED so the eye trusts only the valley).
            When m_gamma is absent this panel shows link_energy with the explicit
            "Coulomb/Higgs unresolved without m_gamma" caveat;
        (c) phase classification -- driven by m_gamma+rho_M (primary) or rho_M-only
            confined-vs-deconfined (fallback, with caveat).
    Returns the saved path or None. The classification thresholds are data-derived and
    are printed on panel (c)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm

    betas = group.beta_grid(); kappas = group.kappa_grid()
    if not betas or not kappas:
        return None
    ridges = detect_boundaries(group)
    have_rho = group.has_rho()
    have_link = HIGGS_OBS in group.present
    have_mg = group.has_mgamma
    be = _edges(np.array(betas)); ke = _edges(np.array(kappas))

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(18.5, 5.6))

    # --- Panel A: rho_M (confinement order parameter) --------------------------------
    if have_rho:
        _, _, ZA = _rho_grid(group)
        mA = axA.pcolormesh(be, ke, np.ma.masked_invalid(ZA), shading="auto", cmap="magma")
        fig.colorbar(mA, ax=axA).set_label(r"$\langle\rho_M\rangle$")
    else:
        axA.text(0.5, 0.5, "rho_M\nnot in data", ha="center", va="center",
                 transform=axA.transAxes)
    axA.set_title("monopole density $\\rho_M$\n(confinement order parameter)")

    # --- Panel B: m_gamma valley (primary) OR link_energy with caveat (fallback) -----
    if have_mg:
        _, _, Zv, Ze = _mgamma_grid(group)
        mthr = _mgamma_threshold(group)
        mB = axB.pcolormesh(be, ke, np.ma.masked_invalid(Zv), shading="auto",
                            cmap="viridis_r")  # reversed: SMALL (Coulomb) = bright
        fig.colorbar(mB, ax=axB).set_label(r"$m_\gamma=\log(C_0/|C_1|)$  (small=Coulomb)")
        # Hatch the heavy/ill-resolved cells: large m_gamma OR error not a small fraction
        # of the value. These are confined/Higgs (heavy photon) -- not the Coulomb valley.
        for i, k in enumerate(kappas):
            for j2, b in enumerate(betas):
                mg = Zv[i, j2]; mge = Ze[i, j2]
                if math.isnan(mg):
                    continue
                if not _is_coulomb_cell(mg, mge, mthr):
                    axB.add_patch(plt.Rectangle(
                        (be[j2], ke[i]), be[j2 + 1] - be[j2], ke[i + 1] - ke[i],
                        fill=False, hatch="xxx", edgecolor="0.3", lw=0.0))
        _overlay_mgamma_lines(axB, mgamma_lines(group))
        axB.set_title("photon mass $m_\\gamma$  (Coulomb discriminator)\n"
                      "bright valley = Coulomb; hatched = heavy/ill-resolved")
    elif have_link:
        _, _, ZB = _mean_grid(group, HIGGS_OBS)
        mB = axB.pcolormesh(be, ke, np.ma.masked_invalid(ZB), shading="auto", cmap="viridis")
        fig.colorbar(mB, ax=axB).set_label(r"$\langle E_{\mathrm{link}}\rangle$")
        _overlay_higgs(axB, ridges)
        axB.set_title("link energy (Higgs proxy)\n"
                      "Coulomb/Higgs UNRESOLVED without $m_\\gamma$")
        axB.text(0.02, 0.02,
                 "No $m_\\gamma$: cannot separate\nCoulomb from Higgs (link energy\n"
                 "is high in BOTH ordered-scalar phases)",
                 transform=axB.transAxes, ha="left", va="bottom", fontsize=7,
                 bbox=dict(boxstyle="round", fc="#ffe9e9", alpha=0.9))
    else:
        axB.text(0.5, 0.5, "neither m_gamma nor link_energy\nin data",
                 ha="center", va="center", transform=axB.transAxes)
        axB.set_title("photon mass / link energy\n(not in data)")

    # --- Panel C: phase classification ----------------------------------------------
    if have_mg:
        ml = mgamma_lines(group)
        P = ml["P"]; j = ml["junction"]; is_junc = _is_junction(j)
        mthr = ml["mgamma_threshold"]; rthr = ml["rho_threshold"]
        cmap = ListedColormap(["#3b4cc0", "#f2f0d8", "#b40426"])  # Confined, Coulomb, Higgs
        norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
        mC = axC.pcolormesh(be, ke, np.ma.masked_invalid(P), shading="auto",
                            cmap=cmap, norm=norm)
        cbar = fig.colorbar(mC, ax=axC, ticks=[0, 1, 2])
        cbar.ax.set_yticklabels(["Confined", "Coulomb", "Higgs"])
        _overlay_mgamma_lines(axC, ml)
        thr_txt = ("$m_\\gamma$ cut = %.3g (max-gap)\n$\\rho_M$ cut = %s" % (
            mthr, ("%.3g" % rthr) if rthr is not None else "n/a"))
        if is_junc:
            axC.plot([j["beta_t"]], [j["kappa_t"]], marker="*", ms=20, color="gold",
                     mec="black", mew=1.3, lw=0, label="triple point")
            axC.set_title("phase classification ($m_\\gamma+\\rho_M$)\n"
                          "(triple point detected)")
        else:
            axC.set_title("phase classification ($m_\\gamma+\\rho_M$)\n"
                          "(Coulomb wedge persists to $\\kappa_{\\max}$)")
        axC.text(0.02, 0.02, thr_txt, transform=axC.transAxes, ha="left", va="bottom",
                 fontsize=7, bbox=dict(boxstyle="round", fc="white", alpha=0.8))
        axC.legend(loc="upper left", fontsize=7, framealpha=0.85)
    else:
        # FALLBACK: confined vs deconfined only; Coulomb/Higgs are NOT separated.
        P, finfo = classify_fallback(group)
        rthr = finfo["rho_threshold"]
        cmap = ListedColormap(["#3b4cc0", "#bdbdbd"])  # Confined, Deconfined (gray)
        # P uses 0 (confined) and 3 (deconfined). Remap 3->1 for the 2-color map.
        Pm = np.where(P == 3.0, 1.0, P)
        norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)
        mC = axC.pcolormesh(be, ke, np.ma.masked_invalid(Pm), shading="auto",
                            cmap=cmap, norm=norm)
        cbar = fig.colorbar(mC, ax=axC, ticks=[0, 1])
        cbar.ax.set_yticklabels(["Confined", "Deconfined\n(Coulomb/Higgs\nUNRESOLVED)"])
        # Show the link-energy crossover line (the only deconfined-internal feature we
        # can honestly draw) but DO NOT label any region Coulomb or Higgs.
        _overlay_higgs(axC, ridges)
        axC.set_title("phase classification ($\\rho_M$ only)\n"
                      "Coulomb vs Higgs UNRESOLVED (no $m_\\gamma$)")
        axC.text(0.02, 0.02,
                 ("$\\rho_M$ cut = %s\n%s" % (
                     ("%.3g" % rthr) if rthr is not None else "n/a",
                     "No Coulomb phase / wedge claimed.")),
                 transform=axC.transAxes, ha="left", va="bottom", fontsize=7,
                 bbox=dict(boxstyle="round", fc="#ffe9e9", alpha=0.9))
        axC.legend(loc="upper left", fontsize=7, framealpha=0.85)

    for ax in (axA, axB, axC):
        ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$")
    sup = "U(1) charge q=%d, L=%d" % (group.q, group.L)
    if not have_mg:
        sup += "  [NO m_gamma: Coulomb/Higgs unresolved]"
    if len(betas) <= 8 or len(kappas) <= 8:
        sup += "  (coarse grid)"
    fig.suptitle(sup)
    fig.tight_layout()
    os.makedirs(figdir, exist_ok=True)
    out = os.path.join(figdir, "phase_q%d_L%d.png" % (group.q, group.L))
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_trajectory(groups, figdir):
    """Multi-q summary: overlay every group's ridges on one (beta,kappa) plot, and
    plot the junction (beta_t,kappa_t) vs q. Returns saved path or None."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 5.5))
    cmap = plt.get_cmap("tab10")
    junctions = []   # (q, L, beta_t, kappa_t)  -- GENUINE triple points only
    wedges = []      # (q, L, kappa_max)         -- m_gamma groups with a persisting wedge
    unresolved = []  # (q, L)                    -- fallback groups (no m_gamma): no claim

    kappa_max_seen = 0.0
    for idx, g in enumerate(sorted(groups, key=lambda gg: (gg.q, gg.L))):
        col = cmap(idx % 10)
        lab = "q=%d L=%d" % (g.q, g.L)
        if g.kappas:
            kappa_max_seen = max(kappa_max_seen, max(g.kappas))
        if g.has_mgamma:
            # PRIMARY: drive the boundaries + junction off the m_gamma classification.
            ml = mgamma_lines(g)
            cc = sorted(ml["coulomb_confine"]); ch = sorted(ml["coulomb_higgs"])
            if cc:
                ax0.plot([p[0] for p in cc], [p[1] for p in cc], "o-", color=col,
                         ms=5, lw=1.6, label="%s Coul-conf" % lab)
            if ch:
                ax0.plot([p[0] for p in ch], [p[1] for p in ch], "s--", color=col,
                         ms=5, lw=1.6, label="%s Coul-Higgs" % lab)
            j = ml["junction"]
            if _is_junction(j):
                ax0.plot([j["beta_t"]], [j["kappa_t"]], marker="*", ms=16, color=col,
                         mec="black")
                junctions.append((g.q, g.L, j["beta_t"], j["kappa_t"]))
            else:
                # No triple point: Coulomb wedge persists up to the scanned kappa_max
                # (only meaningful if there IS a Coulomb region at all).
                if ml["coulomb_confine"] or ml["coulomb_higgs"] or \
                        np.any(ml["P"] == 1.0):
                    wedges.append((g.q, g.L, max(g.kappas) if g.kappas else float("nan")))
        else:
            # FALLBACK: only legacy chi ridges + an HONEST "unresolved" flag. We do NOT
            # plot a Coulomb wedge here -- without m_gamma there is no Coulomb claim.
            ridges = detect_boundaries(g)
            if "higgs" in ridges:
                pts = sorted(ridges["higgs"]["points"], key=lambda p: p[0])
                ax0.plot([p[0] for p in pts], [p[1] for p in pts], "s--", color=col,
                         ms=5, lw=1.6, label="%s link-cross" % lab)
            unresolved.append((g.q, g.L))

    ax0.set_xlabel(r"$\beta$"); ax0.set_ylabel(r"$\kappa$")
    ax0.set_title("phase boundaries, all (q,L)")
    ax0.legend(fontsize=7, ncol=2, framealpha=0.85)

    # Right panel: THE headline figure -- (beta_t,kappa_t) vs q. Genuine junctions are
    # filled markers/lines; charges where the wedge persists (NO triple point <= kappa_max)
    # are NOT given a fabricated coordinate -- they are flagged with OPEN markers parked
    # at the top of the kappa axis with an upward arrow (kappa_t has receded past
    # kappa_max), so the figure can show "kappa_t finite for q<=4, absent for q>=5"
    # WITHOUT inventing junctions where the data shows none.
    if junctions:
        js = sorted(junctions, key=lambda t: t[0])
        qs = [t[0] for t in js]
        bt = [t[2] for t in js]
        kt = [t[3] for t in js]
        ax1.plot(qs, bt, "o-", color="tab:blue", label=r"$\beta_t$ (junction)")
        ax1.plot(qs, kt, "s-", color="tab:red", label=r"$\kappa_t$ (junction)")
        ax1.set_xlabel("charge q"); ax1.set_ylabel("triple-point coord")
        ax1.set_title(r"junction $(\beta_t,\kappa_t)$ vs $q$")
    else:
        ax1.set_xlabel("charge q"); ax1.set_ylabel("triple-point coord")
        ax1.set_title(r"junction vs $q$ (none detected)")

    if wedges:
        ws = sorted(wedges, key=lambda t: t[0])
        wq = [t[0] for t in ws]
        # park the open markers a touch above the scanned kappa_max, with up-arrows.
        ytop = (kappa_max_seen if kappa_max_seen > 0 else 1.0)
        ypark = ytop * 1.04 if ytop > 0 else 0.04
        ax1.plot(wq, [ypark] * len(wq), "^", mfc="none", mec="0.3", mew=1.4, ms=10,
                 lw=0, label=r"wedge persists: no triple point $\leq\kappa_{\max}$")
        for q in wq:
            ax1.annotate("", xy=(q, ypark + 0.06 * ytop), xytext=(q, ypark),
                         arrowprops=dict(arrowstyle="->", color="0.3", lw=1.3))
        if kappa_max_seen > 0:
            ax1.axhline(kappa_max_seen, color="0.6", ls=":", lw=1.0)
            ax1.text(0.99, kappa_max_seen, r" $\kappa_{\max}$", color="0.4",
                     ha="right", va="bottom", fontsize=8,
                     transform=ax1.get_yaxis_transform())
    if unresolved:
        # Mark fallback (no-m_gamma) charges honestly: no Coulomb/triple-point claim.
        uq = sorted(set(q for (q, _L) in unresolved))
        ax1.plot(uq, [0.0] * len(uq), "x", color="0.5", ms=9, mew=1.6, lw=0,
                 label="no $m_\\gamma$: Coulomb/Higgs unresolved")
    if junctions or wedges or unresolved:
        ax1.legend(fontsize=8, loc="best")
    else:
        ax1.text(0.5, 0.5, "no classifiable groups", ha="center",
                 va="center", transform=ax1.transAxes)

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
        summary = load_summary_row(p["tsfile"], p["beta"], p["kappa"])
        key = (p["L"], p["q"])
        if key not in groups:
            groups[key] = Group(p["L"], p["q"], ndim)
        groups[key].add_point(p["beta"], p["kappa"], ts, summary)
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
        if g.has_mgamma:
            ml = mgamma_lines(g)
            j = ml["junction"]
            print("#     classifier=m_gamma+rho_M  m_gamma_cut=%.4g (max-gap)  "
                  "rho_M_cut=%s  [rel_err<=%.2g or abs_err<=%.2g for 'resolved']" % (
                      ml["mgamma_threshold"],
                      ("%.4g" % ml["rho_threshold"]) if ml["rho_threshold"] is not None else "n/a",
                      MGAMMA_REL_ERR_MAX, MGAMMA_ABS_ERR_FLOOR))
            if _is_junction(j):
                print("#     triple point ~ (beta=%.3f, kappa=%.3f)" %
                      (j["beta_t"], j["kappa_t"]))
            else:
                print("#     no triple point: %s" % j.get("reason", ""))
        else:
            _, finfo = classify_fallback(g)
            print("#     classifier=rho_M-only (NO m_gamma)  rho_M_cut=%s" % (
                ("%.4g" % finfo["rho_threshold"]) if finfo["rho_threshold"] is not None else "n/a"))
            print("#     CAVEAT: Coulomb vs Higgs UNRESOLVED without m_gamma "
                  "(no Coulomb phase/wedge claimed)")

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
            p = plot_group_panels(g, figdir)
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


def build_synthetic_campaign(root, rng, with_mgamma=True):
    """Construct a campaign whose ts files have KNOWN ridges:
      chi[avg_plaquette] peaks at beta_peak(kappa)  (a ROW argmax over beta)
      chi[link_energy]  peaks at kappa_peak(beta)   (a COLUMN argmax over kappa)
    We bake the target variances into per-point synthetic series.

    When with_mgamma=True we ALSO emit a per-point summary.csv carrying m_gamma,
    m_gamma_err and rho_M laid out to realize a clean THREE-PHASE structure with a
    triple point, so the new m_gamma-based classifier is exercised:
        Coulomb  (small, tight m_gamma)  : beta >= 1.0 AND kappa <= 0.30
        Higgs    (large m_gamma, low rho): kappa == 0.40 (top row caps the Coulomb wedge)
        Confined (large m_gamma, high rho): beta < 1.0 AND kappa <= 0.30
    The Coulomb band at beta=1.0,1.1,1.2 is CAPPED by the Higgs row at kappa=0.40, while
    at beta=0.8,0.9 there is no Coulomb -> the two boundary lines meet near the lower-
    left corner of the capped wedge => triple point ~ (beta=1.0, kappa=0.30).

    Returns (betas, kappas, beta_peak_of_kappa, kappa_peak_of_beta, expected_triple)
    where expected_triple is (beta_t, kappa_t) or None when with_mgamma is False.
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

    # Synthetic m_gamma/rho_M design (see docstring). Returns (mg, mge, rho, label).
    def cell_phase(b, k):
        coulomb = (b >= 1.0 - 1e-9) and (k <= 0.30 + 1e-9)
        higgs = (k >= 0.40 - 1e-9)
        if coulomb and not higgs:
            return 0.15, 0.01, 0.05, "Coulomb"   # small m_gamma, TIGHT err, low rho
        if higgs:
            return 1.20, 0.30, 0.04, "Higgs"      # large m_gamma (noisy), low rho
        return 1.40, 0.45, 8.0, "Confined"        # large m_gamma (noisy), HIGH rho

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

            mg, mge, rho, _lab = cell_phase(b, k)
            plaq = _synth_series(rng, 0.4 + 0.05 * b, var_plaq, npts)
            link = _synth_series(rng, 0.01 + 0.02 * k, var_link, npts)
            higgs = _synth_series(rng, 0.8, base_var, npts)
            A = _synth_series(rng, 1.0e4, 1.0, npts)
            B = _synth_series(rng, 10.0 * k, 1.0, npts)
            mono = _synth_series(rng, rho, base_var, npts)  # monopole_density (ts fallback)

            with open(fn, "w") as f:
                f.write("# synthetic U(1)+Higgs ts. L=%d q=%d beta=%.6f kappa=%.6f\n" % (L, q, b, k))
                f.write("# columns: traj  A  B  avg_plaquette  higgs_length  link_energy  monopole_density\n")
                for t in range(npts):
                    f.write("%d %.10g %.10g %.10g %.10g %.10g %.10g\n" %
                            (t, A[t], B[t], plaq[t], higgs[t], link[t], mono[t]))

            # per-point summary.csv (NEW photon-enabled layout when with_mgamma).
            with open(os.path.join(d, "summary.csv"), "w", newline="") as sf:
                sf.write("# synthetic per-point summary\n")
                if with_mgamma:
                    sf.write("beta,kappa,L,q,plaq,plaq_err,Llink,Llink_err,phi2,phi2_err,"
                             "acceptance,exp_mdH,rho_M,rho_M_err,m_gamma,m_gamma_err,"
                             "m_gamma_cosh,m_gamma_cosh_err\n")
                    sf.write("%.6f,%.6f,%d,%d,%.6f,0.001,%.6f,0.001,0.8,0.001,0.95,1.0,"
                             "%.6f,0.05,%.6f,%.6f,0.0,0.0\n" % (
                                 b, k, L, q, 0.4 + 0.05 * b, 0.01 + 0.02 * k, rho, mg, mge))
                else:
                    sf.write("beta,kappa,L,q,plaq,plaq_err,Llink,Llink_err,phi2,phi2_err,"
                             "acceptance,exp_mdH,rho_M,rho_M_err\n")
                    sf.write("%.6f,%.6f,%d,%d,%.6f,0.001,%.6f,0.001,0.8,0.001,0.95,1.0,"
                             "%.6f,0.05\n" % (
                                 b, k, L, q, 0.4 + 0.05 * b, 0.01 + 0.02 * k, rho))

            manifest_rows.append((L, q, b, k, d))

    # manifest with outdir paths relative to the campaign's parent (exercise resolution)
    with open(os.path.join(root, "manifest.txt"), "w") as f:
        f.write("# synthetic manifest (NDIM=4)\n")
        f.write("# columns: L q beta kappa base_seed outdir\n")
        seed = 1
        for (L_, q_, b_, k_, d_) in manifest_rows:
            f.write("%d %d %.6f %.6f %d %s\n" % (L_, q_, b_, k_, seed, d_))
            seed += 1

    expected_triple = (1.00, 0.30) if with_mgamma else None
    return betas, kappas, beta_peak_of_kappa, kappa_peak_of_beta, expected_triple


def selftest(ndim=4):
    """Build a synthetic campaign with known ridges + a baked-in m_gamma/rho_M three-
    phase structure, run detection, ASSERT: legacy ridges recovered within one grid
    spacing; prefactor robustness; the m_gamma-based classifier finds the expected
    triple point; the m_gamma cell labels match the construction; outputs written. A
    SECOND (no-m_gamma) campaign verifies the HONEST fallback (no Coulomb/wedge claim)."""
    tmp = tempfile.mkdtemp(prefix="phasediag_selftest_")
    try:
        rng = np.random.default_rng(20260529)
        campaign = os.path.join(tmp, "campaign")
        os.makedirs(campaign, exist_ok=True)
        betas, kappas, beta_peak_of_kappa, kappa_peak_of_beta, expected_triple = \
            build_synthetic_campaign(campaign, rng, with_mgamma=True)

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

        # --- PRIMARY (m_gamma) classification: triple point + cell labels ----------
        if not g.has_mgamma:
            print("SELFTEST FAIL: synthetic group did not load m_gamma from summary.csv")
            return 1
        ml = mgamma_lines(g)
        j = ml["junction"]
        if not _is_junction(j):
            print("SELFTEST FAIL: m_gamma junction not detected on triple-point fixture "
                  "(status=%s reason=%s)" % (j.get("status"), j.get("reason")))
            return 1
        jt = (j["beta_t"], j["kappa_t"])
        if not (abs(jt[0] - expected_triple[0]) <= bspace + 1e-9 and
                abs(jt[1] - expected_triple[1]) <= kspace + 1e-9):
            print("SELFTEST FAIL: m_gamma triple point (%.3f,%.3f) != expected %r "
                  "(>1 grid spacing)" % (jt[0], jt[1], expected_triple))
            return 1
        # Cell-label spot checks against the construction.
        P, info = classify_cells(g)
        bi = {round(b, 9): i for i, b in enumerate(sorted(betas))}
        ki = {round(k, 9): i for i, k in enumerate(sorted(kappas))}

        def lab(b, k):
            return P[ki[round(k, 9)], bi[round(b, 9)]]
        if lab(1.20, 0.10) != 1.0:   # high-beta low-kappa -> Coulomb
            print("SELFTEST FAIL: (1.20,0.10) classified %r, expected Coulomb(1)" % lab(1.20, 0.10))
            return 1
        if lab(0.80, 0.10) != 0.0:   # low-beta low-kappa -> Confined (high rho_M)
            print("SELFTEST FAIL: (0.80,0.10) classified %r, expected Confined(0)" % lab(0.80, 0.10))
            return 1
        if lab(1.20, 0.40) != 2.0:   # high-kappa -> Higgs
            print("SELFTEST FAIL: (1.20,0.40) classified %r, expected Higgs(2)" % lab(1.20, 0.40))
            return 1
        mthr = info["mgamma_threshold"]
        if not (0.15 < mthr < 1.20):  # max-gap cut sits between the two m_gamma clusters
            print("SELFTEST FAIL: m_gamma threshold %.3f not between clusters (0.15,1.20)" % mthr)
            return 1

        # transition_strength must return a well-formed dict for each present ridge.
        for name in ("coulomb", "higgs"):
            s = transition_strength(ridges[name])
            if not s or s.get("label") not in ("transition-like", "crossover-like"):
                print("SELFTEST FAIL: transition_strength('%s') malformed: %r"
                      % (name, s))
                return 1

        # --- FALLBACK path (no m_gamma): must be HONEST, no Coulomb/wedge claim -----
        fb_root = os.path.join(tmp, "campaign_nomg")
        os.makedirs(fb_root, exist_ok=True)
        build_synthetic_campaign(fb_root, rng, with_mgamma=False)
        fgroups, _, _ = build_groups(fb_root, ndim)
        if len(fgroups) != 1 or fgroups[0].has_mgamma:
            print("SELFTEST FAIL: fallback campaign should have exactly 1 group w/o m_gamma")
            return 1
        fg = fgroups[0]
        gj = group_junction(fg)  # active junction for a no-m_gamma group
        Pfb, finfo = classify_fallback(fg)
        # The fallback grid must label ONLY Confined(0) / Deconfined(3); never Coulomb(1)
        # or Higgs(2) -- those would fabricate a distinction the data can't support.
        labs = set(np.unique(Pfb[~np.isnan(Pfb)]).tolist())
        if not labs.issubset({0.0, 3.0}):
            print("SELFTEST FAIL: fallback classification emitted non-{Confined,Deconf} "
                  "labels %r (fabricated Coulomb/Higgs)" % labs)
            return 1
        if 0.0 not in labs or 3.0 not in labs:
            print("SELFTEST FAIL: fallback should show BOTH confined and deconfined cells; got %r" % labs)
            return 1
        if "unresolved" not in finfo["caveat"].lower():
            print("SELFTEST FAIL: fallback caveat missing 'unresolved': %r" % finfo["caveat"])
            return 1
        # The fallback must NOT claim a Coulomb wedge: its junction reason must not be a
        # Coulomb-wedge string (legacy junction may say 'missing ridge' etc., which is fine).
        if _is_junction(gj):
            print("SELFTEST FAIL: fallback (no m_gamma) fabricated a triple point")
            return 1
        rc2 = run(fb_root, ndim, "figs")
        if rc2 != 0:
            print("SELFTEST FAIL: run() on fallback campaign returned %d" % rc2)
            return 1

        print("SELFTEST PASS (groups=%d, m_gamma triple=(%.3f,%.3f), m_gamma_cut=%.3f, "
              "rho_cut=%s, coulomb_pts=%d, higgs_pts=%d, figs=%d; fallback honest=%s)" % (
                  len(groups), jt[0], jt[1], mthr,
                  ("%.3f" % info["rho_threshold"]) if info["rho_threshold"] is not None else "n/a",
                  len(ridges["coulomb"]["points"]), len(ridges["higgs"]["points"]),
                  len(figs), str(labs.issubset({0.0, 3.0}) and not _is_junction(gj))))
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
