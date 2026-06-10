#!/usr/bin/env python3
"""
2D-LLR reconstruction + (beta,kappa) reweighting + first-order transition-line / triple-point
location for the compact U(1) + charge-q Higgs model.

Forked from .llr_ref/scripts/llr2d_reconstruct.py (+ the line-termination logic of
llr_transition_line.py). The integration / continuous-tile-moment / equal-height machinery is
model-agnostic and reused; only the model-specific pieces are adapted:

  * window axis-2 is E2 = -B (scalar hopping), conjugate to kappa  (NOT a 2nd gauge action) ->
    order-parameter normalization re-derived: NPLAQ = V*D*(D-1)/2, NLINKS = V*D,
    <plaquette> = 1 - <A>/NPLAQ,  <link_energy> = <B>/(2 V D) = -<E2>/(2 V D).
  * reweight exponent  lnrho - beta*A - kappa*E2 = lnrho - beta*A + kappa*B = -(beta*A - kappa*B)
    matches the physical exp(-S) sans S_pot (S_pot at fixed lambda,q is folded into rho during
    sampling -> lambda,q are NOT reweight axes).
  * --amax (default 8) replaces the hardcoded |a|<5 divergence filter (physical |a1|~beta,
    |a2|~kappa near a transition can legitimately approach the transition couplings).
  * WEIGHTED lstsq: each edge weighted by 1/(sd1^2+sd2^2+eps) from the ANE2 sd columns.
  * CONNECTIVITY check: the reached-cell set must be one connected component (cold+hot must
    meet) else integrate() leaves ungauge-fixed lnrho islands.
  * THREE transitions (confined<->Coulomb, confined<->Higgs, Coulomb<->Higgs): each gives an
    equal-height (Maxwell) line in (beta,kappa); the triple point is their common intersection.

DENSITY-OF-STATES SCOPE: rho(A,B) carries ONLY beta,kappa dependence. S_pot(lambda,q) is folded
in at the run's FIXED (lambda,q); you cannot reweight in lambda/mass with this rho.

Usage:
  u1_llr_reconstruct.py run.out [run2.out ...]
      [--amax 8] [--q 1] [--beta-range 0 4] [--kappa-range -0.5 1.5]
      [--heatmap out.png] [--obs plaq|link|chi_beta|chi_kappa] [--triple] [--plot lnrho.png]
"""
import sys, argparse, re
import numpy as np
from numpy.linalg import lstsq


# ---------------------------------------------------------------------------- logsumexp
def logsumexp(x):
    x = np.asarray(x, float); m = x.max(); return m + np.log(np.sum(np.exp(x - m)))


# ---------------------------------------------------------------------------- parsing
def parse(fname, amax=8.0):
    cells = []; geo = {}; sd = []; rmdiag = {}
    for l in open(fname):
        if l.startswith("LLR2D"):
            for tok in l.replace("|", " ").split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    try: geo[k] = float(v)
                    except ValueError: pass
            m = re.search(r"E1=\[([-\d.eE+]+),([-\d.eE+]+)\]", l)
            if m:
                geo["E1top"] = float(m.group(1)); geo["E1bot"] = float(m.group(2))
        elif l.startswith("RMDIAG:"):
            t = l.split()
            # RMDIAG: A0 E2_0 res_A res_E2 slope_a1 slope_a2 acc_g acc_hmc nwiden
            try:
                key = (round(float(t[1]), 1), round(float(t[2]), 1))
                rmdiag[key] = dict(res_A=float(t[3]), res_E2=float(t[4]),
                                   acc_hmc=float(t[8]), nwiden=int(t[9]))
            except (ValueError, IndexError):
                pass
        elif l.startswith("ANE2:"):
            t = l.split()
            a1, a2 = float(t[3]), float(t[4])
            if not (np.isfinite(a1) and np.isfinite(a2) and abs(a1) < amax and abs(a2) < amax):
                continue
            cells.append([float(t[1]), float(t[2]), a1, a2])
            sd.append([float(t[5]), float(t[6])])
    return np.array(cells), geo, np.array(sd), rmdiag


def parse_multi(files, amax=8.0):
    geo = None; acc = {}; accsd = {}; rmd = {}
    for f in files:
        cells, g, sd, rmdiag = parse(f, amax)
        if geo is None: geo = g
        rmd.update(rmdiag)
        for (E1, E2, a1, a2), (s1, s2) in zip(cells, sd):
            k = (round(E1, 1), round(E2, 1))
            acc.setdefault(k, []).append((a1, a2))
            accsd.setdefault(k, []).append((s1, s2))
    out = []; outsd = []
    for (E1, E2), lst in acc.items():
        arr = np.array(lst); out.append([E1, E2, arr[:, 0].mean(), arr[:, 1].mean()])
        sarr = np.array(accsd[(E1, E2)]); outsd.append([sarr[:, 0].mean(), sarr[:, 1].mean()])
    return np.array(out), geo, np.array(outsd), {k: len(v) for k, v in acc.items()}, rmd


# ---------------------------------------------------------------------------- grid / integrate
def grid_index(cells, geo):
    E1top = geo["E1top"]; step1 = geo["step1"]; nperp = int(geo["nperp"]); step2 = geo["step2"]
    c0 = geo["c0"]; c1 = geo["c1"]; c2 = geo["c2"]
    idx = []
    for E1, E2, a1, a2 in cells:
        i = int(round((E1top - E1) / step1))
        ridge = c0 + c1 * E1 + c2 * E1 * E1
        j = int(round((E2 - ridge) / step2)) + nperp
        idx.append((i, j))
    return idx


def connectivity(idx):
    """Return list of connected components (4-neighbour) of the reached (i,j) set."""
    cellset = set(idx); seen = set(); comps = []
    for start in cellset:
        if start in seen: continue
        stack = [start]; comp = []
        while stack:
            ij = stack.pop()
            if ij in seen: continue
            seen.add(ij); comp.append(ij)
            i, j = ij
            for nb in [(i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)]:
                if nb in cellset and nb not in seen:
                    stack.append(nb)
        comps.append(comp)
    return comps


def integrate(cells, idx, sd=None, e2_weight=1.0):
    """Weighted least-squares lnrho from the (a1,a2) gradient over (i,j)-adjacent edges.
    e2_weight<1 downweights the E2-direction (i,j+1) edges, which carry a2 (the matter
    slope). a2 is noisier than a1 (gauge) when the matter is under-converged in the gap;
    downweighting keeps the A-direction of lnrho governed by the clean a1 (so the gauge
    transition lands at the right beta), at the cost of a softer B-direction."""
    n = len(cells); pos = {ij: k for k, ij in enumerate(idx)}
    rows = []; rhs = []; wts = []
    eps = 1e-6
    for k, (i, j) in enumerate(idx):
        for (ni, nj) in [(i + 1, j), (i, j + 1)]:
            if (ni, nj) in pos:
                k2 = pos[(ni, nj)]
                x1 = cells[k][:2]; x2 = cells[k2][:2]
                g1 = cells[k][2:]; g2 = cells[k2][2:]
                d = 0.5 * (g1 + g2) @ (x2 - x1)
                r = np.zeros(n); r[k2] = 1; r[k] = -1; rows.append(r); rhs.append(d)
                if sd is not None and len(sd):
                    w = 1.0 / (sd[k, 0] ** 2 + sd[k, 1] ** 2 + sd[k2, 0] ** 2 + sd[k2, 1] ** 2 + eps)
                else:
                    w = 1.0
                if (ni, nj) == (i, j + 1):       # E2-direction edge (carries a2) -> downweight
                    w *= e2_weight * e2_weight
                wts.append(np.sqrt(w))
    # gauge fix: lnrho = 0 at the most-disordered reached cell (max E1)
    kref = int(np.argmax(cells[:, 0]))
    r = np.zeros(n); r[kref] = 1; rows.append(r); rhs.append(0.0); wts.append(1.0)
    A = np.array(rows); b = np.array(rhs); w = np.array(wts)
    Aw = A * w[:, None]; bw = b * w
    lnrho, res, *_ = lstsq(Aw, bw, rcond=None)
    return lnrho


def curl_residual(cells, idx):
    """Circulation of the (a1,a2) field around each unit cell -> localizes non-integrable cells."""
    pos = {ij: k for k, ij in enumerate(idx)}
    bad = []
    for k, (i, j) in enumerate(idx):
        sq = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
        if not all(s in pos for s in sq): continue
        circ = 0.0
        for a, b in [(0, 1), (1, 2), (2, 3), (3, 0)]:
            ka, kb = pos[sq[a]], pos[sq[b]]
            ga, gb = cells[ka][2:], cells[kb][2:]
            dx = cells[kb][:2] - cells[ka][:2]
            circ += 0.5 * (ga + gb) @ dx
        if abs(circ) > 1e-6:
            bad.append((cells[k][0], cells[k][1], circ))
    return bad


# ---------------------------------------------------------------------------- free energy / Maxwell
def Fmarg(cells, lnrho, b1, b2, tol=2):
    """F(A) = -b1 A + logsumexp_j[lnrho - b2 E2] over transverse cells at each A level."""
    E1 = cells[:, 0]; E2 = cells[:, 1]
    levels = sorted(set(np.round(E1, tol)))
    Es = []; F = []
    for e in levels:
        m = np.round(E1, tol) == e
        Es.append(e); F.append(-b1 * e + logsumexp(lnrho[m] - b2 * E2[m]))
    return np.array(Es), np.array(F)


def peakdiff(Es, F):
    """height(frozen=small-A=most-negative E1 level... ) - height(disordered); None if not bimodal.
    Es ascending; index 0 = smallest A = ordered/frozen end."""
    if len(F) < 3: return None
    j = 1 + int(np.argmin(F[1:-1]))
    if not (F[j] < F[j - 1] and F[j] < F[j + 1]): return None
    return F[:j + 1].max() - F[j:].max()


def betac(cells, lnrho, fix, val, lo, hi, npts=300):
    """Maxwell equal-height locator: zero crossing of the two-basin free-energy gap."""
    xs = np.linspace(lo, hi, npts); xk = []; dk = []
    for x in xs:
        b1 = val if fix == 'b1' else x
        b2 = val if fix == 'b2' else x
        Es, F = Fmarg(cells, lnrho, b1, b2); d = peakdiff(Es, F)
        if d is not None: xk.append(x); dk.append(d)
    if len(xk) < 3: return None, 0
    sl, ic = np.polyfit(xk, dk, 1)
    return -ic / sl, len(xk)


def S1mean(cells, lnrho, b1, b2):
    E1 = cells[:, 0]; E2 = cells[:, 1]
    lw = lnrho - b1 * E1 - b2 * E2; lw -= lw.max(); w = np.exp(lw)
    return np.sum(E1 * w) / np.sum(w)


def betac_steep(cells, lnrho, fix, val, lo, hi, npts=500):
    """Steepest drop of <A> over the reweight (robust locator; does NOT certify first-order)."""
    xs = np.linspace(lo, hi, npts)
    O = np.array([S1mean(cells, lnrho, val if fix == 'b1' else x, val if fix == 'b2' else x) for x in xs])
    g = np.gradient(O, xs); return xs[int(np.argmax(np.abs(g)))]


# ---------------------------------------------------------------------------- continuous-tile moments
def _j0_ub_vu(p, h):
    p = np.asarray(p, float); s = p * h; a = np.abs(s)
    logJ0 = np.empty_like(p); ub = np.empty_like(p); vu = np.empty_like(p)
    small = a < 1e-3; large = a > 30.0; mid = ~(small | large)
    logJ0[small] = np.log(2.0 * h) + s[small] ** 2 / 6.0
    ub[small] = p[small] * h * h / 3.0
    vu[small] = h * h / 3.0
    sm = s[mid]; pm = p[mid]
    logJ0[mid] = np.log(2.0) + np.log(np.sinh(sm) / pm)
    ub[mid] = h / np.tanh(sm) - 1.0 / pm
    vu[mid] = 1.0 / (pm * pm) - h * h / (np.sinh(sm) ** 2)
    sl = s[large]; pl = p[large]
    logJ0[large] = a[large] - np.log(np.abs(pl))
    ub[large] = np.sign(pl) * h - 1.0 / pl
    vu[large] = 1.0 / (pl * pl)
    return logJ0, ub, vu


def cell_moments(cells, lnrho, h1, h2, b1, b2):
    E1 = cells[:, 0]; E2 = cells[:, 1]; a1 = cells[:, 2]; a2 = cells[:, 3]
    lJp, ub1, vu1 = _j0_ub_vu(a1 - b1, h1)
    lJq, ub2, vu2 = _j0_ub_vu(a2 - b2, h2)
    logW = lnrho - b1 * E1 - b2 * E2 + lJp + lJq
    logW -= logW.max(); w = np.exp(logW); W = w.sum()
    s1 = E1 + ub1; s1sq = E1 * E1 + 2 * E1 * ub1 + (vu1 + ub1 * ub1)
    s2 = E2 + ub2; s2sq = E2 * E2 + 2 * E2 * ub2 + (vu2 + ub2 * ub2)
    S1 = np.sum(w * s1) / W; S2 = np.sum(w * s2) / W
    return S1, np.sum(w * s1sq) / W - S1 * S1, S2, np.sum(w * s2sq) / W - S2 * S2


# ---------------------------------------------------------------------------- normalization
def norms_from_header(files):
    """(NPLAQ, NLINKS, V) from 'U1q<q> D Nt Nx' in the LLR2D header."""
    for f in files:
        for l in open(f):
            if l.startswith("LLR2D"):
                t = l.split()
                try:
                    D = int(t[2]); Nt = int(t[3]); Nx = int(t[4])
                    V = Nt * Nx ** (D - 1)
                    return V * D * (D - 1) // 2, V * D, V, D
                except Exception:
                    pass
            if l.startswith("ANE2:"): break
    return 1536, 4096, 256, 4


# ---------------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", nargs="+")
    ap.add_argument("--amax", type=float, default=8.0)
    ap.add_argument("--q", type=int, default=1)
    ap.add_argument("--beta-range", type=float, nargs=2, default=[0.0, 4.0])
    ap.add_argument("--kappa-range", type=float, nargs=2, default=[-0.5, 1.5])
    ap.add_argument("--triple", action="store_true",
                    help="trace the three equal-height lines and report the triple point")
    ap.add_argument("--heatmap", default=None, help="write a (beta,kappa) order-parameter heatmap")
    ap.add_argument("--obs", choices=["plaq", "link", "chi_beta", "chi_kappa"], default="plaq")
    ap.add_argument("--plot", default=None, help="scatter plot of lnrho(A,E2)")
    args = ap.parse_args()

    if len(args.file) == 1:
        cells, geo, sd, rmd = parse(args.file[0], args.amax)
        print(f"1 run, {len(cells)} cells (amax={args.amax})")
    else:
        cells, geo, sd, nrep, rmd = parse_multi(args.file, args.amax)
        print(f"{len(args.file)} runs combined -> {len(cells)} unique cells")
    if len(cells) < 4:
        sys.exit("need >=4 reached cells to reconstruct")

    idx = grid_index(cells, geo)

    # connectivity (disconnected ribbon -> ungauge-fixed lnrho islands)
    comps = connectivity(idx)
    if len(comps) > 1:
        comps.sort(key=len, reverse=True)
        print(f"WARNING: reached set is {len(comps)} disconnected components "
              f"(sizes {[len(c) for c in comps]}). lnrho gauge-fix only pins the largest; "
              f"islands carry an undetermined additive constant -> ensure cold+hot passes MEET.")

    lnrho = integrate(cells, idx, sd)
    bad = curl_residual(cells, idx)
    print(f"cells: {len(cells)}   A {cells[:,0].min():.0f}..{cells[:,0].max():.0f}   "
          f"E2 {cells[:,1].min():.0f}..{cells[:,1].max():.0f}   curl-flagged cells: {len(bad)}")
    # RMDIAG: report poorly-converged cells
    poor = [k for k, d in rmd.items() if d.get("acc_hmc", 1.0) < 0.3 or d.get("nwiden", 0) >= 3]
    if poor:
        print(f"RMDIAG: {len(poor)} cells with acc_hmc<0.3 or nwiden>=3 (down-weighted by sd; "
              f"check the B-window sizing there)")

    blo, bhi = args.beta_range; klo, khi = args.kappa_range
    # gauge-bulk line (confined<->Coulomb): vary beta at kappa=0 (pure-gauge limit)
    b1c = betac_steep(cells, lnrho, 'b2', 0.0, blo, bhi)
    print(f"kappa=0  -> beta_c (gauge bulk, steepest <A>) = {b1c:.3f}   "
          f"(compact-U(1) bulk ~1.01 Langfeld-Lucini)")
    bm, nb = betac(cells, lnrho, 'b2', 0.0, blo, bhi)
    if bm is not None:
        print(f"           beta_c (equal-height, {nb} resolved) = {bm:.3f}")

    if args.triple:
        trace_lines(cells, lnrho, args)

    if args.heatmap:
        make_heatmap(cells, lnrho, geo, args)

    if args.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        plt.figure(figsize=(7, 5.5))
        sc = plt.scatter(cells[:, 0], cells[:, 1], c=lnrho, s=40, cmap="viridis")
        plt.colorbar(sc, label="ln rho (arb const)")
        plt.xlabel("A = sum_plaq(1-cos)"); plt.ylabel("E2 = -B")
        plt.title("2D-LLR reconstructed ln rho(A, -B)  [U(1)+Higgs]")
        plt.tight_layout(); plt.savefig(args.plot, dpi=120); print("plot ->", args.plot)


def trace_lines(cells, lnrho, args):
    """Trace the equal-height (Maxwell) first-order line beta_c(kappa) and find where it
    terminates (crossover endpoint, per llr_transition_line.py)."""
    klo, khi = args.kappa_range; blo, bhi = args.beta_range
    print("\n# equal-height first-order line beta_c(kappa)  (None = barrier closed = crossover):")
    fo = {}; cross = []
    for kappa in np.linspace(klo, khi, 25):
        bc, nres = betac(cells, lnrho, 'b2', kappa, blo, bhi)
        if bc is not None and blo <= bc <= bhi and nres >= 3:
            fo[round(kappa, 3)] = bc
            print(f"  kappa={kappa:+.3f}  beta_c={bc:.3f}  ({nres} resolved)")
        else:
            cross.append(round(kappa, 3))
    if cross and fo:
        kfo = sorted(fo)
        cpos = [c for c in cross if c > min(kfo)]
        endpoint = (max(kfo) + min(cpos)) / 2 if cpos else None
        if endpoint is not None:
            print(f"# first-order signal terminates near kappa ~ {endpoint:.3f} "
                  f"(beyond: crossover; the line endpoint bounds the triple-point neighbourhood)")
    if not fo:
        print("  (no resolved barrier over the scanned kappa range at this volume; "
              "use betac_steep as a locator hint or widen the band / increase statistics)")


def make_heatmap(cells, lnrho, geo, args):
    NPLAQ, NLINKS, V, D = norms_from_header(args.file)
    h1 = geo["step1"] / 2.0; h2 = geo["step2"] / 2.0
    blo, bhi = args.beta_range; klo, khi = args.kappa_range
    n = 140
    bb = np.linspace(blo, bhi, n); kk = np.linspace(klo, khi, n)
    Amap = np.empty((n, n)); E2map = np.empty((n, n))
    for ii, b in enumerate(bb):
        for jj, k in enumerate(kk):
            A, _, E2, _ = cell_moments(cells, lnrho, h1, h2, b, k)
            Amap[jj, ii] = A; E2map[jj, ii] = E2
    if args.obs == "plaq":
        Z = 1.0 - Amap / NPLAQ                       # 1=ordered(small A), 0=disordered
        label = r"$\langle{\rm plaq}\rangle = 1-\langle A\rangle/N_{\rm plaq}$"
    elif args.obs == "link":
        Z = (-E2map) / (2.0 * NLINKS)                # <link_energy> = <B>/(2 V D)
        label = r"$\langle{\rm link}\rangle = \langle B\rangle/(2VD)$"
    else:
        from scipy.ndimage import gaussian_filter
        if args.obs == "chi_beta":
            sm = gaussian_filter(Amap, sigma=1.5)
            Z = -np.gradient(sm, bb, axis=1) * V / (NPLAQ * NPLAQ)
            label = r"$\chi_\beta = -\partial\langle A\rangle/\partial\beta \cdot V/N_{\rm plaq}^2$"
        else:
            sm = gaussian_filter(-E2map, sigma=1.5)
            Z = np.gradient(sm, kk, axis=0) * V / (NLINKS * NLINKS)
            label = r"$\chi_\kappa = \partial\langle B\rangle/\partial\kappa \cdot V/N_{\rm links}^2$"
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 6))
    cmap = "inferno" if args.obs.startswith("chi") else "viridis"
    pc = plt.pcolormesh(bb, kk, Z, shading="auto", cmap=cmap)
    plt.colorbar(pc, label=label)
    # overlay the steepest-<A>-drop locator line vs kappa
    line = [(betac_steep(cells, lnrho, 'b2', k, blo, bhi), k) for k in np.linspace(klo, khi, 18)]
    xs = [p[0] for p in line]; ys = [p[1] for p in line]
    plt.plot(xs, ys, ".-", color="cyan", lw=1.5, ms=5, label="locator (steepest <A>)")
    plt.legend(loc="upper right", framealpha=0.6)
    plt.xlabel(r"$\beta$"); plt.ylabel(r"$\kappa$")
    plt.title(f"U(1)+q={args.q} Higgs from 2D-LLR rho(A,-B)")
    plt.tight_layout(); plt.savefig(args.heatmap, dpi=120)
    print(f"heatmap -> {args.heatmap}   (NPLAQ={NPLAQ} NLINKS={NLINKS} V={V} D={D})")


if __name__ == "__main__":
    main()
