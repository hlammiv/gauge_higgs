#!/usr/bin/env python3
"""SU(2)->2T (BT) frozen-Higgs (beta,kappa) MASTER phase diagram -- q=n-master LINE format (cf
scripts/u1f_q5_master.py), NOT heatmaps. Boundaries from the (beta,kappa) plane only; no mode_zq.

TOPOLOGY (center-blind Higgs: spin-3 -> 2T keeps the Z_2 center => at least TWO lines):
  (1) GAUGE line = chi_plaq specific-heat peak beta_f(kappa). FSS DECIDES its nature per kappa:
        * peak height L-INDEPENDENT  -> CROSSOVER (dashed/open)   [pure-SU(2) kappa=0 is a crossover by
          asymptotic freedom -- NO 4D bulk transition; and the L8/L12 data show no volume growth for
          kappa<=2 either: ratio12/8 = 0.77/1.14/0.87 ~ 1 => crossover across the reachable range]
        * peak height GROWS ~volume  -> 1st-order TRANSITION (solid/filled). Expected only as kappa->inf,
          where the link freezes to the DISCRETE 2T subgroup -> pure-2T bulk transition beta_f=2.24 (GLL,
          a genuine 1st-order transition because 2T is discrete). A critical endpoint kappa* separates them.
  (2) MATTER (Higgs) line = where L_link (link energy) orders: chi_link(kappa) peak, else L_link=0.5
        crossing, per beta. L_link is beta-INDEPENDENT (~0.13*kappa) in the disordered regime so this line
        is ~horizontal; it descends toward the gauge line and the two meet near a triple region.
  SU(2)->2T is nonabelian => the high-beta/low-kappa region is weak-coupling DECONFINED, NOT a true
  Coulomb phase (FS: no continuous residual symmetry).

Reusable: edit DATADIR / LS / KMAX / FSS_GROW. Reads su2_bt_B/B_L{L}_b*_k*.out.
"""
import re, glob, os, numpy as np, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

DATADIR = "su2_bt_B"; BETA_F_2T = 2.24; LS = [16, 12, 8]; KMAX = 8.5
BMIN_LINE = 1.4    # ignore the strong-coupling chi_plaq hump below this when locating the gauge peak
FSS_GROW = 1.30    # peak-height ratio (largest L / L=8) above which we call it a 1st-order transition

def grab(f, key):
    if not os.path.exists(f): return np.nan
    m = re.search(rf"{key}\s*=\s*([-\d.eE+]+)", open(f).read()); return float(m.group(1)) if m else np.nan

def load(L):
    rows = []
    for f in glob.glob(f"{DATADIR}/B_L{L}_b*_k*.out"):
        m = re.search(r"_b([0-9.]+)_k([0-9.]+)\.out$", f)
        if not m: continue
        cp = grab(f, "chi_plaq")
        if np.isfinite(cp):
            rows.append((float(m.group(1)), float(m.group(2)), cp, grab(f, "L_link"), grab(f, "chi_link")))
    return sorted(rows)

def hysteresis_bands(kaps, betas, gap_thr=0.05):
    """First-order coexistence from hot-vs-cold plaquette. cold = B_L8cold_*, hot = B_L8_* (600-meas).
    Returns {kappa: (beta_lo, beta_hi)} over the betas where cold(ordered)-hot(disordered) > gap_thr,
    and the hot-spinodal beta (hot plaq crosses the 0.78 mid-level) per kappa. NOTE: precise beta_f
    needs LLR; HMC hysteresis only brackets it (cold stays metastable far below beta_f)."""
    bands = {}; hotsp = {}
    for k in kaps:
        kk = int(k) if float(k).is_integer() else k
        loop = []
        hb = []   # (beta, hot_plaq)
        for b in betas:
            ph = grab(f"{DATADIR}/B_L8_b{b}_k{kk}.out", "plaquette")
            pc = grab(f"{DATADIR}/B_L8cold_b{b}_k{kk}.out", "plaquette")
            if np.isfinite(ph) and np.isfinite(pc):
                if (pc - ph) > gap_thr: loop.append(b)
                hb.append((b, ph))
        if loop: bands[k] = (min(loop), max(loop))
        # hot spinodal: lowest beta where hot has ordered (plaq>0.80); else None
        ord_b = [b for b, p in hb if p > 0.80]
        if ord_b: hotsp[k] = min(ord_b)
    return bands, hotsp

data = {L: load(L) for L in LS}
allk = sorted({k for L in LS for (b, k, cp, ll, cl) in data[L]})

def peak_beta(pts):
    """Locate the gauge feature in beta. A 1st-order SPIKE (max chi_plaq>1) -> its argmax location.
    A broad CROSSOVER hump -> chi-weighted centroid of the top (>=70% of max) points, which is stable
    against the argmax wandering on a flat noisy hump. Returns (beta, peak_height)."""
    pts = sorted((b, c) for b, c in pts if b >= BMIN_LINE and np.isfinite(c))
    if len(pts) < 2: return None
    bs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts]); ymax = ys.max()
    if ymax > 1.0:                                   # 1st-order coexistence spike -> argmax
        return float(bs[int(np.argmax(ys))]), float(ymax)
    sel = ys >= 0.7 * ymax                            # crossover -> centroid of the hump top
    return float(np.sum(bs[sel]*ys[sel]) / np.sum(ys[sel])), float(ymax)

def chi_by_k(L, idx):   # idx: 2=chi_plaq col(=cp), here cp is col2; build {k:[(b,val)]}
    byk = collections.defaultdict(list)
    for b, k, cp, ll, cl in data[L]:
        byk[k].append((b, cp))
    return byk

# --- gauge line beta_f(kappa) (finest L) + FSS nature ---
peaks = {L: {} for L in LS}
for L in LS:
    byk = collections.defaultdict(list)
    for b, k, cp, ll, cl in data[L]:
        byk[k].append((b, cp))
    for k, pts in byk.items():
        pk = peak_beta(pts)
        if pk: peaks[L][k] = pk           # (beta_f, chi_peak)

gauge = []   # (kappa, beta_f, is_transition)
for k in allk:
    bf = next((peaks[L][k][0] for L in LS if k in peaks[L]), None)
    if bf is None: continue
    h8 = peaks[8].get(k, (None, None))[1]
    hbig = next((peaks[L][k][1] for L in LS if k in peaks[L]), None)   # finest available
    grew = (h8 and hbig and (hbig / h8) >= FSS_GROW)
    # 1st-order also flagged by a huge single-volume chi_plaq peak (specific heat ~ V at coexistence)
    is_trans = bool(grew) or (hbig is not None and hbig > 1.0)
    gauge.append((k, bf, is_trans))
gauge.sort()

# --- matter ordering: L_link is a CROSSOVER (smooth 0->1, chi_link flat) -> draw the L_link=0.5
#     iso-contour as a (dashed) crossover line, NOT a transition. beta-averaged (L_link is beta-flat). ---
llk = collections.defaultdict(list)
for L in LS:
    for b, k, cp, ll, cl in data[L]:
        if np.isfinite(ll) and cp < 5:   # exclude the 1st-order coexistence spike from the smooth crossover
            llk[k].append(ll)
llk = {k: float(np.mean(v)) for k, v in llk.items() if v}
kk = sorted(llk)
matter_k = None       # kappa where L_link crosses 0.5 (mid-crossover)
for a, b2 in zip(kk, kk[1:]):
    if llk[a] < 0.5 <= llk[b2]:
        matter_k = a + (0.5-llk[a])/(llk[b2]-llk[a])*(b2-a); break

# --- 1st-order coexistence from hot-vs-cold hysteresis (replaces the unreliable single-start spikes) ---
HBETAS = [2.1, 2.2, 2.3, 2.4, 2.5, 2.6]
bands, hotsp = hysteresis_bands([4, 5, 6, 8], HBETAS)
band_k = sorted(bands)                                  # kappas with a measured coexistence loop
kstar = (min(band_k) + max(k for k, b, t in gauge if not t and k < min(band_k)))/2 if band_k else None

# ============================ PLOT (q=n master style) ============================
fig, ax = plt.subplots(figsize=(10, 7.2))

# gauge CROSSOVER line: only kappas with NO 1st-order loop (low kappa). chi-centroid (smooth).
cross = sorted((g[0], g[1]) for g in gauge if g[0] not in bands and (not band_k or g[0] < min(band_k)))
ax.plot([b for k, b in cross], [k for k, b in cross], '--o', color='darkred', lw=1.8, ms=6,
        mfc='white', zorder=4, label='gauge freezing CROSSOVER (chi_plaq centroid; FSS: no vol. growth)')

# 1st-order COEXISTENCE band (hysteresis): shade [beta_lo, beta_hi] per kappa; hot-spinodal envelope
if band_k:
    lo = [bands[k][0] for k in band_k]; hi = [bands[k][1] for k in band_k]
    ax.fill_betweenx(band_k, lo, hi, color='darkred', alpha=0.25, zorder=2,
                     label='1st-order COEXISTENCE (hot/cold hysteresis loop; widens with kappa)')
    mid = [0.5*(bands[k][0]+bands[k][1]) for k in band_k]      # beta_f bracket = band midpoint
    ax.plot(mid, band_k, '--s', color='darkred', lw=2.0, ms=7, zorder=5,
            label='beta_f bracket (band midpoint; LLR for precise value)')
# bridge crossover -> band at kappa*
if cross and band_k:
    ax.plot([cross[-1][1], bands[min(band_k)][0]], [cross[-1][0], min(band_k)],
            ':', color='darkred', lw=1.2, alpha=0.6, zorder=3)

# matter ordering CROSSOVER line (horizontal: L_link beta-flat; dashed = crossover, chi_link has no peak)
if matter_k is not None:
    ax.axhline(matter_k, color='tab:purple', lw=2.0, ls='--', zorder=4,
               label=f'matter ordering CROSSOVER (L_link=0.5, kappa~{matter_k:.1f}; chi_link flat)')

# kappa->inf pure-2T endpoint
ax.axvline(BETA_F_2T, ls=':', color='red', alpha=0.7, label='pure-2T beta_f=2.24 (kappa->inf, 1st order)')

# confined/deconfined boundary across ALL kappa: crossover beta (low k) + band lower edge (high k)
bnd = sorted([(k, b) for k, b in cross] + [(k, bands[k][0]) for k in band_k])
bk = [k for k, b in bnd]; bb = [b for k, b in bnd]
ax.fill_betweenx(bk, 0.3, bb, alpha=0.08, color='tab:green', zorder=0)
ax.fill_betweenx(bk, bb, 3.2, alpha=0.08, color='tab:blue', zorder=0)
ax.text(0.95, KMAX*0.55, "CONFINED\n(area law;\n2T keeps Z2)", fontsize=11, color='tab:green', ha='center', va='center')
ax.text(2.72, 1.7, "weak-coupling\nDECONFINED\n(\"Coulomb-like\";\nno true Coulomb)",
        fontsize=8.5, color='tab:blue', ha='center', va='center')
ax.text(2.78, KMAX*0.86, "HIGGS\n(matter\nordered)", fontsize=9.5, color='tab:purple', ha='center', va='center')
# annotations
ax.annotate("kappa=0: pure SU(2) -> CROSSOVER only\n(no 4D bulk transition)", xy=(cross[0][1] if cross else 1.7, 0.02),
            xytext=(0.85, KMAX*0.17), fontsize=7.5, color='dimgray',
            arrowprops=dict(arrowstyle='->', color='dimgray', lw=0.8))
if kstar:
    ax.annotate(f"critical endpoint\nkappa*~{kstar:.0f} (crossover->1st order)", xy=(bands[min(band_k)][0], min(band_k)),
                xytext=(1.25, min(band_k)+0.6), fontsize=7.5, color='darkred',
                arrowprops=dict(arrowstyle='->', color='darkred', lw=0.8))
ax.text(2.52, 3.25, "precise beta_f -> LLR\n(HMC can't tunnel this line)", fontsize=7, color='darkred', ha='center', style='italic')

ax.set_ylim(-0.1, KMAX); ax.set_xlim(0.4, 3.0)
ax.set_xlabel('beta (gauge)'); ax.set_ylabel('kappa (Higgs hopping)')
ax.set_title('SU(2)->2T (BT) frozen-Higgs phase diagram (spin-3, |phi|=1)\n'
             'gauge freezing crossover -> 1st-order at kappa->inf (beta_f=2.24); matter/Higgs line; no mode_zq')
ax.legend(fontsize=8, loc='upper left'); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig('su2_bt_master.png', dpi=130); print('wrote su2_bt_master.png')

print('gauge line: ' + ', '.join(f'k{k}:b{b:.2f}{"[1st-order]" if t else "[crossover]"}' for k, b, t in gauge))
print('L_link(kappa): ' + ', '.join(f'k{k}:{llk[k]:.3f}' for k in sorted(llk)))
print(f'matter ordering crossover (L_link=0.5) at kappa ~ {matter_k:.2f}' if matter_k else 'matter <0.5 in range')
