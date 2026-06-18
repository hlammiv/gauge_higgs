#!/usr/bin/env python3
"""SU(2)->2T first-order freezing line beta_f(kappa) from the LLR Maxwell construction (su2_bt_B/llr_L8_k*.out),
overlaid on the hot/cold hysteresis coexistence band and the pure-2T endpoint. This is the PRECISE first-order
line that HMC cannot tunnel -- the deliverable of the LLR campaign. Error bar = resid-filter sensitivity
(resid_max in 0.06-0.10); points with a large spread (sparse usable cells, e.g. deep kappa) are drawn open."""
import glob, re, os, numpy as np
import importlib.util as u
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

spec = u.spec_from_file_location('bf', os.path.join(os.path.dirname(__file__), 'su2_llr_betaf.py'))
bf = u.module_from_spec(spec); spec.loader.exec_module(bf)

# --- LLR beta_f(kappa) with filter-sensitivity error ---
pts = []
for f in sorted(glob.glob('su2_bt_B/llr_L8_k*.out')):
    m = re.search(r'llr_L8_k(\d+)\.out', f); k = int(m.group(1))
    vals = []
    for rm in (0.06, 0.08, 0.10):
        e, A, a1 = bf.load(f, resid_max=rm)
        if len(A) >= 4:
            v = bf.beta_f(A, a1)[0]
            if v: vals.append(v)
    if vals:
        pts.append((k, float(np.mean(vals)), (max(vals) - min(vals)) / 2, len(vals)))
pts.sort()

fig, ax = plt.subplots(figsize=(8.5, 6))
# hysteresis coexistence band (from the hot/cold scan: where plaq_cold-plaq_hot>0.05), per kappa
HYST = {4: (2.2, 2.6), 5: (2.2, 2.6), 6: (2.2, 2.4), 8: (2.2, 2.6)}   # measured loops (B_L8_* vs B_L8cold_*)
hk = sorted(HYST); ax.fill_betweenx(hk, [HYST[k][0] for k in hk], [HYST[k][1] for k in hk],
            alpha=0.15, color='darkred', zorder=0, label='hysteresis coexistence band (HMC, L=8)')

# pure-2T endpoint
ax.axvline(2.24, ls=':', color='gray', alpha=0.8); ax.text(2.24, 8.4, 'pure-2T\n2.24', fontsize=8, ha='center', color='gray')

# LLR points: solid if filter-stable (err<0.05), open if uncertain (sparse cells)
for k, b, err, n in pts:
    stable = err < 0.05
    ax.errorbar([b], [k], xerr=[err], fmt='o' if stable else 's', ms=9, color='darkred',
                mfc='darkred' if stable else 'white', mec='darkred', lw=1.5, capsize=4, zorder=5)
    ax.annotate(f'{b:.2f}{"" if stable else "?"}', (b, k), textcoords='offset points', xytext=(8, 6), fontsize=8, color='darkred')
# connect the stable points
stab = [(b, k) for k, b, err, n in pts if err < 0.05]
if len(stab) >= 2:
    o = sorted(stab, key=lambda p: p[1]); ax.plot([p[0] for p in o], [p[1] for p in o], '-', color='darkred', lw=2, zorder=4,
              label='LLR beta_f(kappa) first-order line (Maxwell, L=8)')
ax.plot([], [], 's', mfc='white', mec='darkred', label='LLR (filter-sensitive: sparse cells, needs more stats)')

ax.set_xlabel('beta_f (gauge)'); ax.set_ylabel('kappa (Higgs hopping)')
ax.set_ylim(3, 9); ax.set_xlim(1.6, 3.0)
ax.set_title('SU(2)->2T (BT): precise first-order freezing line beta_f(kappa) via LLR\n'
             '(density-of-states Maxwell construction; HMC cannot tunnel this line)')
ax.legend(fontsize=8, loc='lower right'); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig('su2_llr_line.png', dpi=130); print('wrote su2_llr_line.png')
print('LLR beta_f(kappa): ' + ', '.join(f'k{k}:{b:.3f}+/-{e:.3f}' for k, b, e, n in pts))
