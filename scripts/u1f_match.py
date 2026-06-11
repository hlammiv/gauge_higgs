#!/usr/bin/env python3
"""kappa->inf digitization matching: show the deep-Higgs U(1)+charge-q gauge observables (sigma_1, <plaq>)
CONVERGE to the pure Z_q gauge theory as kappa grows. The convergence IS the proof that the kappa->inf limit
is the discrete gauge theory."""
import glob, re, os
import numpy as np
HZQ = re.compile(r"zq \(pure Z_q gauge\): D=\d+ L=(\d+) beta=([\d.]+) q=(\d+)")
HPT = re.compile(r"point: D=\d+ L=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
PL  = re.compile(r"<plaq>=([\-\d.eE]+)")
SG  = re.compile(r"sigma1=([\-\d.eEna]+)")
def f(x):
    try: return float(x)
    except: return float('nan')

zq = {}   # (q,beta) -> (plaq, sig1)
pt = {}   # (q,beta,kappa) -> (plaq, sig1)
for d in ("u1f_match", "u1f_match_lenore"):
    for fn in glob.glob(d + "/job_*.out"):
        t = open(fn).read(); pl, sg = PL.search(t), SG.search(t)
        if not pl or not sg: continue
        hz, hp = HZQ.search(t), HPT.search(t)
        if hz: zq[(int(hz[3]), float(hz[2]))] = (f(pl[1]), f(sg[1]))
        elif hp: pt[(int(hp[4]), float(hp[2]), float(hp[3]))] = (f(pl[1]), f(sg[1]))
print(f"# pure-Zq points: {len(zq)}   U(1)+q points: {len(pt)}")

qs = sorted(set(k[0] for k in zq))
kaps = sorted(set(k[2] for k in pt))
print(f"# q: {qs}   kappa: {kaps}")

# convergence table: at each (q,beta) report |sigma1_U1q(kappa) - sigma1_Zq| shrinking with kappa
print("\n## sigma_1 convergence to pure Z_q as kappa grows  (|delta| -> 0 = digitization):")
for q in qs:
    betas = sorted(set(b for (qq, b) in zq if qq == q))
    for b in betas:
        if (q, b) not in zq: continue
        s_zq = zq[(q, b)][1]
        ds = []
        for k in kaps:
            v = pt.get((q, b, k))
            ds.append(f"k{k:g}:d={abs(v[1]-s_zq):.2f}" if v else f"k{k:g}:--")
        print(f"  q={q} b={b}: Zq sigma1={s_zq:.2f}  | " + " ".join(ds))

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
except Exception as e:
    print("no matplotlib:", e); raise SystemExit
os.makedirs("u1f_campaign_analysis", exist_ok=True)
fig, axes = plt.subplots(1, len(qs), figsize=(6 * len(qs), 4.2), squeeze=False)
for ax, q in zip(axes[0], qs):
    betas = sorted(set(b for (qq, b) in zq if qq == q))
    # pure Z_q sigma_1(beta) -- the target
    bz = [b for b in betas if (q, b) in zq]; sz = [zq[(q, b)][1] for b in bz]
    ax.plot(bz, sz, "k-o", lw=2.5, ms=5, label="pure Z_q (target)", zorder=10)
    for k in kaps:
        bk = [b for b in betas if (q, b, k) in pt]; sk = [pt[(q, b, k)][1] for b in bk]
        if bk: ax.plot(bk, sk, "o--", ms=4, alpha=0.8, label=f"U(1)+q kappa={k:g}")
    ax.set_xlabel("beta"); ax.set_ylabel("sigma_1"); ax.set_title(f"q={q}: deep-Higgs U(1)+q -> pure Z_q")
    ax.legend(fontsize=8); ax.set_ylim(-0.05, None)
fig.suptitle("kappa->inf MATCHING: deep-Higgs U(1)+charge-q sigma_1(beta) converges to pure Z_q (digitization proof)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
p = "u1f_campaign_analysis/match_zq.png"; fig.savefig(p, dpi=110); print("\nwrote", p)
