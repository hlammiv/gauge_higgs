#!/usr/bin/env python3
"""Parse the frozen-U(1) campaign point-job outputs and build phase-diagram heatmaps + validation tables.
Only COMPLETE jobs (with a '# <A>=' observable line) are used; running jobs (header only) are skipped."""
import sys, os, glob, re
import numpy as np

HDR = re.compile(r"L=(\d+)\s+beta=([\d.]+)\s+kappa=([\d.]+)\s+q=(\d+)")
OBS = re.compile(r"<A>=([\-\d.]+)\s+<B>=([\-\d.]+)\s+<plaq>=([\-\d.]+)\s+<cos>_link=([\-\d.]+)")

def parse(dirs):
    recs = []
    for d in dirs:
        for f in glob.glob(os.path.join(d, "job_point_*.out")):
            txt = open(f).read()
            h, o = HDR.search(txt), OBS.search(txt)
            if not h or not o:
                continue
            recs.append(dict(L=int(h[1]), beta=float(h[2]), kappa=float(h[3]), q=int(h[4]),
                             A=float(o[1]), B=float(o[2]), plaq=float(o[3]), cos=float(o[4])))
    return recs

def main():
    dirs = sys.argv[1:] or ["u1f_campaign", "u1f_campaign_lenore"]
    recs = parse(dirs)
    print(f"# complete point jobs: {len(recs)}")
    qs = sorted(set(r['q'] for r in recs))
    Ls = sorted(set(r['L'] for r in recs))
    print(f"# q values: {qs}   L values: {Ls}")
    for q in qs:
        for L in Ls:
            n = sum(1 for r in recs if r['q'] == q and r['L'] == L)
            if n: print(f"#   q={q} L={L}: {n} points")

    # ---- VALIDATION 1: pure-gauge beta_c (kappa=0, q=2): plaq(beta) per L ----
    print("\n## VALIDATION: pure-gauge transition (kappa=0, q=2) -- compact-U(1) bulk ~1.01")
    for L in Ls:
        pts = sorted([r for r in recs if r['q'] == 2 and r['kappa'] == 0.0 and r['L'] == L], key=lambda r: r['beta'])
        if pts:
            print(f"  L={L:2d}: " + "  ".join(f"b={r['beta']:.2f}:plaq={r['plaq']:.3f}" for r in pts))

    # ---- VALIDATION 2: q=2 deep-kappa gauge transition (-> Z_2 self-dual ~0.4407) ----
    print("\n## VALIDATION: q=2 deep-kappa gauge transition (Bowler: -> Z_2 self-dual beta_eff~0.4407)")
    for k in sorted(set(r['kappa'] for r in recs if r['q'] == 2 and r['kappa'] >= 4.0)):
        pts = sorted([r for r in recs if r['q'] == 2 and r['kappa'] == k], key=lambda r: r['beta'])
        if pts:
            print(f"  kappa={k}: " + "  ".join(f"b={r['beta']:.2f}:plaq={r['plaq']:.3f},cos={r['cos']:.3f}" for r in pts))

    # ---- heatmaps per q (L=8): <cos>_link and <plaq> over (beta,kappa) ----
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print("no matplotlib:", e); return
    os.makedirs("u1f_campaign_analysis", exist_ok=True)
    L = 8
    for q in qs:
        pts = [r for r in recs if r['q'] == q and r['L'] == L]
        if len(pts) < 6: continue
        bs = sorted(set(r['beta'] for r in pts)); ks = sorted(set(r['kappa'] for r in pts))
        cosg = np.full((len(ks), len(bs)), np.nan); plg = np.full((len(ks), len(bs)), np.nan)
        for r in pts:
            i, j = ks.index(r['kappa']), bs.index(r['beta'])
            cosg[i, j] = r['cos']; plg[i, j] = r['plaq']
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for a, g, t in [(ax[0], cosg, "<cos>_link (matter order)"), (ax[1], plg, "<plaq> (gauge order)")]:
            im = a.imshow(g, origin="lower", aspect="auto", cmap="viridis",
                          extent=[min(bs), max(bs), min(ks), max(ks)])
            a.set_xlabel("beta"); a.set_ylabel("kappa"); a.set_title(t); fig.colorbar(im, ax=a)
        fig.suptitle(f"frozen U(1)+charge-{q} Higgs, L={L}")
        fig.tight_layout(); p = f"u1f_campaign_analysis/phase_q{q}_L{L}.png"; fig.savefig(p, dpi=110); plt.close(fig)
        print(f"  wrote {p}  ({len(pts)} pts)")

if __name__ == "__main__":
    main()
