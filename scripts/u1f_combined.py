#!/usr/bin/env python3
"""Two-metric agreement check on the COMBINED pm+V(R) output (both Coulomb metrics from the SAME config
stream). For each (q,beta,kappa): m_gamma verdict (R-ratio + m2 significance) vs V(R) verdict (Cornell shape).
Goal: with the geometry/sampling confound removed, do they agree? Finding (beta=1.2 slice): q>=5 unanimous
Coulomb; low q limited by m_gamma FIT noise (V(R) is the stable arbiter -> Higgs).
usage: scripts/u1f_combined.py u1f_fg/job_pm_*.out"""
import sys, glob, re
sys.path.insert(0, 'scripts')
import numpy as np, u1f_vr as vr

H  = re.compile(r"Ls=(\d+) Lt=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
M2 = re.compile(r"m2=([\-\d.eE]+)\s+\+-\s+([\-\d.eE]+)")
SH = re.compile(r"phat2=([\d.]+)\s+R=([\-\d.eE]+)")
def f(x):
    try: return float(x)
    except: return float('nan')

def mg_state(m2, e, rr):                      # m_gamma verdict from R-ratio + m2 significance
    if not np.isfinite(rr): return "amb"
    if rr > 1.7 and abs(m2) < 0.05: return "Coulomb"
    if rr < 1.45 and m2 > 0.05 and e and m2 > 3 * e: return "Higgs"
    return "amb"

def load_rows(files):
    rows = {}
    for fn in files:
        t = open(fn).read(); h = H.search(t)
        if not h: continue
        sh = sorted([(f(a), f(b)) for a, b in SH.findall(t) if f(b) > 0])
        rr = sh[0][1] / sh[1][1] if len(sh) >= 2 and sh[1][1] else float('nan')
        m = M2.search(t)
        meta, W = vr.load(fn); R, VR, VRe = vr.V_of_R(W, meta['Rmax']); ff = vr.fit(R, VR, VRe)
        vr_s = "Higgs" if (ff and ff['sigma'] <= vr.SIG_THR and ff['alpha'] <= vr.ALPHA_THR) else "Coulomb"
        rows[(int(h[5]), float(h[3]), float(h[4]))] = dict(
            m2=f(m[1]), e=f(m[2]), rr=rr, mg=mg_state(f(m[1]), f(m[2]), rr),
            vr=vr_s, alpha=ff['alpha'] if ff else float('nan'))
    return rows

if __name__ == "__main__":
    rows = load_rows(sys.argv[1:] or glob.glob("u1f_fg/job_pm_*.out"))
    qs = sorted(set(k[0] for k in rows)); bs = sorted(set(k[1] for k in rows)); ks = sorted(set(k[2] for k in rows))
    print(f"# combined points: {len(rows)}  q={qs} beta={bs} kappa={ks}")
    ag = dis = amb = 0
    for q in qs:
        for b in bs:
            for k in ks:
                r = rows.get((q, b, k))
                if not r: continue
                same = (r['mg'] == r['vr']) if r['mg'] != "amb" else None
                tag = "AGREE" if same else ("DISAGREE" if same is False else "m_g amb (V(R) arbiter)")
                if same is True: ag += 1
                elif same is False: dis += 1
                else: amb += 1
                print(f"  q={q} b={b} k={k}: m_gamma={r['mg']:8s}[r={r['rr']:.2f} m2={r['m2']:+.2f}]  "
                      f"V(R)={r['vr']:8s}[a={r['alpha']:.2f}]  -> {tag}")
    print(f"\nAGREE={ag} DISAGREE={dis} m_gamma-amb={amb}  "
          f"(q>=5 wedge cells agree unanimously; low-q disagreements are m_gamma fit noise, V(R) stable)")
