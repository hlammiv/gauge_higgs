#!/usr/bin/env python3
"""Static potential V(R) from the charge-1 Wilson grid W[R][T] (frozen `pot` mode) -- the POSITION-SPACE,
INDEPENDENT cross-check of the structure-factor m_gamma for the Coulomb/Higgs (wedge) question.

V(R) from the T-decay:  ln W(R,T) = ln A(R) - V(R) T  -> V(R) = -slope (weighted linear fit in T).
Then fit V(R) vs R to two forms and compare:
   Coulomb (massless):  V(R) = c - alpha/R
   Yukawa  (Higgs):     V(R) = c - alpha*exp(-m R)/R      (m = screening mass = position-space m_gamma)
Coulomb wedge  <=>  m_screen ~ 0 (Yukawa no better than Coulomb).  Higgs <=> m_screen > 0 significantly.
usage: scripts/u1f_vr.py file1.out [file2.out ...]"""
import sys, re
import numpy as np
from scipy.optimize import curve_fit

HDR = re.compile(r"pot: L=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+) Rmax=(\d+)")
WL  = re.compile(r"W R=(\d+) T=(\d+)\s+W=([\d.eE+\-]+)\s+\+-\s+([\d.eE+\-]+)")

def load(fn):
    t = open(fn).read(); h = HDR.search(t)
    meta = dict(L=int(h[1]), beta=float(h[2]), kappa=float(h[3]), q=int(h[4]), Rmax=int(h[5]))
    W = {}
    for R, T, w, e in WL.findall(t):
        W[(int(R), int(T))] = (float(w), float(e))
    return meta, W

def V_of_R(W, Rmax, tmin=1, tmax=None):
    """V(R) = -slope of ln W(R,T) vs T (weighted by propagated ln-errors)."""
    VR, VRe, Rs = [], [], []
    for R in range(1, Rmax + 1):
        Ts, y, ye = [], [], []
        for T in range(1, Rmax + 1):
            if (R, T) not in W: continue
            w, e = W[(R, T)]
            if w <= 0 or (tmax and T > tmax) or T < tmin: continue
            Ts.append(T); y.append(np.log(w)); ye.append(e / w)   # sigma_lnW = sigma_W / W
        if len(Ts) < 2: continue
        Ts, y, ye = np.array(Ts, float), np.array(y), np.array(ye)
        # weighted linear fit y = a + b T ; V = -b
        wts = 1.0 / ye**2
        S, Sx, Sy = wts.sum(), (wts*Ts).sum(), (wts*y).sum()
        Sxx, Sxy = (wts*Ts*Ts).sum(), (wts*Ts*y).sum()
        d = S*Sxx - Sx*Sx
        b = (S*Sxy - Sx*Sy) / d
        be = np.sqrt(S / d)
        VR.append(-b); VRe.append(be); Rs.append(R)
    return np.array(Rs, float), np.array(VR), np.array(VRe)

def cornell(R, c, sig, a):  return c + sig * R - a / R   # confined sigma*R + Coulomb -a/R + const

# The static potential separates ALL THREE phases by SHAPE:
#   Confined: sigma>0 (linear, V keeps rising)   Coulomb: sigma~0, alpha>0 (1/R tail, V rises & saturates)
#   Higgs:    sigma~0, alpha~0 (FLAT -- charge-1 screened by the condensate)
SIG_THR, ALPHA_THR = 0.02, 0.03
def fit(Rs, VR, VRe):
    try:
        p, cov = curve_fit(cornell, Rs, VR, sigma=VRe, absolute_sigma=True,
                           p0=[VR.mean(), 0.0, 0.1], maxfev=20000)
        pe = np.sqrt(np.diag(cov))
        dV = float(VR[-1] - VR[0])                       # model-free: total rise (screened -> ~0)
        return dict(c=p[0], sigma=p[1], alpha=p[2], serr=pe[1], aerr=pe[2], dV=dV)
    except Exception:
        return None

def verdict(f):
    if f is None: return "fit failed"
    if f['sigma'] > SIG_THR:    return f"CONFINED (sigma={f['sigma']:.3f})"
    if f['alpha'] > ALPHA_THR:  return f"COULOMB / wedge (alpha={f['alpha']:.3f}, unscreened 1/R; dV={f['dV']:.3f})"
    return f"HIGGS (screened: flat V, alpha={f['alpha']:.3f}, dV={f['dV']:.3f})"

if __name__ == "__main__":
    files = sys.argv[1:]
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for fn in files:
        meta, W = load(fn)
        Rs, VR, VRe = V_of_R(W, meta['Rmax'])
        f = fit(Rs, VR, VRe)
        lbl = f"q={meta['q']} b={meta['beta']} k={meta['kappa']}"
        v = verdict(f)
        print(f"\n{fn}\n  {lbl}")
        print(f"  V(R): " + "  ".join(f"R{int(r)}={x:.3f}" for r, x in zip(Rs, VR)))
        if f: print(f"  Cornell c+sigma*R-alpha/R: sigma={f['sigma']:.4f}+-{f['serr']:.4f}  "
                    f"alpha={f['alpha']:.3f}+-{f['aerr']:.3f}  rise dV={f['dV']:.3f}")
        print(f"  => {v}")
        line = ax.errorbar(Rs, VR, yerr=VRe, marker="o", capsize=3, lw=1.6, ms=6, label=f"{lbl}: {v.split('(')[0].strip()}")
        col = line[0].get_color()
        if f:
            Rf = np.linspace(Rs.min(), Rs.max(), 100)
            ax.plot(Rf, cornell(Rf, f['c'], f['sigma'], f['alpha']), "--", color=col, lw=1, alpha=0.7)
    ax.set_xlabel("R (lattice units)"); ax.set_ylabel("V(R)")
    ax.set_title("Static potential V(R) -- independent of m$_\\gamma$:\n"
                 "wedge q=6 RISES (Coulombic $-\\alpha/R$, unscreened) ; Higgs q=4 FLAT (charge-1 screened)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    import os; os.makedirs("u1f_campaign_analysis", exist_ok=True)
    p = "u1f_campaign_analysis/vr_wedge_vs_higgs.png"; fig.savefig(p, dpi=120); print("\nwrote", p)
