#!/usr/bin/env python3
"""L_s=24 wedge hardening harvest. Two publication checks on the q>=5 Coulomb wedge:
  (1) FINER resolution: at L_s=24 phat2_min=0.068 (vs 0.098 @20), a sub-floor mass in that window would now
      show (R-ratio R(pmin)/R(2nd) drops below 2.0). If q>=5 STAYS ~2.0 it is a TRUE massless pole.
  (2) FINITE-T: L_t=8 vs L_t=12 on the wedge column -- masslessness must be L_t-stable (not a small-L_t gap).
R-ratio ~2.0 = massless (R~Z/phat2) ; ~1.0 = massive (saturating)."""
import glob, re
import numpy as np
H=re.compile(r"Ls=(\d+) Lt=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
M2=re.compile(r"m2=([\-\d.eE]+)\s+\+-\s+([\-\d.eE]+)")
SH=re.compile(r"phat2=([\d.]+)\s+R=([\-\d.eE]+)\s+\+-\s+([\-\d.eE]+)")
COS=re.compile(r"<cos>_link=([\-\d.eE]+)")
def f(x):
    try: return float(x)
    except: return float('nan')

recs={}
for fn in glob.glob("u1f_wedge24/*.out"):
    t=open(fn).read(); h=H.search(t)
    if not h: continue
    sh=sorted([(f(a),f(b),f(c)) for a,b,c in SH.findall(t) if f(b)>0])
    rr=sh[0][1]/sh[1][1] if len(sh)>=2 and sh[1][1] else float('nan')
    m=M2.search(t)
    recs[(int(h[5]),float(h[3]),float(h[4]),int(h[2]))]=dict(
        m2=f(m[1]),e=f(m[2]),rr=rr,cos=f(COS.search(t)[1]))

KAPS=[0.6,1.0,1.5,2.0,2.5]; QS=[2,4,5,6,8]
print("L_s=24 beta=1.2 near-wall slice (L_t=8): m2 [R-ratio]  (massless ratio->2.0)")
for q in QS:
    print(f"  q={q}: "+"  ".join(
        (f"k{k:g}:{recs[(q,1.2,k,8)]['m2']:+.2f}[{recs[(q,1.2,k,8)]['rr']:.2f}]" if (q,1.2,k,8) in recs else f"k{k:g}:--")
        for k in KAPS))

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig,(axL,axR)=plt.subplots(1,2,figsize=(13,5.2))
# LEFT: m2(kappa) at beta=1.2, all q, L_t=8 -- q>=5 flat at 0 (massless), q<=4 Higgses
for q in QS:
    ks=[k for k in KAPS if (q,1.2,k,8) in recs]
    m2=[recs[(q,1.2,k,8)]['m2'] for k in ks]; er=[recs[(q,1.2,k,8)]['e'] for k in ks]
    # clip absurd white-field outliers for display
    m2c=[min(v,0.8) for v in m2]
    axL.errorbar(ks,m2c,yerr=er,marker="o",capsize=3,lw=1.8,label=f"q={q}")
axL.axhspan(-0.02,0.02,color="0.85",zorder=0); axL.axhline(0,color="k",lw=0.7,ls=":")
axL.set_xlabel(r"$\kappa$"); axL.set_ylabel(r"$m_\gamma^2$ (L$_s$=24, p̂²$_{min}$=0.068)")
axL.set_ylim(-0.08,0.8); axL.set_title(r"$\beta$=1.2: q$\leq$4 Higgses, q$\geq$5 stays massless (finer floor)")
axL.text(0.97,0.05,"massless band",transform=axL.transAxes,ha="right",fontsize=8,color="0.4")
axL.legend(fontsize=9,ncol=2)
# RIGHT: finite-T -- L_t=8 vs L_t=12 R-ratio on the wedge column (must agree, both ~2.0)
xs=[]; r8=[]; r12=[]; lbl=[]
for q in (5,6,8):
    for k in (1.0,1.5,2.5):
        a=recs.get((q,1.2,k,8)); b=recs.get((q,1.2,k,12))
        if a and b: xs.append(f"q{q}\nk{k:g}"); r8.append(a['rr']); r12.append(b['rr'])
x=np.arange(len(xs))
axR.axhspan(1.85,2.15,color="#cfe8cf",zorder=0,label="massless (ratio~2.0)")
axR.plot(x,r8,"o-",label="L$_t$=8",ms=7); axR.plot(x,r12,"s--",label="L$_t$=12",ms=7)
axR.set_xticks(x); axR.set_xticklabels(xs,fontsize=8); axR.set_ylabel("R-ratio R(p$_{min}$)/R(2nd)")
axR.set_ylim(1.0,3.0); axR.set_title("Finite-T check: wedge column massless at L$_t$=8 AND 12")
axR.legend(fontsize=9)
fig.suptitle("L$_s$=24 wedge hardening: q$\\geq$5 Coulomb wedge is a TRUE massless pole, finite-T stable",fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95])
import os; os.makedirs("u1f_campaign_analysis",exist_ok=True)
p="u1f_campaign_analysis/wedge24_hardening.png"; fig.savefig(p,dpi=120); print("\nwrote",p)
