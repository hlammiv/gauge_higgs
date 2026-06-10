#!/usr/bin/env python3
"""U(1)+charge-q phase diagram from the lucia tempered (beta,kappa) scan.

Honest two-axis reading with the gauge-sector bistability shown as physics:

  kappa <= 0.5  : matter uncondensed, hot==cold (CURED). Gauge axis decides:
                  Confined (small beta, high rho_M, charge-1 confined)
                  vs Coulomb (large beta, low rho_M, charge-1 free).
  kappa >= 0.75 : matter CONDENSED (L_link agrees hot==cold, Higgs in matter),
                  but the GAUGE sector is bistable -- hot trapped low-plaq
                  (confined-like), cold high-plaq (deconfined). This is the
                  Fradkin-Shenker first-order Higgs/confinement line; HMC cannot
                  tunnel it at L=8. The two branches BRACKET the transition.
                  For q>=4 the hot (confined) branch keeps charge-1 CONFINED
                  (residual Z_q, sigma1 up, |P1| down); q=2 does not.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, "u1_l8_temper_lucia")
C = dict(beta=0, kappa=1, plaq=4, Llink=6, rhoM=12, sig1=18, absP1=22)

def load(q, beta, start):
    f = os.path.join(SRC, f"q{q}_b{beta}_{start}", "summary.csv")
    rows = []
    for ln in open(f):
        if ln.startswith("#") or ln.startswith("beta"):
            continue
        p = ln.strip().split(",")
        if len(p) < 23:
            continue
        rows.append([float(p[i]) for i in (C["beta"],C["kappa"],C["plaq"],
                      C["Llink"],C["rhoM"],C["sig1"],C["absP1"])])
    return np.array(rows)  # beta,kappa,plaq,Llink,rhoM,sig1,absP1

BETAS  = ["0.30","0.70","1.10","1.60","2.20"]
QS     = [2,4,8]
LLINK_HIGGS = 1.5     # matter condensed
RHOM_DECONF = 3.0
SIG1_CONF   = 0.15

CONF="#c0392b"; COUL="#2980b9"; HIGGS="#27ae60"; BIST="#7f8c8d"

fig, axes = plt.subplots(1, 3, figsize=(16,5.4), sharey=True)
for ax,q in zip(axes,QS):
    for b in BETAS:
        h=load(q,b,"hot"); c=load(q,b,"cold")
        if h.size==0: continue
        bb=float(b)
        for i in range(len(h)):
            k=h[i,1]; j=np.argmin(np.abs(c[:,1]-k))
            plaq_h,Ll_h,rho_h,sig_h,P_h = h[i,2:7]
            plaq_c,Ll_c,rho_c,sig_c,P_c = c[j,2:7]
            gap=abs(plaq_h-plaq_c)
            condensed = (Ll_h>LLINK_HIGGS) and (Ll_c>LLINK_HIGGS)
            if not condensed and gap<0.10:
                # CURED single-phase point: gauge axis
                if rho_h>RHOM_DECONF and sig_h>SIG1_CONF:
                    col=CONF
                else:
                    col=COUL
                ax.scatter(bb,k,c=col,s=300,marker="s",ec="k",lw=0.5,zorder=2)
            elif condensed and gap>0.10:
                # gauge-bistable Higgs (FS first-order coexistence)
                ax.scatter(bb,k,c=HIGGS,s=300,marker="s",ec="k",lw=0.5,zorder=2)
                ax.scatter(bb,k,marker="|",c=BIST,s=300,lw=2.4,zorder=3)  # bistable bar
                # q>=4 hot-branch charge-1 confinement
                if sig_h>SIG1_CONF and P_h<0.35:
                    ax.scatter(bb,k,facecolors="none",ec="yellow",s=360,
                               marker="s",lw=2.4,zorder=4)
            else:
                # condensed & cured (rare) -> plain Higgs
                ax.scatter(bb,k,c=HIGGS,s=300,marker="s",ec="k",lw=0.5,zorder=2)
    ttl=f"q={q}  " + ("(q≤4: Z_q has NO Coulomb)" if q<=4 else "(q≥5: Z_q Coulomb)")
    ax.set_title(ttl,fontsize=11)
    ax.set_xlabel(r"$\beta=4/g^2$  (gauge coupling)")
    ax.set_xlim(0.05,2.45); ax.set_ylim(-0.18,2.7)
    ax.axhline(0.625,color="k",ls=":",lw=1,alpha=0.6)
    ax.text(2.4,0.66,"matter condenses",ha="right",va="bottom",fontsize=8,alpha=0.7)
    ax.grid(alpha=0.25)
axes[0].set_ylabel(r"$\kappa$  (hopping / Higgs)")

handles=[Patch(fc=CONF,ec="k",label="Confined (κ≤0.5, cured)"),
         Patch(fc=COUL,ec="k",label="Coulomb (κ≤0.5, cured)"),
         Patch(fc=HIGGS,ec="k",label="Higgs in matter (κ≥0.75)"),
         Line2D([],[],marker="|",color=BIST,ls="",ms=14,mew=2.4,
                label="gauge bistable (FS 1st-order, HMC can't tunnel)"),
         Line2D([],[],marker="s",mfc="none",mec="yellow",mew=2.4,ls="",ms=13,
                label="hot branch: charge-1 confined (residual Z_q)")]
fig.legend(handles=handles,loc="upper center",ncol=3,fontsize=9.5,
           bbox_to_anchor=(0.5,1.02),frameon=False)
fig.suptitle("U(1)+charge-q Higgs (L=8⁴, tempered κ-ladder, hot+cold). "
             "Lower κ cured; upper κ = gauge-bistable Higgs/confinement coexistence.",
             y=1.10,fontsize=11.5)
fig.tight_layout(rect=[0,0,1,0.95])
out=os.path.join(ROOT,"u1_phase_diagram_q248.png")
fig.savefig(out,dpi=130,bbox_inches="tight")
print("wrote",out)
PY = None
