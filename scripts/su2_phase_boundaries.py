#!/usr/bin/env python3
"""SU(2)->H 3-phase diagram: draw BOTH transition boundaries on the (beta,kappa) plane.

Two transitions live in two different sectors, located by two different observables and
scan directions (locators only -- inflection/peak at finite volume; no FSS/order claim):

  * HIGGS boundary  kappa_c(beta): matter sector. Scan in kappa at fixed beta; the Higgs
    onset is the steepest rise of <L_link> (frozen gauge-invariant hopping). Separates
    confining/Coulomb (low kappa) from Higgs (high kappa).
  * CONFINING<->COULOMB boundary  beta_c(kappa): gauge sector. Scan in beta at fixed kappa;
    the gauge crossover is the peak of the gauge specific heat chi_plaq = V*Var(plaq).
    Separates confining (low beta) from Coulomb/deconfined (high beta); at small kappa this
    is the SU(2) crossover (~may be only a crossover), deep in Higgs it -> discrete freezing.

Usage: python3 scripts/su2_phase_boundaries.py [scan_dir]   (default su2scan_grid)"""
import sys, os, csv, glob, re
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "su2scan_grid"
V = 16
ORDER = [("adj","adjoint -> U(1)"),("Q8soft","spin-2 -> Q8"),
         ("2T","2T (spin-3)"),("2O","2O (spin-4)"),("2I","2I (spin-6)")]
betaf = {"2T":2.24,"2O":3.26,"2I":5.82}   # literature deep-Higgs freezing beta_f

def chi_plaq_map(rep):
    out={}
    for f in glob.glob(os.path.join(D,"ts",f"{rep}_*.dat")):
        try:
            hdr=open(f).readline()
            mb=re.search(r"beta=([\d.]+)",hdr); mk=re.search(r"kappa=([\d.]+)",hdr)
            d=np.loadtxt(f,comments="#")
            if mb and mk and d.ndim==2 and len(d)>=10:
                out[(round(float(mb.group(1)),3),round(float(mk.group(1)),3))]=V*np.var(d[:,1])
        except Exception: pass
    return out

def load(rep):
    p=os.path.join(D,f"{rep}.csv")
    if not os.path.exists(p): return None
    L={}
    for r in csv.DictReader(open(p)):
        L[(round(float(r["beta"]),3),round(float(r["kappa"]),3))]=float(r["Llink"])
    return L

def higgs_boundary(L):
    """kappa_c(beta) = midpoint of the kappa-interval with the steepest rise in <L_link>."""
    bs=sorted({b for b,k in L}); out=[]
    for b in bs:
        ks=sorted(k for bb,k in L if bb==b)
        vals=[L[(b,k)] for k in ks]
        if len(ks)<2: continue
        slopes=[(vals[i+1]-vals[i])/(ks[i+1]-ks[i]) for i in range(len(ks)-1)]
        i=int(np.argmax(slopes))
        out.append((b, 0.5*(ks[i]+ks[i+1])))
    return np.array(out)

def coulomb_boundary(cp, L):
    """beta_c(kappa) = beta of the chi_plaq peak (gauge specific heat) per kappa."""
    ks=sorted({k for b,k in L}); out=[]
    for k in ks:
        prof=[(b,cp.get((b,k),np.nan)) for b in sorted({bb for bb,kk in L if kk==k})]
        prof=[(b,v) for b,v in prof if not np.isnan(v)]
        if len(prof)<3: continue
        b,v=zip(*prof); v=np.array(v)
        if v.max()<=0: continue
        out.append((k, b[int(np.argmax(v))]))
    return np.array(out)

present=[(l,p,load(l)) for l,p in ORDER]; present=[(l,p,L) for l,p,L in present if L]
n=len(present)
fig,axes=plt.subplots(1,n,figsize=(4.6*n,4.6),squeeze=False)
for j,(lab,pretty,L) in enumerate(present):
    ax=axes[0][j]; cp=chi_plaq_map(lab)
    bs=np.array(sorted({b for b,k in L})); ks=np.array(sorted({k for b,k in L}))
    # faint <L_link> background
    M=np.full((len(ks),len(bs)),np.nan)
    for (b,k),v in L.items(): M[list(ks).index(k),list(bs).index(b)]=v
    ax.imshow(M,origin="lower",aspect="auto",cmap="Greys",alpha=0.35,
              extent=[bs.min()-.5,bs.max()+.5,ks.min()-.5,ks.max()+.5])
    hb=higgs_boundary(L); cb=coulomb_boundary(cp,L)
    if len(hb): ax.plot(hb[:,0],hb[:,1],"o-",color="crimson",lw=2,ms=5,label="Higgs κ_c(β)  [L_link rise]")
    if len(cb): ax.plot(cb[:,1],cb[:,0],"s--",color="navy",lw=2,ms=5,label="conf↔Coulomb β_c(κ)  [χ_plaq peak]")
    if lab in betaf: ax.axvline(betaf[lab],ls=":",color="gray",alpha=.7,label=f"lit. freeze β_f={betaf[lab]}")
    # region labels
    ax.text(0.6*bs.max(),0.85*ks.max(),"HIGGS",fontsize=11,ha="center",color="darkred",weight="bold")
    ax.text(0.18*bs.max(),0.18*ks.max(),"conf.",fontsize=9,ha="center",color="dimgray")
    ax.text(0.8*bs.max(),0.18*ks.max(),"Coulomb",fontsize=9,ha="center",color="navy")
    ax.set(title=f"{pretty}",xlabel="β",ylabel="κ",xlim=(bs.min()-.5,bs.max()+.5),ylim=(ks.min()-.5,ks.max()+.5))
    ax.legend(fontsize=7,loc="upper left"); ax.grid(alpha=0.2)
fig.suptitle(f"SU(2)→H phase diagram (L=2, locators only): Higgs κ_c(β) [matter] + confining↔Coulomb β_c(κ) [gauge]  — {D}")
fig.tight_layout(rect=[0,0,1,0.95])
out=os.path.join(D,"su2_phase_boundaries.png"); fig.savefig(out,dpi=120); print("wrote",out)
for lab,pretty,L in present:
    cp=chi_plaq_map(lab); hb=higgs_boundary(L); cb=coulomb_boundary(cp,L)
    print(f"\n{lab}: Higgs κ_c(β)= "+" ".join(f"{b:.0f}:{k:.1f}" for b,k in hb))
    print(f"      conf↔Coulomb β_c(κ)= "+" ".join(f"{k:.0f}:{b:.0f}" for k,b in cb))
