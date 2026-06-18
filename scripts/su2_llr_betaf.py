#!/usr/bin/env python3
"""Maxwell (equal-area) construction on an LLR a1(e) curve -> first-order beta_f.

Input: a su2_llr output file (columns: e0 A0 a1 acc_window acc_hmc resid_e avg_plaq L_link).
a1(A0) = d ln rho / dA  (the microcanonical beta as a function of the extensive gauge energy A).
The canonical weight is ln P_beta(A) = ln rho(A) - beta*A,  ln rho(A) = cumulative integral of a1 dA.
At a 1st-order transition ln P_beta is double-peaked; beta_f = the beta where the two peaks have EQUAL
height (equal weight) -- equivalently the Maxwell equal-area level of the back-bending a1(A) curve. We
locate it as the beta at which the GLOBAL maximum of ln P jumps between the ordered (low-A) and the
disordered (high-A) basin. Latent heat = the A-gap between the two coexisting peaks (-> 0 at kappa*).

  usage: su2_llr_betaf.py <llr_out_file> [--plot out.png]
"""
import sys, re, numpy as np

def load(fn, resid_max=0.04):
    """columns: e0 A0 a1 acc_hmc resid_e meanE avg_plaq L_link [side]. Keep only well-sampled cells
    (resid_e <= resid_max); on duplicate e (two-sided overlap) keep the smaller-resid one."""
    best = {}                                    # e0 -> (resid, A, a1)
    for ln in open(fn):
        if ln.startswith('#') or not ln.strip(): continue
        p = ln.split()
        if len(p) < 3: continue
        try: ee, AA, aa = float(p[0]), float(p[1]), float(p[2])
        except ValueError: continue
        if 'DRIVE_FAIL' in ln: continue
        resid = float(p[4]) if len(p) > 4 else 0.0
        if resid > resid_max: continue           # poorly-sampled cell (failed to confine)
        key = round(ee, 5)
        if key not in best or resid < best[key][0]: best[key] = (resid, AA, aa)
    e = np.array(sorted(best))
    A = np.array([best[k][1] for k in e]); a1 = np.array([best[k][2] for k in e])
    return e, A, a1

def beta_f(A, a1):
    """ln rho via trapezoid cumulative integral of a1 over A; beta_f = global-max-basin crossover."""
    lnrho = np.concatenate([[0.0], np.cumsum(0.5*(a1[1:]+a1[:-1])*np.diff(A))])
    nplaq_guess = None
    betas = np.linspace(float(min(a1)), float(max(a1)), 4000)
    Amid = 0.5*(A[0]+A[-1])
    locs = np.array([A[np.argmax(lnrho - b*A)] for b in betas])
    side = locs > Amid                             # True = disordered basin wins
    # beta_f = where the winning basin flips from disordered (low beta) to ordered (high beta)
    flip = np.where(side[:-1] & ~side[1:])[0]
    if len(flip) == 0:
        return None, lnrho, None, None
    bi = flip[0]
    bf = 0.5*(betas[bi]+betas[bi+1])
    # coexisting A's at beta_f: the two peaks of lnP
    lnP = lnrho - bf*A
    lo = A[:np.searchsorted(A, Amid)]; hi = A[np.searchsorted(A, Amid):]
    A1 = A[np.argmax((lnrho-bf*A)[:len(lo)])] if len(lo) else None
    A2 = A[len(lo)+np.argmax((lnrho-bf*A)[len(lo):])] if len(hi) else None
    return bf, lnrho, A1, A2

def main():
    if len(sys.argv) < 2: print(__doc__); return
    fn = sys.argv[1]
    e, A, a1 = load(fn)
    if len(A) < 4: print(f"only {len(A)} usable cells; need >=4"); return
    # recover kappa + n_plaq from the header for reporting
    txt = open(fn).read()
    mk = re.search(r'kappa=([\d.]+)', txt); kap = float(mk.group(1)) if mk else float('nan')
    mn = re.search(r'n_plaq=([\d.]+)', txt); nplaq = float(mn.group(1)) if mn else 1.0
    bf, lnrho, A1, A2 = beta_f(A, a1)
    if bf is None:
        print(f"kappa={kap}: NO equal-weight crossover found (a1(e) monotonic -> crossover, not 1st-order). "
              f"a1 range [{a1.min():.3f},{a1.max():.3f}]")
        return
    lat = (A2-A1)/nplaq if (A1 is not None and A2 is not None) else float('nan')
    print(f"kappa={kap:.3f}:  beta_f = {bf:.4f}   (latent heat Delta e = {lat:.4f}; "
          f"coexist e1={A1/nplaq:.3f}, e2={A2/nplaq:.3f})")
    if '--plot' in sys.argv:
        out = sys.argv[sys.argv.index('--plot')+1]
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        ax[0].plot(A/nplaq, a1, '-o', ms=4); ax[0].axhline(bf, ls='--', color='r', label=f'beta_f={bf:.3f}')
        ax[0].set_xlabel('e = A/n_plaq = 1-avg_plaq'); ax[0].set_ylabel('a1 = d ln rho/dA (microcanonical beta)')
        ax[0].set_title(f'LLR a1(e), kappa={kap:.2f} (Maxwell level = beta_f)'); ax[0].legend(); ax[0].grid(alpha=.3)
        lnP = lnrho - bf*A
        ax[1].plot(A/nplaq, lnP-lnP.max(), '-o', ms=4); ax[1].set_xlabel('e'); ax[1].set_ylabel('ln P_betaf(A) (shifted)')
        ax[1].set_title('canonical ln P at beta_f (equal-height double peak)'); ax[1].grid(alpha=.3)
        if A1: ax[1].axvline(A1/nplaq, ls=':', color='g');
        if A2: ax[1].axvline(A2/nplaq, ls=':', color='purple')
        plt.tight_layout(); plt.savefig(out, dpi=120); print(f"wrote {out}")

if __name__ == '__main__': main()
