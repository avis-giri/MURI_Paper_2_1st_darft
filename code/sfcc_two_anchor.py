# Two-anchor SFCC reconstruction test.
# QUESTION: if the only constraints are the two KSSL retention pressures mapped
# through Clapeyron (33 kPa -> -0.027 C, 1500 kPa -> -1.227 C), how well is the
# full measured freezing curve recovered, and how large is the extrapolation
# penalty below -1.23 C?
# Method: for each QC-passing measured curve, (i) read the anchor values by
# monotone interpolation, (ii) solve a van Genuchten-form SFCC exactly through
# the two anchors (theta_r fixed at 0, m = 1-1/n), (iii) score RMSE against the
# measured points in 0..-1.23 (interpolation zone) and -1.23..-10 C
# (extrapolation zone), normalized by theta at the first anchor.
# Benchmark: free 3-parameter least-squares fit to ALL points (best case).
# Unit = curve; report medians and IQRs, stratified organic vs mineral.
import csv, math
import numpy as np
from collections import defaultdict
from scipy.optimize import brentq, least_squares

B="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/anchors/sfcc"
COEF=1.223e3   # kPa per K
T1,T2=-0.027,-1.227
def num(x):
    try: return float(x)
    except: return None
meta={}
with open(B+"/SFCCMetadata.csv", encoding="utf-8", errors="replace") as f:
    for i,row in enumerate(csv.DictReader(f), start=1): meta[str(i)]=row
curves=defaultdict(list)
with open(B+"/SFCCData.csv", encoding="utf-8", errors="replace") as f:
    for row in csv.DictReader(f):
        t=num(row["temperature"]); v=num(row["volumetric_water_content"])
        if t is not None and v is not None and t<=0.05: curves[row["index"]].append((t,v))

def vg(T, ws, alpha, n):
    psi=COEF*np.abs(T)
    m=1-1/n
    return ws*(1+(alpha*psi)**n)**(-m)

def interp_at(pts, T):
    """monotone linear interpolation of w at temperature T (pts warm->cold)"""
    for (ta,wa),(tb,wb) in zip(pts,pts[1:]):
        if tb<=T<=ta:
            if ta==tb: return wa
            return wa+(wb-wa)*(T-ta)/(tb-ta)
    return None

res=[]
for k,v in curves.items():
    pts=sorted(set(v), key=lambda p:-p[0])
    if len(pts)<8: continue
    if pts[0][0] < -0.2: continue
    if min(t for t,_ in pts) > -3.0: continue
    w1=interp_at(pts,T1); w2=interp_at(pts,T2)
    if w1 is None or w2 is None or w1<=0 or w2<=0 or w2>=w1: continue
    # exact two-anchor solve: ratio r = w2/w1 determines n given alpha... two
    # unknowns (alpha,n), two equations. Solve nested: for given n, alpha from
    # anchor 2; require anchor 1 satisfied.
    def anchor_gap(n):
        m=1-1/n
        # from anchor2: (1+(a*psi2)^n)^-m = w2/w1s  where w1s = ws; take ws from
        # anchor 1 relation simultaneously -> reduce: define x_i=(a*psi_i)^n
        # w_i = ws (1+x_i)^-m ; ratio: (1+x1)^-m / (1+x2)^-m = w1/w2
        # and x2/x1 = (psi2/psi1)^n. One unknown x1.
        R=(psi2/psi1)**n
        target=w1/w2
        def g(x1): return ((1+x1*R)/(1+x1))**(1-1/n) - target
        try:
            x1=brentq(g, 1e-12, 1e12, maxiter=200)
        except ValueError:
            return None
        ws=w1*(1+x1)**(1-1/n)
        a=x1**(1/n)/psi1
        return ws,a
    psi1,psi2=COEF*abs(T1),COEF*abs(T2)
    # choose n by matching curvature is impossible with 2 anchors: n is NOT
    # identified. Honest approach: sweep plausible n in [1.1, 3.0] and report
    # the BEST and WORST reconstruction, bracketing the identifiability gap;
    # plus a fixed literature default n=1.6.
    Ts=np.array([t for t,_ in pts]); Ws=np.array([w for _,w in pts])
    zi=(Ts<=T1)&(Ts>=T2); ze=(Ts<T2)&(Ts>=-10.0)
    if zi.sum()<2 or ze.sum()<2: continue
    def rmse_for(n):
        sol=anchor_gap(n)
        if sol is None: return None
        ws,a=sol
        pred=vg(Ts,ws,a,n)
        ri=float(np.sqrt(np.mean((pred[zi]-Ws[zi])**2)))/w1
        re=float(np.sqrt(np.mean((pred[ze]-Ws[ze])**2)))/w1
        return ri,re
    grid=[rmse_for(n) for n in np.arange(1.15,3.01,0.05)]
    grid=[g for g in grid if g]
    if not grid: continue
    fixed=rmse_for(1.6)
    # free fit benchmark on all points (3 params)
    def resid(p):
        ws,la,ln=p
        return vg(Ts,ws,math.exp(la),1+math.exp(ln))-Ws
    try:
        fit=least_squares(resid,[w1,math.log(1e-3),math.log(0.6)],max_nfev=2000)
        pred=vg(Ts,fit.x[0],math.exp(fit.x[1]),1+math.exp(fit.x[2]))
        free_i=float(np.sqrt(np.mean((pred[zi]-Ws[zi])**2)))/w1
        free_e=float(np.sqrt(np.mean((pred[ze]-Ws[ze])**2)))/w1
    except Exception:
        free_i=free_e=None
    r=meta.get(k,{})
    org=num(r.get("organic"))
    res.append(dict(k=k, org=(org or 0)>0,
        best_i=min(g[0] for g in grid), best_e=min(g[1] for g in grid),
        worst_e=max(g[1] for g in grid),
        fix_i=fixed[0] if fixed else None, fix_e=fixed[1] if fixed else None,
        free_i=free_i, free_e=free_e))

import statistics as st
def med(xs): 
    xs=[x for x in xs if x is not None]
    return (st.median(xs), np.percentile(xs,25), np.percentile(xs,75), len(xs)) if xs else (None,)*4
print("Curves scored: %d (organic: %d)" % (len(res), sum(1 for r in res if r["org"])))
print()
print("Normalized RMSE (fraction of water content at the -0.03 C anchor):")
for lab,key_i,key_e in [("two-anchor, n fixed 1.6","fix_i","fix_e"),
                        ("two-anchor, best n in sweep","best_i","best_e"),
                        ("free 3-par fit to ALL points","free_i","free_e")]:
    mi=med([r[key_i] for r in res]); me=med([r[key_e] for r in res])
    print("  %-30s  0..-1.23C: median %.3f (IQR %.3f-%.3f)   -1.23..-10C: median %.3f (IQR %.3f-%.3f)"
          % (lab, mi[0],mi[1],mi[2], me[0],me[1],me[2]))
we=med([r["worst_e"] for r in res])
print("  n-identifiability bracket, extrapolation zone worst case: median %.3f (IQR %.3f-%.3f)" % (we[0],we[1],we[2]))
print()
for sub,lab in [([r for r in res if r["org"]],"organic"),([r for r in res if not r["org"]],"mineral/unknown")]:
    if len(sub)>=5:
        me=med([r["fix_e"] for r in sub])
        print("  %-16s n=%3d  extrapolation nRMSE (n=1.6): median %.3f (IQR %.3f-%.3f)" % (lab,len(sub),me[0],me[1],me[2]))
