# Phase-partition index test on real PFRR spectra.
# Index: I = ln[R(1500)/R(1450)]  (albedo-invariant; physics predicts
# I more negative = more ice-shifted, from Warren-Brandt / Segelstein constants).
# Groups: (a) in-situ UNFROZEN horizons (Aug 2025, soil T +1.2 to +25.6 C);
#         (b) in-situ FROZEN faces (depth classes 45+/60+/permafrost, thaw front);
#         (c) legacy frozen cores at logged sub-zero temperatures.
# Unit of independence = site/depth group (pre-registration Part A3 discipline).
# QC: reject any spectrum with a value <=0.05% or >100% inside 1330-1650 nm.
import csv, math, statistics as st
from collections import defaultdict
BASE="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/pilot_analysis/data"
def num(x):
    try: return float(str(x).strip())
    except: return None
def load(p):
    with open(p, encoding="utf-8", errors="replace") as fh: return list(csv.DictReader(fh))

def band_cols(rows):
    bands=[c for c in rows[0] if c.startswith("W_")]
    wl=[int(c[2:]) for c in bands]
    return bands, {w:i for i,w in enumerate(wl)}

def qc_ok(s, idx, lo=1330, hi=1650):
    for w in range(lo,hi+1):
        v=s[idx[w]]
        if v is None or v<=0.05 or v>100: return False
    return True

def index_of(s, idx):
    # mean over +/-5 nm to reduce single-channel noise
    r1450=st.mean(s[idx[w]] for w in range(1445,1456))
    r1500=st.mean(s[idx[w]] for w in range(1495,1506))
    return math.log(r1500/r1450)

# ---------- in-situ ----------
ins=load(BASE+"/insitu_pffr_psr_data.csv")
bands,idx=band_cols(ins)
groups=defaultdict(list)
for r in ins:
    d=str(r.get("Depth (cm)","")).strip()
    soil=str(r.get("Soil Type","")).strip().lower()
    loc=str(r.get("Sample Location","")).strip()
    if "calibration" in loc.lower(): continue
    s=[num(r[c]) for c in bands]
    if not qc_ok(s,idx): continue
    frozen = d.endswith("+") or "permafrost" in soil or "frozen" in soil
    groups[(loc,d,"frozen" if frozen else "unfrozen")].append(index_of(s,idx))

# ---------- legacy cores ----------
old=load(BASE+"/Old_Samples_Go.csv")
TC=[c for c in old[0] if c.strip().lower().startswith("temp")][0]
tmap=defaultdict(list)
for r in old:
    t=num(r[TC]); d=num(r["Depth (cm)"]); l=r["Sample Location"].strip()
    if t is not None and d is not None: tmap[(l,d)].append(t)
leg=load(BASE+"/go_old_samples_spectrally_merged.csv")
lb,lidx=band_cols(leg)
leg_groups=defaultdict(list); leg_T={}
for r in leg:
    l=r["Sample Location"].strip(); d=num(r["Depth (cm)"])
    ts=tmap.get((l,d))
    if not ts: continue
    s=[num(r[c]) for c in lb]
    if not qc_ok(s,lidx): continue
    leg_groups[(l,d)].append(index_of(s,lidx))
    leg_T[(l,d)]=st.mean(ts)

# ---------- report at GROUP level ----------
unf=[st.median(v) for (l,d,f),v in groups.items() if f=="unfrozen"]
fro=[st.median(v) for (l,d,f),v in groups.items() if f=="frozen"]
legm=[(leg_T[k], st.median(v), k, len(v)) for k,v in leg_groups.items()]
nu=sum(len(v) for (l,d,f),v in groups.items() if f=="unfrozen")
nf=sum(len(v) for (l,d,f),v in groups.items() if f=="frozen")

print("Phase-partition index I = ln[R(1500)/R(1450)] on REAL PFRR spectra")
print("(physics prediction: liquid-dominated I > ice-dominated I)\n")
print("(a) In-situ UNFROZEN horizons : %2d groups (%3d spectra)  median I = %+.4f  IQR %+.4f..%+.4f"
      % (len(unf), nu, st.median(unf), sorted(unf)[len(unf)//4], sorted(unf)[3*len(unf)//4]))
print("(b) In-situ FROZEN faces      : %2d groups (%3d spectra)  median I = %+.4f  IQR %+.4f..%+.4f"
      % (len(fro), nf, st.median(fro), sorted(fro)[len(fro)//4] if len(fro)>=4 else min(fro), sorted(fro)[3*len(fro)//4] if len(fro)>=4 else max(fro)))
print("(c) Legacy cores (sub-zero)   : %2d groups (%3d spectra)  median I = %+.4f"
      % (len(legm), sum(n for _,_,_,n in legm), st.median([m for _,m,_,_ in legm])))

# exact rank test (group level, a vs b)
import itertools
def mann_whitney_exact(x,y):
    nx,ny=len(x),len(y)
    allv=[(v,0) for v in x]+[(v,1) for v in y]
    allv.sort()
    # U statistic
    ranks={}
    u=0
    for xi in x:
        u+=sum(1 for yi in y if yi<xi)+0.5*sum(1 for yi in y if yi==xi)
    # permutation p (exact if small, else monte carlo with fixed seed)
    import random
    random.seed(20260812)
    comb=[v for v in x+y]
    obs=st.median(x)-st.median(y)
    B=20000; cnt=0
    for _ in range(B):
        random.shuffle(comb)
        px=comb[:nx]; py=comb[nx:]
        if st.median(px)-st.median(py) >= obs: cnt+=1
    return u/(nx*ny), (cnt+1)/(B+1)

if len(fro)>=3:
    auc,p=mann_whitney_exact(unf,fro)
    print("\nGroup-level test, unfrozen vs frozen faces: rank AUC = %.2f, one-sided permutation p = %.4f (20000 perms, seed fixed)" % (auc,p))

print("\nLegacy core groups, coldest to warmest:")
for T,m,k,n in sorted(legm):
    print("   %6.1f C  I = %+.4f   %s depth %.0f cm  (n=%d)" % (T,m,k[0],k[1],n))
# correlation with temperature across legacy groups
if len(legm)>=5:
    xs=[t for t,_,_,_ in legm]; ys=[m for _,m,_,_ in legm]
    mx,my=st.mean(xs),st.mean(ys)
    r=sum((a-mx)*(b-my) for a,b in zip(xs,ys))/ (math.sqrt(sum((a-mx)**2 for a in xs))*math.sqrt(sum((b-my)**2 for b in ys)))
    print("\nAcross legacy groups: Pearson r(T, I) = %+.2f (n=%d groups; one T per specimen, confounded with lithology)" % (r,len(legm)))
