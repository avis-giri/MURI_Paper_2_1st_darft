# THE LADDER DESIGN: derive the laboratory temperature ladder from the
# temperature distribution the field actually occupies during freeze-up.
# Rationale: a ladder with equal steps wastes measurements on states the soil
# rarely occupies. Weighting by residence time makes every laboratory point
# carry field-relevant information.
import csv, statistics as st
from collections import defaultdict
import numpy as np
P="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/multitemporal"
rows=list(csv.DictReader(open(P+"/snotel_pfrr.csv")))
def C(f): return (float(f)-32)*5/9
sto=defaultdict(list)
for r in rows:
    if r["element"]!="STO": continue
    if r["unit"]!="degF": continue
    try: t=C(r["value"])
    except: continue
    sto[(r["triplet"],r["depth_in"])].append((r["date"],t))
print("Soil temperature series (converted to C):")
for k,v in sorted(sto.items()):
    ts=[t for _,t in v]
    print("   %-14s depth %4s in  n=%5d  min %6.1f  max %6.1f" % (k[0],k[1],len(v),min(ts),max(ts)))

# Residence-time distribution in the SUB-ZERO band, shallow sensor (-2 in = 5 cm)
KEY=("947:AK:SNTL","-2")
series=sorted(sto[KEY])
ts=np.array([t for _,t in series])
sub=ts[(ts<0.0)&(ts>-25)]
print()
print("Little Chena Ridge, 5 cm depth: %d daily values, %d of them sub-zero (%.0f%% of the year)"
      % (len(ts), len(sub), 100*len(sub)/len(ts)))
edges=[0,-0.5,-1,-2,-3,-5,-7,-10,-15,-25]
print()
print("RESIDENCE TIME by temperature band (days per year, 5 cm):")
nyears=len(ts)/365.25
cum=0
for a,b in zip(edges,edges[1:]):
    n=int(((ts<=a)&(ts>b)).sum()); cum+=n
    print("   %6.1f to %6.1f C : %5d days total = %5.1f d/yr  (%4.1f%% of sub-zero time)"
          % (a,b,n,n/nyears,100*n/max(len(sub),1)))
print()
# The Clapeyron-bracketed window vs what the field occupies
inwin=int(((ts<0)&(ts>=-1.227)).sum())
print("Fraction of sub-zero time inside the KSSL-bracketed window (0 to -1.23 C): %.0f%%" % (100*inwin/len(sub)))
print("Fraction of sub-zero time warmer than -5 C: %.0f%%" % (100*int(((ts<0)&(ts>=-5)).sum())/len(sub)))
print()
# Proposed ladder: equal-residence-time quantiles of the sub-zero distribution
K=12
qs=np.quantile(sub,np.linspace(0,1,K+1))
print("PROPOSED LADDER (%d equal-residence-time steps, 5 cm freeze-up climatology):" % K)
print("   " + "  ".join("%.2f" % q for q in qs))
print()
# Freeze-up timing: first autumn crossing below 0 C, per year
byyear=defaultdict(list)
for d,t in series: byyear[d[:4]].append((d,t))
print("FREEZE-UP DATE (first autumn day the 5 cm soil goes below 0 C):")
fu=[]
for y in sorted(byyear):
    v=[(d,t) for d,t in sorted(byyear[y]) if d[5:7] in ("09","10","11","12")]
    first=next((d for d,t in v if t<0), None)
    if first: fu.append(first); print("   %s  %s" % (y,first))
if fu:
    doys=[int(d[5:7])*30+int(d[8:10]) for d in fu]
    print("   -> median freeze-up around day-of-year band %d (Sept=270, Oct=300)" % int(np.median(doys)))
# Soil moisture at freeze-up = the total water that must partition
sms=defaultdict(list)
for r in rows:
    if r["element"]=="SMS" and r["triplet"]==KEY[0] and r["depth_in"]==KEY[1]:
        try: sms[r["date"]]=float(r["value"])
        except: pass
vals=[sms[d] for d,_ in series if d in sms and sms[d]>0]
if vals:
    print()
    print("Volumetric soil moisture at 5 cm (%%): median %.1f, IQR %.1f-%.1f  -> the water available to freeze"
          % (np.median(vals), np.percentile(vals,25), np.percentile(vals,75)))
