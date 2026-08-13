# Tier-2 index for saturated organic soils: I2 = ln[R(1250)/R(1200)].
# Physics: alpha(1/mm) liquid 0.126 vs ice 0.070 at 1200; liquid 0.111 vs ice
# 0.130 at 1250 -> ratio flips, same construction as 1450/1500 but ~20x weaker
# absorption, so it stays unsaturated at mm-scale films in wet organics.
import csv, math, statistics as st
from collections import defaultdict
src=open(__file__.replace("index_tier2","index_test")).read()
exec(src.split("def index_of")[0])
def index2(s, idx):
    r1200=st.mean(s[idx[w]] for w in range(1195,1206))
    r1250=st.mean(s[idx[w]] for w in range(1245,1256))
    return math.log(r1250/r1200)
def qc2(s, idx):
    for w in range(1150,1311):
        v=s[idx[w]]
        if v is None or v<=0.05 or v>100: return False
    return True

ins=load(BASE+"/insitu_pffr_psr_data.csv")
bands,idx=band_cols(ins)
groups=defaultdict(list)
for r in ins:
    d=str(r.get("Depth (cm)","")).strip(); soil=str(r.get("Soil Type","")).strip().lower()
    loc=str(r.get("Sample Location","")).strip()
    if "calibration" in loc.lower(): continue
    s=[num(r[c]) for c in bands]
    if not qc2(s,idx): continue
    frozen = d.endswith("+") or "permafrost" in soil or "frozen" in soil
    if not frozen: groups[(loc,d)].append(index2(s,idx))
unf=[st.median(v) for v in groups.values()]

old=load(BASE+"/Old_Samples_Go.csv")
TC=[c for c in old[0] if c.strip().lower().startswith("temp")][0]
tmap=defaultdict(list)
for r in old:
    t=num(r[TC]); dd=num(r["Depth (cm)"]); l=r["Sample Location"].strip()
    if t is not None and dd is not None: tmap[(l,dd)].append(t)
leg=load(BASE+"/go_old_samples_spectrally_merged.csv")
lb,lidx=band_cols(leg)
legg=defaultdict(list)
for r in leg:
    l=r["Sample Location"].strip(); dd=num(r["Depth (cm)"])
    if (l,dd) not in tmap: continue
    s=[num(r[c]) for c in lb]
    if not qc2(s,lidx): continue
    legg[(l,dd,round(st.mean(tmap[(l,dd)]),1))].append(index2(s,lidx))

print("Tier-2 index I2 = ln[R(1250)/R(1200)]  (weak bands, for saturated organics)")
print("Unfrozen in-situ groups: n=%d  median I2 = %+.4f  IQR %+.4f..%+.4f" % (
    len(unf), st.median(unf), sorted(unf)[len(unf)//4], sorted(unf)[3*len(unf)//4]))
print("Legacy core groups (QC on 1150-1310 only, so more survive):")
for (l,dd,t),v in sorted(legg.items(), key=lambda kv: kv[0][2]):
    lab=""
    if l=="PF-Q2" and dd==0: lab="  <-- saturated peat, tier-1 index FAILED here"
    print("   %6.1f C  I2 = %+.4f   %-6s depth %3.0f cm (n=%d)%s" % (t, st.median(v), l, dd, len(v), lab))
cold2=[st.median(v) for (l,dd,t),v in legg.items() if t<=-2]
if cold2:
    import random
    random.seed(20260812)
    obs=st.median(unf)-st.median(cold2)
    comb=unf+cold2; nx=len(unf); cnt=0
    for _ in range(100000):
        random.shuffle(comb)
        if st.median(comb[:nx])-st.median(comb[nx:])>=obs: cnt+=1
    print("Unfrozen vs cold cores on I2: diff=%+.4f, one-sided p=%.5f" % (obs,(cnt+1)/100001))
