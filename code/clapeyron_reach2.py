import csv, statistics as st
from collections import defaultdict
def num(x):
    try: return float(x)
    except: return None

# Clapeyron: dP/dT = rho_w * L_f / T_0
rho, Lf, T0 = 1000.0, 334000.0, 273.15
coef = rho*Lf/T0/1e6
print("Clapeyron coefficient = %.3f MPa/K" % coef)
for p,lab in [(0.033,"33 kPa (field capacity)"),(1.5,"1500 kPa (wilting point, KSSL coldest)"),
              (6.11,"6.11 MPa"),(12.2,"12.2 MPa")]:
    print("   %-38s <-> %+.2f C" % (lab, -p/coef))
CL = -1.5/coef
print()

meta={}
with open("SFCCMetadata.csv", encoding="utf-8", errors="replace") as f:
    for i,row in enumerate(csv.DictReader(f), start=1): meta[str(i)]=row
curves=defaultdict(list)
with open("SFCCData.csv", encoding="utf-8", errors="replace") as f:
    for row in csv.DictReader(f):
        t=num(row["temperature"]); v=num(row["volumetric_water_content"])
        if t is not None and v is not None and t<=0.05: curves[row["index"]].append((t,v))

fracs=[]; org=[]; fine=[]
for k,v in curves.items():
    v=sorted(set(v), key=lambda p:-p[0])
    if len(v)<8: continue
    if v[0][0] < -0.2: continue          # must actually start at the freezing point
    if min(t for t,_ in v) > -3.0: continue
    w0=v[0][1]; wend=v[-1][1]; total=w0-wend
    if total<=1e-9: continue
    below=[p for p in v if p[0]<=CL]
    if not below: continue
    fr=100*(w0-below[0][1])/total
    if not (0 <= fr <= 100): continue
    fracs.append(fr)
    r=meta.get(k,{})
    o=num(r.get("organic")); s=num(r.get("silt")); c=num(r.get("clay"))
    if o is not None and o>0: org.append(fr)
    if (s or 0)>=40 or (c or 0)>=30: fine.append(fr)

fracs.sort()
print("Curves passing strict QC (>=8 pts, start warmer than -0.2 C, reach colder than -3 C): %d" % len(fracs))
print()
print("Percent of the TOTAL unfrozen-water decline that occurs between 0 and %.2f C" % CL)
print("(i.e. the part KSSL's two retention points already bracket):")
print("   median %.0f%%   mean %.0f%%   IQR %.0f-%.0f%%" % (
    st.median(fracs), st.mean(fracs), fracs[len(fracs)//4], fracs[3*len(fracs)//4]))
for thr in (50,70,80,90):
    n=sum(1 for f in fracs if f>thr)
    print("   >%d%% inside: %3d of %d (%.0f%%)" % (thr,n,len(fracs),100*n/len(fracs)))
if org: print("\n   organic soils (organic>0), n=%d: median %.0f%%" % (len(org), st.median(org)))
if fine: print("   fine-textured (silt>=40 or clay>=30), n=%d: median %.0f%%" % (len(fine), st.median(fine)))
