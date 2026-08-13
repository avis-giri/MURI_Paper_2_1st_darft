# Fair comparison with the SNOTEL station: restrict to the SHALLOWEST NEON
# sensor level (verticalPosition 501, nominally about 2 cm) and aggregate to
# DAILY MEANS, matching the daily 5 cm station record.
import csv, os, re, statistics as st
from collections import defaultdict
import numpy as np
OUT="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/neon"
files=os.listdir(OUT+"/raw")
pat=re.compile(r"\.([A-Z]{4})\.DP1\.00041\.001\.(\d{3})\.(\d{3})\.030\.ST_30_minute\.(\d{4}-\d{2})\.")
daily=defaultdict(lambda: defaultdict(list))   # site -> date -> values
nfiles=0
for fn in files:
    m=pat.search(fn)
    if not m: continue
    site,hpos,vpos,month=m.groups()
    if vpos!="501": continue          # shallowest level only
    nfiles+=1
    with open(os.path.join(OUT,"raw",fn), newline="") as fh:
        for rec in csv.DictReader(fh):
            v=rec.get("soilTempMean",""); q=rec.get("finalQF","0")
            if not v or q not in ("0","0.0"): continue
            d=rec["startDateTime"].strip("\"")[:10]
            try: daily[site][d].append(float(v))
            except: pass
print("shallowest-level files used:", nfiles)
print()
print("NEON Alaska sites, shallowest sensor (~2 cm), DAILY MEANS, Sep-Nov, QF=0")
print("%-6s %8s %8s %10s %14s %14s" % ("site","n days","min C","sub-zero","of sub-zero:","" ))
print("%-6s %8s %8s %10s %14s %14s" % ("","","","% of days",">-1.23 C",">-5 C"))
res={}
for site in ["BONA","DEJU","HEAL"]:
    dd=daily[site]
    if not dd: print("  %s no data"%site); continue
    a=np.array([st.mean(v) for v in dd.values()])
    sub=a[a<0]
    r=(a.size,a.min(),100*sub.size/a.size,
       100*(sub>-1.227).sum()/sub.size if sub.size else float("nan"),
       100*(sub>-5).sum()/sub.size if sub.size else float("nan"))
    res[site]=r
    print("%-6s %8d %8.1f %9.1f%% %13.1f%% %13.1f%%" % (site,*r))
print()
print("Reference, SNOTEL Little Chena Ridge, 5 cm, daily, full year: min -5.6 C,")
print("70%% of sub-zero time warmer than -1.23 C, 100%% warmer than -5 C.")
print()
print("NOTE: the NEON window here is Sep-Nov only (freeze-up), the station figure")
print("is full-year. Deep winter is excluded from the NEON numbers.")
