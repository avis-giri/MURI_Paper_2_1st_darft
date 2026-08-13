# Download the shoulder-season 30-minute soil temperature files and test
# whether the residence-time result from the single monitoring station holds
# at the three NEON Alaska sites.
import csv, urllib.request, os, io, statistics as st
from collections import defaultdict
import numpy as np
OUT="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/neon"
TOK=open("/u/agiri1/.secrets/neon_token").read().strip()
rows=list(csv.DictReader(open(OUT+"/neon_file_inventory.csv")))
temp=[r for r in rows if r["product"]=="DP1.00041.001"]
print("soil temperature files to pull:", len(temp))
os.makedirs(OUT+"/raw", exist_ok=True)
vals=defaultdict(list)
for i,r in enumerate(temp):
    fp=os.path.join(OUT,"raw",r["name"])
    if not os.path.exists(fp):
        try:
            req=urllib.request.Request(r["url"], headers={"X-API-Token":TOK})
            with urllib.request.urlopen(req,timeout=180) as h, open(fp,"wb") as o: o.write(h.read())
        except Exception as e:
            print("  fail",r["name"][:50],e); continue
    with open(fp, newline="") as fh:
        for rec in csv.DictReader(fh):
            v=rec.get("soilTempMean","")
            q=rec.get("finalQF","0")
            if v and q in ("0","0.0"):
                try: vals[r["site"]].append(float(v))
                except: pass
    if (i+1)%40==0: print("  ...",i+1,"files")
print()
print("NEON 30-minute soil temperature, shoulder season (Sep-Nov), QF=0 only")
print("(all sensor depths pooled; shallowest NEON depth is about 2 cm)")
for site in ["BONA","DEJU","HEAL"]:
    a=np.array(vals[site])
    if a.size==0: print("  %s: no data" % site); continue
    sub=a[a<0]
    print("  %-5s n=%7d  min %6.1f C  sub-zero %5.1f%%" % (site,a.size,a.min(),100*sub.size/a.size), end="")
    if sub.size:
        print("  | of sub-zero time: %4.1f%% warmer than -1.23C, %4.1f%% warmer than -5C"
              % (100*(sub>-1.227).sum()/sub.size, 100*(sub>-5).sum()/sub.size))
    else: print()
np.savez_compressed(OUT+"/neon_soiltemp_shoulder.npz", **{k:np.array(v) for k,v in vals.items()})
print("saved neon_soiltemp_shoulder.npz")
