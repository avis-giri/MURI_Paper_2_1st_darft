# NEON soil temperature and water content at the three Alaska sites, to test
# whether the freeze-up residence-time result generalises beyond one station.
import json, urllib.request, os, csv, time
TOK=open("/u/agiri1/.secrets/neon_token").read().strip()
OUT="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/neon"
API="https://data.neonscience.org/api/v0"
def get(u):
    r=urllib.request.Request(u, headers={"X-API-Token":TOK})
    with urllib.request.urlopen(r, timeout=120) as h: return json.load(h)
SITES=["BONA","DEJU","HEAL"]
PRODUCTS={"DP1.00041.001":"soil temperature","DP1.00094.001":"soil water content"}
inventory=[]
for pid,label in PRODUCTS.items():
    d=get(f"{API}/products/{pid}")["data"]
    for s in d.get("siteCodes",[]):
        if s["siteCode"] in SITES:
            months=sorted(s.get("availableMonths",[]))
            inventory.append((pid,label,s["siteCode"],len(months),months[0],months[-1]))
            print("%-16s %-20s %-5s %3d months  %s to %s" % (pid,label,s["siteCode"],len(months),months[0],months[-1]))
# Pull a freeze-up month per site per year and record file inventory (not bulk download yet)
rows=[]
for pid,label in PRODUCTS.items():
    for site in SITES:
        for yr in range(2018,2026):
            for mo in ("09","10","11"):
                try:
                    d=get(f"{API}/data/{pid}/{site}/{yr}-{mo}")["data"]
                except Exception:
                    continue
                fs=[f for f in d.get("files",[]) if f["name"].endswith(".csv")]
                # 30-minute aggregation files
                pick=[f for f in fs if "30_minute" in f["name"] or "30min" in f["name"]]
                for f in pick[:2]:
                    rows.append((pid,label,site,f"{yr}-{mo}",f["name"],f["size"],f["url"]))
                time.sleep(0.1)
print("\ncandidate 30-minute files found:", len(rows))
with open(OUT+"/neon_file_inventory.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["product","label","site","month","name","size","url"]); w.writerows(rows)
tot=sum(r[5] for r in rows)/1e6
print("total size if all pulled: %.0f MB" % tot)
for r in rows[:6]: print("   ",r[2],r[3],r[4][:66],"%.1f MB"%(r[5]/1e6))
