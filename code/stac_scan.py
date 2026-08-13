import json, urllib.request, urllib.parse
from collections import Counter
BBOX="-147.58,65.02,-147.28,65.22"
ES="https://earth-search.aws.element84.com/v1/search"
def search(coll,d0,d1):
    out=[]; page=1
    q={"collections":coll,"bbox":BBOX,"datetime":f"{d0}T00:00:00Z/{d1}T23:59:59Z","limit":100}
    url=ES+"?"+urllib.parse.urlencode(q)
    while url and page<=5:
        try:
            with urllib.request.urlopen(url,timeout=180) as r: d=json.load(r)
        except Exception as e:
            print("  ERR",coll,e); break
        for f in d.get("features",[]):
            p=f["properties"]
            out.append((p.get("datetime","")[:10], p.get("eo:cloud_cover"), f.get("id","")))
        nxt=[l for l in d.get("links",[]) if l.get("rel")=="next"]
        url=nxt[0]["href"] if nxt else None; page+=1
    return out
for coll in ["sentinel-2-l2a","landsat-c2-l2"]:
    print("==",coll)
    allr=[]
    for yr in range(2019,2026):
        r=search(coll,f"{yr}-09-01",f"{yr}-11-15")
        allr+=r
        clear=[x for x in r if x[1] is not None and x[1]<40]
        print("   %d  Sep1-Nov15: %3d scenes | %2d with cloud<40%%" % (yr,len(r),len(clear)))
    print("   by month:", dict(sorted(Counter(d[5:7] for d,_,_ in allr).items())))
    clear=[x for x in allr if x[1] is not None and x[1]<40]
    print("   TOTAL %d scenes, %d clear (<40%%)" % (len(allr),len(clear)))
    print("   clearest:")
    for d,c,i in sorted(clear,key=lambda x:x[1])[:8]:
        print("      %s  cloud %4.1f%%  %s" % (d,c,i[:46]))
    print()
# also: how late in the year does usable imagery persist?
print("LATE-SEASON limit check (Sentinel-2, cloud<40, by day-of-year):")
r=[]
for yr in range(2019,2026): r+=search("sentinel-2-l2a",f"{yr}-09-01",f"{yr}-12-31")
cl=[x for x in r if x[1] is not None and x[1]<40]
if cl:
    latest=sorted(cl,key=lambda x:x[0][5:])[-6:]
    for d,c,i in latest: print("   %s cloud %4.1f%%" % (d,c))
