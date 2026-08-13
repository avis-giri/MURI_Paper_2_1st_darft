# Soil temperature and moisture time series near Poker Flat, from NRCS SNOTEL.
# Purpose: measure the DISTRIBUTION of near-surface soil temperature the field
# actually occupies during freeze-up, so the laboratory temperature ladder is
# sampled from reality rather than chosen arbitrarily.
import json, urllib.request, urllib.parse, csv, os
OUT="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/multitemporal"
BASE="https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/data"
STATIONS={"947:AK:SNTL":"Little Chena Ridge (32 km from PFRR)",
          "1090:AK:SNTL":"Upper Nome Creek (48 km)",
          "948:AK:SNTL":"Mt. Ryan (61 km)",
          "1302:AK:SNTL":"Creamers Field (32 km, Fairbanks)"}
ELEMS=["STO:-2","STO:-8","STO:-20","SMS:-2","SMS:-8","SMS:-20","TOBS","SNWD"]
rows=[]
for trip,name in STATIONS.items():
    for el in ELEMS:
        q={"stationTriplets":trip,"elements":el,"duration":"DAILY",
           "beginDate":"2003-10-01","endDate":"2026-08-01"}
        url=BASE+"?"+urllib.parse.urlencode(q)
        try:
            with urllib.request.urlopen(url,timeout=300) as r: d=json.load(r)
        except Exception as e:
            print("ERR",trip,el,e); continue
        n=0
        for rec in d:
            for de in rec.get("data",[]):
                se=de.get("stationElement",{})
                unit=se.get("storedUnitCode")
                for v in de.get("values",[]):
                    if v.get("value") is None: continue
                    rows.append((trip,name,se.get("elementCode"),se.get("heightDepth"),unit,v["date"],v["value"]))
                    n+=1
        print("%-14s %-8s %6d values" % (trip,el,n))
with open(OUT+"/snotel_pfrr.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["triplet","name","element","depth_in","unit","date","value"])
    w.writerows(rows)
print("total rows:",len(rows))
