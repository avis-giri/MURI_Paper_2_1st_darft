# Knowledge-graph / structural-hole analysis of the field, from OpenAlex.
# Design: a SENSING MODALITY x SOIL-STATE co-occurrence matrix over titles+
# abstracts. A genuine research gap appears as an empty cell whose row and
# column margins are both large (i.e. both topics are active, but never met).
import json, urllib.request, urllib.parse, time, re, itertools, os
from collections import Counter, defaultdict
MAILTO="aviskar.giri@taylorgeospatial.org"
OUT="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/kg"

def fetch(query, per=200, maxpages=6):
    """search OpenAlex works; return list of dicts with title+abstract text"""
    out=[]; cursor="*"
    for _ in range(maxpages):
        params={"search":query,"per-page":str(per),"cursor":cursor,"mailto":MAILTO,
                "select":"id,title,publication_year,cited_by_count,abstract_inverted_index,primary_location"}
        url="https://api.openalex.org/works?"+urllib.parse.urlencode(params)
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                d=json.load(r)
        except Exception as e:
            print("  ERR", query[:40], e); break
        for w in d.get("results",[]):
            ab=w.get("abstract_inverted_index")
            txt=""
            if ab:
                pos={}
                for term,idxs in ab.items():
                    for i in idxs: pos[i]=term
                txt=" ".join(pos[k] for k in sorted(pos))
            ven=""
            pl=w.get("primary_location") or {}
            src=pl.get("source") or {}
            ven=src.get("display_name") or ""
            out.append(dict(id=w["id"], title=w.get("title") or "", year=w.get("publication_year"),
                            cites=w.get("cited_by_count",0), abstract=txt, venue=ven))
        cursor=(d.get("meta") or {}).get("next_cursor")
        if not cursor: break
        time.sleep(0.3)
    return out

QUERIES=[
 "soil reflectance spectroscopy visible near infrared prediction",
 "soil mid infrared spectroscopy prediction properties library",
 "soil freezing characteristic curve unfrozen water content",
 "frozen soil dielectric unfrozen water microwave remote sensing",
 "permafrost remote sensing active layer thickness",
 "soil moisture optical reflectance radiative transfer model",
 "imaging spectroscopy soil properties mapping",
 "ice water absorption spectroscopy near infrared phase",
 "spectral induced polarization frozen ground permafrost",
 "soil water retention curve pedotransfer prediction",
]
corpus={}
for q in QUERIES:
    r=fetch(q)
    print("%-58s %5d works" % (q[:58], len(r)))
    for w in r: corpus[w["id"]]=w
print("unique works:", len(corpus))
json.dump(list(corpus.values()), open(OUT+"/corpus.json","w"))

# ---- lexicons ----
MODALITY={
 "optical/reflectance (VNIR-SWIR)": r"\b(reflectance|vis-?nir|visnir|visible.{0,12}near.?infrared|shortwave infrared|swir|hyperspectral|imaging spectroscop|spectroradiomet)\w*",
 "mid-infrared (MIR/FTIR)": r"\b(mid.?infrared|\bmir\b|ftir|diffuse reflectance infrared|drift)\w*",
 "dielectric/TDR/capacitance": r"\b(dielectric|permittivity|time.?domain reflectometry|\btdr\b|capacitance probe)\w*",
 "microwave/radar/SAR": r"\b(microwave|radar|\bsar\b|backscatter|brightness temperature|radiometer|l-?band|p-?band)\w*",
 "NMR/calorimetry/DSC": r"\b(nuclear magnetic resonance|\bnmr\b|calorimet|differential scanning|\bdsc\b)\w*",
 "electrical resistivity/SIP": r"\b(resistivity|induced polarization|electrical tomography|\bert\b)\w*",
}
STATE={
 "unfrozen water / freezing curve": r"\b(unfrozen water|freezing characteristic|soil freezing curve|liquid water content.{0,20}frozen|freezing point depression|premelt)\w*",
 "frozen ground / permafrost": r"\b(permafrost|frozen ground|frozen soil|active layer|freeze.?thaw)\w*",
 "soil moisture (unfrozen)": r"\b(soil moisture|water content|volumetric water)\w*",
 "water retention / matric potential": r"\b(water retention|matric potential|soil water characteristic|pedotransfer)\w*",
 "soil carbon / organic matter": r"\b(organic carbon|organic matter|\bsoc\b)\w*",
 "texture / mineralogy": r"\b(clay content|texture|mineralog|particle size)\w*",
}
MOD={k:re.compile(v,re.I) for k,v in MODALITY.items()}
STA={k:re.compile(v,re.I) for k,v in STATE.items()}

cell=defaultdict(list)
mrow=Counter(); scol=Counter()
for w in corpus.values():
    t=(w["title"] or "")+" "+(w["abstract"] or "")
    ms=[k for k,p in MOD.items() if p.search(t)]
    ss=[k for k,p in STA.items() if p.search(t)]
    for m in ms: mrow[m]+=1
    for s in ss: scol[s]+=1
    for m in ms:
        for s in ss: cell[(m,s)].append(w)

mods=list(MODALITY); stas=list(STATE)
print("\nCO-OCCURRENCE MATRIX (works mentioning both), corpus n=%d" % len(corpus))
w0=32
print(" "*w0 + "".join("%9s" % s.split()[0][:8] for s in stas))
for m in mods:
    print("%-32s" % m[:32] + "".join("%9d" % len(cell[(m,s)]) for s in stas))
print("\nrow margins:", {k:v for k,v in mrow.most_common()})
print("col margins:", {k:v for k,v in scol.most_common()})

# structural hole score: expected under independence vs observed
N=len(corpus)
print("\nSTRUCTURAL HOLES (observed vs expected under independence; both margins large):")
rows=[]
for m in mods:
    for s in stas:
        o=len(cell[(m,s)]); e=mrow[m]*scol[s]/N if N else 0
        if mrow[m]>=25 and scol[s]>=25:
            rows.append((o-e, o, e, m, s))
rows.sort()
for d,o,e,m,s in rows[:8]:
    print("  obs %3d  exp %6.1f  deficit %7.1f   %-30s x %s" % (o,e,d,m[:30],s))
print("\nTop over-represented (the crowded cells):")
for d,o,e,m,s in rows[-5:]:
    print("  obs %3d  exp %6.1f  surplus %+6.1f   %-30s x %s" % (o,e,d,m[:30],s))

# who occupies the target cell, if anyone
tgt=cell[("optical/reflectance (VNIR-SWIR)","unfrozen water / freezing curve")]
print("\nTARGET CELL: optical/reflectance x unfrozen water/freezing curve -> %d works" % len(tgt))
for w in sorted(tgt, key=lambda x:-x["cites"])[:12]:
    print("   %4s c=%-5d %-52s | %s" % (w["year"], w["cites"], (w["title"] or "")[:52], w["venue"][:34]))
json.dump({ "matrix": {f"{m} || {s}": len(cell[(m,s)]) for m in mods for s in stas},
            "row_margins": dict(mrow), "col_margins": dict(scol), "n": N},
          open(OUT+"/kg_matrix.json","w"), indent=1)
print("\nsaved corpus.json + kg_matrix.json")
