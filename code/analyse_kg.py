# Second pass: (a) is the gap closing over time? (b) which properties does
# optical spectroscopy already own vs not? (c) which venues publish the
# adjacent crowded cells -> venue targeting evidence.
import json, re
from collections import Counter, defaultdict
OUT="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/kg"
corpus=json.load(open(OUT+"/corpus.json"))
def txt(w): return ((w["title"] or "")+" "+(w["abstract"] or ""))
OPT=re.compile(r"\b(reflectance|vis-?nir|visnir|shortwave infrared|swir|hyperspectral|imaging spectroscop|spectroradiomet)\w*",re.I)
MIR=re.compile(r"\b(mid.?infrared|\bmir\b|ftir|diffuse reflectance infrared)\w*",re.I)
UFW=re.compile(r"\b(unfrozen water|freezing characteristic|soil freezing curve|freezing point depression|premelt)\w*",re.I)
FRZ=re.compile(r"\b(permafrost|frozen ground|frozen soil|active layer|freeze.?thaw)\w*",re.I)
DIE=re.compile(r"\b(dielectric|permittivity|time.?domain reflectometry|\btdr\b)\w*",re.I)

# (a) temporal
print("(a) TEMPORAL: works per 5-year bin")
bins=defaultdict(Counter)
for w in corpus:
    y=w.get("year") or 0
    if y<1990: continue
    b=(y//5)*5
    t=txt(w)
    if OPT.search(t) or MIR.search(t):
        if UFW.search(t): bins[b]["optical x unfrozen-water"]+=1
        if FRZ.search(t): bins[b]["optical x frozen-ground"]+=1
    if DIE.search(t) and UFW.search(t): bins[b]["dielectric x unfrozen-water"]+=1
    if UFW.search(t): bins[b]["unfrozen-water (any method)"]+=1
keys=["unfrozen-water (any method)","dielectric x unfrozen-water","optical x frozen-ground","optical x unfrozen-water"]
print("%-6s"%"bin" + "".join("%28s"%k for k in keys))
for b in sorted(bins):
    print("%-6d"%b + "".join("%28d"%bins[b][k] for k in keys))

# (b) property coverage by optical
print("\n(b) WHICH PROPERTIES OPTICAL SPECTROSCOPY ALREADY OWNS (works with optical+property)")
PROPS={
 "organic carbon": r"\b(organic carbon|organic matter|\bsoc\b)",
 "clay / texture": r"\b(clay content|texture|particle size)",
 "pH": r"\bph\b",
 "electrical conductivity/salinity": r"\b(electrical conductivity|salinity|\bec\b)",
 "CEC": r"\b(cation exchange|\bcec\b)",
 "moisture content": r"\b(soil moisture|water content)",
 "bulk density": r"\bbulk density\b",
 "water retention": r"\b(water retention|matric potential|field capacity|wilting point)",
 "carbonate": r"\b(carbonate|caco3)",
 "ice content": r"\bice content\b",
 "unfrozen water": r"\bunfrozen water\b",
 "thermal conductivity": r"\bthermal conductivity\b",
 "Atterberg / plasticity": r"\b(atterberg|liquid limit|plastic limit|plasticity index)",
}
rows=[]
for name,pat in PROPS.items():
    p=re.compile(pat,re.I)
    n_any=sum(1 for w in corpus if p.search(txt(w)))
    n_opt=sum(1 for w in corpus if p.search(txt(w)) and (OPT.search(txt(w)) or MIR.search(txt(w))))
    rows.append((n_opt/max(n_any,1), n_opt, n_any, name))
rows.sort(reverse=True)
print("%-34s %7s %7s %8s" % ("property","optical","any","share"))
for share,no,na,name in rows:
    print("%-34s %7d %7d %7.0f%%" % (name,no,na,100*share))
print("\n-> low share + high any = the properties optical has NOT claimed = where to lead")

# (c) venues for the adjacent cells
print("\n(c) VENUES publishing optical x frozen-ground and dielectric x unfrozen-water")
for lab,cond in [("optical x frozen-ground", lambda t:(OPT.search(t) or MIR.search(t)) and FRZ.search(t)),
                 ("unfrozen-water (any)", lambda t: UFW.search(t))]:
    c=Counter(w["venue"] for w in corpus if cond(txt(w)) and w["venue"])
    print("  %s:" % lab)
    for v,n in c.most_common(8): print("     %3d  %s" % (n,v[:62]))
# IEEE presence
ieee=[w for w in corpus if "IEEE" in (w["venue"] or "")]
print("\nIEEE works in corpus: %d" % len(ieee))
c=Counter(w["venue"] for w in ieee)
for v,n in c.most_common(6): print("   %4d  %s" % (n,v[:62]))
iu=[w for w in ieee if UFW.search(txt(w))]
print("IEEE works touching unfrozen water/freezing curve: %d" % len(iu))
for w in sorted(iu,key=lambda x:-x["cites"])[:8]:
    print("   %4s c=%-4d %-56s" % (w["year"],w["cites"],(w["title"] or "")[:56]))
