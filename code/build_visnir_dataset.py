# Paper-1 analog dataset at national scale: KSSL VisNIR spectra (350-2500 nm,
# the same range as the PSR+ used at PFRR) paired with the Paper-1-relatable
# property suite and pedon group keys. md5-gated inputs.
import gzip, csv, hashlib
import numpy as np
B="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP"
MAN={}
for line in open(B+"/manifests/ossl_md5_manifest.txt"):
    h,f=line.split(); MAN[f.split("/")[-1]]=h
def gate(p):
    assert MAN[p.split("/")[-1]]==hashlib.md5(open(p,"rb").read()).hexdigest(), "md5 fail "+p
    print("md5 OK", p.split("/")[-1])
for f in ["ossl_visnir_L0_v1.2.csv.gz","ossl_soillab_L1_v1.2.csv.gz","ossl_soilsite_L0_v1.2.csv.gz"]:
    gate(B+"/ossl/"+f)

site={}; taxa={}; depth={}
with gzip.open(B+"/ossl/ossl_soilsite_L0_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r); I={c:i for i,c in enumerate(h)}
    for row in r:
        u=row[I["id.layer_uuid_txt"]]
        site[u]=(row[I["dataset.code_ascii_txt"]], row[I["id.dataset.site_ascii_txt"]])
        taxa[u]=row[I["pedon.taxa_usda_txt"]]
        depth[u]=(row[I["layer.upper.depth_usda_cm"]], row[I["layer.lower.depth_usda_cm"]])

# Paper-1-relatable properties: pH (2 methods), EC, OC, clay/sand/silt, CaCO3,
# CEC, retention endpoints (moisture analogs), bulk density
PROPS=["ph.h2o_usda.a268_index","ph.cacl2_usda.a481_index","ec_usda.a364_ds.m",
       "oc_usda.c729_w.pct","clay.tot_usda.a334_w.pct","sand.tot_usda.c60_w.pct",
       "silt.tot_usda.c62_w.pct","caco3_usda.a54_w.pct","cec_usda.a723_cmolc.kg",
       "wr.33kPa_usda.a415_w.pct","wr.1500kPa_usda.a417_w.pct","bd_usda.a4_g.cm3"]
lab={}
with gzip.open(B+"/ossl/ossl_soillab_L1_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r); I={c:i for i,c in enumerate(h)}
    for row in r:
        u=row[I["id.layer_uuid_txt"]]
        lab[u]=[float(row[I[c]]) if row[I[c]].strip() else np.nan for c in PROPS]

X=[]; y=[]; g=[]; ds=[]; uu=[]; tx=[]; dep=[]
with gzip.open(B+"/ossl/ossl_visnir_L0_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r)
    ui=h.index("id.layer_uuid_txt")
    si=[i for i,c in enumerate(h) if c.startswith("scan_visnir.")]
    wls=[int(h[i].split(".")[1].split("_")[0]) for i in si]
    probe=h.index("scan_visnir.1426_ref")
    for row in r:
        if not row[probe].strip(): continue
        u=row[ui]
        if u not in lab or u not in site: continue
        try: spec=np.array([float(row[i]) for i in si], dtype=np.float32)
        except ValueError: continue
        X.append(spec); y.append(lab[u]); g.append(site[u][1]); ds.append(site[u][0])
        uu.append(u); tx.append(taxa.get(u,"")); dep.append(depth.get(u,("","")))
X=np.vstack(X); y=np.array(y,dtype=np.float32)
print("VisNIR paired dataset:", X.shape, "| wl %d-%d nm (%d bands)" % (min(wls),max(wls),len(wls)))
import collections
c=collections.Counter(ds)
print("by dataset:", dict(c))
np.savez_compressed(B+"/pilot_analysis/visnir_paper1_analog.npz",
    X=X, y=y, groups=np.array(g), dataset=np.array(ds), uuids=np.array(uu),
    taxa=np.array(tx), depths=np.array(dep), props=np.array(PROPS), wl=np.array(wls))
print("property coverage (non-NaN):")
for j,p in enumerate(PROPS):
    n=int(np.sum(~np.isnan(y[:,j])))
    print("   %-38s %6d" % (p,n))
print("saved visnir_paper1_analog.npz")
