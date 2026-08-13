# Build the pre-registered Part-B dataset: KSSL layers with MIR spectra + both
# retention endpoints (+ clay, OC), with pedon group keys. md5-gated inputs.
import gzip, csv, hashlib, sys
import numpy as np
B="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP"
MAN={}
for line in open(B+"/manifests/ossl_md5_manifest.txt"):
    h,f=line.split(); MAN[f.split("/")[-1]]=h
def gate(path):
    h=hashlib.md5(open(path,"rb").read()).hexdigest()
    name=path.split("/")[-1]
    assert MAN[name]==h, f"md5 mismatch for {name}"
    print("md5 OK", name)
for f in ["ossl_mir_L0_v1.2.csv.gz","ossl_soillab_L1_v1.2.csv.gz","ossl_soilsite_L0_v1.2.csv.gz"]:
    gate(B+"/ossl/"+f)

# soilsite: uuid -> (dataset, pedon site key)
site={}
with gzip.open(B+"/ossl/ossl_soilsite_L0_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r); I={c:i for i,c in enumerate(h)}
    for row in r:
        if row[I["dataset.code_ascii_txt"]]=="KSSL.SSL":
            site[row[I["id.layer_uuid_txt"]]]=row[I["id.dataset.site_ascii_txt"]]
# soillab: retention + clay + oc
lab={}
with gzip.open(B+"/ossl/ossl_soillab_L1_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r); I={c:i for i,c in enumerate(h)}
    cols=["wr.33kPa_usda.a415_w.pct","wr.1500kPa_usda.a417_w.pct","clay.tot_usda.a334_w.pct","oc_usda.c729_w.pct"]
    for row in r:
        u=row[I["id.layer_uuid_txt"]]
        if u not in site: continue
        vals=[]
        ok=True
        for c in cols[:2]:
            v=row[I[c]].strip()
            if not v: ok=False; break
            vals.append(float(v))
        if not ok: continue
        for c in cols[2:]:
            v=row[I[c]].strip()
            vals.append(float(v) if v else np.nan)
        lab[u]=vals
print("layers with both retention endpoints:", len(lab))
# MIR spectra for those layers
X=[]; y=[]; g=[]; uu=[]
with gzip.open(B+"/ossl/ossl_mir_L0_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r)
    ui=h.index("id.layer_uuid_txt")
    si=[i for i,c in enumerate(h) if c.startswith("scan_mir.")]
    probe=h.index("scan_mir.2300_abs")
    for row in r:
        u=row[ui]
        if u not in lab: continue
        if not row[probe].strip(): continue
        try:
            spec=np.array([float(row[i]) for i in si], dtype=np.float32)
        except ValueError:
            continue
        X.append(spec); y.append(lab[u]); g.append(site[u]); uu.append(u)
X=np.vstack(X); y=np.array(y,dtype=np.float32)
groups=np.array(g); uuids=np.array(uu)
print("dataset:", X.shape, "| pedons:", len(set(g)))
np.savez_compressed(B+"/pilot_analysis/retention_mir.npz", X=X, y=y, groups=groups, uuids=uuids)
print("saved retention_mir.npz")
