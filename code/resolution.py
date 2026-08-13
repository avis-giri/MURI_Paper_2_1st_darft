# SFCC resolution: minimum detectable change in liquid fraction f at 3-sigma,
# per layer, from the tier-1 index slope dI/df and the contact-probe noise.
import gzip, csv, math, re
import numpy as np
B="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP"
d=np.load(B+"/pilot_analysis/retention_mir.npz", allow_pickle=True)
y=d["y"]; uuids=d["uuids"]
w33,w1500,clay,oc=y[:,0],y[:,1],y[:,2],y[:,3]
bd={}
with gzip.open(B+"/ossl/ossl_soillab_L1_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r); I={c:i for i,c in enumerate(h)}
    need=set(uuids.tolist())
    for row in r:
        u=row[I["id.layer_uuid_txt"]]
        if u in need and row[I["bd_usda.a4_g.cm3"]].strip(): bd[u]=float(row[I["bd_usda.a4_g.cm3"]])
pat=re.compile(r"(turbel|orthel|histel|gelisol)", re.I)
taxa={}
with gzip.open(B+"/ossl/ossl_soilsite_L0_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r); I={c:i for i,c in enumerate(h)}
    for row in r:
        u=row[I["id.layer_uuid_txt"]]
        if u in set(uuids.tolist()): taxa[u]=row[I["pedon.taxa_usda_txt"]]
BD=np.array([bd.get(u,1.35) for u in uuids])
gel=np.array([bool(pat.search(taxa.get(u,""))) for u in uuids])
theta33=np.clip(w33/100.0*BD,0.01,0.65)
med=np.median(theta33[~gel&(oc<3)])
L=0.1*theta33/med
AW1450,AI1450,AW1500,AI1500=3.1491,1.8311,1.8833,4.5497
# dI/df is exactly 2L[(AI1500-AW1500)-(AI1450-AW1450)] (log-ratio is linear in f)
slope=2*L*((AI1500-AW1500)-(AI1450-AW1450))
sigma_I=math.sqrt(2)*0.01
df_min=3*sigma_I/np.abs(slope)
# translate to volumetric water resolution: dtheta = df * theta33
dth=df_min*theta33
for name,m in [("ALL",np.ones(len(L),bool)),("Gelisols",gel),("organic OC>=8",oc>=8),("mineral",oc<8)]:
    print("%-14s min detectable df: median %.3f (IQR %.3f-%.3f) | as volumetric water: %.3f cm3/cm3"
          % (name, np.median(df_min[m]), np.percentile(df_min[m],25), np.percentile(df_min[m],75),
             np.median(dth[m])))
print()
print("Reference-method context: TDR/dielectric unfrozen-water uncertainty is typically ~0.02-0.03 cm3/cm3.")
