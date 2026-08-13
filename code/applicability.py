# National applicability analysis: for each of the 17,612 KSSL layers with MIR
# + both retention endpoints, predict the tier-1 phase-partition index swing
# over the Clapeyron-bracketed freeze-out (0 to -1.23 C) and test it against
# measurement noise. The freeze-out is exactly f: 1 -> w1500/w33 (the liquid
# fraction remaining at the window edge), so no extrapolation is involved.
# Optical constants: Warren-Brandt 2008 (ice), Segelstein 1981 (liquid).
# Two-pass film model; equivalent film thickness L = k * theta_33 with k
# calibrated so the median mineral soil sits at L = 0.1 mm; k swept 2x both ways.
import gzip, csv, math, re
import numpy as np
B="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP"
d=np.load(B+"/pilot_analysis/retention_mir.npz", allow_pickle=True)
y=d["y"]; groups=d["groups"]; uuids=d["uuids"]
w33, w1500, clay, oc = y[:,0], y[:,1], y[:,2], y[:,3]

# BD + taxa lookups for these uuids
need=set(uuids.tolist())
bd={}; taxa={}
with gzip.open(B+"/ossl/ossl_soillab_L1_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r); I={c:i for i,c in enumerate(h)}
    for row in r:
        u=row[I["id.layer_uuid_txt"]]
        if u in need and row[I["bd_usda.a4_g.cm3"]].strip():
            bd[u]=float(row[I["bd_usda.a4_g.cm3"]])
pat=re.compile(r"(turbel|orthel|histel|gelisol)", re.I)
with gzip.open(B+"/ossl/ossl_soilsite_L0_v1.2.csv.gz","rt") as f:
    r=csv.reader(f); h=next(r); I={c:i for i,c in enumerate(h)}
    for row in r:
        u=row[I["id.layer_uuid_txt"]]
        if u in need: taxa[u]=row[I["pedon.taxa_usda_txt"]]
BD=np.array([bd.get(u,np.nan) for u in uuids])
gel=np.array([bool(pat.search(taxa.get(u,""))) for u in uuids])
print("layers: %d | with BD: %d | Gelisol layers: %d" % (len(uuids), np.isfinite(BD).sum(), gel.sum()))
BDf=np.where(np.isfinite(BD), BD, 1.35)   # class default where missing, counted above

# volumetric water at field capacity; liquid fraction remaining at -1.23 C
theta33=np.clip(w33/100.0*BDf, 0.01, 0.65)
f_end=np.clip(w1500/np.maximum(w33,1e-3), 0.02, 0.999)

# optical constants at the tier-1 pair (1/mm, from the verified tables)
AW1450, AI1450 = 3.1491, 1.8311
AW1500, AI1500 = 1.8833, 4.5497
RD=0.30
def index_swing(L):
    def R(nm_aw, nm_ai, f, LL): return RD*np.exp(-2*(f*nm_aw+(1-f)*nm_ai)*LL)
    I1=np.log(R(AW1500,AI1500,1.0,L)/R(AW1450,AI1450,1.0,L))
    I2=np.log(R(AW1500,AI1500,f_end,L)/R(AW1450,AI1450,f_end,L))
    # dark-floor check at the end state (reflectance must stay above 1.5%)
    Rmin=np.minimum(R(AW1500,AI1500,f_end,L), R(AW1450,AI1450,f_end,L))
    return I1-I2, Rmin

med_theta=np.median(theta33[~gel & (oc<3)])
for K in [0.2, 0.4, 0.8]:
    L=K*theta33/med_theta*0.1    # median mineral soil -> L=0.1mm at K=0.4... normalize
L0=0.1*theta33/med_theta         # base calibration
print("median mineral theta33 = %.3f -> L calibrated so that soil = 0.1 mm" % med_theta)
print()
print("Noise floors: contact probe CV 1%% -> sigma_I = %.3f; PFRR standoff CV 18%% -> sigma_I = %.3f"
      % (math.sqrt(2)*0.01, math.sqrt(2)*0.18))
print()
hdr="%-34s %8s %8s %8s %8s"
print(hdr % ("population (k=L multiplier)","median", "dI>3s", "dI>3s", "tier-2"))
print(hdr % ("", "dI", "contact", "standoff", "needed"))
for K in [0.5, 1.0, 2.0]:
    dI, Rmin = index_swing(K*L0)
    s_c, s_s = math.sqrt(2)*0.01, math.sqrt(2)*0.18
    dark = Rmin < 0.015
    det_c = (dI > 3*s_c) & ~dark
    det_s = (dI > 3*s_s) & ~dark
    for name, m in [("ALL layers (k=%.1f)"%K, np.ones(len(dI),bool)),
                    ("Gelisols (k=%.1f)"%K, gel),
                    ("organic OC>=8%% (k=%.1f)"%K, oc>=8),
                    ("mineral OC<8%% (k=%.1f)"%K, oc<8)]:
        mm=m & np.isfinite(dI)
        if mm.sum()==0: continue
        print("%-34s %8.3f %7.0f%% %7.0f%% %7.0f%%" % (name, np.median(dI[mm]),
            100*det_c[mm].mean(), 100*det_s[mm].mean(), 100*dark[mm].mean()))
    print()
# how much of the water actually freezes in the window
frz=1-f_end
print("Freeze-out magnitude inside the KSSL-bracketed window (1 - w1500/w33):")
for name,m in [("ALL",np.ones(len(frz),bool)),("Gelisols",gel),("organic OC>=8",oc>=8)]:
    print("   %-14s median %.0f%% of field-capacity water freezes by -1.23 C" % (name, 100*np.median(frz[m])))
