# Signal-size estimate for the optical freezing curve.
# Sources: ice = Warren & Brandt 2008 compilation (atmos.uw.edu);
#          liquid water = Segelstein 1981 thesis table.
# Model: MARMIT-style two-pass film attenuation R = R_dry * exp(-2*alpha_eff*L),
# alpha_eff = f*alpha_liquid + (1-f)*alpha_ice, f = liquid fraction of pore water
# along the photon path, L = equivalent water film thickness (MARMIT fits:
# ~0.05-0.5 mm for wet soils).
import math

def load_ice(p):
    d={}
    for line in open(p):
        t=line.split()
        if len(t)>=3:
            try: wl=float(t[0]); k=float(t[2])
            except: continue
            if 0.3<=wl<=3.0: d[round(wl*1000)]=k
    return d
def load_water(p):
    d={}
    for line in open(p):
        t=line.split()
        if len(t)>=3:
            try: wl=float(t[0]); k=float(t[2])
            except: continue
            if 0.3<=wl<=3.0: d[round(wl*1000)]=k
    return d
ice=load_ice("ice_wb2008.dat"); wat=load_water("water_segelstein.txt")
def interp(d, nm):
    ks=sorted(d)
    lo=max(k for k in ks if k<=nm); hi=min(k for k in ks if k>=nm)
    if lo==hi: return d[lo]
    return d[lo]+(d[hi]-d[lo])*(nm-lo)/(hi-lo)
def alpha(d, nm):   # absorption coefficient in 1/mm
    k=interp(d,nm)
    return 4*math.pi*k/(nm*1e-6)   # lambda in mm

print("Absorption coefficients (1/mm):")
print(f"  {'nm':>5s} {'liquid':>10s} {'ice':>10s} {'ice/liq':>8s}")
for nm in [970,1030,1200,1250,1300,1400,1450,1500,1600,1850,1940,2020,2100,2250]:
    aw,ai=alpha(wat,nm),alpha(ice,nm)
    print(f"  {nm:5d} {aw:10.4f} {ai:10.4f} {ai/aw:8.2f}")

R_DRY=0.30
def refl(nm,f,L):
    a=f*alpha(wat,nm)+(1-f)*alpha(ice,nm)
    return R_DRY*math.exp(-2*a*L)

print("\nScenario: silt, fixed total pore water on path, freeze-out from f=1.00 (all liquid,")
print("-0.1 C) to f=0.25 (-5 C typical for silt). Two-pass film model, R_dry=0.30.")
for L in [0.05,0.1,0.2,0.5]:
    print(f"\n  L={L} mm equivalent water film:")
    print(f"  {'nm':>5s} {'R(f=1)':>8s} {'R(f=.25)':>9s} {'dR abs':>8s} {'dR rel%':>8s}")
    for nm in [970,1030,1450,1500,1940,2020]:
        r1,r2=refl(nm,1.0,L),refl(nm,0.25,L)
        print(f"  {nm:5d} {r1:8.4f} {r2:9.4f} {r2-r1:8.4f} {100*(r2-r1)/r1:8.1f}")

print("\nPhase-partition INDEX (dual band-depth ratio), robust to albedo and L:")
print("  ratio = ln[R(1030)/R(970)] type differential absorption.")
for L in [0.05,0.1,0.2,0.5]:
    row=[]
    for f in [1.0,0.75,0.5,0.25,0.0]:
        idx=math.log(refl(1030,f,L)/refl(970,f,L))
        idx2=math.log(refl(1500,f,L)/refl(1450,f,L))
        row.append((f,idx,idx2))
    print(f"  L={L} mm: " + "  ".join(f"f={f:.2f}: I970/1030={a:+.4f} I1450/1500={b:+.4f}" for f,a,b in row[::2]))
