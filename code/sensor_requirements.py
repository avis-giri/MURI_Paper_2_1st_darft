# Instrument-requirements analysis: which sensors can measure the phase-
# partition index, and which index tier survives the atmosphere from orbit.
# Index tiers (from the forward model, absorption coefficients in 1/mm):
#   tier 1  1450 / 1500   liquid 3.15 / 1.88   ice 1.83 / 4.55   STRONGEST
#   tier 2  1200 / 1250   liquid 0.126/0.111   ice 0.070/0.130   weak, clean
#   tier 3  1940 / 2020   liquid 12.4 / 5.42   ice 6.90 / 9.74   saturates
import numpy as np
PAIRS={"tier 1 (1450/1500)":(1450,1500),"tier 2 (1200/1250)":(1200,1250),"tier 3 (1940/2020)":(1940,2020)}
# Atmospheric water-vapour absorption regions (opaque from orbit) and windows.
OPAQUE=[(1340,1465),(1800,1960)]        # strong H2O bands
WINDOWS=[(400,1330),(1500,1780),(2000,2350)]
def atm(nm):
    for a,b in OPAQUE:
        if a<=nm<=b: return "OPAQUE"
    for a,b in WINDOWS:
        if a<=nm<=b: return "window"
    return "edge"
print("ATMOSPHERIC VIABILITY OF EACH INDEX TIER (spaceborne/airborne above the water column)")
for name,(a,b) in PAIRS.items():
    print("   %-20s %4d nm = %-7s | %4d nm = %-7s -> %s" % (
        name,a,atm(a),b,atm(b),
        "USABLE from orbit" if atm(a)!="OPAQUE" and atm(b)!="OPAQUE" else "LAB / CONTACT ONLY"))
print()
# Sensor inventory: (name, ranges as list of (lo,hi), sampling nm, platform)
SENSORS=[
 ("Sentinel-2 MSI",[(443,443),(490,490),(560,560),(665,665),(704,704),(740,740),(783,783),(842,842),(865,865),(945,945),(1374,1374),(1610,1610),(2190,2190)],None,"satellite, multispectral"),
 ("Landsat 8/9 OLI",[(440,440),(480,480),(560,560),(650,650),(870,870),(1610,1610),(2200,2200)],None,"satellite, multispectral"),
 ("EnMAP",[(420,2450)],10,"satellite, imaging spectrometer"),
 ("PRISMA",[(400,2500)],12,"satellite, imaging spectrometer"),
 ("EMIT",[(380,2500)],7.4,"ISS, +/-52 deg latitude only"),
 ("Tanager-1",[(400,2500)],5,"satellite, imaging spectrometer"),
 ("DESIS",[(400,1000)],2.55,"ISS, VNIR only"),
 ("AVIRIS-NG",[(380,2510)],5,"airborne"),
 ("NEON AOP",[(380,2510)],5,"airborne"),
 ("HySpex VNIR+SWIR (indoor)",[(400,1000),(960,2500)],[1.2,1.7],"laboratory benchtop"),
 ("PSR+ 3500 (field/contact)",[(350,2500)],[2.8,8.0],"field spectroradiometer"),
]
def covers(ranges, nm):
    for lo,hi in ranges:
        if lo<=nm<=hi: return True
    return False
def multispectral(ranges):
    return all(lo==hi for lo,hi in ranges)
print("%-28s %-32s %-9s %-9s %-9s" % ("sensor","platform","tier1","tier2","tier3"))
for name,rng,samp,plat in SENSORS:
    cells=[]
    for _,(a,b) in PAIRS.items():
        if multispectral(rng):
            # a discrete band must sit within ~15 nm of the required centre
            ok=all(any(abs(lo-t)<=15 for lo,hi in rng) for t in (a,b))
        else:
            ok=covers(rng,a) and covers(rng,b)
        cells.append("yes" if ok else "NO")
    print("%-28s %-32s %-9s %-9s %-9s" % (name,plat,*cells))
print()
print("Sentinel-2 gap check: no band between 945 and 1374 nm, none between 1374 and 1610 nm.")
print("Landsat OLI gap check: no band between 870 and 1610 nm.")
print()
print("CONCLUSION")
print(" - The workhorse multispectral satellites cannot sample any index pair.")
print("   The measurement requires imaging spectroscopy. This is a plausible reason")
print("   the gap persisted: the sensors most people use are blind to it.")
print(" - Tier 1 is the strongest signal but 1450 nm sits inside a water-vapour")
print("   absorption band, so it is a laboratory and contact measurement.")
print(" - Tier 2 is weaker but both bands fall in a clean atmospheric window,")
print("   making it the spaceborne-viable form of the index.")
