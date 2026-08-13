# Group-level hypothesis tests for the phase-partition index I = ln[R(1500)/R(1450)].
# H1: unfrozen in-situ groups have HIGHER I than cold (T <= -2 C) legacy core groups.
# H2 (confound demonstration): field-exposed frozen faces do NOT differ from unfrozen.
# Unit of independence = site/depth group. Permutation tests, fixed seed.
import csv, math, statistics as st, random
from collections import defaultdict
src=open(__file__.replace("index_test2","index_test")).read()
exec(src.split("# ---------- report")[0])

unf=[st.median(v) for (l,d,f),v in groups.items() if f=="unfrozen"]
fro=[st.median(v) for (l,d,f),v in groups.items() if f=="frozen"]
cold=[(leg_T[k],st.median(v)) for k,v in leg_groups.items() if leg_T[k]<=-2.0]
warm=[(leg_T[k],st.median(v)) for k,v in leg_groups.items() if leg_T[k]>-2.0]

def perm_test(x, y, B=100000, seed=20260812):
    random.seed(seed)
    obs=st.median(x)-st.median(y)
    comb=list(x)+list(y); nx=len(x); cnt=0
    for _ in range(B):
        random.shuffle(comb)
        if st.median(comb[:nx])-st.median(comb[nx:]) >= obs: cnt+=1
    return obs,(cnt+1)/(B+1)

cx=[m for _,m in cold]
print("Groups: unfrozen in-situ n=%d | legacy cores T<=-2C n=%d | legacy T>-2C n=%d | exposed faces n=%d" %
     (len(unf), len(cx), len(warm), len(fro)))
print()
obs,p=perm_test(unf,cx)
print("H1  unfrozen (%+.3f) vs frozen cores (%+.3f): diff=%+.3f, one-sided p=%.5f" % (st.median(unf),st.median(cx),obs,p))
obs2,p2=perm_test(fro,cx)
print("H1b exposed faces (%+.3f) vs frozen cores: diff=%+.3f, one-sided p=%.5f" % (st.median(fro),obs2,p2))
obs3,p3=perm_test(unf,fro)
print("H2  unfrozen vs exposed faces: diff=%+.3f, one-sided p=%.3f (expected null: surface melt)" % (obs3,p3))
print()
print("Direction inside legacy set: warm (T>-2C) median I=%+.3f  vs  cold median I=%+.3f" %
     (st.median([m for _,m in warm]), st.median(cx)))
print()
print("Every group, sorted by I:")
allg=[("unfrozen-insitu",m) for m in unf]+[("exposed-face",m) for m in fro]+\
     [("core %+.1fC"%t,m) for t,m in cold+warm]
for lab,m in sorted(allg,key=lambda x:x[1]):
    print("   %+.4f  %s" % (m,lab))
