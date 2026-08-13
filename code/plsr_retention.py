# Pre-registered Part-B primary: pedon-grouped 10-fold PLSR for water retention
# endpoints from KSSL MIR. Nested group 5-fold selects n_components from the
# fixed grid {5,10,15,20,25,30} inside training folds only (prereg C2).
# Also reports the OPTIMISTIC random row split once, labeled as such (B3),
# and the leakage magnitude (grouped minus random).
import numpy as np, json, platform, sklearn
from sklearn.model_selection import GroupKFold, KFold
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
SEED=20260812
B="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/pilot_analysis"
d=np.load(B+"/retention_mir.npz", allow_pickle=True)
X,y,groups=d["X"],d["y"],d["groups"]
# SNV per spectrum (pre-specified single pipeline)
X=(X-X.mean(axis=1,keepdims=True))/(X.std(axis=1,keepdims=True)+1e-9)
GRID=[5,10,15,20,25,30]
TARGETS=[("wr33",0),("wr1500",1)]
rng=np.random.RandomState(SEED)
# shuffle pedon order deterministically for fold assignment
order=rng.permutation(len(X))
X,y,groups=X[order],y[order],groups[order]

def run_cv(splitter, grouped, label):
    out={}
    for name,col in TARGETS:
        yy=y[:,col]
        r2s=[]; rmses=[]; chosen=[]
        splits=splitter.split(X,yy,groups) if grouped else splitter.split(X,yy)
        for tr,te in splits:
            # inner selection on training pedons only
            best=(None,-1e9)
            inner=GroupKFold(n_splits=5) if grouped else KFold(n_splits=5,shuffle=True,random_state=SEED)
            for nc in GRID:
                scores=[]
                isplits=inner.split(X[tr],yy[tr],groups[tr]) if grouped else inner.split(X[tr],yy[tr])
                for itr,ite in isplits:
                    m=PLSRegression(n_components=nc,scale=False)
                    m.fit(X[tr][itr],yy[tr][itr])
                    p=m.predict(X[tr][ite]).ravel()
                    ss=1-np.sum((p-yy[tr][ite])**2)/np.sum((yy[tr][ite]-yy[tr][ite].mean())**2)
                    scores.append(ss)
                sc=np.mean(scores)
                if sc>best[1]: best=(nc,sc)
            nc=best[0]; chosen.append(nc)
            m=PLSRegression(n_components=nc,scale=False)
            m.fit(X[tr],yy[tr])
            p=m.predict(X[te]).ravel()
            r2=1-np.sum((p-yy[te])**2)/np.sum((yy[te]-yy[te].mean())**2)
            r2s.append(float(r2)); rmses.append(float(np.sqrt(np.mean((p-yy[te])**2))))
        out[name]=dict(r2_mean=float(np.mean(r2s)), r2_sd=float(np.std(r2s)),
                       rmse_mean=float(np.mean(rmses)), folds=r2s, ncomp=chosen)
        print("%s %-7s R2 = %.2f +/- %.2f   RMSE = %.2f   ncomp=%s" %
              (label,name,np.mean(r2s),np.std(r2s),np.mean(rmses),chosen))
    return out

print("n=%d layers, %d pedons | sklearn %s" % (len(X), len(set(groups.tolist())), sklearn.__version__))
res={}
res["grouped"]=run_cv(GroupKFold(n_splits=10), True,  "PEDON-GROUPED")
res["random"] =run_cv(KFold(n_splits=10,shuffle=True,random_state=SEED), False, "RANDOM(optim.)")
for name,_ in TARGETS:
    gap=res["random"][name]["r2_mean"]-res["grouped"][name]["r2_mean"]
    print("LEAKAGE MAGNITUDE %-7s: random - grouped = %+.2f R2" % (name,gap))
json.dump(res, open(B+"/plsr_retention_results.json","w"), indent=1)
print("saved plsr_retention_results.json | python", platform.python_version())
