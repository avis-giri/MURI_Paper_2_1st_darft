# Figure: three independent constraints converge on the same temperature window.
import csv
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
M="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/multitemporal"
FIG="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/figures"
rows=list(csv.DictReader(open(M+"/snotel_pfrr.csv")))
def C(f): return (float(f)-32)*5/9
ts=[C(r["value"]) for r in rows if r["element"]=="STO" and r["triplet"]=="947:AK:SNTL"
    and r["depth_in"]=="-2" and r["unit"]=="degF"]
ts=np.array(ts); sub=ts[(ts<0)&(ts>-25)]
AW1450,AI1450,AW1500,AI1500=3.1491,1.8311,1.8833,4.5497
L=0.1; RD=0.30
f=np.linspace(0,1,200)
I=np.log((RD*np.exp(-2*(f*AW1500+(1-f)*AI1500)*L))/(RD*np.exp(-2*(f*AW1450+(1-f)*AI1450)*L)))
plt.rcParams.update({"font.family":"serif","font.size":9})
fig,axes=plt.subplots(1,3,figsize=(7.2,2.5))
ax=axes[0]
ax.hist(sub,bins=np.arange(-6,0.25,0.25),color="#3b6ea5",edgecolor="white",linewidth=0.3)
ax.axvspan(-1.227,0,color="#c8a45c",alpha=0.35,zorder=0)
ax.set_xlabel("soil temperature at 5 cm (C)"); ax.set_ylabel("days")
ax.set_title("(a) where the ground actually is",fontsize=8.5)
ax.set_xlim(-6,0)
ax=axes[1]
ax.plot(f,I,color="#8c3a2b",lw=1.8)
ax.axhline(0,color="#888888",lw=0.6)
ax.set_xlabel("liquid fraction of pore water"); ax.set_ylabel(r"index  ln[$R_{1500}/R_{1450}$]")
ax.set_title("(b) what the physics predicts",fontsize=8.5)
ax=axes[2]
labels=["0 to -1.23 C\n(retention bracket)","-1.23 to -5 C","colder than -5 C"]
field=[100*np.mean((ts<0)&(ts>=-1.227))/np.mean(ts<0),
       100*np.mean((ts<-1.227)&(ts>=-5))/np.mean(ts<0),
       100*np.mean(ts<-5)/np.mean(ts<0)]
ax.barh(range(3),field,color=["#c8a45c","#3b6ea5","#bbbbbb"])
ax.set_yticks(range(3)); ax.set_yticklabels(labels,fontsize=7)
ax.set_xlabel("percent of sub-zero time")
ax.set_title("(c) overlap of the three",fontsize=8.5)
ax.invert_yaxis()
for i,v in enumerate(field): ax.text(v+1.5,i,"%.0f%%"%v,va="center",fontsize=7.5)
ax.set_xlim(0,100)
fig.tight_layout()
fig.savefig(FIG+"/fig_convergence.pdf",bbox_inches="tight")
fig.savefig(FIG+"/fig_convergence.png",dpi=200,bbox_inches="tight")
print("wrote fig_convergence.pdf/.png")
print("field percentages:", [round(x,1) for x in field])
