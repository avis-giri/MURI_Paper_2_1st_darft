# Figure: literature structure as a modality x soil-state matrix, plotted as
# the ratio of observed to expected co-occurrence. Values below one mark
# combinations the field has not made.
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
KG="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/kg"
FIG="/projects/benq/agiri1/KSSL_NEON_FLAGSHIP/figures"
d=json.load(open(KG+"/kg_matrix.json"))
N=d["n"]; rm=d["row_margins"]; cm=d["col_margins"]
mods=["optical/reflectance (VNIR-SWIR)","mid-infrared (MIR/FTIR)","dielectric/TDR/capacitance",
      "microwave/radar/SAR","NMR/calorimetry/DSC","electrical resistivity/SIP"]
stas=["unfrozen water / freezing curve","frozen ground / permafrost","soil moisture (unfrozen)",
      "water retention / matric potential","soil carbon / organic matter","texture / mineralogy"]
short_m=["Optical VNIR-SWIR","Mid-infrared","Dielectric / TDR","Microwave / SAR","NMR / calorimetry","Resistivity / SIP"]
short_s=["Unfrozen water\n& freezing curve","Frozen ground\n& permafrost","Soil moisture\n(unfrozen)",
         "Water retention\n& potential","Soil carbon\n& organic matter","Texture\n& mineralogy"]
O=np.zeros((len(mods),len(stas))); E=np.zeros_like(O)
for i,m in enumerate(mods):
    for j,s in enumerate(stas):
        O[i,j]=d["matrix"][f"{m} || {s}"]
        E[i,j]=rm[m]*cm[s]/N
R=np.log2((O+0.5)/(E+0.5))
plt.rcParams.update({"font.family":"serif","font.size":9})
fig,ax=plt.subplots(figsize=(7.0,4.3))
norm=TwoSlopeNorm(vmin=min(R.min(),-4), vcenter=0, vmax=max(R.max(),2))
im=ax.imshow(R,cmap="RdBu_r",norm=norm,aspect="auto")
for i in range(len(mods)):
    for j in range(len(stas)):
        ax.text(j,i,"%d\n(%.0f)"%(O[i,j],E[i,j]),ha="center",va="center",fontsize=7.5,
                color="white" if abs(R[i,j])>2.2 else "black")
# highlight the target cell
ax.add_patch(plt.Rectangle((-0.5,-0.5),1,1,fill=False,edgecolor="#111111",lw=2.4))
ax.add_patch(plt.Rectangle((-0.5,0.5),1,1,fill=False,edgecolor="#111111",lw=2.4,ls=":"))
ax.set_xticks(range(len(stas))); ax.set_xticklabels(short_s,fontsize=7.5)
ax.set_yticks(range(len(mods))); ax.set_yticklabels(short_m,fontsize=8)
ax.set_xlabel("soil state or property",fontsize=9)
ax.set_ylabel("sensing modality",fontsize=9)
cb=fig.colorbar(im,ax=ax,fraction=0.035,pad=0.02)
cb.set_label("log$_2$(observed / expected)",fontsize=8)
cb.ax.tick_params(labelsize=7)
ax.set_title("Co-occurrence of sensing modality and soil state in %d works\n"
             "cell shows observed count (expected under independence)"%N,fontsize=9,pad=8)
fig.tight_layout()
fig.savefig(FIG+"/fig_kg_matrix.pdf",bbox_inches="tight")
fig.savefig(FIG+"/fig_kg_matrix.png",dpi=200,bbox_inches="tight")
print("wrote fig_kg_matrix.pdf/.png")
print("target cell observed %d expected %.1f  log2ratio %.2f" % (O[0,0],E[0,0],R[0,0]))
