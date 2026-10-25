"""Redraw Figures 1 and 2 from deposited-value lesion-level scores (no new analysis)."""
import pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
SRC="/workspace/gibson_research_group_2026-10-11/10_Psoriasis_Dual_Pathway_ExpDermatol/analysis/outputs_reference/out_deposited/lesion_level_scores.tsv"
d=pd.read_csv(SRC,sep="\t")
order=["GSE13355","GSE30999","GSE14905","GSE117468","GSE121212","GSE54456"]
plat={"GSE13355":"array","GSE30999":"array","GSE14905":"array","GSE117468":"array","GSE121212":"RNA-seq","GSE54456":"RNA-seq"}
cls=["co-dominant","IL-17-dominant","IL-23-dominant","neither"]
col={"co-dominant":"#0072B2","IL-17-dominant":"#E69F00","IL-23-dominant":"#CC79A7","neither":"#999999"}
mk={"co-dominant":"o","IL-17-dominant":"s","IL-23-dominant":"^","neither":"x"}
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8})
fig,axs=plt.subplots(2,3,figsize=(7.0,4.8))
for ax,g in zip(axs.flat,order):
    s=d[d.dataset==g]
    for c in cls:
        t=s[s.class_primary==c]
        if len(t): ax.scatter(t.il23_mod,t.il17_resp_mod,s=10,c=col[c],marker=mk[c],label=c,linewidths=0.6)
    ax.axvline(2,ls="--",c="k",lw=0.7); ax.axhline(2,ls="--",c="k",lw=0.7)
    rho=spearmanr(s.il23_mod,s.il17_resp_mod)[0]
    ax.set_title(f"{g} ({plat[g]}, n = {len(s)})\nSpearman rho {rho:.2f}",fontsize=8)
    ax.set_xlabel("IL-23 axis module Z"); ax.set_ylabel("IL-17 response module Z")
h,l=[],[]
for a in axs.flat:
    for hh,ll in zip(*a.get_legend_handles_labels()):
        if ll not in l: h.append(hh); l.append(ll)
fig.legend(h,l,loc="lower center",ncol=4,frameon=False)
fig.tight_layout(rect=(0,0.05,1,1)); fig.savefig("Figure1_module_scatter.png",dpi=600); fig.savefig("Figure1_module_scatter.tif",dpi=600)
tab=pd.crosstab(d.dataset,d.class_primary).reindex(order).fillna(0)
for c in cls:
    if c not in tab: tab[c]=0
pct=tab[cls].div(tab.sum(1),axis=0)*100
fig,ax=plt.subplots(figsize=(7.0,3.0)); left=[0]*len(order)
lab=[f"{g}\n({plat[g]}, n = {int(tab.loc[g].sum())})" for g in order]
for c in cls:
    ax.barh(lab,pct[c],left=left,color=col[c],label=f"{c} (total {int(tab[c].sum())})",edgecolor="white",hatch={"co-dominant":"","IL-17-dominant":"//","IL-23-dominant":"..","neither":"xx"}[c])
    left=[a+b for a,b in zip(left,pct[c])]
ax.invert_yaxis(); ax.set_xlim(0,100); ax.set_xlabel("Lesional biopsies (%)")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.22),ncol=4,frameon=False,fontsize=7)
fig.tight_layout(); fig.savefig("Figure2_class_proportions.png",dpi=600); fig.savefig("Figure2_class_proportions.tif",dpi=600)
tab.to_csv("figure2_counts.csv"); print(tab)
