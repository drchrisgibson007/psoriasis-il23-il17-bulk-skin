import pandas as pd, numpy as np, json
from scipy import stats
L=pd.read_csv('../out_rma/lesion_level_scores.tsv',sep='\t')
def wilson(k,n,z=1.96):
    p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [round(100*(c-h),1),round(100*(c+h),1)]
out={}
for ds,s in L.groupby('dataset',sort=False):
    n=len(s); r={}
    for thr in [1.645,2,3]:
        a=s.il23_mod>thr;b=s.il17_resp_mod>thr
        r[f'thr{thr}']=dict(co=int((a&b).sum()),il17=int((b&~a).sum()),il23=int((a&~b).sum()),neither=int((~a&~b).sum()))
    k=int((s.class_primary=='co-dominant').sum()); k17=int((s.class_primary=='IL-17-dominant').sum()); k23=int((s.class_primary=='IL-23-dominant').sum())
    r['n']=n; r['co_pct']=round(100*k/n,1); r['co_CI']=wilson(k,n); r['il17_pct']=round(100*k17/n,1); r['il17_CI']=wilson(k17,n)
    r['il23dom_pct']=round(100*k23/n,1); r['il23dom_CI_upper']=wilson(k23,n)[1]
    r['il17_resp_elev']=int((s.il17_resp_mod>2).sum())
    out[ds]=r
T=pd.read_csv('../out_rma/GSE117468_baseline_response.tsv',sep='\t')
fis={}
for a,s in T.groupby('arm'):
    x=s[s.cls.isin(['co-dominant','IL-17-dominant'])]
    tab=pd.crosstab(x.cls,x.pasi75).reindex(index=['co-dominant','IL-17-dominant'],columns=[True,False],fill_value=0)
    fis[a]=dict(table=tab.values.tolist(),p=round(stats.fisher_exact(tab.values)[1],4))
s=T[T.arm.str.startswith('Brod')]; x=s[s.cls.isin(['co-dominant','IL-17-dominant'])]
tab=pd.crosstab(x.cls,x.pasi75).reindex(index=['co-dominant','IL-17-dominant'],columns=[True,False],fill_value=0)
fis['Brodalumab pooled']=dict(table=tab.values.tolist(),p=round(stats.fisher_exact(tab.values)[1],4))
out['GSE117468_fisher_PASI75_codom_vs_il17dom']=fis
json.dump(out,open('../out_rma/sensitivity.json','w'),indent=1)
print(json.dumps(out,indent=1))
