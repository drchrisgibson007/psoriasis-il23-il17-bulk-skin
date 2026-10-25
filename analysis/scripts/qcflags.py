import pandas as pd, numpy as np, json, sys
from scipy.stats import fisher_exact, chi2_contingency
G=pd.read_csv('groups.tsv',sep='\t')
out={}
for g in sys.argv[1:]:
    q=pd.read_csv(f'{g}_qc.tsv',sep='\t').merge(G[G.cohort==g],on='gsm')
    q['f_nuse']=q.NUSE_median>1.05
    q['f_rle_med']=q.RLE_median.abs()>0.1
    q3,q1=q.RLE_IQR.quantile(.75),q.RLE_IQR.quantile(.25)
    q['f_rle_iqr']=q.RLE_IQR>q3+1.5*(q3-q1)
    q['flag_any']=q.f_nuse|q.f_rle_med|q.f_rle_iqr
    q['n_flags']=q[['f_nuse','f_rle_med','f_rle_iqr']].sum(axis=1)
    t=pd.crosstab(q.grp,q.flag_any)
    p=chi2_contingency(t)[1] if t.shape==(t.shape[0],2) and t.shape[1]==2 else None
    out[g]=dict(n=len(q),nuse_gt_1_05=int(q.f_nuse.sum()),abs_rle_median_gt_0_1=int(q.f_rle_med.sum()),rle_iqr_outlier=int(q.f_rle_iqr.sum()),
        any_flag=int(q.flag_any.sum()),any_flag_pct=round(100*q.flag_any.mean(),1),two_or_more=int((q.n_flags>=2).sum()),
        by_group={k:f"{int(v.sum())}/{len(v)}" for k,v in q.groupby('grp').flag_any},group_chi2_p=None if p is None else round(p,4),
        nuse_median_range=[round(q.NUSE_median.min(),3),round(q.NUSE_median.max(),3)])
    q.to_csv(f'{g}_qcflags.tsv',sep='\t',index=False)
print(json.dumps(out,indent=1)); json.dump(out,open('qcflags_'+'_'.join(sys.argv[1:])+'.json','w'),indent=1)
