import pandas as pd, numpy as np, gzip, json
from scipy import stats
from load import series_matrix as _sm, D
import os as _o; RMA=_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..','rma')+'/'
def series_matrix(gse):
    X,M=_sm(gse)
    if gse in ('GSE13355','GSE30999','GSE14905','GSE117468'):
        R=pd.read_csv(RMA+f'{gse}_rma.tsv.gz',sep='\t',index_col=0)
        keep=[c for c in M.index if c in R.columns]
        import os
        if os.environ.get('QCDROP'):
            q=pd.read_csv(RMA+f'{gse}_qcflags.tsv',sep='\t'); bad=set(q.gsm[q.flag_any])
            keep=[c for c in keep if c not in bad]
        return R[keep],M
    return X,M
import os
OUT='../out_rma_qc/' if os.environ.get('QCDROP') else '../out_rma/'
IL23=['IL23A','IL12B','IL23R']
IL23X=IL23+['IL22','IL26']
IL17R=['CCL20','DEFB4','LCN2','S100A7','S100A8','S100A9','PI3','IL36G','CXCL1','CXCL8']
IL17L=['IL17A','IL17F']
ALL=sorted(set(IL23X+IL17R+IL17L))
# Ensembl 116 eligible probe sets (from _sdata_build Table 3 / tier_u133.tsv)
PROBES={'CCL20':['205476_at'],'DEFB4':['207356_at'],'LCN2':['212531_at'],'S100A7':['205916_at'],'S100A8':['202917_s_at'],
'S100A9':['203535_at'],'PI3':['203691_at','41469_at'],'IL36G':['220322_at'],'CXCL1':['204470_at'],'CXCL8':['202859_x_at','211506_s_at'],
'IL23A':['220054_at'],'IL12B':['207901_at'],'IL23R':['1552912_a_at','1561853_a_at'],'IL22':['221165_s_at','222974_at'],'IL26':['221111_at'],
'IL17A':['208402_at','216876_s_at'],'IL17F':['234408_at'],
'XIST':['214218_s_at','221728_x_at','224588_at','224589_at','224590_at','227671_at','243712_at'],'RPS4Y1':['201909_at'],'DDX3Y':['205000_at'],'KDM5D':['206700_s_at']}
YPROBES=['201909_at','205000_at','206700_s_at']
def array_genes(X):
    G={};chosen={}
    for g,ps in PROBES.items():
        ps=[p for p in ps if p in X.index]
        if not ps: continue
        best=X.loc[ps].mean(axis=1).idxmax(); G[g]=X.loc[best]; chosen[g]=best
    return pd.DataFrame(G).T, chosen
def infer_sex(G):
    x=G.loc['XIST']; y=G.loc['RPS4Y1']
    tx=(x.quantile(.05)+x.quantile(.95))/2; ty=(y.quantile(.05)+y.quantile(.95))/2
    return np.where((x>tx)&(y<ty),'F',np.where((x<tx)&(y>ty),'M','indeterminate'))
def rnaseq_genes(E,sym):
    alias={'DEFB4':['DEFB4A'],'CXCL8':['CXCL8','IL8'],'IL36G':['IL36G','IL1F9'],'S100A7':['S100A7']}
    G={}
    for g in list(PROBES):
        for a in alias.get(g,[g]):
            if a in E.index: G[g]=E.loc[a]; break
    return pd.DataFrame(G).T
datasets={}
# ---- microarrays
X,M=series_matrix('GSE13355'); G,ch=array_genes(X)
grp=M['char0'].map({'involved skin from cases':'LS','uninvolved skin from cases':'NL','normal skin from controls':'HC'})
subj=M['title'].str.extract(r'Individual_(\d+)_')[0]
datasets['GSE13355']=dict(G=G,grp=grp,subj=subj,platform='GPL570 array',values='RMA from CEL files (this study)',X=X,chosen=ch)
X,M=series_matrix('GSE30999'); G,ch=array_genes(X)
grp=M['biopsy type'].map({'non-lesion (NL)':'NL','psoriasis lesion (LS)':'LS'})
datasets['GSE30999']=dict(G=G,grp=grp,subj=M['subject'],platform='GPL570 array',values='RMA from CEL files (this study)',X=X,chosen=ch)
X,M=series_matrix('GSE14905'); G,ch=array_genes(X)
grp=M['source_name_ch1'].map({'Skin biopsy, lesional':'LS','Skin biopsy, non-lesional':'NL','Skin biopsy, normal':'HC'})
subj=M['title'].str.extract(r'(Patient-\d+|Normal-\d+)')[0]
datasets['GSE14905']=dict(G=G,grp=grp,subj=subj,platform='GPL570 array',values='RMA from CEL files (this study)',X=X,chosen=ch)
X,M=series_matrix('GSE117468'); G,ch=array_genes(X)
bl=M['visit']=='BL'
grp=pd.Series(np.where(~bl,'OT',M['tissue'].map({'lesional skin':'LS','non-lesional skin':'NL'})),index=M.index)
datasets['GSE117468']=dict(G=G,grp=grp,subj=M['patientid'],platform='GPL570 array',values='RMA from CEL files (this study)',X=X,chosen=ch,M=M)
# ---- RNA-seq
C=pd.read_csv(D+'GSE121212_readcount.txt.gz',sep='\t',index_col=0)
C=C.groupby(level=0).sum()
lib=C.sum(); E=np.log2(C/lib*1e6+0.5)
keep=[c for c in C.columns if c.startswith('PSO') or c.startswith('CTRL')]
E=E[keep]; G=rnaseq_genes(E,None)
grp=pd.Series(['LS' if c.endswith('_lesional') else 'NL' if c.endswith('non-lesional') else 'HC' for c in keep],index=keep)
subj=pd.Series([c.rsplit('_',1)[0] for c in keep],index=keep)
datasets['GSE121212']=dict(G=G,grp=grp,subj=subj,platform='RNA-seq (Illumina HiSeq 2500), read counts',values='log2 counts per million from deposited read counts',CPM=(C[keep]/lib[keep]*1e6))
R=pd.read_csv(D+'GSE54456_RPKM_samples.txt.gz',sep='\t',index_col=0)
_,M=series_matrix('GSE54456')
t2g=dict(zip(M['title'],M['tissue type']))
R=R.groupby(level=0).mean()
E=np.log2(R+1); G=rnaseq_genes(E,None)
grp=pd.Series([ 'LS' if 'lesional' in t2g[c] else 'HC' for c in R.columns],index=R.columns)
datasets['GSE54456']=dict(G=G,grp=grp,subj=pd.Series(R.columns,index=R.columns),platform='RNA-seq (Illumina), RPKM',values='log2(RPKM+1) from deposited RPKM',RPKM=R)

res={}; lesion_tables=[]
THR=2.0
for name,d in datasets.items():
    G=d['G']; grp=d['grp'].reindex(G.columns)
    ref='NL' if (grp=='NL').sum()>=10 else 'HC'
    ls=grp[grp=='LS'].index; rf=grp[grp==ref].index
    # detection
    det={}
    if 'X' in d:
        sex=infer_sex(G); fem=G.columns[sex=='F']
        Y=d['X'].loc[[p for p in YPROBES if p in d['X'].index],fem]
        floor=float(np.nanpercentile(Y.values,95)) if len(fem)>=10 else float(np.nanpercentile(d['X'].mean(axis=1).nsmallest(len(d['X'])//10),95))
        for g in ALL:
            if g in G.index: det[g]=float((G.loc[g,ls]>floor).mean())
        dmeta=dict(rule='array noise floor = 95th percentile of Y-chromosome probe sets in expression-inferred female arrays',floor=round(floor,3),n_female_inferred=int(len(fem)))
    elif 'CPM' in d:
        cpm=d['CPM']
        al={'DEFB4':'DEFB4A','CXCL8':'CXCL8','IL36G':'IL36G'}
        for g in ALL:
            k=al.get(g,g)
            if k in cpm.index: det[g]=float((cpm.loc[k,ls]>=1).mean())
        dmeta=dict(rule='CPM >= 1')
    elif 'RPKM' in d:
        r=d['RPKM']; al={'DEFB4':'DEFB4A','CXCL8':'IL8','IL36G':'IL1F9'}
        for g in ALL:
            k=al.get(g,g)
            if k in r.index: det[g]=float((r.loc[k,ls]>=1).mean())
        dmeta=dict(rule='RPKM >= 1')
    else:
        dmeta=dict(rule='not computable: deposited values are centred residuals'); det={}
    # gene z vs reference
    mu=G[rf].mean(axis=1); sd=G[rf].std(axis=1).replace(0,np.nan)
    Z=G.sub(mu,axis=0).div(sd,axis=0)
    def mod(genes):
        gg=[g for g in genes if g in Z.index]
        s=Z.loc[gg].mean()
        return (s-s[rf].mean())/s[rf].std(), gg
    m23,g23=mod(IL23); m17,g17=mod(IL17R); m23x,_=mod(IL23X); m17l,_=mod(IL17L)
    def classify(a,b):
        return pd.Series(np.select([(a>THR)&(b>THR),(b>THR)&~(a>THR),(a>THR)&~(b>THR)],['co-dominant','IL-17-dominant','IL-23-dominant'],'neither'),index=a.index)
    cls=classify(m23[ls],m17[ls]); clsx=classify(m23x[ls],m17[ls]); clsl=classify(m23[ls],m17l[ls])
    # single-gene sensitivity: IL23A alone
    s=Z.loc['IL23A']; a=(s-s[rf].mean())/s[rf].std()
    cls_a=classify(a[ls],m17[ls])
    # paired lesional - reference within subject for paired sets
    rho,p=stats.spearmanr(m23[ls],m17[ls])
    # lesional vs ref single gene log2FC and paired p
    fc={}
    for g in ALL:
        if g in G.index:
            fc[g]=dict(log2FC_LS_vs_ref=round(float(G.loc[g,ls].mean()-G.loc[g,rf].mean()),2),
                       p_MWU=float(stats.mannwhitneyu(G.loc[g,ls],G.loc[g,rf]).pvalue))
    # residual approach: IL-17 output not explained by IL-23 module
    sl,ic,_,_,_=stats.linregress(m23[ls],m17[ls]); resid=m17[ls]-(ic+sl*m23[ls])
    res[name]=dict(platform=d['platform'],values=d['values'],reference=ref,n_LS=int(len(ls)),n_ref=int(len(rf)),
      il23_genes=g23,il17_genes=g17,detection=dmeta,detect_frac_LS={k:round(v,3) for k,v in det.items()},
      class_primary=cls.value_counts().to_dict(),class_ext23=clsx.value_counts().to_dict(),class_ligand17=clsl.value_counts().to_dict(),class_IL23A_only=cls_a.value_counts().to_dict(),
      il23_mod_median=round(float(m23[ls].median()),2),il23_mod_IQR=[round(float(x),2) for x in m23[ls].quantile([.25,.75])],
      il17_mod_median=round(float(m17[ls].median()),2),il17_mod_IQR=[round(float(x),2) for x in m17[ls].quantile([.25,.75])],
      il17lig_mod_median=round(float(m17l[ls].median()),2),
      spearman_il23_il17=[round(float(rho),3),float(p)],fc=fc,chosen=d.get('chosen'))
    t=pd.DataFrame(dict(dataset=name,sample=ls,subject=d['subj'].reindex(ls).values,il23_mod=m23[ls].values,il23ext_mod=m23x[ls].values,il17_resp_mod=m17[ls].values,il17_lig_mod=m17l[ls].values,class_primary=cls.values,class_ext=clsx.values,class_lig=clsl.values,class_IL23A=cls_a.values,il17_resid=resid.values))
    lesion_tables.append(t)
    d['m23']=m23;d['m17']=m17;d['m17l']=m17l;d['ls']=ls;d['ref']=rf
L=pd.concat(lesion_tables); L.to_csv(OUT+'lesion_level_scores.tsv',sep='\t',index=False)
json.dump(res,open(OUT+'summary.json','w'),indent=1,default=str)
# ---- GSE117468 exploratory response analysis
d=datasets['GSE117468']; M=d['M']; G=d['G']
M=M.copy(); M['pasi']=pd.to_numeric(M['pasi'],errors='coerce')
blp=M[(M.visit=='BL')&(M.tissue=='lesional skin')].groupby('patientid')['pasi'].first()
w12=M[(M.visit=='W12')].groupby('patientid')['pasi'].first()
arm=M.groupby('patientid')['treatment'].first()
lsb=M[(M.visit=='BL')&(M.tissue=='lesional skin')]
T=pd.DataFrame(dict(pid=lsb['patientid'],arm=lsb['treatment'],m23=d['m23'][lsb.index],m17=d['m17'][lsb.index]))
T['cls']=L.set_index('sample').loc[lsb.index,'class_primary'].values
T['pasi0']=T.pid.map(blp); T['pasi12']=T.pid.map(w12)
T['impr']=(T.pasi0-T.pasi12)/T.pasi0; T['pasi75']=T.impr>=0.75
T=T.dropna(subset=['pasi12'])
T.to_csv(OUT+'GSE117468_baseline_response.tsv',sep='\t',index=False)
out={}
for a,s in T.groupby('arm'):
    r=stats.spearmanr(s.m23,s.impr); r2=stats.spearmanr(s.m17,s.impr)
    out[a]=dict(n=len(s),pasi75=int(s.pasi75.sum()),class_counts=s.cls.value_counts().to_dict(),
       rho_il23_impr=[round(r.statistic,3),round(r.pvalue,4)],rho_il17_impr=[round(r2.statistic,3),round(r2.pvalue,4)],
       pasi75_by_class={k:f"{int(v.sum())}/{len(v)}" for k,v in s.groupby('cls').pasi75})
# interaction: active IL-17RA (brodalumab pooled) vs ustekinumab, continuous m23 on percent improvement (linear)
import numpy.linalg as la
s=T[T.arm!='Placebo'].copy(); s['brod']=s.arm.str.startswith('Brod').astype(float)
Xd=np.column_stack([np.ones(len(s)),s.brod,s.m23,s.brod*s.m23]); y=s.impr.values
b,res_,_,_=la.lstsq(Xd,y,rcond=None); e=y-Xd@b; sig2=e@e/(len(y)-4); cov=sig2*la.inv(Xd.T@Xd); se=np.sqrt(np.diag(cov))
tval=b[3]/se[3]; pint=2*stats.t.sf(abs(tval),len(y)-4)
out['interaction_brod_vs_ust_x_il23mod_on_pct_improvement']=dict(n=len(s),beta=round(b[3],4),se=round(se[3],4),p=round(pint,4))
json.dump(out,open(OUT+'GSE117468_response.json','w'),indent=1,default=str)
print(json.dumps({k:{kk:v[kk] for kk in ['n_LS','n_ref','reference','class_primary','class_ext23','class_ligand17','class_IL23A_only','il23_mod_median','il17_mod_median','spearman_il23_il17','detect_frac_LS','detection']} for k,v in res.items()},indent=1,default=str))
print(json.dumps(out,indent=1,default=str))
