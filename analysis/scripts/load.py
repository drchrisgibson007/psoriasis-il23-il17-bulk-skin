import gzip, pandas as pd, numpy as np, io, re
import os as _o; D=_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..','data')+'/'
def series_matrix(gse):
    meta={}; rows=[]; intable=False; chars=[]
    with gzip.open(D+f'{gse}_series_matrix.txt.gz','rt',errors='replace') as f:
        lines=[]
        for line in f:
            if line.startswith('!series_matrix_table_begin'): intable=True; continue
            if line.startswith('!series_matrix_table_end'): break
            if intable: lines.append(line)
            elif line.startswith('!Sample_'):
                p=line.rstrip('\n').split('\t'); k=p[0][8:]; v=[x.strip('"') for x in p[1:]]
                if k=='characteristics_ch1': chars.append(v)
                else: meta.setdefault(k,v)
    X=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',index_col=0)
    X.columns=[c.strip('"') for c in X.columns]; X.index=[str(i).strip('"') for i in X.index]
    M=pd.DataFrame({k:v for k,v in meta.items() if k in('geo_accession','title','source_name_ch1')})
    for i,c in enumerate(chars):
        key=c[0].split(':')[0] if ':' in c[0] else f'char{i}'
        M[key]=[x.split(':',1)[1].strip() if ':' in x else x for x in c]
    M.index=M['geo_accession']
    return X.apply(pd.to_numeric,errors='coerce'),M
