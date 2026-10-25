#!/usr/bin/env python3
"""Check for donor overlap between cohorts using GEO metadata only.
Reads the sample headers of the series matrices in ../data/ and the files in ../data/raw/overlap_check/.
Writes ../outputs/overlap_check.txt."""
import gzip,re,os,html,collections
B=os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
def meta(g):
    d=collections.defaultdict(list)
    for l in gzip.open(f'{B}/data/{g}_series_matrix.txt.gz','rt',errors='ignore'):
        if l.startswith('!series_matrix_table_begin'): break
        if l.startswith('!'):
            k,*v=l.rstrip('\n').split('\t'); d[k].append([x.strip('"') for x in v])
    return d
G=['GSE13355','GSE30999','GSE14905','GSE117468','GSE121212','GSE54456']
M={g:meta(g) for g in G}
out=[]
for g in G:
    d=M[g]
    out.append(f"{g}: submitted {d['!Series_submission_date'][0][0]}; contact institute {d['!Series_contact_institute'][0][0]}; contributors {', '.join(c[0].replace(',,',' ').replace(',',' ') for c in d['!Series_contributor'])}; sample title example {d['!Sample_title'][0][0]}")
t13=set(re.findall(r'Individual_(\d+)',' '.join(M['GSE13355']['!Sample_title'][0])))
t54=set(x[1:] for x in M['GSE54456']['!Sample_title'][0] if re.fullmatch(r'M\d+',x))
ov=set(l.strip().lstrip('M') for l in gzip.open(f'{B}/data/raw/overlap_check/GSE54456_MAoverlappedsamples.txt.gz','rt') if l.strip())
t121=sorted(set(re.sub(r'_(lesional|non-lesional|chronic_lesion)$','',x) for x in M['GSE121212']['!Sample_title'][0]))
out+=['',f'GSE13355 individual numbers: {len(t13)}; GSE54456 M-numbers: {len(t54)}; shared numbers: {len(t13&t54)}',
f'GSE54456 MAoverlappedsamples.txt lists {len(ov)} samples; {len(ov&t13)} match GSE13355 individual numbers and {len(ov&t54)} match GSE54456 titles',
f'GSE121212 donor labels: {len(t121)} (examples {", ".join(t121[:3])} ... {t121[-1]}); numeric labels shared with GSE13355 or GSE54456 numbers: none possible (3-digit cohort codes vs 4-6 digit Michigan individual numbers)']
x=re.sub(r'\s+',' ',html.unescape(re.sub(r'<[^>]+>',' ',open(f'{B}/data/raw/overlap_check/Tsoi2019_PMC6711380_fulltext.xml',errors='ignore').read())))
m=re.search(r'Informed written consent was obtained[^.]*\.',x)
out+=['','Tsoi 2019 (GSE121212 source paper), Methods: '+(m.group(0) if m else 'NOT FOUND')]
for a,b in [('GSE30999','GSE117468')]:
    s1=set(r.split(':')[1].strip() for r in sum(M[a]['!Sample_characteristics_ch1'],[]) if r.startswith('subject'))
    s2=set(r.split(':')[1].strip() for r in sum(M[b]['!Sample_characteristics_ch1'],[]) if r.startswith('patientid'))
    out.append(f'{a} subject IDs {len(s1)}, {b} patient IDs {len(s2)}, shared {len(s1&s2)}')
os.makedirs(f'{B}/outputs',exist_ok=True)
open(f'{B}/outputs/overlap_check.txt','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
