#!/usr/bin/env python3
"""Clinically corrected summary analysis and QC for the measles TCR-beta project.

Uses explicit clinical T0/T1/T2/T3 values from config/sample_manifest.csv instead
of assuming that raw suffix .1/.2/.3 always denotes the same biological stage.
The old analysis scripts remain in the repository for audit/history; the figures
and tables produced here are the authoritative biological summary.
"""
from pathlib import Path
import re, warnings
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results';CONFIG=ROOT/'config'/'sample_manifest.csv';PATIENT_META=ROOT/'config'/'clinical_metadata_sequenced.csv';OUT=BASE/'corrected';FIG=ROOT/'figures'/'corrected';PATFIG=FIG/'patients';TAB=OUT
for p in [FIG,PATFIG,TAB]:p.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(1031)
manifest=pd.read_csv(CONFIG);fastq=pd.read_csv(BASE/'fastq_audit.csv');qc=pd.read_csv(BASE/'sample_qc.csv');patient_meta=pd.read_csv(PATIENT_META)
for c in ['include','is_control']:
    if c in manifest:manifest[c]=manifest[c].astype(str).str.lower().isin(['true','1','yes'])
manifest['analysis_timepoint']=pd.to_numeric(manifest.analysis_timepoint,errors='coerce')
patient_lookup=patient_meta.set_index('patient_id').to_dict('index')

def num(raw):
    m=re.search(r'\.(\d+)',str(raw));return int(m.group(1)) if m else np.nan
rows=[]
for r in manifest.itertuples(index=False):
    n=num(r.raw_sample_id);tp=str(r.clinical_timepoint).strip() if pd.notna(r.clinical_timepoint) else None
    if tp in ['', 'nan']:tp=None
    if r.batch=='legacy':tp='CONTROL' if bool(r.is_control) else 'LEGACY_CASE'
    pm=patient_lookup.get(r.patient_id,{})
    age_text=pm.get('age_text','') if pd.notna(pm.get('age_text','')) else ''
    src=pm.get('metadata_source','') if pd.notna(pm.get('metadata_source','')) else ''
    note=str(r.exclusion_reason).strip() if pd.notna(r.exclusion_reason) else ''
    rows.append({'library_id':r.library_id,'sequencing_index':r.sequencing_index,'batch':r.batch,'raw_sample_id':r.raw_sample_id,'patient_id':r.patient_id,'source_numeric_timepoint':n,'clinical_timepoint':tp,'age_group':r.age_group,'age_text':age_text,'disease_status':r.disease_status,'is_control':bool(r.is_control),'include_config':bool(r.include),'exclusion_reason':note,'clinical_source':src,'metadata_note':note})
meta=pd.DataFrame(rows).merge(fastq,on=['library_id','sequencing_index','batch'],how='left').merge(qc[['library_id','productive_trb_reads','productive_trb_clonotypes','largest_clone_fraction']],on='library_id',how='left')
meta['fastq_ready']=meta['fastq_ready'].astype(str).str.lower().isin(['true','1','yes'])
meta['productive_trb_reads']=pd.to_numeric(meta.productive_trb_reads,errors='coerce').fillna(0).astype(int);meta['productive_trb_clonotypes']=pd.to_numeric(meta.productive_trb_clonotypes,errors='coerce').fillna(0).astype(int);meta['productive_reads_millions']=meta.productive_trb_reads/1e6;meta['largest_clone_percent']=100*pd.to_numeric(meta.largest_clone_fraction,errors='coerce')
meta['qc_status']=np.where(meta.productive_trb_reads>0,'PASS','FAIL_EMPTY_TRB');meta['identity_qc']='PASS';meta.loc[meta.library_id=='new_39','identity_qc']='SUSPECT';meta.loc[meta.library_id=='new_39','metadata_note']='ММ1_033.2 identity-QC outlier; exclude pending sample identity check.'
meta['analysis_include']=meta.include_config & meta.fastq_ready & meta.qc_status.eq('PASS') & meta.identity_qc.eq('PASS') & meta.clinical_timepoint.notna()
raw_depth_path=BASE/'raw_fastq_depth.csv'
if raw_depth_path.exists():
    raw_depth=pd.read_csv(raw_depth_path)
    meta=meta.merge(raw_depth[['library_id','raw_read_pairs','raw_reads_total','compressed_bytes_all','status']],on='library_id',how='left',suffixes=('','_raw'))
    meta['raw_fastq_total_read_pairs']=meta['raw_read_pairs']
    meta['raw_input_note']='Raw FASTQ read counts collected from R1 files; raw_reads_total counts R1+R2 reads.'
else:
    meta['raw_fastq_total_read_pairs']=np.nan;meta['raw_reads_total']=np.nan;meta['compressed_bytes_all']=np.nan
    meta['raw_input_note']='Run scripts/collect_fastq_depth.py to populate total raw FASTQ read counts; productive TRB depth is available.'
inc=meta.loc[meta.analysis_include,'productive_trb_reads'];q25,q75=inc.quantile(.25),inc.quantile(.75);meta['productive_depth_flag']=meta.productive_trb_reads.map(lambda x:'empty' if x<=0 else 'lower_quartile' if x<q25 else 'upper_quartile' if x>q75 else 'middle_50%')

D=pd.read_csv(BASE/'diversity_metrics.csv');D['analysis_timepoint']=pd.to_numeric(D.analysis_timepoint,errors='coerce');D['is_control']=D.is_control.astype(str).str.lower().isin(['true','1','yes']);D['batch']=np.where(D['patient_id'].astype(str).str.startswith('UDP'),'legacy','new')
phase_map={(r.patient_id,int(r.source_numeric_timepoint)):r.clinical_timepoint for r in meta.itertuples(index=False) if r.batch=='new' and pd.notna(r.source_numeric_timepoint)}
inc_keys={(r.patient_id,int(r.source_numeric_timepoint)) for r in meta[meta.analysis_include].itertuples(index=False) if r.batch=='new' and pd.notna(r.source_numeric_timepoint)};inc_legacy=set(meta.loc[(meta.batch=='legacy')&meta.analysis_include,'patient_id'])
def phase(r):
    if r.is_control:return 'CONTROL'
    if str(r.patient_id).startswith('UDP'):return 'LEGACY_CASE'
    return phase_map.get((r.patient_id,int(r.analysis_timepoint))) if pd.notna(r.analysis_timepoint) else None
def keep(r):
    if str(r.patient_id).startswith('UDP'):return r.patient_id in inc_legacy
    return (r.patient_id,int(r.analysis_timepoint)) in inc_keys if pd.notna(r.analysis_timepoint) else False
D['clinical_timepoint']=D.apply(phase,axis=1);D['analysis_include']=D.apply(keep,axis=1);D=D[D.analysis_include&D.clinical_timepoint.notna()].copy();D.to_csv(TAB/'corrected_diversity_metrics.csv',index=False)
bio=D[['patient_id','analysis_timepoint','clinical_timepoint','raw_richness','rarefaction_depth','included_in_rarefaction','rarefied_richness','shannon','hill_q1','top10_fraction','rarefied_top10_fraction']]
meta2=meta.merge(bio,how='left',left_on=['patient_id','source_numeric_timepoint','clinical_timepoint'],right_on=['patient_id','analysis_timepoint','clinical_timepoint']).drop(columns=['analysis_timepoint'],errors='ignore');meta2['top10_percent']=100*meta2.top10_fraction;meta2['rarefied_top10_percent']=100*meta2.rarefied_top10_fraction;meta2.to_csv(TAB/'sample_metadata_enriched_with_analysis.csv',index=False)
ps=meta2.groupby('patient_id',as_index=False).agg(age_group=('age_group','first'),disease_status=('disease_status','first'),batch=('batch','first'),n_libraries=('library_id','nunique'),n_analyzable_libraries=('analysis_include','sum'),total_productive_trb_reads=('productive_trb_reads','sum'),min_productive_trb_reads=('productive_trb_reads','min'),max_productive_trb_reads=('productive_trb_reads','max'));pt=meta2[meta2.analysis_include].groupby('patient_id').clinical_timepoint.apply(lambda x:';'.join(sorted(set(x),key=lambda z:{'CONTROL':0,'LEGACY_CASE':1,'T0':2,'T1':3,'T2':4,'T3':5}.get(z,99))));ps=ps.merge(pt.rename('analyzable_clinical_timepoints'),on='patient_id',how='left');ps.to_csv(TAB/'patient_summary.csv',index=False)

def cliffs(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);return np.nan if len(x)==0 or len(y)==0 else (np.sum(x[:,None]>y)-np.sum(x[:,None]<y))/(len(x)*len(y))
def fmt(x):return 'NA' if not np.isfinite(x) else f'{x:.3g}'
def first_case(batch,age):
    z=D[(D.batch==batch)&(D.age_group==age)&(~D.is_control)].copy();z['ord']=z.clinical_timepoint.map({'T0':0,'T1':1,'T2':2,'T3':3,'LEGACY_CASE':0});return z.sort_values(['patient_id','ord']).groupby('patient_id',as_index=False).first()
def controls(batch,age):return D[(D.batch==batch)&(D.age_group==age)&D.is_control].groupby('patient_id',as_index=False).first()
coh=[('Adults, new batch','new','adult'),('Adults, legacy batch','legacy','adult'),('Children, new batch','new','child')];sr=[]
for label,b,a in coh:
    ca=first_case(b,a);co=controls(b,a)
    for met,sc in [('rarefied_top10_fraction',100),('rarefied_richness',1),('shannon',1)]:
        x=pd.to_numeric(ca[met],errors='coerce').dropna().to_numpy()*sc;y=pd.to_numeric(co[met],errors='coerce').dropna().to_numpy()*sc;p=mannwhitneyu(x,y).pvalue if len(x)>=2 and len(y)>=2 else np.nan;sr.append([label,met,len(x),len(y),np.median(x) if len(x) else np.nan,np.median(y) if len(y) else np.nan,np.median(x)-np.median(y) if len(x) and len(y) else np.nan,cliffs(x,y),p,np.nan])
S=pd.DataFrame(sr,columns=['cohort','metric','n_measles','n_healthy','median_measles','median_healthy','median_difference','cliffs_delta','p_value','fdr'])
for z in S.cohort.unique():
    idx=S.index[(S.cohort==z)&S.p_value.notna()]
    if len(idx):S.loc[idx,'fdr']=multipletests(S.loc[idx,'p_value'],method='fdr_bh')[1]
S.to_csv(TAB/'cross_sectional_statistics.csv',index=False)
newcase=D[(D.batch=='new')&(~D.is_control)&D.clinical_timepoint.isin(['T0','T1','T2','T3'])];pr=[]
for age in ['adult','child']:
    for a,b in [('T0','T1'),('T1','T2'),('T2','T3')]:
        da=newcase[(newcase.age_group==age)&(newcase.clinical_timepoint==a)].set_index('patient_id');db=newcase[(newcase.age_group==age)&(newcase.clinical_timepoint==b)].set_index('patient_id');pp=sorted(set(da.index)&set(db.index))
        for met,sc in [('rarefied_top10_fraction',100),('rarefied_richness',1),('shannon',1)]:
            xa=pd.to_numeric(da.loc[pp,met],errors='coerce') if pp else pd.Series(dtype=float);xb=pd.to_numeric(db.loc[pp,met],errors='coerce') if pp else pd.Series(dtype=float);mask=xa.notna()&xb.notna();x=xa[mask].to_numpy()*sc;y=xb[mask].to_numpy()*sc;p=wilcoxon(y,x).pvalue if len(x)>=3 and np.any(y!=x) else np.nan;pr.append([age,f'{a}->{b}',met,len(x),np.median(x) if len(x) else np.nan,np.median(y) if len(y) else np.nan,np.median(y-x) if len(x) else np.nan,p,np.nan])
P=pd.DataFrame(pr,columns=['age_group','transition','metric','n_pairs','median_before','median_after','median_change','p_value','fdr'])
for age in P.age_group.unique():
    idx=P.index[(P.age_group==age)&P.p_value.notna()]
    if len(idx):P.loc[idx,'fdr']=multipletests(P.loc[idx,'p_value'],method='fdr_bh')[1]
P.to_csv(TAB/'longitudinal_statistics.csv',index=False)

phase_libs=meta2[meta2.analysis_include].groupby(['patient_id','clinical_timepoint'])['library_id'].apply(list).to_dict()
def clean_gene(s):
    if pd.isna(s):return ''
    return str(s).split(',')[0].split('(')[0].split('*')[0].strip()
cache={}
def load_phase(pid,tp):
    key=(pid,tp)
    if key in cache:return cache[key]
    libs=phase_libs.get(key,[]);parts=[]
    for lib in libs:
        f=BASE/'clonotypes'/f'{lib}.tsv'
        if not f.exists():continue
        d=pd.read_csv(f,sep='\t',usecols=['readCount','aaSeqCDR3','allVHitsWithScore','allJHitsWithScore'])
        if d.empty:continue
        d=d[pd.to_numeric(d.readCount,errors='coerce').fillna(0)>0].copy();d['count']=pd.to_numeric(d.readCount);d['v_gene']=d.allVHitsWithScore.map(clean_gene);d['j_gene']=d.allJHitsWithScore.map(clean_gene);parts.append(d[['aaSeqCDR3','count','v_gene','j_gene']])
    if not parts:return None
    d=pd.concat(parts,ignore_index=True);best=d.sort_values('count',ascending=False).drop_duplicates('aaSeqCDR3')[['aaSeqCDR3','v_gene','j_gene']];cnt=d.groupby('aaSeqCDR3',as_index=False)['count'].sum();d=cnt.merge(best,on='aaSeqCDR3');d['aav']=list(zip(d.aaSeqCDR3.astype(str),d.v_gene.astype(str)));cache[key]=d;return d

def load_library(lib):
    f=BASE/'clonotypes'/f'{lib}.tsv'
    if not f.exists():return None
    d=pd.read_csv(f,sep='\t',usecols=['readCount','aaSeqCDR3','allVHitsWithScore','allJHitsWithScore'])
    if d.empty:return None
    d=d[pd.to_numeric(d.readCount,errors='coerce').fillna(0)>0].copy()
    if d.empty:return None
    d['count']=pd.to_numeric(d.readCount);d['v_gene']=d.allVHitsWithScore.map(clean_gene);d['j_gene']=d.allJHitsWithScore.map(clean_gene)
    best=d.sort_values('count',ascending=False).drop_duplicates('aaSeqCDR3')[['aaSeqCDR3','v_gene','j_gene']];cnt=d.groupby('aaSeqCDR3',as_index=False)['count'].sum();d=cnt.merge(best,on='aaSeqCDR3');d['aav']=list(zip(d.aaSeqCDR3.astype(str),d.v_gene.astype(str)));return d

def similarity(x,y):
    sx,sy=set(x.aav),set(y.aav);union=sx|sy;j=len(sx&sy)/len(union) if union else np.nan
    vx=x.groupby('aav')['count'].sum();vy=y.groupby('aav')['count'].sum();idx=vx.index.union(vy.index);xx=vx.reindex(idx,fill_value=0).to_numpy(float);yy=vy.reindex(idx,fill_value=0).to_numpy(float);px=xx/xx.sum();py=yy/yy.sum();mh=2*np.sum(px*py)/(np.sum(px*px)+np.sum(py*py));return j,mh

pairrows=[];trbv_diffs={'adult':{},'child':{}}
for age in ['adult','child']:
    a=newcase[(newcase.age_group==age)&(newcase.clinical_timepoint=='T1')].set_index('patient_id');b=newcase[(newcase.age_group==age)&(newcase.clinical_timepoint=='T2')].set_index('patient_id')
    for pid in sorted(set(a.index)&set(b.index)):
        x=load_phase(pid,'T1');y=load_phase(pid,'T2')
        if x is None or y is None:continue
        sx,sy=set(x.aav),set(y.aav);j=len(sx&sy)/len(sx|sy);vx=x.groupby('aav')['count'].sum();vy=y.groupby('aav')['count'].sum();idx=vx.index.union(vy.index);xx=vx.reindex(idx,fill_value=0).to_numpy(float);yy=vy.reindex(idx,fill_value=0).to_numpy(float);px=xx/xx.sum();py=yy/yy.sum();mh=2*np.sum(px*py)/(np.sum(px*px)+np.sum(py*py))
        xv=x.groupby('v_gene')['count'].sum();yv=y.groupby('v_gene')['count'].sum();genes=xv.index.union(yv.index);xp=xv.reindex(genes,fill_value=0)/xv.sum();yp=yv.reindex(genes,fill_value=0)/yv.sum();tv=.5*np.abs(xp-yp).sum();pairrows.append([pid,age,j,mh,tv])
        for g in genes:trbv_diffs[age].setdefault(g,[]).append((yp.get(g,0)-xp.get(g,0))*100)
Q=pd.DataFrame(pairrows,columns=['patient_id','age_group','jaccard','morisita_horn','trbv_total_variation']);Q.to_csv(TAB/'T1_T2_repertoire_similarity.csv',index=False)
tr=[]
for age,gd in trbv_diffs.items():
    for g,dif in gd.items():
        dif=np.asarray(dif);p=wilcoxon(dif).pvalue if len(dif)>=3 and np.any(dif!=0) else np.nan;tr.append([age,g,len(dif),np.median(dif),p,np.nan])
T=pd.DataFrame(tr,columns=['age_group','trbv','n_pairs','median_change_pp','p_value','fdr'])
for age in T.age_group.unique():
    idx=T.index[(T.age_group==age)&T.p_value.notna()]
    if len(idx):T.loc[idx,'fdr']=multipletests(T.loc[idx,'p_value'],method='fdr_bh')[1]
T.to_csv(TAB/'TRBV_T1_T2_paired_statistics.csv',index=False)

usage=[]
for (pid,tp),libs in phase_libs.items():
    d=load_phase(pid,tp)
    if d is None:continue
    mr=meta2[(meta2.patient_id==pid)&(meta2.clinical_timepoint==tp)].iloc[0]
    for gene_col,label in [('v_gene','TRBV'),('j_gene','TRBJ')]:
        vc=d[d[gene_col]!=''].drop_duplicates('aaSeqCDR3').groupby(gene_col).size();den=vc.sum()
        for g,n in vc.items():usage.append([pid,tp,mr.age_group,bool(mr.is_control),mr.batch,label,g,n/den])
usage=pd.DataFrame(usage,columns=['patient_id','clinical_timepoint','age_group','is_control','batch','gene_type','gene','fraction']);usage.to_csv(TAB/'corrected_vj_usage_by_sample.csv',index=False)

ref=pd.read_csv(BASE/'antigen_reference_used.csv');am=pd.read_csv(BASE/'antigen_specific_clonotypes.csv');paths=['Measles','Hepatitis B','Poliovirus','Tetanus'];cov=pd.DataFrame({'pathogen':paths,'unique_reference_CDR3':[ref.loc[ref.pathogen==p,'cdr3_trimmed'].nunique() for p in paths],'exact_matches_detected':[am.loc[am.pathogen==p,'cdr3aa'].nunique() if not am.empty else 0 for p in paths]});cov.to_csv(TAB/'specificity_reference_coverage.csv',index=False)

plt.rcParams.update({'font.size':11,'axes.titlesize':15,'axes.labelsize':12})
def save(fig,n):fig.tight_layout();fig.savefig(FIG/f'{n}.png',dpi=180,bbox_inches='tight');fig.savefig(FIG/f'{n}.pdf',bbox_inches='tight');plt.close(fig)
def pts(ax,pos,v):v=np.asarray(v,float);ax.scatter(np.full(len(v),pos)+rng.normal(0,.05,len(v)),v,s=30,alpha=.8,zorder=3)
cc=[]
for label,b,a in coh:cc+=[(label+'\nHealthy',len(controls(b,a))),(label+'\nMeasles',len(first_case(b,a)))]
fig,ax=plt.subplots(figsize=(11,5));xx=np.arange(len(cc));vv=[x[1] for x in cc];ax.bar(xx,vv);ax.set_xticks(xx,[x[0] for x in cc],rotation=18,ha='right');ax.set_ylabel('Patients');ax.set_title('Analyzable cohort after metadata and QC corrections');[ax.text(i,n+.2,str(n),ha='center') for i,n in enumerate(vv)];save(fig,'00_cohort_after_QC')
for met,yl,sc,nm in [('rarefied_top10_fraction','Reads represented by top 10 TRB clonotypes (%)',100,'01_oligoclonality_case_control'),('rarefied_richness','Unique CDR3β clonotypes after common-depth downsampling',1,'02_rarefied_richness_case_control'),('shannon','Shannon diversity index',1,'03_shannon_case_control')]:
    fig,ax=plt.subplots(figsize=(13,6.5));pos=[1,2,4,5,7,8];labs=[]
    for i,(label,b,a) in enumerate(coh):
        ca=first_case(b,a);co=controls(b,a);x=pd.to_numeric(co[met],errors='coerce').dropna().to_numpy()*sc;y=pd.to_numeric(ca[met],errors='coerce').dropna().to_numpy()*sc;p1,p2=pos[2*i:2*i+2];ax.boxplot([x,y],positions=[p1,p2],widths=.55,showfliers=False);pts(ax,p1,x);pts(ax,p2,y);labs += [f'Healthy\nn={len(x)}',f'Measles\nn={len(y)}'];st=S[(S.cohort==label)&(S.metric==met)].iloc[0];yy=max(np.r_[x,y]);ax.text((p1+p2)/2,yy*1.08,f'{label}\np={fmt(st.p_value)}, q={fmt(st.fdr)}, δ={st.cliffs_delta:.2f}',ha='center',fontsize=9)
    ax.set_xticks(pos,labs);ax.set_ylabel(yl);ax.set_title('Measles vs healthy controls, separated by age and sequencing batch');ax.grid(axis='y',alpha=.2);save(fig,nm)
for met,yl,sc,nm in [('rarefied_top10_fraction','Top-10 read share (%)',100,'04_oligoclonality_T1_T2'),('rarefied_richness','Rarefied CDR3β richness',1,'05_richness_T1_T2'),('shannon','Shannon diversity',1,'06_shannon_T1_T2')]:
    fig,ax=plt.subplots(figsize=(10,6));ticks=[];labs=[]
    for gi,age in enumerate(['adult','child']):
        a=newcase[(newcase.age_group==age)&(newcase.clinical_timepoint=='T1')].set_index('patient_id');b=newcase[(newcase.age_group==age)&(newcase.clinical_timepoint=='T2')].set_index('patient_id');pp=sorted(set(a.index)&set(b.index));x1=1+3*gi;x2=2+3*gi;vals=[]
        for pid in pp:
            v1=a.loc[pid,met];v2=b.loc[pid,met]
            if np.isfinite(v1) and np.isfinite(v2):v1*=sc;v2*=sc;ax.plot([x1,x2],[v1,v2],marker='o',alpha=.55,linewidth=1.2);vals += [v1,v2]
        st=P[(P.age_group==age)&(P.transition=='T1->T2')&(P.metric==met)].iloc[0]
        if vals:ax.text((x1+x2)/2,max(vals)*1.08,f'{age.capitalize()}: n={int(st.n_pairs)}\np={fmt(st.p_value)}, q={fmt(st.fdr)}',ha='center',fontsize=9)
        ticks += [x1,x2];labs += [f'{age.capitalize()} T1',f'{age.capitalize()} T2']
    ax.set_xticks(ticks,labs);ax.set_ylabel(yl);ax.set_title('Within-patient change from clinical T1 to T2');ax.grid(axis='y',alpha=.2);save(fig,nm)
for col,yl,sc,nm in [('jaccard','aaV clonotype-set overlap, Jaccard (%)',100,'07_aav_overlap_T1_T2'),('morisita_horn','Abundance-weighted repertoire similarity (Morisita–Horn)',1,'08_weighted_similarity_T1_T2'),('trbv_total_variation','TRBV composition redistributed between T1 and T2 (%)',100,'09_TRBV_composition_shift_T1_T2')]:
    fig,ax=plt.subplots(figsize=(7,6))
    for po,age in [(1,'adult'),(2,'child')]:
        z=Q.loc[Q.age_group==age,col].to_numpy()*sc;pts(ax,po,z);ax.hlines(np.median(z),po-.18,po+.18,linewidth=2);ax.text(po,max(z)*1.05,f'n={len(z)}\nmedian={np.median(z):.2f}',ha='center',fontsize=9)
    a=Q.loc[Q.age_group=='adult',col].to_numpy()*sc;c=Q.loc[Q.age_group=='child',col].to_numpy()*sc;p=mannwhitneyu(a,c).pvalue if len(a)>=2 and len(c)>=2 else np.nan;ax.set_xticks([1,2],['Adults','Children']);ax.set_ylabel(yl);ax.set_title(f'T1→T2 repertoire change by age\np={fmt(p)}, Cliff δ={cliffs(a,c):.2f}');ax.grid(axis='y',alpha=.2);save(fig,nm)
for age in ['adult','child']:
    z=T[T.age_group==age].copy();z['abs']=z.median_change_pp.abs();z=z.nlargest(12,'abs').sort_values('median_change_pp');fig,ax=plt.subplots(figsize=(8,6));y=np.arange(len(z));ax.scatter(z.median_change_pp,y,s=55);ax.axvline(0,linewidth=1);ax.set_yticks(y,z.trbv);ax.set_xlabel('Median T1→T2 change in TRBV read share (percentage points)');ax.set_title(f'{age.capitalize()}: largest paired TRBV shifts');[ax.text(r.median_change_pp,i+.18,f'q={fmt(r.fdr)}',ha='center',fontsize=8) for i,r in enumerate(z.itertuples(index=False))];ax.grid(axis='x',alpha=.2);save(fig,f'10_TRBV_shifts_{age}')
for age in ['adult','child']:
    for gt,nm in [('TRBV',f'11_TRBV_usage_{age}'),('TRBJ',f'12_TRBJ_usage_{age}')]:
        z=usage[(usage.age_group==age)&(usage.gene_type==gt)].copy();z['group']=np.where(z.is_control,'Healthy',z.clinical_timepoint);mean=z.groupby(['group','gene']).fraction.mean().reset_index();top=mean.groupby('gene').fraction.mean().nlargest(12 if gt=='TRBV' else 8).index;order=[x for x in ['Healthy','LEGACY_CASE','T0','T1','T2','T3'] if x in mean.group.unique()];mat=mean[mean.gene.isin(top)].pivot(index='gene',columns='group',values='fraction').fillna(0).reindex(columns=order);fig,ax=plt.subplots(figsize=(8,6));im=ax.imshow(mat.to_numpy(),aspect='auto');ax.set_xticks(np.arange(len(mat.columns)),mat.columns,rotation=30,ha='right');ax.set_yticks(np.arange(len(mat.index)),mat.index);ax.set_title(f'{age.capitalize()}: mean unique-clonotype {gt} usage');ax.set_xlabel('Clinical group');ax.set_ylabel(gt);fig.colorbar(im,ax=ax,label='Mean fraction of unique CDR3β clonotypes');save(fig,nm)
fig,ax=plt.subplots(figsize=(8,5));x=np.arange(len(cov));ax.bar(x,cov.unique_reference_CDR3);ax.set_xticks(x,cov.pathogen);ax.set_ylabel('Unique CDR3β sequences in exact-match reference');ax.set_title('Specificity reference coverage: current limitation');[ax.text(i,r.unique_reference_CDR3+max(cov.unique_reference_CDR3.max()*.025,1),f'ref={r.unique_reference_CDR3}\nfound={r.exact_matches_detected}',ha='center',fontsize=9) for i,r in enumerate(cov.itertuples(index=False))];save(fig,'13_specificity_reference_QC')
ca=meta2[meta2.patient_id=='СА1_013'].sort_values('source_numeric_timepoint');fig,ax=plt.subplots(figsize=(8,5));x=np.arange(len(ca));ax.bar(x,ca.productive_trb_clonotypes);ax.set_xticks(x,[f'.{int(r.source_numeric_timepoint)} / {r.clinical_timepoint}' for r in ca.itertuples(index=False)]);ax.set_ylabel('Productive TRB clonotypes');ax.set_title('СА1_013: T1 exists but has no productive TRB after MiXCR');[ax.text(i,r.productive_trb_clonotypes+max(ca.productive_trb_clonotypes.max()*.03,1),f'{r.productive_trb_reads:,} reads\n{r.productive_trb_clonotypes:,} clones',ha='center',fontsize=9) for i,r in enumerate(ca.itertuples(index=False))];save(fig,'14_CA1_013_QC')
mm_specs=[('.1 vs .2','new_38','new_39'),('.2 vs .3','new_39','new_43'),('.1 vs .3','new_38','new_43')];mm_rows=[]
for label,la,lb in mm_specs:
    xa,xb=load_library(la),load_library(lb)
    if xa is None or xb is None:continue
    j,mh=similarity(xa,xb);mm_rows.append([label,la,lb,j,mh])
mm=pd.DataFrame(mm_rows,columns=['pair','library_a','library_b','jaccard','morisita_horn']);mm.to_csv(TAB/'MM1_033_identity_QC.csv',index=False)
if not mm.empty:
    fig,ax=plt.subplots(figsize=(7,5));x=np.arange(len(mm));ax.bar(x,mm.morisita_horn);ax.set_xticks(x,mm.pair);ax.set_ylim(0,1);ax.set_ylabel('Morisita–Horn similarity');ax.set_title('ММ1_033: .2 is an identity-QC outlier');[ax.text(i,v+.03,f'{v:.4f}',ha='center') for i,v in enumerate(mm.morisita_horn)];save(fig,'15_MM1_033_identity_QC')
phaseord={'T0':0,'T1':1,'T2':2,'T3':3}
for pid,g in meta2[meta2.analysis_include&meta2.clinical_timepoint.isin(phaseord)].groupby('patient_id'):
    phases=sorted(set(g.clinical_timepoint),key=lambda x:phaseord[x])
    if len(phases)<2:continue
    sets={}
    for tp in phases:
        d=load_phase(pid,tp);sets[tp]=set(d.aav) if d is not None else set()
    fig,ax=plt.subplots(figsize=(7,6));v=venn2([sets[phases[0]],sets[phases[1]]],set_labels=phases,ax=ax) if len(phases)==2 else venn3([sets[p] for p in phases[:3]],set_labels=phases[:3],ax=ax)
    if v:
        for patch in v.patches:
            if patch is not None:patch.set_alpha(.5);patch.set_edgecolor('black')
    ax.set_title(f'{pid}: TRB aaV overlap by clinical time point');fig.tight_layout();fig.savefig(PATFIG/f'{pid}_clinical_aav_venn.png',dpi=180,bbox_inches='tight');plt.close(fig)
(OUT/'README.txt').write_text('Corrected clinical-timepoint analysis. See sample_metadata_enriched_with_analysis.csv for per-library QC and analysis inclusion. Run scripts/collect_fastq_depth.py before this script to populate raw FASTQ read counts.\n',encoding='utf-8')
print('Corrected analysis complete.')
print('Corrected figures:', FIG)
print('Corrected tables:', TAB)
