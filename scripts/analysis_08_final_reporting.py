#!/usr/bin/env python3
"""Authoritative final reporting for the measles TRB project.

Key rules:
- child: age < 6 years; adult: age >= 6 years;
- clinical T0/T1/T2/T3 are taken from curated metadata, never inferred from .1/.2/.3;
- biological analyses require >= 5,000 productive TRB clonotypes per library;
- failed/suspect libraries are shown in QC but excluded from biological comparisons;
- patient overlap figures explicitly show missing/failed clinical points (e.g. CA1_013 T1).
"""
from pathlib import Path
import re, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results'
CORR = BASE / 'corrected'
OUT = BASE / 'current'
FIG = ROOT / 'figures' / 'current'
PATFIG = FIG / 'patients'
MIN_PRODUCTIVE_CLONOTYPES = 5000
AGE_SPLIT_YEARS = 6.0
RAREFACTION_DEPTH = 154184

for p in [OUT, FIG, PATFIG]:
    p.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(1031)
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 14, 'axes.labelsize': 12})

meta = pd.read_csv(CORR / 'sample_metadata_enriched_with_analysis.csv')
clono = pd.read_csv(BASE / 'library_trb_clonotypes.csv')


def parse_age_years(x):
    if pd.isna(x) or str(x).strip() == '':
        return np.nan
    s = str(x).lower().replace(' ', '')
    y = 0.0
    my = re.search(r'(\d+)(?:г|лет|год)', s)
    if my:
        y += int(my.group(1))
    mm = re.search(r'(\d+)мес', s)
    if mm:
        y += int(mm.group(1))/12.0
    return y if y > 0 else np.nan


meta['age_years'] = meta['age_text'].map(parse_age_years)
meta['age_group_previous'] = meta['age_group']
meta['age_group'] = np.where(meta['age_years'].notna() & (meta['age_years'] < AGE_SPLIT_YEARS), 'child', 'adult')
meta.loc[meta['batch'].eq('legacy'), 'age_group'] = 'adult'
meta.loc[meta['age_years'].isna() & ~meta['batch'].eq('legacy'), 'age_group'] = meta.loc[
    meta['age_years'].isna() & ~meta['batch'].eq('legacy'), 'age_group_previous'
]
for c in ['include_config','fastq_ready']:
    if c in meta:
        meta[c] = meta[c].astype(str).str.lower().isin(['true','1','yes'])
meta['productive_trb_reads'] = pd.to_numeric(meta['productive_trb_reads'], errors='coerce').fillna(0).astype(int)
meta['productive_trb_clonotypes'] = pd.to_numeric(meta['productive_trb_clonotypes'], errors='coerce').fillna(0).astype(int)
meta['clone_count_qc'] = np.where(meta['productive_trb_clonotypes'] >= MIN_PRODUCTIVE_CLONOTYPES, 'PASS', 'FAIL_LOW_CLONOTYPES')
meta.loc[meta['productive_trb_clonotypes'].eq(0), 'clone_count_qc'] = 'FAIL_EMPTY_TRB'
meta['quality_grade'] = np.select(
    [
        meta['identity_qc'].ne('PASS'),
        meta['productive_trb_clonotypes'].eq(0),
        meta['productive_trb_clonotypes'].lt(MIN_PRODUCTIVE_CLONOTYPES),
        meta['productive_trb_clonotypes'].ge(20000) & meta['productive_trb_reads'].ge(1000000),
        meta['productive_trb_clonotypes'].ge(10000) & meta['productive_trb_reads'].ge(500000),
    ],
    ['IDENTITY_SUSPECT','FAILED','LOW','VERY_HIGH','GOOD'],
    default='ADEQUATE'
)
base_ok = (
    meta['include_config']
    & meta['fastq_ready']
    & meta['qc_status'].eq('PASS')
    & meta['identity_qc'].eq('PASS')
    & meta['clinical_timepoint'].notna()
)
meta['analysis_include_final'] = base_ok & meta['productive_trb_clonotypes'].ge(MIN_PRODUCTIVE_CLONOTYPES)
meta['final_exclusion_reason'] = ''
meta.loc[~meta['include_config'], 'final_exclusion_reason'] = meta.loc[~meta['include_config'], 'exclusion_reason'].fillna('excluded by curated metadata')
meta.loc[meta['qc_status'].ne('PASS'), 'final_exclusion_reason'] = 'No productive TRB repertoire after MiXCR'
meta.loc[meta['identity_qc'].ne('PASS'), 'final_exclusion_reason'] = 'Repertoire identity QC suspect; excluded pending sample identity check'
meta.loc[base_ok & meta['productive_trb_clonotypes'].lt(MIN_PRODUCTIVE_CLONOTYPES), 'final_exclusion_reason'] = 'Below final QC threshold: <5,000 productive TRB clonotypes'
meta['group_simple'] = np.where(
    meta['disease_status'].eq('healthy') | meta['is_control'].astype(str).str.lower().isin(['true','1','yes']),
    'Healthy', meta['clinical_timepoint']
)

keep_cols = [
    'library_id','raw_sample_id','patient_id','clinical_timepoint','source_numeric_timepoint','age_text','age_years','age_group',
    'disease_status','is_control','batch','raw_fastq_total_read_pairs','raw_reads_total','productive_trb_reads',
    'productive_trb_clonotypes','largest_clone_percent','quality_grade','clone_count_qc','identity_qc',
    'analysis_include_final','final_exclusion_reason','rarefaction_depth','included_in_rarefaction','rarefied_richness',
    'shannon','hill_q1','top10_percent','rarefied_top10_percent'
]
keep_cols = [c for c in keep_cols if c in meta.columns]
meta[keep_cols].sort_values(['patient_id','source_numeric_timepoint']).to_csv(OUT / 'sample_metadata_current.csv', index=False)

valid = meta[meta['analysis_include_final']].copy()
metric_cols = ['rarefied_richness','shannon','hill_q1','rarefied_top10_percent','top10_percent']
for c in metric_cols:
    valid[c] = pd.to_numeric(valid[c], errors='coerce')
phase_metrics = valid.groupby(
    ['patient_id','age_group','disease_status','batch','group_simple','clinical_timepoint'], as_index=False
).agg(
    productive_trb_reads=('productive_trb_reads','sum'),
    productive_trb_clonotypes=('productive_trb_clonotypes','sum'),
    rarefied_richness=('rarefied_richness','mean'),
    shannon=('shannon','mean'),
    hill_q1=('hill_q1','mean'),
    rarefied_top10_percent=('rarefied_top10_percent','mean'),
    top10_percent=('top10_percent','mean')
)
phase_metrics.to_csv(OUT / 'patient_phase_metrics_current.csv', index=False)


def cliffs_delta(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    x=x[np.isfinite(x)]; y=y[np.isfinite(y)]
    if len(x)==0 or len(y)==0: return np.nan
    return (np.sum(x[:,None]>y)-np.sum(x[:,None]<y))/(len(x)*len(y))


def fnum(x):
    return 'NA' if not np.isfinite(x) else f'{x:.3g}'


def save(fig,name):
    fig.tight_layout()
    fig.savefig(FIG/f'{name}.png', dpi=200, bbox_inches='tight')
    fig.savefig(FIG/f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


fig, ax = plt.subplots(figsize=(10,5.5))
qc_counts = meta['quality_grade'].value_counts().reindex(['VERY_HIGH','GOOD','ADEQUATE','LOW','FAILED','IDENTITY_SUSPECT']).fillna(0)
ax.bar(np.arange(len(qc_counts)), qc_counts.values)
ax.set_xticks(np.arange(len(qc_counts)), ['Very high','Good','Adequate','Low\n<5k clones','Failed','Identity\nsuspect'])
ax.set_ylabel('Libraries')
ax.set_title(f'Library QC overview\nFinal biological threshold: ≥ {MIN_PRODUCTIVE_CLONOTYPES:,} productive TRB clonotypes')
for i,v in enumerate(qc_counts.values): ax.text(i,v+0.4,str(int(v)),ha='center')
save(fig,'00_library_QC_overview')

phase_order={'T0':0,'T1':1,'T2':2,'T3':3,'LEGACY_CASE':0}
def first_case(batch, age):
    z=phase_metrics[(phase_metrics['batch']==batch)&(phase_metrics['age_group']==age)&phase_metrics['disease_status'].eq('measles')].copy()
    z['ord']=z['clinical_timepoint'].map(phase_order).fillna(99)
    return z.sort_values(['patient_id','ord']).groupby('patient_id',as_index=False).first()
def controls(batch, age):
    z=phase_metrics[(phase_metrics['batch']==batch)&(phase_metrics['age_group']==age)&phase_metrics['disease_status'].eq('healthy')]
    return z.groupby('patient_id',as_index=False).first()

cohorts=[('Adults, new','new','adult'),('Adults, legacy','legacy','adult'),('Children, new','new','child')]
stats=[]
for label,batch,age in cohorts:
    ca=first_case(batch,age); co=controls(batch,age)
    for metric in ['rarefied_top10_percent','rarefied_richness','shannon']:
        x=pd.to_numeric(ca[metric],errors='coerce').dropna().to_numpy(); y=pd.to_numeric(co[metric],errors='coerce').dropna().to_numpy()
        p=mannwhitneyu(x,y,alternative='two-sided').pvalue if len(x)>=2 and len(y)>=2 else np.nan
        stats.append([label,metric,len(x),len(y),np.nanmedian(x) if len(x) else np.nan,np.nanmedian(y) if len(y) else np.nan,cliffs_delta(x,y),p,np.nan])
S=pd.DataFrame(stats,columns=['cohort','metric','n_measles','n_healthy','median_measles','median_healthy','cliffs_delta','p_value','fdr'])
for coh in S['cohort'].unique():
    idx=S.index[(S['cohort']==coh)&S['p_value'].notna()]
    if len(idx): S.loc[idx,'fdr']=multipletests(S.loc[idx,'p_value'],method='fdr_bh')[1]
S.to_csv(OUT/'cross_sectional_statistics_current.csv',index=False)

for metric,ylabel,name in [
    ('rarefied_top10_percent','Reads represented by top 10 TRB clones (%)','01_oligoclonality_case_control'),
    ('rarefied_richness','Unique TRB clonotypes after common-depth downsampling','02_richness_case_control'),
    ('shannon','Shannon diversity index','03_shannon_case_control')]:
    fig,ax=plt.subplots(figsize=(13,6)); positions=[]; labels=[]; pos=1
    for label,batch,age in cohorts:
        ca=first_case(batch,age); co=controls(batch,age)
        a=pd.to_numeric(co[metric],errors='coerce').dropna().to_numpy(); b=pd.to_numeric(ca[metric],errors='coerce').dropna().to_numpy()
        if len(a)+len(b)==0: continue
        ax.boxplot([a,b],positions=[pos,pos+1],widths=.55,showfliers=False)
        for xx,v in [(pos,a),(pos+1,b)]: ax.scatter(np.full(len(v),xx)+rng.normal(0,.05,len(v)),v,s=28,alpha=.75)
        st=S[(S['cohort']==label)&(S['metric']==metric)].iloc[0]
        ymax=np.nanmax(np.r_[a,b]); ymin=np.nanmin(np.r_[a,b])
        ax.text(pos+.5,ymax+(ymax-ymin+1)*.06,f"p={fnum(st.p_value)}, q={fnum(st.fdr)}, δ={fnum(st.cliffs_delta)}",ha='center',fontsize=8)
        positions += [pos,pos+1]; labels += [f'{label}\nHealthy\nn={len(a)}',f'{label}\nMeasles\nn={len(b)}']; pos+=3
    ax.set_xticks(positions,labels); ax.set_ylabel(ylabel)
    ax.set_title(f'Final QC-filtered analysis (≥{MIN_PRODUCTIVE_CLONOTYPES:,} productive clonotypes/library)')
    ax.grid(axis='y',alpha=.2); save(fig,name)

paired_stats=[]
for age in ['adult','child']:
    a=phase_metrics[(phase_metrics.age_group==age)&(phase_metrics.group_simple=='T1')].set_index('patient_id')
    b=phase_metrics[(phase_metrics.age_group==age)&(phase_metrics.group_simple=='T2')].set_index('patient_id')
    pats=sorted(set(a.index)&set(b.index))
    for metric in ['rarefied_top10_percent','rarefied_richness','shannon']:
        x=pd.to_numeric(a.loc[pats,metric],errors='coerce') if pats else pd.Series(dtype=float)
        y=pd.to_numeric(b.loc[pats,metric],errors='coerce') if pats else pd.Series(dtype=float)
        mask=x.notna()&y.notna(); x=x[mask].to_numpy();y=y[mask].to_numpy()
        p=wilcoxon(y,x).pvalue if len(x)>=3 and np.any(y!=x) else np.nan
        paired_stats.append([age,metric,len(x),np.nanmedian(x) if len(x) else np.nan,np.nanmedian(y) if len(y) else np.nan,np.nanmedian(y-x) if len(x) else np.nan,p,np.nan])
P=pd.DataFrame(paired_stats,columns=['age_group','metric','n_pairs','median_T1','median_T2','median_change','p_value','fdr'])
for age in P.age_group.unique():
    idx=P.index[(P.age_group==age)&P.p_value.notna()]
    if len(idx): P.loc[idx,'fdr']=multipletests(P.loc[idx,'p_value'],method='fdr_bh')[1]
P.to_csv(OUT/'T1_T2_paired_statistics_current.csv',index=False)

for metric,ylabel,name in [
    ('rarefied_top10_percent','Top-10 read share (%)','04_T1_T2_oligoclonality'),
    ('rarefied_richness','Rarefied TRB richness','05_T1_T2_richness'),
    ('shannon','Shannon diversity','06_T1_T2_shannon')]:
    fig,ax=plt.subplots(figsize=(10,6)); xpos={'adult':(1,2),'child':(4,5)}
    for age in ['adult','child']:
        a=phase_metrics[(phase_metrics.age_group==age)&(phase_metrics.group_simple=='T1')].set_index('patient_id')
        b=phase_metrics[(phase_metrics.age_group==age)&(phase_metrics.group_simple=='T2')].set_index('patient_id')
        pats=sorted(set(a.index)&set(b.index)); used=0
        for pid in pats:
            v1=pd.to_numeric(pd.Series([a.loc[pid,metric]]),errors='coerce').iloc[0]; v2=pd.to_numeric(pd.Series([b.loc[pid,metric]]),errors='coerce').iloc[0]
            if not(np.isfinite(v1)&np.isfinite(v2)): continue
            ax.plot(xpos[age],[v1,v2],marker='o',alpha=.55,linewidth=1.2); used+=1
        st=P[(P.age_group==age)&(P.metric==metric)].iloc[0]
        ax.text(np.mean(xpos[age]),.98,f'{age.capitalize()}: n={used}, p={fnum(st.p_value)}, q={fnum(st.fdr)}',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=9)
    ax.set_xticks([1,2,4,5],['Adult T1','Adult T2','Child T1','Child T2']); ax.set_ylabel(ylabel)
    ax.set_title('Within-patient clinical T1 → T2'); ax.grid(axis='y',alpha=.2); save(fig,name)

cl = clono.merge(valid[['library_id','patient_id','age_group','disease_status','group_simple']],on='library_id',how='inner',suffixes=('','_meta'))
if 'patient_id_meta' in cl: cl['patient_id']=cl['patient_id_meta']
cl['count']=pd.to_numeric(cl['count'],errors='coerce').fillna(0)
for c in ['cdr3aa','v_gene','j_gene']: cl[c]=cl[c].fillna('').astype(str)
cl['aav']=list(zip(cl['cdr3aa'],cl['v_gene']))
phase_sets={}; phase_counts={}
for (pid,tp),g in cl.groupby(['patient_id','group_simple']):
    phase_sets[(pid,tp)]=set(g['aav']); phase_counts[(pid,tp)]=g.groupby('aav')['count'].sum()

def jaccard(a,b): return len(a&b)/len(a|b) if (a|b) else np.nan
def morisita(pid,tp1,tp2):
    x=phase_counts.get((pid,tp1));y=phase_counts.get((pid,tp2))
    if x is None or y is None:return np.nan
    idx=x.index.union(y.index);a=x.reindex(idx,fill_value=0).to_numpy(float);b=y.reindex(idx,fill_value=0).to_numpy(float)
    a/=a.sum();b/=b.sum();den=(a*a).sum()+(b*b).sum();return 2*(a*b).sum()/den if den else np.nan

healthy_refs={}
for age,g in valid[valid.disease_status=='healthy'].groupby('age_group'):
    healthy_refs[age]=[]
    for pid in g.patient_id.unique():
        s=phase_sets.get((pid,'Healthy'))
        if s: healthy_refs[age].append((pid,s))
vr=[]
for pid,g in valid[valid.disease_status=='measles'].groupby('patient_id'):
    age=g.age_group.iloc[0];s1=phase_sets.get((pid,'T1'));s2=phase_sets.get((pid,'T2'))
    if not s1 or not s2: continue
    hj=[jaccard(s1,hs) for _,hs in healthy_refs.get(age,[]) if hs]
    vr.append({'patient_id':pid,'age_group':age,'T1_clonotypes':len(s1),'T2_clonotypes':len(s2),'shared':len(s1&s2),'jaccard_T1_T2':jaccard(s1,s2),'morisita_horn_T1_T2':morisita(pid,'T1','T2'),'mean_jaccard_healthy_vs_T1':np.nanmean(hj) if hj else np.nan,'n_healthy_refs':len(hj)})
V=pd.DataFrame(vr);V.to_csv(OUT/'venn_T1_T2_vs_healthy_current.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(14,6));ax=axes[0]
for age,off in [('adult',0),('child',3)]:
    z=V[V.age_group==age]
    if z.empty:continue
    for r in z.itertuples(index=False):ax.plot([1+off,2+off],[100*r.mean_jaccard_healthy_vs_T1,100*r.jaccard_T1_T2],marker='o',alpha=.65)
    ax.hlines(100*np.nanmedian(z.mean_jaccard_healthy_vs_T1),.8+off,1.2+off,linewidth=3);ax.hlines(100*np.nanmedian(z.jaccard_T1_T2),1.8+off,2.2+off,linewidth=3)
    ax.text(1.5+off,.97,f'{age.capitalize()} n={len(z)}',transform=ax.get_xaxis_transform(),ha='center')
ax.set_xticks([1,2,4,5],['Healthy vs T1','T1 vs T2','Healthy vs T1','T1 vs T2']);ax.set_ylabel('Jaccard overlap of aaV clonotypes (%)');ax.set_title('All valid patient Venn comparisons');ax.grid(axis='y',alpha=.2)
ax=axes[1]
for i,age in enumerate(['adult','child']):
    z=V[V.age_group==age]
    if z.empty:continue
    vals=100*z.morisita_horn_T1_T2.dropna().to_numpy();ax.scatter(np.full(len(vals),i+1)+rng.normal(0,.04,len(vals)),vals,s=35,alpha=.75);ax.hlines(np.median(vals),i+.8,i+1.2,linewidth=3);ax.text(i+1,max(vals)+3,f'n={len(vals)}\nmedian={np.median(vals):.1f}%',ha='center')
ax.set_xticks([1,2],['Adults','Children']);ax.set_ylabel('Morisita-Horn similarity (%)');ax.set_title('Abundance-weighted T1 → T2 repertoire similarity');ax.grid(axis='y',alpha=.2);save(fig,'07_venn_and_weighted_similarity_summary')


def gene_heatmap(gene_col,gene_label,age,top_n,name):
    x=cl.groupby(['patient_id','group_simple',gene_col],as_index=False)['count'].sum();totals=cl.groupby(['patient_id','group_simple'],as_index=False)['count'].sum().rename(columns={'count':'total'});x=x.merge(totals,on=['patient_id','group_simple']);x['fraction']=x['count']/x['total']
    info=valid[valid.age_group==age][['patient_id','group_simple']].drop_duplicates();groups=[g for g in ['Healthy','T0','T1','T2','T3'] if g in set(info.group_simple)]
    mean_gene=x[x.patient_id.isin(info.patient_id)].groupby(gene_col)['fraction'].mean().sort_values(ascending=False);genes=[g for g in mean_gene.index if g][:top_n]
    mat=np.zeros((len(genes),len(groups)));ann=[['' for _ in groups] for _ in genes]
    for i,gene in enumerate(genes):
        for j,grp in enumerate(groups):
            patients=info[info.group_simple==grp].patient_id.unique();nall=len(patients);sub=x[(x.patient_id.isin(patients))&(x.group_simple==grp)&(x[gene_col]==gene)];vals=sub.set_index('patient_id')['fraction'].reindex(patients,fill_value=0)
            mat[i,j]=100*vals.mean() if nall else 0;present=int((vals>0).sum());reads=int(sub['count'].sum());ann[i][j]=f'{present}/{nall}\n{reads/1000:.0f}k' if reads>=1000 else f'{present}/{nall}\n{reads}'
    fig,ax=plt.subplots(figsize=(max(8,1.2*len(groups)+5),max(5.5,.55*len(genes)+2)));im=ax.imshow(mat,aspect='auto');ax.set_xticks(range(len(groups)),groups);ax.set_yticks(range(len(genes)),genes)
    ax.set_title(f'{age.capitalize()}: {gene_label} usage\nColor = mean read share; cell = patients with gene / all patients, total reads');ax.set_xlabel('Clinical group');ax.set_ylabel(gene_label);vmax=np.nanmax(mat) if mat.size else 0
    for i in range(len(genes)):
        for j in range(len(groups)):ax.text(j,i,ann[i][j],ha='center',va='center',fontsize=8,color='white' if vmax and mat[i,j]>.55*vmax else 'black')
    fig.colorbar(im,ax=ax,label='Mean productive-read share (%)');save(fig,name)
for age in ['adult','child']:
    if (valid.age_group==age).any():
        gene_heatmap('v_gene','TRBV',age,12,f'08_TRBV_usage_{age}');gene_heatmap('j_gene','TRBJ',age,8,f'09_TRBJ_usage_{age}')

phase_order2={'T0':0,'T1':1,'T2':2,'T3':3};patient_qc=[]
for pid,g_all in meta[(meta.batch=='new') & meta.disease_status.eq('measles')].groupby('patient_id'):
    clinical_rows=g_all[g_all.clinical_timepoint.isin(['T0','T1','T2','T3'])].sort_values('clinical_timepoint',key=lambda s:s.map(phase_order2))
    if clinical_rows['clinical_timepoint'].nunique()<2: continue
    for r in clinical_rows.itertuples(index=False):patient_qc.append({'patient_id':pid,'library_id':r.library_id,'raw_sample_id':r.raw_sample_id,'clinical_timepoint':r.clinical_timepoint,'productive_trb_reads':r.productive_trb_reads,'productive_trb_clonotypes':r.productive_trb_clonotypes,'quality_grade':r.quality_grade,'identity_qc':r.identity_qc,'analysis_include_final':r.analysis_include_final,'final_exclusion_reason':r.final_exclusion_reason})
    analyzable=clinical_rows[clinical_rows.analysis_include_final];phase_names=sorted(analyzable.clinical_timepoint.unique(),key=lambda z:phase_order2[z])
    fig=plt.figure(figsize=(11,8));gs=fig.add_gridspec(2,2,height_ratios=[1,2.3]);axq=fig.add_subplot(gs[0,:]);x=np.arange(len(clinical_rows));vals=clinical_rows.productive_trb_clonotypes.to_numpy();axq.bar(x,vals);axq.axhline(MIN_PRODUCTIVE_CLONOTYPES,linestyle='--',linewidth=1,label=f'QC threshold = {MIN_PRODUCTIVE_CLONOTYPES:,} clonotypes')
    labels=[]
    for i,r in enumerate(clinical_rows.itertuples(index=False)):
        labels.append(f'{r.clinical_timepoint}\n{r.raw_sample_id}');text=f'{int(r.productive_trb_clonotypes):,} clones\n{int(r.productive_trb_reads):,} reads';text += '\nEXCLUDED' if not r.analysis_include_final else '';axq.text(i,vals[i]+max(max(vals),MIN_PRODUCTIVE_CLONOTYPES)*.03,text,ha='center',fontsize=8)
    axq.set_xticks(x,labels);axq.set_ylabel('Productive TRB clonotypes');axq.set_title(f'{pid}: sample QC by clinical time point');axq.legend(fontsize=8);axq.grid(axis='y',alpha=.2)
    axv=fig.add_subplot(gs[1,0]);axt=fig.add_subplot(gs[1,1]);axt.axis('off');venn_done=False
    if len(phase_names)>=2:
        if len(phase_names)==2:
            p1,p2=phase_names[:2];s1=phase_sets.get((pid,p1),set());s2=phase_sets.get((pid,p2),set())
            if s1 and s2:venn2([s1,s2],set_labels=(p1,p2),ax=axv);venn_done=True;axv.set_title('aaV overlap: CDR3β amino acid + TRBV')
        else:
            p1,p2,p3=phase_names[:3];sets=[phase_sets.get((pid,p),set()) for p in [p1,p2,p3]]
            if all(sets):venn3(sets,set_labels=(p1,p2,p3),ax=axv);venn_done=True;axv.set_title('aaV overlap: CDR3β amino acid + TRBV')
    if not venn_done:axv.axis('off');axv.text(.5,.5,'Fewer than 2 QC-passing clinical phases\nNo biological Venn comparison',ha='center',va='center',fontsize=12)
    txt=['Clinical interpretation:']
    for r in clinical_rows.itertuples(index=False):txt.append(f'• {r.clinical_timepoint}: included ({int(r.productive_trb_clonotypes):,} clonotypes)' if r.analysis_include_final else f'• {r.clinical_timepoint}: excluded — {r.final_exclusion_reason}')
    txt += ['',f'Final minimum = {MIN_PRODUCTIVE_CLONOTYPES:,} productive TRB clonotypes/library.','Only QC-passing phases enter overlap and group statistics.'];axt.text(0,.98,'\n'.join(txt),va='top',fontsize=10,wrap=True)
    fig.suptitle(f'{pid}: clinical repertoire QC and overlap',fontsize=16);fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(PATFIG/f'{pid}_repertoire_QC_and_overlap.png',dpi=200,bbox_inches='tight');fig.savefig(PATFIG/f'{pid}_repertoire_QC_and_overlap.pdf',bbox_inches='tight');plt.close(fig)
pd.DataFrame(patient_qc).to_csv(OUT/'patient_sample_QC_current.csv',index=False)

ref_path=BASE/'antigen_reference_used.csv';match_path=BASE/'antigen_specific_clonotypes.csv'
if ref_path.exists():
    ref=pd.read_csv(ref_path);paths=['Measles','Hepatitis B','Poliovirus','Tetanus'];rr=[];matches=pd.read_csv(match_path) if match_path.exists() else pd.DataFrame()
    for p in paths:
        nref=ref.loc[ref.pathogen.eq(p),'cdr3_trimmed'].nunique() if 'pathogen' in ref else 0;nmatch=matches.loc[matches.pathogen.eq(p),'cdr3aa'].nunique() if (not matches.empty and 'pathogen' in matches) else 0;rr.append([p,nref,nmatch])
    R=pd.DataFrame(rr,columns=['pathogen','reference_CDR3','exact_matches_detected']);R.to_csv(OUT/'specificity_reference_coverage_current.csv',index=False);fig,ax=plt.subplots(figsize=(8,5));x=np.arange(len(R));ax.bar(x,R.reference_CDR3);ax.set_xticks(x,R.pathogen);ax.set_ylabel('Unique reference CDR3β sequences');ax.set_title('Exact-match specificity reference coverage (QC, not a biological endpoint)')
    for i,r in enumerate(R.itertuples(index=False)):ax.text(i,r.reference_CDR3+max(R.reference_CDR3.max()*.02,1),f'ref={r.reference_CDR3}\nmatches={r.exact_matches_detected}',ha='center',fontsize=9)
    save(fig,'10_specificity_reference_coverage')

summary=[
    'CURRENT authoritative outputs',f'Generated with final QC threshold: >= {MIN_PRODUCTIVE_CLONOTYPES:,} productive TRB clonotypes/library.','Age rule: child < 6 years; adult >= 6 years.','Clinical timepoints come from curated metadata.','CA1_013 T1 is shown in patient QC but excluded because it has 0 productive TRB clonotypes.','MM1_033.2 remains excluded because repertoire identity QC is suspect.','',f'Libraries passing final biological QC: {int(meta.analysis_include_final.sum())}/{len(meta)}',f'Patient QC figures: {len(list(PATFIG.glob("*.png")))}',
]
(OUT/'README_CURRENT.txt').write_text('\n'.join(summary),encoding='utf-8');(FIG/'GENERATED_CURRENT.txt').write_text('\n'.join(summary),encoding='utf-8');print('\n'.join(summary))
