#!/usr/bin/env python3
"""Additional reporting with the project age split fixed at 6 years.

Age rule used everywhere in this script:
  child: age < 6 years
  adult: age >= 6 years

This script runs after analysis_06_clinical_summary.py and adds:
- annotated TRBV/TRBJ heatmaps with patient prevalence and total reads per cell;
- TRBV priority scores combining healthy-vs-T1 and paired T1-to-T2 shifts;
- one summary of all valid patient T1/T2 Venn comparisons against age-matched healthy repertoires;
- audit tables for age-group assignments and clinical timepoints.
"""
from pathlib import Path
import re, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results'
CORR = BASE / 'corrected'
OUT = BASE / 'corrected_age6'
FIG = ROOT / 'figures' / 'corrected_age6'
TAB = OUT
for p in [OUT, FIG]:
    p.mkdir(parents=True, exist_ok=True)

AGE_SPLIT_YEARS = 6.0
rng = np.random.default_rng(1031)
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 14, 'axes.labelsize': 12})

meta = pd.read_csv(CORR / 'sample_metadata_enriched_with_analysis.csv')
clono = pd.read_csv(BASE / 'library_trb_clonotypes.csv')
manifest = pd.read_csv(ROOT / 'config' / 'sample_manifest.csv')
clinical = pd.read_csv(ROOT / 'config' / 'clinical_metadata_sequenced.csv')


def parse_age_years(x):
    if pd.isna(x):
        return np.nan
    s = str(x).lower().strip().replace(' ', '')
    years = 0.0
    my = re.search(r'(\d+)(?:г|лет|год)', s)
    if my:
        years += int(my.group(1))
    mm = re.search(r'(\d+)мес', s)
    if mm:
        years += int(mm.group(1)) / 12.0
    return years if years > 0 else np.nan


meta['age_years'] = meta['age_text'].map(parse_age_years)
meta['age_group_from_age6_rule'] = np.where(
    meta['age_years'].notna() & (meta['age_years'] < AGE_SPLIT_YEARS), 'child', 'adult'
)
# Legacy patients have no age_text in the clinical table but are known adults.
meta.loc[meta['batch'].eq('legacy'), 'age_group_from_age6_rule'] = 'adult'
# Unresolved sample with no age remains unknown and is excluded from biology anyway.
meta.loc[meta['age_years'].isna() & ~meta['batch'].eq('legacy'), 'age_group_from_age6_rule'] = meta.loc[
    meta['age_years'].isna() & ~meta['batch'].eq('legacy'), 'age_group'
]
meta['age_group_manifest'] = meta['age_group']
meta['age_group'] = meta['age_group_from_age6_rule']
meta['age_group_changed_vs_previous'] = meta['age_group_manifest'] != meta['age_group']
meta['age_group_rule'] = 'child if age < 6 years; adult if age >= 6 years'
meta['analysis_include'] = meta['analysis_include'].astype(bool)
meta['group_simple'] = np.where(
    meta['disease_status'].eq('healthy') | meta['is_control'].fillna(False),
    'Healthy', meta['clinical_timepoint']
)

# Fail loudly if config and parsed age rule disagree for included patients with known age.
known = meta[meta['analysis_include'] & meta['age_years'].notna()].copy()
wrong = known[known['age_group_manifest'] != known['age_group_from_age6_rule']]
if not wrong.empty:
    cols = ['library_id', 'patient_id', 'age_text', 'age_group_manifest', 'age_group_from_age6_rule']
    raise RuntimeError('sample_manifest age_group is inconsistent with the 6-year rule:\n' + wrong[cols].to_string(index=False))

inc_meta = meta[meta['analysis_include']].copy()
inc_meta.to_csv(TAB / 'sample_metadata_age6.csv', index=False)

# Audit table: one row per library, preserving source suffix and clinical timepoint.
audit_cols = [
    'library_id','raw_sample_id','patient_id','source_numeric_timepoint','clinical_timepoint',
    'age_text','age_years','age_group','disease_status','analysis_include','productive_trb_reads',
    'productive_trb_clonotypes','identity_qc','metadata_note'
]
meta[audit_cols].sort_values(['patient_id','source_numeric_timepoint']).to_csv(TAB / 'timepoint_age_audit.csv', index=False)

# Analysis cohort counts.
summary = []
for age in ['adult','child']:
    z = inc_meta[inc_meta['age_group'] == age]
    summary.append({
        'age_group': age,
        'n_samples_included': len(z),
        'n_patients_included': z['patient_id'].nunique(),
        'n_healthy_samples': int(z['group_simple'].eq('Healthy').sum()),
        'n_T0_samples': int(z['group_simple'].eq('T0').sum()),
        'n_T1_samples': int(z['group_simple'].eq('T1').sum()),
        'n_T2_samples': int(z['group_simple'].eq('T2').sum()),
        'n_T3_samples': int(z['group_simple'].eq('T3').sum()),
    })
pd.DataFrame(summary).to_csv(TAB / 'age6_group_summary.csv', index=False)

# Join clonotypes to authoritative corrected metadata.
join_cols = ['library_id','patient_id','age_group','disease_status','group_simple','clinical_timepoint']
cl = clono.merge(inc_meta[join_cols], on='library_id', how='inner', suffixes=('','_meta'))
if 'patient_id_meta' in cl:
    cl['patient_id'] = cl['patient_id_meta']
    cl = cl.drop(columns=[c for c in cl.columns if c.endswith('_meta')])
for c in ['v_gene','j_gene','cdr3aa']:
    cl[c] = cl[c].fillna('').astype(str)
cl['count'] = pd.to_numeric(cl['count'], errors='coerce').fillna(0)

sample_info = inc_meta[['library_id','patient_id','age_group','disease_status','group_simple','clinical_timepoint']].copy()


def build_gene_matrices(gene_col):
    x = cl.groupby(['library_id', gene_col], as_index=False)['count'].sum()
    totals = cl.groupby('library_id', as_index=False)['count'].sum().rename(columns={'count':'sample_total'})
    x = x.merge(totals, on='library_id', how='left')
    x['fraction'] = x['count'] / x['sample_total']
    frac = x.pivot_table(index='library_id', columns=gene_col, values='fraction', aggfunc='sum', fill_value=0)
    cnt = x.pivot_table(index='library_id', columns=gene_col, values='count', aggfunc='sum', fill_value=0)
    libs = sample_info[['library_id']]
    frac = libs.merge(frac, left_on='library_id', right_index=True, how='left').fillna(0).set_index('library_id')
    cnt = libs.merge(cnt, left_on='library_id', right_index=True, how='left').fillna(0).set_index('library_id')
    return frac, cnt


vfrac, vcounts = build_gene_matrices('v_gene')
jfrac, jcounts = build_gene_matrices('j_gene')


def fmt_reads(x):
    if x >= 1e6:
        return f'{x/1e6:.1f}M'
    if x >= 1e3:
        return f'{x/1e3:.0f}k'
    return str(int(round(x)))


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f'{name}.png', dpi=190, bbox_inches='tight')
    fig.savefig(FIG / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


# TRBV priority score: difference Healthy->T1 plus paired T1->T2 shift.
priority_parts = []
for age in ['adult','child']:
    info = sample_info[sample_info['age_group'] == age]
    healthy = info.loc[info['group_simple'].eq('Healthy'), 'library_id'].tolist()
    t1 = info.loc[info['group_simple'].eq('T1') & info['disease_status'].eq('measles'), 'library_id'].tolist()
    pair_tab = info[info['group_simple'].isin(['T1','T2'])].pivot_table(
        index='patient_id', columns='group_simple', values='library_id', aggfunc='first'
    )
    pair_tab = pair_tab.dropna(subset=['T1','T2']) if {'T1','T2'}.issubset(pair_tab.columns) else pd.DataFrame()
    p1libs = pair_tab['T1'].tolist() if not pair_tab.empty else []
    p2libs = pair_tab['T2'].tolist() if not pair_tab.empty else []
    if not healthy or not t1:
        continue
    rec = []
    for gene in [g for g in vfrac.columns if g]:
        h = vfrac.loc[healthy, gene].to_numpy() * 100
        a = vfrac.loc[t1, gene].to_numpy() * 100
        d_ht1 = float(np.mean(a) - np.mean(h))
        p_ht1 = mannwhitneyu(a, h, alternative='two-sided').pvalue if len(a)>=2 and len(h)>=2 else np.nan
        if p1libs:
            d = (vfrac.loc[p2libs, gene].to_numpy() - vfrac.loc[p1libs, gene].to_numpy()) * 100
            d_t12 = float(np.mean(d))
            p_t12 = wilcoxon(d).pvalue if len(d)>=3 and np.any(np.abs(d)>1e-12) else np.nan
        else:
            d = np.array([]); d_t12 = np.nan; p_t12 = np.nan
        rec.append({
            'age_group': age, 'gene': gene,
            'n_healthy': len(healthy), 'n_T1': len(t1), 'n_pairs_T1_T2': len(d),
            'patients_present_healthy': int((vfrac.loc[healthy, gene] > 0).sum()),
            'patients_present_T1': int((vfrac.loc[t1, gene] > 0).sum()),
            'delta_healthy_to_T1_pp': d_ht1, 'p_healthy_to_T1': p_ht1,
            'delta_T1_to_T2_pp': d_t12, 'p_T1_to_T2': p_t12,
        })
    z = pd.DataFrame(rec)
    for pcol, qcol in [('p_healthy_to_T1','q_healthy_to_T1'),('p_T1_to_T2','q_T1_to_T2')]:
        z[qcol] = np.nan
        idx = z[pcol].notna()
        if idx.any():
            z.loc[idx, qcol] = multipletests(z.loc[idx, pcol], method='fdr_bh')[1]
    s1 = z['delta_healthy_to_T1_pp'].std(ddof=0) or 1.0
    s2 = z['delta_T1_to_T2_pp'].fillna(0).std(ddof=0) or 1.0
    z['score_healthy_to_T1'] = z['delta_healthy_to_T1_pp'].abs() / s1
    z['score_T1_to_T2'] = z['delta_T1_to_T2_pp'].abs().fillna(0) / s2
    z['returns_toward_healthy'] = np.sign(z['delta_healthy_to_T1_pp']).replace(0,np.nan) != np.sign(z['delta_T1_to_T2_pp']).replace(0,np.nan)
    z['priority_score'] = (z['score_healthy_to_T1'] + z['score_T1_to_T2']) / 2
    z.loc[z['returns_toward_healthy'].fillna(False), 'priority_score'] *= 1.25
    priority_parts.append(z)
priority = pd.concat(priority_parts, ignore_index=True) if priority_parts else pd.DataFrame()
priority.to_csv(TAB / 'TRBV_priority_scores.csv', index=False)


def make_usage_heatmap(age, gene_type, priority_mode=False):
    frac, cnt = (vfrac, vcounts) if gene_type == 'TRBV' else (jfrac, jcounts)
    info = sample_info[sample_info['age_group'] == age]
    groups = [g for g in ['Healthy','T0','T1','T2','T3'] if g in set(info['group_simple'])]
    totals = info.groupby('group_simple')['library_id'].nunique().to_dict()
    if priority_mode and gene_type == 'TRBV' and not priority.empty:
        genes = priority[priority['age_group'].eq(age)].nlargest(10, 'priority_score')['gene'].tolist()
        title_tail = 'most biologically shifted'
        name = f'TRBV_priority_{age}_annotated'
    else:
        rows = []
        for gene in [g for g in frac.columns if g]:
            means = []
            for group in groups:
                libs = info.loc[info['group_simple'].eq(group), 'library_id']
                means.append(frac.loc[libs, gene].mean() if len(libs) else 0)
            rows.append((gene, np.mean(means)))
        top_n = 12 if gene_type == 'TRBV' else 8
        genes = [x[0] for x in sorted(rows, key=lambda x:x[1], reverse=True)[:top_n]]
        title_tail = 'most abundant'
        name = f'{gene_type}_usage_{age}_annotated'
    mat = np.zeros((len(genes), len(groups)))
    annotations = [['' for _ in groups] for _ in genes]
    for i,gene in enumerate(genes):
        for j,group in enumerate(groups):
            libs = info.loc[info['group_simple'].eq(group), 'library_id'].tolist()
            vals = frac.loc[libs, gene] if libs else pd.Series(dtype=float)
            mat[i,j] = 100*vals.mean() if len(vals) else 0
            present = int((vals > 0).sum()) if len(vals) else 0
            reads = cnt.loc[libs, gene].sum() if libs else 0
            annotations[i][j] = f'{present}/{totals.get(group,0)}\n{fmt_reads(reads)}'
    fig,ax = plt.subplots(figsize=(max(8,1.2*len(groups)+5), max(5.5,.55*len(genes)+2)))
    im = ax.imshow(mat, aspect='auto')
    ax.set_xticks(np.arange(len(groups)), groups)
    ax.set_yticks(np.arange(len(genes)), genes)
    ax.set_xlabel('Clinical group')
    ax.set_ylabel(gene_type)
    ax.set_title(
        f'{age.capitalize()}: {gene_type} usage, {title_tail}\n'
        'Age split: child < 6 years, adult >= 6 years\n'
        'Color: mean read share (%) | Cell: patients with gene / all patients, total reads'
    )
    vmax = float(np.nanmax(mat)) if mat.size else 0
    for i in range(len(genes)):
        for j in range(len(groups)):
            color = 'white' if vmax and mat[i,j] > .55*vmax else 'black'
            ax.text(j,i,annotations[i][j],ha='center',va='center',fontsize=8,color=color)
    fig.colorbar(im, ax=ax, label='Mean share of productive TRB reads (%)')
    save(fig, name)


for age in ['adult','child']:
    if (sample_info['age_group']==age).any():
        make_usage_heatmap(age, 'TRBV')
        make_usage_heatmap(age, 'TRBJ')
        if not priority.empty and age in set(priority['age_group']):
            make_usage_heatmap(age, 'TRBV', priority_mode=True)
            z = priority[priority['age_group'].eq(age)].nlargest(10,'priority_score').sort_values('priority_score')
            fig,ax = plt.subplots(figsize=(10,6.5)); y=np.arange(len(z)); ax.barh(y,z['priority_score']); ax.set_yticks(y,z['gene'])
            ax.set_xlabel('Combined priority score')
            ax.set_title(f'{age.capitalize()}: TRBV genes prominent in both contrasts\nHealthy vs T1 and paired T1 to T2')
            mx = z['priority_score'].max() if len(z) else 1
            for i,r in enumerate(z.itertuples(index=False)):
                ax.text(r.priority_score+.02*mx,i,f'H→T1 {r.delta_healthy_to_T1_pp:+.2f} pp | T1→T2 {r.delta_T1_to_T2_pp:+.2f} pp',va='center',fontsize=8)
            save(fig, f'TRBV_priority_scores_{age}')

# Venn-summary logic using CDR3beta sets.
pt_sets = {(pid,tp):set(g['cdr3aa']) for (pid,tp),g in cl.groupby(['patient_id','group_simple']) if tp in {'Healthy','T0','T1','T2','T3'}}
healthy_sets = {}
for age,info in sample_info[sample_info['group_simple'].eq('Healthy')].groupby('age_group'):
    healthy_sets[age] = [(lib, set(cl.loc[cl['library_id'].eq(lib),'cdr3aa'])) for lib in info['library_id']]

def jaccard(a,b):
    u = a | b
    return len(a & b)/len(u) if u else np.nan

venn_rows = []
for pid,g in inc_meta.groupby('patient_id'):
    gps = set(g['group_simple'])
    if not {'T1','T2'}.issubset(gps):
        continue
    age = g['age_group'].iloc[0]
    s1,s2 = pt_sets.get((pid,'T1'),set()), pt_sets.get((pid,'T2'),set())
    if not s1 or not s2:
        continue
    hj = [jaccard(s1,h) for _,h in healthy_sets.get(age,[]) if h]
    shared = len(s1 & s2)
    venn_rows.append({
        'patient_id':pid,'age_group':age,
        't1_only':len(s1-s2),'shared':shared,'t2_only':len(s2-s1),
        'jaccard_T1_T2':jaccard(s1,s2),
        'mean_jaccard_healthy_vs_T1':np.nanmean(hj) if hj else np.nan,
        'retained_from_T1_pct':100*shared/len(s1),
        'new_in_T2_pct':100*len(s2-s1)/len(s2),
        'n_healthy_refs':len(hj),
    })
venn = pd.DataFrame(venn_rows)
venn.to_csv(TAB / 'venn_T1_T2_vs_healthy_summary.csv', index=False)

fig,axes = plt.subplots(1,2,figsize=(14,6.5)); ax=axes[0]
for age,offset in [('adult',0),('child',3)]:
    z=venn[venn['age_group'].eq(age)] if not venn.empty else pd.DataFrame()
    if z.empty: continue
    x1,x2=1+offset,2+offset
    for r in z.itertuples(index=False):
        ax.plot([x1,x2],[100*r.mean_jaccard_healthy_vs_T1,100*r.jaccard_T1_T2],marker='o',alpha=.6)
    ax.hlines(100*np.nanmedian(z['mean_jaccard_healthy_vs_T1']),x1-.18,x1+.18,linewidth=3)
    ax.hlines(100*np.nanmedian(z['jaccard_T1_T2']),x2-.18,x2+.18,linewidth=3)
    ymax=np.nanmax(np.r_[100*z['mean_jaccard_healthy_vs_T1'],100*z['jaccard_T1_T2']]); ax.text((x1+x2)/2,ymax+.25,f'{age.capitalize()}\nn={len(z)}',ha='center')
ax.set_xticks([1,2,4,5],['Healthy vs T1','T1 vs T2','Healthy vs T1','T1 vs T2']); ax.set_ylabel('Jaccard overlap of unique CDR3β clonotypes (%)')
ax.set_title('Summary of patient Venn diagrams\nAge split: child < 6 years, adult >= 6 years'); ax.grid(axis='y',alpha=.2)
ax=axes[1]
if not venn.empty:
    z=venn.groupby('age_group')[['t1_only','shared','t2_only']].median().reindex(['adult','child']).dropna(how='all'); x=np.arange(len(z)); bottom=np.zeros(len(z))
    for col in ['t1_only','shared','t2_only']:
        vals=z[col].fillna(0).to_numpy(); ax.bar(x,vals,bottom=bottom,label=col.replace('_',' '))
        for i,v in enumerate(vals):
            if v>0: ax.text(i,bottom[i]+v/2,str(int(round(v))),ha='center',va='center',fontsize=9)
        bottom += vals
    ns=venn.groupby('age_group').size().reindex(z.index); ax.set_xticks(x,[f'{a.capitalize()}\nn={ns[a]}' for a in z.index]); ax.set_ylabel('Median number of unique clonotypes'); ax.set_title('Typical T1/T2 Venn composition'); ax.legend(title='Venn region')
save(fig,'venn_summary_T1_T2_vs_healthy')

# Three-or-more distinct clinical timepoints after QC.
three=[]
for pid,g in inc_meta.groupby('patient_id'):
    tps=[x for x in ['T0','T1','T2','T3'] if x in set(g['group_simple'])]
    if len(tps)>=3:
        three.append({'patient_id':pid,'age_group':g['age_group'].iloc[0],'timepoints_present':'; '.join(tps)})
pd.DataFrame(three).to_csv(TAB/'patients_with_3plus_timepoints.csv',index=False)

(OUT/'README.txt').write_text(
    'Age-6 reporting outputs. Age rule: child < 6 years; adult >= 6 years.\n'
    'Annotated V/J cells show: patients with gene / all patients in that group, then total reads.\n',
    encoding='utf-8'
)
print('Age-6 reporting complete:', OUT)
