#!/usr/bin/env python3
"""Additive nucleotide-level clone fate analysis for the 22 paired T1-T2 patients.

Existing aaV analyses are preserved. This module adds:
- CDR3 nucleotide + TRBV + TRBJ tracking of every T1 clonotype
- CDR3 amino acid + TRBV + TRBJ functional/public key
- 0.10 percent primary T1-expanded candidate pool
- sensitivity at 0.05, 0.10, 0.20, 0.50, 0.75 and 1.00 percent
- analytical detection probabilities with no pseudocount for zero T2
- repeated depth-controlled resampling
- public/private and control prevalence annotations
- simple sequence-neighbour clusters
- T3 follow-up and repertoire identity QC

Sampling controls operate on productive read counts and do not correct PCR bias because
these bulk genomic DNA libraries do not have UMI information.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import math
import os
import re
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FIG = ROOT / "figures" / "current"
CURRENT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

PRIMARY = 0.001
THRESHOLDS = [0.0005, 0.001, 0.002, 0.005, 0.0075, 0.01]
REQUESTED_DEPTH = int(os.environ.get("MEASLES_RAREFACTION_DEPTH", "50000"))
N_ITER = int(os.environ.get("MEASLES_RAREFACTION_ITERATIONS", "500"))
SEED = int(os.environ.get("MEASLES_RAREFACTION_SEED", "1031"))
ROBUST_P = 0.99
MAX_EDIT = 2
KEY_NT = ["cdr3nt", "cdr3aa", "v_gene", "j_gene"]


def pick_col(df, names):
    for name in names:
        if name in df.columns:
            return name
    return None


def clean_gene(value):
    if pd.isna(value):
        return ""
    return str(value).split(",")[0].split("(")[0].split("*")[0].strip()


def age_years(text):
    s = str(text).lower()
    y = 0.0
    m = re.search(r"(\d+)\s*(?:год|года|лет|г\b)", s)
    if m:
        y += int(m.group(1))
    m = re.search(r"(\d+)\s*мес", s)
    if m:
        y += int(m.group(1)) / 12
    m = re.search(r"(\d+)\s*(?:дн|день|дня|дней)", s)
    if m:
        y += int(m.group(1)) / 365.25
    return y if y > 0 else np.nan


def fdr(values):
    values = pd.to_numeric(values, errors="coerce")
    out = np.full(len(values), np.nan)
    ok = values.notna().to_numpy()
    if ok.any():
        out[ok] = multipletests(values[ok], method="fdr_bh")[1]
    return out


def ensure_inputs():
    if not (CURRENT / "sample_metadata_current.csv").exists():
        subprocess.run([sys.executable, str(ROOT / "scripts" / "analysis_01_prepare.py")], check=True)


ensure_inputs()
meta = pd.read_csv(CURRENT / "sample_metadata_current.csv")
for col in ["analysis_include_final", "is_control"]:
    if col in meta.columns:
        meta[col] = meta[col].astype(str).str.lower().isin(["true", "1", "yes"])
if "age_years" not in meta.columns:
    meta["age_years"] = meta["age_text"].map(age_years)

# Read nucleotide CDR3 directly from existing MiXCR export tables. No FASTQ rerun is
# required when nSeqCDR3 is already present in the export.
frames = []
missing_nt = []
for r in meta[meta.analysis_include_final].itertuples(index=False):
    path = RESULTS / "clonotypes" / f"{r.library_id}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing accepted clonotype table {path}")
    df = pd.read_csv(path, sep="\t")
    nt = pick_col(df, ["nSeqCDR3", "cdr3nt", "cdr3_nt", "nucleotideSeqCDR3"])
    aa = pick_col(df, ["aaSeqCDR3", "cdr3aa", "cdr3_aa"])
    count = pick_col(df, ["cloneCount", "readCount", "count"])
    v = pick_col(df, ["allVHitsWithScore", "bestVHit", "vGene", "v_call"])
    j = pick_col(df, ["allJHitsWithScore", "bestJHit", "jGene", "j_call"])
    if nt is None:
        missing_nt.append(str(r.library_id))
        continue
    if aa is None or count is None:
        raise ValueError(f"Missing CDR3 amino acid or count column in {path.name}")
    z = pd.DataFrame({
        "library_id": str(r.library_id),
        "patient_id": str(r.patient_id),
        "clinical_timepoint": str(r.clinical_timepoint),
        "disease_status": str(r.disease_status),
        "age_group": str(r.age_group),
        "age_years": float(r.age_years) if pd.notna(r.age_years) else np.nan,
        "batch": str(r.batch),
        "cdr3nt": df[nt].astype(str),
        "cdr3aa": df[aa].astype(str),
        "v_gene": df[v].map(clean_gene) if v else "",
        "j_gene": df[j].map(clean_gene) if j else "",
        "count": pd.to_numeric(df[count], errors="coerce").fillna(0.0),
    })
    z = z[z.cdr3nt.ne("nan") & z.cdr3nt.ne("") & z.cdr3aa.ne("nan") & z.cdr3aa.ne("") & z["count"].gt(0)]
    frames.append(z)

if missing_nt:
    raise RuntimeError(
        "Nucleotide CDR3 is absent from accepted MiXCR exports: " + ", ".join(missing_nt[:20]) +
        ". Re-export the existing .clns files with nSeqCDR3; FASTQ analysis does not need to be repeated."
    )

lib = pd.concat(frames, ignore_index=True)
group_cols = ["library_id", "patient_id", "clinical_timepoint", "disease_status", "age_group", "age_years", "batch"] + KEY_NT
lib = lib.groupby(group_cols, as_index=False, dropna=False)["count"].sum()
totals = lib.groupby("library_id")["count"].sum().rename("library_total")
lib = lib.merge(totals, on="library_id", how="left")
lib["library_fraction"] = lib["count"] / lib["library_total"]
lib.to_csv(RESULTS / "library_trb_clonotypes_nt.csv", index=False)

case = meta[meta.analysis_include_final & meta.disease_status.eq("measles")].copy()
paired = []
for patient, z in case.groupby("patient_id"):
    if {"T1", "T2"}.issubset(set(z.clinical_timepoint.astype(str))):
        paired.append(str(patient))
paired = sorted(paired)
if len(paired) != 22:
    raise RuntimeError(f"Expected 22 T1-T2 patients, found {len(paired)}")


def libs_for(patient, point):
    return sorted(case.loc[case.patient_id.eq(patient) & case.clinical_timepoint.eq(point), "library_id"].astype(str))


def point_table(patient, point):
    libs = libs_for(patient, point)
    if not libs:
        return pd.DataFrame(columns=KEY_NT + ["fraction", "raw_count_sum", "read_equivalents", "total_detection_reads"])
    z = lib[lib.library_id.isin(libs)]
    a = z.groupby(KEY_NT, as_index=False, dropna=False).agg(sum_fraction=("library_fraction", "sum"), raw_count_sum=("count", "sum"))
    a["fraction"] = a.sum_fraction / len(libs)
    norm = a.fraction.sum()
    if norm > 0:
        a["fraction"] /= norm
    mean_reads = float(totals.reindex(libs).mean())
    a["read_equivalents"] = a.fraction * mean_reads
    a["total_detection_reads"] = float(totals.reindex(libs).sum())
    a["n_libraries"] = len(libs)
    return a


points = {(p, t): point_table(p, t) for p in paired for t in ["T1", "T2", "T3"]}

# Add all accepted controls and unpaired observations for public-clone prevalence.
all_points = dict(points)
for r in meta[meta.analysis_include_final].drop_duplicates(["patient_id", "clinical_timepoint"]).itertuples(index=False):
    key = (str(r.patient_id), str(r.clinical_timepoint))
    if key not in all_points:
        # point_table uses case metadata, so build directly for non-case observations.
        ll = sorted(meta.loc[meta.analysis_include_final & meta.patient_id.eq(key[0]) & meta.clinical_timepoint.eq(key[1]), "library_id"].astype(str))
        if ll:
            z = lib[lib.library_id.isin(ll)]
            a = z.groupby(KEY_NT, as_index=False, dropna=False).agg(sum_fraction=("library_fraction", "sum"), raw_count_sum=("count", "sum"))
            a["fraction"] = a.sum_fraction / len(ll)
            n = a.fraction.sum()
            if n > 0:
                a["fraction"] /= n
            a["read_equivalents"] = a.fraction * float(totals.reindex(ll).mean())
            a["total_detection_reads"] = float(totals.reindex(ll).sum())
            a["n_libraries"] = len(ll)
            all_points[key] = a

# Amino-acid public/private prevalence.
aa_people = defaultdict(set)
aa_controls = defaultdict(set)
healthy = set(meta.loc[meta.analysis_include_final & meta.disease_status.eq("healthy"), "patient_id"].astype(str))
for (patient, _), z in all_points.items():
    if z.empty:
        continue
    for rr in z[["cdr3aa", "v_gene", "j_gene"]].drop_duplicates().itertuples(index=False):
        key = (str(rr.cdr3aa), str(rr.v_gene), str(rr.j_gene))
        aa_people[key].add(patient)
        if patient in healthy:
            aa_controls[key].add(patient)

# Existing exact antigen annotations, when present.
antigen = {}
antigen_path = RESULTS / "antigen_specific_clonotypes.csv"
if antigen_path.exists() and antigen_path.stat().st_size > 0:
    a = pd.read_csv(antigen_path)
    if {"cdr3aa", "pathogen"}.issubset(a.columns):
        for aa, z in a.dropna(subset=["cdr3aa", "pathogen"]).groupby("cdr3aa"):
            antigen[str(aa)] = "; ".join(sorted(set(z.pathogen.astype(str))))

rows = []
for patient in paired:
    t1 = points[(patient, "T1")]
    t2 = points[(patient, "T2")]
    t3 = points[(patient, "T3")]
    t2i = t2.set_index(KEY_NT, drop=False)
    t3i = t3.set_index(KEY_NT, drop=False) if not t3.empty else None
    pm = case[case.patient_id.eq(patient)]
    age_group = str(pm.age_group.iloc[0])
    age = float(pm.age_years.dropna().iloc[0]) if pm.age_years.notna().any() else np.nan
    b1 = sorted(set(pm.loc[pm.clinical_timepoint.eq("T1"), "batch"].astype(str)))
    b2 = sorted(set(pm.loc[pm.clinical_timepoint.eq("T2"), "batch"].astype(str)))
    t1_mass = float(t1.loc[t1.fraction.ge(PRIMARY), "fraction"].sum())
    N2 = float(t2.total_detection_reads.iloc[0])
    zero_upper = 1 - 0.05 ** (1 / max(N2, 1))

    for r in t1.itertuples(index=False):
        key = (str(r.cdr3nt), str(r.cdr3aa), str(r.v_gene), str(r.j_gene))
        f1 = float(r.fraction)
        if key in t2i.index:
            q = t2i.loc[key]
            if isinstance(q, pd.DataFrame): q = q.iloc[0]
            f2, c2, re2 = float(q.fraction), float(q.raw_count_sum), float(q.read_equivalents)
        else:
            f2 = c2 = re2 = 0.0
        if f1 >= PRIMARY:
            fate = "retained expanded" if f2 >= PRIMARY else ("detected below 0.10 percent" if f2 > 0 else "not detected")
        else:
            fate = "T1 below primary candidate threshold"

        pdet = 1 - (1 - min(f1, 1 - 1e-15)) ** int(round(N2)) if f1 > 0 else 0.0
        plower = float(binom.cdf(int(round(c2)), int(round(N2)), min(max(f1, 0), 1)))
        if f2 > 0:
            log2ratio = math.log2(f2 / f1)
            fold = f1 / f2
            fold_lower = np.nan
        else:
            log2ratio = fold = np.nan
            fold_lower = f1 / zero_upper if zero_upper > 0 else np.inf

        aa_key = (str(r.cdr3aa), str(r.v_gene), str(r.j_gene))
        if t3i is None:
            f3, fate3 = np.nan, "not available"
        elif key in t3i.index:
            q3 = t3i.loc[key]
            if isinstance(q3, pd.DataFrame): q3 = q3.iloc[0]
            f3 = float(q3.fraction)
            fate3 = "expanded at T3" if f3 >= PRIMARY else "detected below 0.10 percent at T3"
        else:
            f3, fate3 = 0.0, "not detected at T3"

        rows.append({
            "patient": patient,
            "age_group": age_group,
            "age_years": age,
            "clone_id_nt": "|".join([str(r.cdr3nt), str(r.v_gene), str(r.j_gene)]),
            "clone_id_aa": "|".join([str(r.cdr3aa), str(r.v_gene), str(r.j_gene)]),
            "cdr3nt": str(r.cdr3nt), "cdr3aa": str(r.cdr3aa), "TRBV": str(r.v_gene), "TRBJ": str(r.j_gene),
            "T1_raw_count_sum": float(r.raw_count_sum), "T2_raw_count_sum": c2,
            "T1_read_equivalents": float(r.read_equivalents), "T2_read_equivalents": re2,
            "T1_fraction": f1, "T2_fraction": f2, "T1_percent": 100*f1, "T2_percent": 100*f2,
            "expanded_T1_0_10_percent": bool(f1 >= PRIMARY), "expanded_T2_0_10_percent": bool(f2 >= PRIMARY),
            "fate_0_10_percent": fate, "absolute_delta_fraction": f2-f1,
            "absolute_delta_percentage_points": 100*(f2-f1), "log2_T2_over_T1": log2ratio,
            "fold_contraction_T1_over_T2": fold, "fold_contraction_lower_bound_if_T2_zero": fold_lower,
            "T1_expanded_mass_share": f1/t1_mass if f1 >= PRIMARY and t1_mass > 0 else np.nan,
            "T2_zero_detection_upper_95_fraction": zero_upper if f2 == 0 else np.nan,
            "T2_zero_detection_upper_95_percent": 100*zero_upper if f2 == 0 else np.nan,
            "P_detection_if_unchanged": pdet, "P_observed_or_lower_if_unchanged": plower,
            "sampling_robust_not_detected": bool(f2 == 0 and pdet >= ROBUST_P),
            "public_private": "public" if len(aa_people[aa_key]) >= 2 else "private",
            "n_patients_with_aa_clone": len(aa_people[aa_key]),
            "control_prevalence": len(aa_controls[aa_key]) / len(healthy) if healthy else np.nan,
            "known_antigen": antigen.get(str(r.cdr3aa), ""), "HLA_compatible": "not available",
            "T3_fraction": f3, "T3_fate": fate3,
            "T1_batches": ";".join(b1), "T2_batches": ";".join(b2), "cross_batch_pair": b1 != b2,
        })

master = pd.DataFrame(rows)
master["sampling_contraction_fdr"] = fdr(master.P_observed_or_lower_if_unchanged)

# Threshold sensitivity.
sens_rows = []
for patient in paired:
    t1 = points[(patient, "T1")].set_index(KEY_NT, drop=False)
    t2 = points[(patient, "T2")].set_index(KEY_NT, drop=False)
    age_group = str(case.loc[case.patient_id.eq(patient), "age_group"].iloc[0])
    for threshold in THRESHOLDS:
        e1 = set(t1.index[t1.fraction.ge(threshold)])
        e2 = set(t2.index[t2.fraction.ge(threshold)])
        retained = e1 & e2
        below = {k for k in e1-e2 if k in t2.index}
        lost = {k for k in e1 if k not in t2.index}
        mass = float(t1.loc[list(e1), "fraction"].sum()) if e1 else 0.0
        def mass_pct(keys):
            return 100*float(t1.loc[list(keys), "fraction"].sum())/mass if keys and mass else 0.0
        union = e1 | e2
        sens_rows.append({
            "patient": patient, "age_group": age_group, "threshold_fraction": threshold, "threshold_percent": 100*threshold,
            "T1_expanded_count": len(e1), "T2_expanded_count": len(e2), "retained_count": len(retained),
            "below_threshold_count": len(below), "not_detected_count": len(lost),
            "expanded_jaccard": len(retained)/len(union) if union else np.nan,
            "retained_T1_mass_percent": mass_pct(retained), "below_T1_mass_percent": mass_pct(below),
            "lost_T1_mass_percent": mass_pct(lost),
        })
sensitivity = pd.DataFrame(sens_rows)
sensitivity.to_csv(CURRENT / "clone_fate_threshold_sensitivity.csv", index=False)

# Repeated equal-depth resampling. Clone-level detection probability above remains analytical.
paired_libs = case[case.patient_id.isin(paired) & case.clinical_timepoint.isin(["T1","T2"])].library_id.astype(str)
min_depth = int(totals.reindex(paired_libs).min())
depth = min(REQUESTED_DEPTH, min_depth)
rng = np.random.default_rng(SEED)
rare_rows = []
for patient in paired:
    a, b = points[(patient,"T1")], points[(patient,"T2")]
    keys = sorted(set(map(tuple,a[KEY_NT].to_numpy())) | set(map(tuple,b[KEY_NT].to_numpy())))
    amap = {tuple(x[k] for k in KEY_NT): float(x.fraction) for x in a.itertuples(index=False) for k in []}
    # Explicit dict comprehensions avoid dependence on namedtuple column aliases.
    amap = {(str(r.cdr3nt),str(r.cdr3aa),str(r.v_gene),str(r.j_gene)): float(r.fraction) for r in a.itertuples(index=False)}
    bmap = {(str(r.cdr3nt),str(r.cdr3aa),str(r.v_gene),str(r.j_gene)): float(r.fraction) for r in b.itertuples(index=False)}
    pa = np.array([amap.get(tuple(k),0.0) for k in keys], float); pa /= pa.sum()
    pb = np.array([bmap.get(tuple(k),0.0) for k in keys], float); pb /= pb.sum()
    age_group = str(case.loc[case.patient_id.eq(patient),"age_group"].iloc[0])
    store = defaultdict(list)
    for _ in range(N_ITER):
        ca = rng.multinomial(depth, pa); cb = rng.multinomial(depth, pb)
        fa = ca/depth; fb = cb/depth
        p1 = fa[fa>0]; p2 = fb[fb>0]
        sh1 = float(-(p1*np.log(p1)).sum()); sh2 = float(-(p2*np.log(p2)).sum())
        store["T1_richness"].append(len(p1)); store["T2_richness"].append(len(p2))
        store["T1_shannon"].append(sh1); store["T2_shannon"].append(sh2)
        store["T1_clonality"].append(1-sh1/np.log(len(p1)) if len(p1)>1 else np.nan)
        store["T2_clonality"].append(1-sh2/np.log(len(p2)) if len(p2)>1 else np.nan)
        all1=set(np.flatnonzero(ca>0)); all2=set(np.flatnonzero(cb>0)); u=all1|all2
        store["complete_jaccard"].append(len(all1&all2)/len(u) if u else np.nan)
        for threshold in [PRIMARY,0.01]:
            e1=set(np.flatnonzero(fa>=threshold)); e2=set(np.flatnonzero(fb>=threshold)); eu=e1|e2
            store[f"expanded_jaccard_{100*threshold:g}pct"].append(len(e1&e2)/len(eu) if eu else np.nan)
    rec={"patient":patient,"age_group":age_group,"rarefaction_depth":depth,"iterations":N_ITER}
    for name, values in store.items():
        v=np.asarray(values,float); v=v[np.isfinite(v)]
        rec[name+"_median"]=float(np.median(v)) if len(v) else np.nan
        rec[name+"_q025"]=float(np.quantile(v,0.025)) if len(v) else np.nan
        rec[name+"_q975"]=float(np.quantile(v,0.975)) if len(v) else np.nan
    rare_rows.append(rec)
rarefaction = pd.DataFrame(rare_rows)
rarefaction.to_csv(CURRENT / "rarefaction_summary.csv", index=False)

# Simple amino-acid sequence-neighbour clusters for primary candidates. This is not TCRdist.
def edit_distance(a,b,limit=2):
    if abs(len(a)-len(b))>limit: return limit+1
    prev=list(range(len(b)+1))
    for i,ca in enumerate(a,1):
        cur=[i]
        for j,cb in enumerate(b,1):
            cur.append(min(cur[-1]+1,prev[j]+1,prev[j-1]+(ca!=cb)))
        if min(cur)>limit: return limit+1
        prev=cur
    return prev[-1]

cand = master[master.expanded_T1_0_10_percent].copy().reset_index(drop=True)
parent=list(range(len(cand)))
def find(i):
    while parent[i]!=i:
        parent[i]=parent[parent[i]]; i=parent[i]
    return i
def union(i,j):
    a,b=find(i),find(j)
    if a!=b: parent[b]=a
for _, ids in cand.groupby("TRBV").groups.items():
    ids=list(ids)
    for ii in range(len(ids)):
        for jj in range(ii+1,len(ids)):
            i,j=ids[ii],ids[jj]; a=str(cand.at[i,"cdr3aa"]); b=str(cand.at[j,"cdr3aa"])
            if abs(len(a)-len(b))<=1 and edit_distance(a,b,MAX_EDIT)<=MAX_EDIT: union(i,j)
cluster_map={}; cluster=[]
for i in range(len(cand)):
    root=find(i)
    if root not in cluster_map: cluster_map[root]=f"SEQ{len(cluster_map)+1:04d}"
    cluster.append(cluster_map[root])
cand["motif_cluster"]=cluster

# Exact-clone disappearance with a nearby T2 amino-acid sequence is flagged as possible clonotypic replacement.
neighbour=[]
for r in cand.itertuples(index=False):
    t2=points[(str(r.patient),"T2")]
    z=t2[t2.v_gene.astype(str).eq(str(r.TRBV))]
    found=False
    for aa in z.cdr3aa.astype(str):
        if aa!=str(r.cdr3aa) and abs(len(aa)-len(str(r.cdr3aa)))<=1 and edit_distance(aa,str(r.cdr3aa),MAX_EDIT)<=MAX_EDIT:
            found=True; break
    neighbour.append(found)
cand["T2_sequence_neighbour_present"]=neighbour
master=master.merge(cand[["patient","clone_id_nt","motif_cluster","T2_sequence_neighbour_present"]],on=["patient","clone_id_nt"],how="left")
master["motif_cluster"]=master.motif_cluster.fillna("")
master["T2_sequence_neighbour_present"]=master.T2_sequence_neighbour_present.fillna(False)

def candidate_class(r):
    if not r.expanded_T1_0_10_percent: return "outside primary T1-expanded pool"
    if r.fate_0_10_percent=="not detected" and r.sampling_robust_not_detected: return "sampling-robust loss"
    if pd.notna(r.fold_contraction_T1_over_T2) and r.fold_contraction_T1_over_T2>=10: return "strong contraction"
    if r.fate_0_10_percent=="detected below 0.10 percent": return "contracted below expansion threshold"
    if r.fate_0_10_percent=="retained expanded": return "retained expanded"
    return "other"
master["candidate_class"]=master.apply(candidate_class,axis=1)
master.to_csv(CURRENT / "paired_T1_clone_master_nt.csv",index=False)

clusters=cand.merge(master[["patient","clone_id_nt","candidate_class"]],on=["patient","clone_id_nt"],how="left").groupby("motif_cluster",as_index=False).agg(n_clones=("clone_id_nt","nunique"),n_patients=("patient","nunique"),n_robust_loss=("candidate_class",lambda x:int((x=="sampling-robust loss").sum())),n_strong_contraction=("candidate_class",lambda x:int((x=="strong contraction").sum())))
clusters["loss_or_strong_fraction"]=(clusters.n_robust_loss+clusters.n_strong_contraction)/clusters.n_clones
clusters.sort_values(["n_patients","n_clones"],ascending=False).to_csv(CURRENT / "sequence_neighbour_cluster_summary.csv",index=False)

# T3 validation of T1-expanded clones.
t3=[]
for patient,z in master[master.expanded_T1_0_10_percent].groupby("patient"):
    if z.T3_fate.eq("not available").all(): continue
    lost=z.fate_0_10_percent.eq("not detected")
    t3.append({"patient":patient,"n_T1_expanded":len(z),"n_not_detected_T2":int(lost.sum()),"n_T2_lost_reappeared_T3":int((lost & z.T3_fraction.fillna(0).gt(0)).sum()),"n_T2_lost_still_absent_T3":int((lost & z.T3_fraction.fillna(0).eq(0)).sum())})
t3=pd.DataFrame(t3); t3.to_csv(CURRENT / "T3_clone_fate_validation.csv",index=False)

# Morisita-Horn identity fingerprint for paired observations.
def aa_freq(z):
    if z.empty: return {}
    a=z.groupby(["cdr3aa","v_gene","j_gene"],as_index=False).fraction.sum()
    return {(str(r.cdr3aa),str(r.v_gene),str(r.j_gene)):float(r.fraction) for r in a.itertuples(index=False)}
def mh(p,q):
    keys=set(p)|set(q)
    x=np.array([p.get(k,0.0) for k in keys]); y=np.array([q.get(k,0.0) for k in keys])
    d=(x*x).sum()+(y*y).sum()
    return float(2*np.dot(x,y)/d) if d>0 else np.nan
obs=[]
for p in paired:
    for t in ["T1","T2","T3"]:
        f=aa_freq(points[(p,t)])
        if f: obs.append((f"{p} {t}",p,t,f))
labels=[x[0] for x in obs]; mat=np.zeros((len(obs),len(obs)))
for i,x in enumerate(obs):
    for j,y in enumerate(obs): mat[i,j]=mh(x[3],y[3])
pd.DataFrame(mat,index=labels,columns=labels).to_csv(CURRENT / "repertoire_identity_morisita_matrix.csv")
ids=[]; t2idx=[i for i,x in enumerate(obs) if x[2]=="T2"]
for i,x in enumerate(obs):
    if x[2]!="T1": continue
    own=[j for j in t2idx if obs[j][1]==x[1]]
    if not own: continue
    ownsim=mat[i,own[0]]; rank=1+sum(mat[i,j]>ownsim for j in t2idx if obs[j][1]!=x[1])
    ids.append({"patient":x[1],"own_T1_T2_morisita":ownsim,"own_T2_rank_among_T2":rank,"identity_pass_rank1":rank==1})
pd.DataFrame(ids).to_csv(CURRENT / "repertoire_identity_qc.csv",index=False)

# Explicit cross-batch sensitivity for GA1_003.
primary=sensitivity[sensitivity.threshold_fraction.eq(PRIMARY)]
ga=[]
for age in ["child","adult"]:
    for name,z in [("all paired patients",primary[primary.age_group.eq(age)]),("excluding GA1_003",primary[primary.age_group.eq(age)&~primary.patient.eq("ГА1_003")])]:
        ga.append({"age_group":age,"analysis_set":name,"n_patients":z.patient.nunique(),"mean_lost_T1_mass_percent":z.lost_T1_mass_percent.mean(),"median_lost_T1_mass_percent":z.lost_T1_mass_percent.median(),"mean_expanded_jaccard":z.expanded_jaccard.mean()})
pd.DataFrame(ga).to_csv(CURRENT / "GA1_003_cross_batch_sensitivity.csv",index=False)

# Compact overview figure. Existing figures 01-13 and patient panels remain untouched.
fig,axs=plt.subplots(2,2,figsize=(14,10),constrained_layout=True)
for age,label in [("child","Children younger than 18 years"),("adult","Adults 18 years and older")]:
    g=sensitivity[sensitivity.age_group.eq(age)].groupby("threshold_percent").lost_T1_mass_percent.mean()
    axs[0,0].plot(g.index,g.values,marker="o",label=label)
axs[0,0].set(title="Apparent loss across expansion thresholds",xlabel="Expansion threshold, percent",ylabel="Mean percent of T1 expanded mass not detected at T2"); axs[0,0].legend(fontsize=8); axs[0,0].grid(alpha=.2)
rob=master[master.expanded_T1_0_10_percent].groupby(["patient","age_group"],as_index=False).agg(loss=("sampling_robust_not_detected","sum"),n=("clone_id_nt","size")); rob["pct"]=100*rob.loss/rob.n
for x,age in enumerate(["child","adult"]):
    v=rob.loc[rob.age_group.eq(age),"pct"].to_numpy(float); axs[0,1].scatter(np.full(len(v),x)+np.linspace(-.08,.08,len(v)),v,alpha=.75); axs[0,1].hlines(np.median(v),x-.18,x+.18,lw=2)
axs[0,1].set_xticks([0,1],["Children","Adults"]); axs[0,1].set(title="Sampling-robust disappearance",ylabel="Percent of T1 expanded clones"); axs[0,1].grid(axis="y",alpha=.2)
raw=primary[["patient","age_group","expanded_jaccard"]].merge(rarefaction[["patient","expanded_jaccard_0.1pct_median"]],on="patient")
for age,label in [("child","Children"),("adult","Adults")]:
    z=raw[raw.age_group.eq(age)]; axs[1,0].scatter(z.expanded_jaccard,z["expanded_jaccard_0.1pct_median"],label=label,alpha=.75)
mx=max(.05,np.nanmax(raw[["expanded_jaccard","expanded_jaccard_0.1pct_median"]].to_numpy())); axs[1,0].plot([0,mx],[0,mx],ls="--",lw=1); axs[1,0].set(title=f"Expanded overlap after depth control to {depth:,} reads",xlabel="Observed Jaccard",ylabel="Median rarefied Jaccard"); axs[1,0].legend(fontsize=8); axs[1,0].grid(alpha=.2)
rank=master[master.expanded_T1_0_10_percent].copy(); rank["score"]=rank.T1_fraction*np.where(rank.candidate_class.eq("sampling-robust loss"),2,1); rank=rank.sort_values("score",ascending=False).head(12).iloc[::-1]
axs[1,1].barh(np.arange(len(rank)),100*rank.T1_fraction); axs[1,1].set_yticks(np.arange(len(rank)),[f"{r.patient} {r.cdr3aa[:14]} {r.TRBV}" for r in rank.itertuples(index=False)],fontsize=7); axs[1,1].set(title="Largest T1 expanded candidates",xlabel="T1 repertoire fraction, percent"); axs[1,1].grid(axis="x",alpha=.2)
fig.suptitle("Nucleotide-level T1 clone fate with sampling-depth controls",fontweight="bold",fontsize=15)
fig.savefig(FIG/"14_clone_fate_v2_overview.png",dpi=220,bbox_inches="tight"); fig.savefig(FIG/"14_clone_fate_v2_overview.pdf",bbox_inches="tight"); plt.close(fig)

if len(labels):
    fig,ax=plt.subplots(figsize=(12,11)); im=ax.imshow(mat,vmin=0,vmax=1,aspect="auto"); ax.set_xticks(np.arange(len(labels)),labels,rotation=90,fontsize=5); ax.set_yticks(np.arange(len(labels)),labels,fontsize=5); ax.set_title("Repertoire identity fingerprint across paired observations\nMorisita-Horn similarity on CDR3aa plus TRBV plus TRBJ frequencies"); fig.colorbar(im,ax=ax,label="Morisita-Horn similarity"); fig.tight_layout(); fig.savefig(FIG/"17_repertoire_identity_fingerprint.png",dpi=220,bbox_inches="tight"); fig.savefig(FIG/"17_repertoire_identity_fingerprint.pdf",bbox_inches="tight"); plt.close(fig)

if not t3.empty:
    fig,ax=plt.subplots(figsize=(8,5.5)); x=np.arange(len(t3)); ax.bar(x,t3.n_T2_lost_still_absent_T3,label="Still absent at T3"); ax.bar(x,t3.n_T2_lost_reappeared_T3,bottom=t3.n_T2_lost_still_absent_T3,label="Reappeared at T3"); ax.set_xticks(x,t3.patient,rotation=30,ha="right"); ax.set_ylabel("T1 expanded clones not detected at T2"); ax.set_title("T3 follow-up of clones apparently lost at T2"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(FIG/"18_T3_clone_fate_validation.png",dpi=220,bbox_inches="tight"); fig.savefig(FIG/"18_T3_clone_fate_validation.pdf",bbox_inches="tight"); plt.close(fig)

(CURRENT/"clone_fate_v2_README.txt").write_text("\n".join([
    f"Paired T1-T2 patients: {len(paired)}",
    "Physical tracking key: CDR3 nucleotide + TRBV + TRBJ",
    "Functional/public key: CDR3 amino acid + TRBV + TRBJ",
    "Primary expanded threshold: 0.10 percent",
    "Sensitivity thresholds percent: 0.05, 0.10, 0.20, 0.50, 0.75, 1.00",
    f"Repeated depth-control: {N_ITER} multinomial resamples at {depth} productive reads",
    "Clone-level P_detection_if_unchanged is analytical, not estimated by brute-force iterations",
    "T2 zero values use a 95 percent detection upper bound and no arbitrary pseudocount",
    "Sequence-neighbour clusters use amino-acid edit distance and are not claimed to be TCRdist clusters",
    "HLA compatibility is not available without patient HLA typing",
    "Sampling controls do not remove PCR amplification bias because there are no UMI",
    "Existing aaV analyses were preserved",
])+"\n",encoding="utf-8")

print("Clone fate v2 completed")
print("Paired patients",len(paired))
print("T1 nucleotide clonotypes",len(master))
print("Primary T1 expanded candidates",int(master.expanded_T1_0_10_percent.sum()))
print("Rarefaction depth",depth,"iterations",N_ITER)
