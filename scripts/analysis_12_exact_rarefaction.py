#!/usr/bin/env python3
"""Exact read downsampling for the 22 paired T1-T2 patients.

This module is additive. It performs true without-replacement rarefaction of each
sequencing library to a common productive-read depth. Biological replicate libraries
are rarefied separately and only then averaged at the patient time-point level.

Outputs include Shannon diversity, the established entropy-based clonality, Simpson D,
sqrt Simpson D, normalized richness, complete-repertoire Jaccard and expanded-repertoire
Jaccard at 0.10 and 1.00 percent.

Read downsampling controls sequencing depth only. It cannot remove PCR amplification
bias because these bulk genomic DNA libraries do not contain UMI information.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("MEASLES_REPO", Path.cwd())).resolve()
RESULTS = ROOT / "results"
CURRENT = RESULTS / "current"
FIG = ROOT / "figures" / "current"
CURRENT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

REQUESTED_DEPTH = int(os.environ.get("MEASLES_RAREFACTION_DEPTH", "50000"))
N_ITER = int(os.environ.get("MEASLES_RAREFACTION_ITERATIONS", "500"))
SEED = int(os.environ.get("MEASLES_RAREFACTION_SEED", "1031"))
THRESHOLDS = [0.001, 0.01]
KEY = ["cdr3nt", "cdr3aa", "v_gene", "j_gene"]

nt_path = RESULTS / "library_trb_clonotypes_nt.csv"
meta_path = CURRENT / "sample_metadata_current.csv"
if not nt_path.exists():
    subprocess.run([sys.executable, str(ROOT / "scripts" / "analysis_10_clone_fate_v2.py")], check=True)
lib = pd.read_csv(nt_path)
meta = pd.read_csv(meta_path)
meta["analysis_include_final"] = meta.analysis_include_final.astype(str).str.lower().isin(["true","1","yes"])
lib["count"] = pd.to_numeric(lib["count"], errors="coerce").fillna(0).round().astype(np.int64)

case = meta[meta.analysis_include_final & meta.disease_status.eq("measles")].copy()
paired = []
for patient, z in case.groupby("patient_id"):
    if {"T1","T2"}.issubset(set(z.clinical_timepoint.astype(str))):
        paired.append(str(patient))
paired = sorted(paired)
if len(paired) != 22:
    raise RuntimeError(f"Expected 22 T1-T2 patients, found {len(paired)}")

paired_libs = case[case.patient_id.astype(str).isin(paired) & case.clinical_timepoint.isin(["T1","T2"])].library_id.astype(str)
lib_totals = lib.groupby("library_id")["count"].sum()
minimum_observed_depth = int(lib_totals.reindex(paired_libs).min())
depth = min(REQUESTED_DEPTH, minimum_observed_depth)
if depth <= 0:
    raise RuntimeError("Common rarefaction depth is not positive")

# Precompute integer count vectors once. A global key space is used only inside each
# patient, so the sampling work remains compact.
patient_data = {}
for patient in paired:
    pmeta = case[case.patient_id.astype(str).eq(patient) & case.clinical_timepoint.isin(["T1","T2"])]
    library_ids = pmeta.library_id.astype(str).tolist()
    z = lib[lib.library_id.astype(str).isin(library_ids)].copy()
    keys = sorted(set(map(tuple, z[KEY].astype(str).to_numpy())))
    key_to_i = {k:i for i,k in enumerate(keys)}
    library_vectors = {}
    for library_id, zz in z.groupby("library_id"):
        vec = np.zeros(len(keys), dtype=np.int64)
        for r in zz.itertuples(index=False):
            k = (str(r.cdr3nt),str(r.cdr3aa),str(r.v_gene),str(r.j_gene))
            vec[key_to_i[k]] += int(r.count)
        library_vectors[str(library_id)] = vec
    point_libraries = {
        point: pmeta.loc[pmeta.clinical_timepoint.eq(point),"library_id"].astype(str).tolist()
        for point in ["T1","T2"]
    }
    age_group = str(pmeta.age_group.iloc[0])
    patient_data[patient] = (keys, library_vectors, point_libraries, age_group)

rng = np.random.default_rng(SEED)


def rarefy(vec):
    total = int(vec.sum())
    if total < depth:
        raise RuntimeError(f"Library depth {total} is below common rarefaction depth {depth}")
    if total == depth:
        return vec.copy()
    # NumPy's multivariate hypergeometric samples observed reads without replacement.
    return rng.multivariate_hypergeometric(vec, depth).astype(np.int64)


def point_frequency(libraries, vectors):
    fractions = []
    for lib_id in libraries:
        sampled = rarefy(vectors[lib_id])
        fractions.append(sampled / float(depth))
    if not fractions:
        raise RuntimeError("A paired patient time point has no library")
    f = np.mean(np.stack(fractions, axis=0), axis=0)
    s = float(f.sum())
    return f/s if s > 0 else f


def metrics(f1, f2):
    out = {}
    for name, f in [("T1",f1),("T2",f2)]:
        p = f[f>0]
        richness = int(len(p))
        shannon = float(-(p*np.log(p)).sum()) if richness else np.nan
        out[f"{name}_richness"] = richness
        out[f"{name}_normalized_richness"] = richness / float(depth)
        out[f"{name}_shannon"] = shannon
        out[f"{name}_clonality"] = float(1-shannon/np.log(richness)) if richness>1 else np.nan
        simpson_d = float(np.sum(p*p)) if richness else np.nan
        out[f"{name}_simpson_D"] = simpson_d
        out[f"{name}_sqrt_simpson_D"] = float(np.sqrt(simpson_d)) if np.isfinite(simpson_d) else np.nan
    s1=set(np.flatnonzero(f1>0)); s2=set(np.flatnonzero(f2>0)); u=s1|s2
    out["complete_jaccard"] = len(s1&s2)/len(u) if u else np.nan
    for th in THRESHOLDS:
        e1=set(np.flatnonzero(f1>=th)); e2=set(np.flatnonzero(f2>=th)); eu=e1|e2
        label=str(100*th).rstrip("0").rstrip(".").replace(".","_")
        out[f"expanded_jaccard_{label}_percent"] = len(e1&e2)/len(eu) if eu else np.nan
        out[f"T1_expanded_count_{label}_percent"] = len(e1)
        out[f"T2_expanded_count_{label}_percent"] = len(e2)
    return out

raw_rows=[]
summary_rows=[]
for patient in paired:
    keys, vectors, point_libraries, age_group = patient_data[patient]
    store=defaultdict(list)
    for iteration in range(N_ITER):
        f1=point_frequency(point_libraries["T1"],vectors)
        f2=point_frequency(point_libraries["T2"],vectors)
        rec=metrics(f1,f2)
        for k,v in rec.items(): store[k].append(v)
        if iteration < 5:
            # Small audit sample only. The full iteration matrix would be unnecessarily large.
            raw_rows.append({"patient":patient,"age_group":age_group,"iteration":iteration+1,"depth":depth,**rec})
    row={"patient":patient,"age_group":age_group,"rarefaction_depth":depth,"iterations":N_ITER,"n_T1_libraries":len(point_libraries["T1"]),"n_T2_libraries":len(point_libraries["T2"])}
    for name,values in store.items():
        a=np.asarray(values,float); a=a[np.isfinite(a)]
        row[name+"_median"]=float(np.median(a)) if len(a) else np.nan
        row[name+"_q025"]=float(np.quantile(a,.025)) if len(a) else np.nan
        row[name+"_q975"]=float(np.quantile(a,.975)) if len(a) else np.nan
    summary_rows.append(row)

summary=pd.DataFrame(summary_rows)
summary.to_csv(CURRENT/"exact_rarefaction_summary.csv",index=False)
pd.DataFrame(raw_rows).to_csv(CURRENT/"exact_rarefaction_audit_first5.csv",index=False)

# Compact figure: each patient remains one biological replicate. Lines connect paired T1/T2.
fig,axs=plt.subplots(2,2,figsize=(13,10),constrained_layout=True)
metric_specs=[
    ("clonality","Entropy based clonality"),
    ("shannon","Shannon diversity"),
    ("simpson_D","Simpson D"),
    ("normalized_richness","Normalized richness")
]
for ax,(metric,title) in zip(axs.flat,metric_specs):
    for age,x0 in [("child",0),("adult",3)]:
        z=summary[summary.age_group.eq(age)].copy()
        t1=z[f"T1_{metric}_median"].to_numpy(float); t2=z[f"T2_{metric}_median"].to_numpy(float)
        for i in range(len(z)):
            ax.plot([x0,x0+1],[t1[i],t2[i]],alpha=.35,lw=1)
            ax.scatter([x0,x0+1],[t1[i],t2[i]],s=22)
    ax.set_xticks([0,1,3,4],["Child T1","Child T2","Adult T1","Adult T2"],rotation=15)
    ax.set_title(title)
    ax.grid(axis="y",alpha=.2)
fig.suptitle(f"Paired repertoire metrics after exact rarefaction to {depth:,} productive reads\nEach library sampled without replacement, {N_ITER} iterations",fontweight="bold",fontsize=14)
fig.savefig(FIG/"19_exact_rarefaction_repertoire_metrics.png",dpi=220,bbox_inches="tight")
fig.savefig(FIG/"19_exact_rarefaction_repertoire_metrics.pdf",bbox_inches="tight")
plt.close(fig)

# Expanded overlap at 0.10 and 1.00 percent after exact depth control.
fig,axs=plt.subplots(1,2,figsize=(12,5.5),constrained_layout=True)
for ax,(col,title) in zip(axs,[
    ("expanded_jaccard_0_1_percent_median","Expanded Jaccard at 0.10 percent"),
    ("expanded_jaccard_1_percent_median","Expanded Jaccard at 1.00 percent"),
]):
    for x,age in enumerate(["child","adult"]):
        vals=summary.loc[summary.age_group.eq(age),col].dropna().to_numpy(float)
        jitter=np.linspace(-.07,.07,len(vals)) if len(vals)>1 else np.array([0.0])
        ax.scatter(np.full(len(vals),x)+jitter,vals,alpha=.75)
        if len(vals): ax.hlines(np.median(vals),x-.18,x+.18,lw=2)
    ax.set_xticks([0,1],["Children","Adults"]); ax.set_ylim(0,1); ax.set_ylabel("Jaccard similarity"); ax.set_title(title); ax.grid(axis="y",alpha=.2)
fig.suptitle("Expanded repertoire similarity after exact read-depth control",fontweight="bold",fontsize=14)
fig.savefig(FIG/"20_exact_rarefaction_expanded_overlap.png",dpi=220,bbox_inches="tight")
fig.savefig(FIG/"20_exact_rarefaction_expanded_overlap.pdf",bbox_inches="tight")
plt.close(fig)

(CURRENT/"exact_rarefaction_README.txt").write_text("\n".join([
    f"Common productive-read depth: {depth}",
    f"Requested depth: {REQUESTED_DEPTH}",
    f"Iterations: {N_ITER}",
    f"Seed: {SEED}",
    "Sampling method: NumPy multivariate hypergeometric without replacement",
    "Each sequencing library is rarefied separately",
    "Biological repeat libraries are converted to frequencies and averaged only after rarefaction",
    "Metrics: Shannon diversity, entropy-based clonality, Simpson D, sqrt Simpson D, normalized richness, complete Jaccard, expanded Jaccard at 0.10 and 1.00 percent",
    "This controls sequencing depth but not PCR amplification bias because UMI are absent",
])+"\n",encoding="utf-8")

print("Exact rarefaction completed")
print("Paired patients",len(summary))
print("Common depth",depth)
print("Iterations",N_ITER)
