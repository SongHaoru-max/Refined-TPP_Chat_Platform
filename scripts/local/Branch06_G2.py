# -*- coding: utf-8 -*-
"""
Branch06_G2 - Local reference implementation
-----------------------------------------------
Full TMT proteome and N-glycosite quantification QC / DeltaRm pipeline.

This reference version is intentionally kept as a directly executable local
Python script. It does NOT depend on an Agent or JSON input. The main goals of
this revision are:
1. keep the original scientific workflow recognizable;
2. make parameter semantics and statistical rules explicit;
3. standardize output organization and naming;
4. make later JSON / Agent integration straightforward without changing the
   scientific core.

Important definitions in this reference version
------------------------------------------------
- N-glycopeptide detection: an N residue is considered deamidated when one of
  its attached modification annotations contains a mass shift around +0.98 Da.
  Other modifications attached to the same N are allowed.
- Rm normalization: replicate-specific multiplicative median normalization.
  The three control-proteome Rm medians are aligned to a common reference
  center. The same correction factors are applied to protein and glycosite Rm.
- CV QC: CV is calculated from RAW Rm values, before normalization.
- DeltaRm effect-size threshold: applied to MEAN DeltaRm across biological
  replicates, not to every replicate separately.
- BH correction: user-selectable. If disabled, raw p-values are used for
  significance decisions.
- Paired tests compare normalized protein Rm and normalized glycosite Rm.
"""

import os
import re
import json
from datetime import datetime
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.interpolate import make_interp_spline
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
from matplotlib_venn import venn2, venn3
from upsetplot import from_contents, UpSet
import statsmodels.api as sm


plt.rcParams.update({'figure.dpi':300})

PIPELINE_NAME = "Branch06_G2"
PIPELINE_VERSION = "reference-1.0"
DELTA_RM_THRESHOLD_DEFAULT = 0.10
NORMALITY_ALPHA_DEFAULT = 0.05


def stage(message):
    """Print a consistent stage message for local execution."""
    print(f"\n[Branch06_G2] {message}")


def prompt_yes_no(prompt, default=True):
    """Read a YES/NO option while preserving simple local-script usage."""
    raw = input(prompt).strip().upper()
    if raw == "":
        return bool(default)
    if raw in {"Y", "YES"}:
        return True
    if raw in {"N", "NO"}:
        return False
    raise ValueError("Please answer YES or NO.")


def build_output_dirs(work_dir):
    """Create a stable output layout instead of scattering files in cwd."""
    root = os.path.join(work_dir, "Branch06_G2_output")
    tables = os.path.join(root, "tables")
    figures = os.path.join(root, "figures")
    metadata = os.path.join(root, "metadata")
    intermediate = os.path.join(root, "intermediate")
    for path in (root, tables, figures, metadata, intermediate):
        os.makedirs(path, exist_ok=True)
    return root, tables, figures, metadata, intermediate


def write_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def contains_target_mass_shift(modifications, target=0.98, tolerance=0.01):
    """
    Return True when any modification annotation contains a numeric mass shift
    close to +0.98 Da. This is deliberately more tolerant than exact string
    matching, so forms such as '(+42.01)(+0.98)' or 'Deamidated(+0.984)' can be
    recognized when the +0.98-like value is attached to the N residue.
    """
    for mod in modifications:
        for token in re.findall(r"[+-]?\d+(?:\.\d+)?", str(mod)):
            try:
                value = float(token)
            except ValueError:
                continue
            if abs(value - target) <= tolerance:
                return True
    return False


def compute_replicate_median_normalization(df, rm_cols):
    """
    Compute multiplicative factors that align replicate-specific Rm medians.

    For replicate i:
        M_i = median(Rm_i across control proteins)
        reference = exp(median(log(M_i)))
        CF_i = reference / M_i
        Rm_i_normalized = Rm_i * CF_i

    With three positive replicate medians, the log-median reference is
    numerically the middle replicate median, but the log-space formulation makes
    the multiplicative interpretation explicit.
    """
    medians = df[rm_cols].median().astype(float)
    if (medians <= 0).any() or (~np.isfinite(medians)).any():
        raise ValueError(f"Invalid Rm medians for normalization: {medians.to_dict()}")
    reference = float(np.exp(np.median(np.log(medians.values))))
    factors = reference / medians
    return medians, reference, factors


def coefficient_of_variation_from_raw_rm(df, rm_cols):
    """CV of raw Rm across biological replicates; normalization is not used."""
    mean_rm = df[rm_cols].mean(axis=1)
    sd_rm = df[rm_cols].std(axis=1, ddof=1)
    return sd_rm / mean_rm


def adjust_pvalues(pvalues, apply_bh):
    """Return the significance metric and its label (BH-FDR or raw p-value)."""
    pvalues = np.asarray(pvalues, dtype=float)
    if apply_bh:
        _, adjusted, _, _ = multipletests(pvalues, method="fdr_bh")
        return adjusted, "fdr_bh"
    return pvalues.copy(), "p_value"


PROTEIN_EXPORT_RENAME = {
    "Rm1": "protein_rm_rep1_raw",
    "Rm2": "protein_rm_rep2_raw",
    "Rm3": "protein_rm_rep3_raw",
    "Rm1_normalized": "protein_rm_rep1_norm",
    "Rm2_normalized": "protein_rm_rep2_norm",
    "Rm3_normalized": "protein_rm_rep3_norm",
    "adjusted_Rm_unmod": "protein_rm_norm_mean",
    "CV_Rm_filtered": "protein_rm_raw_cv",
    "CV_Rm": "protein_rm_raw_cv",
}

SITE_EXPORT_RENAME = {
    "Rm_glyco1": "glycosite_rm_rep1_raw",
    "Rm_glyco2": "glycosite_rm_rep2_raw",
    "Rm_glyco3": "glycosite_rm_rep3_raw",
    "CV_Rm_glyco_complete": "glycosite_rm_raw_cv",
    "CV_Rm_glyco_filtered": "glycosite_rm_raw_cv",
    "Rm_glyco1_norm": "glycosite_rm_rep1_norm",
    "Rm_glyco2_norm": "glycosite_rm_rep2_norm",
    "Rm_glyco3_norm": "glycosite_rm_rep3_norm",
    "adjusted_Rm_glyco": "glycosite_rm_norm_mean",
    "adjusted_Rm_unmod": "protein_rm_norm_mean",
    "Rm_unmod1_norm": "protein_rm_rep1_norm",
    "Rm_unmod2_norm": "protein_rm_rep2_norm",
    "Rm_unmod3_norm": "protein_rm_rep3_norm",
    "ΔRm_rep1": "delta_rm_rep1",
    "ΔRm_rep2": "delta_rm_rep2",
    "ΔRm_rep3": "delta_rm_rep3",
    "ΔRm": "delta_rm_mean",
    "Protein_CV": "protein_rm_raw_cv",
    "Rep_CV": "representative_raw_rm_cv",
}


def export_csv(df, path, rename_map=None, columns=None):
    """
    Export a table with stable machine-friendly derived-column names while
    keeping the internal dataframe close to the original implementation.
    """
    out = df.copy() if columns is None else df.loc[:, columns].copy()
    if rename_map:
        out = out.rename(columns=rename_map)
    duplicated = out.columns[out.columns.duplicated()].tolist()
    if duplicated:
        raise ValueError(f"Duplicate export columns after renaming: {duplicated}")
    out.to_csv(path, index=False)


# ------------------------------
# Helper functions
# ------------------------------

def normalize_path(path):
    return os.path.normpath(os.path.abspath(path))

def canonical_protein_id(value):
    if pd.isna(value):
        return ''
    text = str(value).strip()
    if not text:
        return ''

    text = text.split(';')[0].strip()
    if '|' in text:
        parts = [p.strip() for p in text.split('|') if p.strip()]
        if len(parts) >= 2 and parts[0].lower() in {'sp', 'tr', 'up', 'ref', 'gb', 'emb', 'dbj', 'gi'}:
            text = parts[1]
        elif parts:
            text = parts[0]

    text = text.split()[0].strip()
    return text

def detect_sample_blocks(df):
    """Detect Sample columns grouped by Sample N"""
    pattern = re.compile(r"Sample\s*(\d+)", flags=re.I)
    sample_cols = defaultdict(list)
    sample_order = []
    for col in df.columns:
        m = pattern.search(col)
        if m:
            idx = int(m.group(1))
            sample_cols[idx].append(col)
            if idx not in sample_order:
                sample_order.append(idx)
    sample_order.sort()
    return sample_order, sample_cols

def extract_reporter_fragment(colname):
    """Extract TMT reporter ion channel from column name"""
    m = re.search(r"(126|127|128|129|130|131)", colname)
    if m:
        return f"TMT-{m.group(1)}"
    return colname

def parse_peptide(peptide):
    n_positions = []
    bare_seq = ''
    i = 0
    while i < len(peptide):
        if peptide[i] == 'N':
            mods = []
            j = i + 1
            while j < len(peptide) and peptide[j] == '(':
                k = peptide.find(')', j)
                if k == -1:
                    break
                mods.append(peptide[j+1:k])
                j = k + 1
            if contains_target_mass_shift(mods, target=0.98, tolerance=0.01):
                n_positions.append(len(bare_seq))
            bare_seq += 'N'
            i = j
        elif peptide[i] == '(':
            k = peptide.find(')', i)
            i = k + 1
        else:
            bare_seq += peptide[i]
            i += 1
    return bare_seq, n_positions

# 5) Function to compute robust left/right sigma from percentiles
def compute_sigma(vals):
    """Compute p50, left_sigma, right_sigma for a given array."""
    if len(vals) == 0:
        return np.nan, np.nan, np.nan
    p15 = np.percentile(vals, 15.87)
    p50 = np.percentile(vals, 50.0)
    p84 = np.percentile(vals, 84.13)
    left_sigma = p50 - p15
    right_sigma = p84 - p50
    left_sigma = left_sigma if left_sigma > 0 else np.nan
    right_sigma = right_sigma if right_sigma > 0 else np.nan
    return p50, left_sigma, right_sigma


# 6) Function to compute p-values for one bin (one-tailed)
def compute_bin_pvals(vals, side="right"):
    """
    Robust one-tailed p-values within one bin.
    Use robust piecewise sigma, but avoid assigning all opposite-side values to p=1,
    which causes severe p-value ties and discrete FDR artifacts.
    """
    vals = np.asarray(vals, dtype=float)
    pvals = np.ones(vals.shape[0], dtype=float)

    valid_vals = vals[~np.isnan(vals)]
    if valid_vals.size == 0:
        return pvals

    p50, left_sigma, right_sigma = compute_sigma(valid_vals)

    # Fallback sigma for small or degenerate bins.
    fallback_sigma = np.nanstd(valid_vals, ddof=1)
    if not np.isfinite(fallback_sigma) or fallback_sigma <= 0:
        fallback_sigma = np.nanstd(valid_vals, ddof=0)
    if not np.isfinite(fallback_sigma) or fallback_sigma <= 0:
        fallback_sigma = 1.0

    for i, val in enumerate(vals):
        if np.isnan(val):
            continue

        if val >= p50:
            sigma = right_sigma if (np.isfinite(right_sigma) and right_sigma > 0) else fallback_sigma
        else:
            sigma = left_sigma if (np.isfinite(left_sigma) and left_sigma > 0) else fallback_sigma

        z = (val - p50) / sigma
        if side == "right":
            pv = norm.sf(z)
        elif side == "left":
            pv = norm.cdf(z)
        else:
            raise ValueError("side must be 'right' or 'left'")

        pvals[i] = float(np.clip(pv, 0.0, 1.0))

    return pvals

# 7) Function to compute bin indices given bin count
def get_bins(n_sites, bin_count):
    """Return list of index arrays for each bin. Last bin may merge if too small."""
    bin_size = int(round(n_sites / bin_count))
    bins = []
    i = 0
    while i < n_sites:
        j = i + bin_size
        if j >= n_sites:
            if len(bins) > 0 and (n_sites - i) < bin_size:
                bins[-1] = np.concatenate([bins[-1], np.arange(i, n_sites)])
            else:
                bins.append(np.arange(i, n_sites))
            break
        else:
            bins.append(np.arange(i, j))
        i = j
    return bins

# 8) Function to compute p-value using percentile test (Non-parametric, single-tail, default: right)
def percentile_pvals(vals, side="right"):
    """
    Non-parametric bin-wise p-values based on percentiles.
    One-tailed: right or left.
    """
    n = len(vals)
    if n == 0:
        return np.ones(0)

    ranks = stats.rankdata(vals, method='average')  # rank within bin
    F = ranks / (n + 1e-12)

    if side == "right":
        # p = P(X >= x) = 1 - F
        pvals = 1 - F

    elif side == "left":
        # p = P(X <= x) = F
        pvals = F

    else:
        raise ValueError("side must be 'right' or 'left'")

    return pvals


def paired_one_tailed_pvalue(protein_vals, glycosite_vals, method="wilcoxon", side="right"):
    """Compute a one-tailed paired p-value from normalized protein/glycosite Rm pairs."""
    protein_vals = np.asarray(protein_vals, dtype=float)
    glycosite_vals = np.asarray(glycosite_vals, dtype=float)

    valid_mask = (~np.isnan(protein_vals)) & (~np.isnan(glycosite_vals))
    protein_valid = protein_vals[valid_mask]
    glycosite_valid = glycosite_vals[valid_mask]
    n_pairs = len(protein_valid)

    if n_pairs < 2:
        return 1.0, n_pairs

    diffs = glycosite_valid - protein_valid
    if np.allclose(diffs, 0):
        return 1.0, n_pairs

    if method == "paired_t":
        t_stat, p_two = stats.ttest_rel(glycosite_valid, protein_valid, nan_policy="omit")
        if not np.isfinite(t_stat) or not np.isfinite(p_two):
            return 1.0, n_pairs

        if side == "right":
            p_one = (p_two / 2.0) if (t_stat >= 0) else (1.0 - p_two / 2.0)
        elif side == "left":
            p_one = (p_two / 2.0) if (t_stat <= 0) else (1.0 - p_two / 2.0)
        else:
            raise ValueError("side must be 'right' or 'left'")

        return float(np.clip(p_one, 0.0, 1.0)), n_pairs

    if method == "wilcoxon":
        alt = "greater" if side == "right" else "less"
        try:
            _, p_one = stats.wilcoxon(
                glycosite_valid,
                protein_valid,
                alternative=alt,
                zero_method="wilcox",
                method="auto"
            )
        except TypeError:
            _, p_one = stats.wilcoxon(
                glycosite_valid,
                protein_valid,
                alternative=alt,
                zero_method="wilcox",
                mode="auto"
            )
        except Exception:
            p_one = 1.0

        return float(np.clip(p_one, 0.0, 1.0)), n_pairs

    raise ValueError("method must be 'wilcoxon' or 'paired_t'")


def assess_distribution_normality(values, alpha=0.05):
    """Assess normality with a sample-size-adaptive strategy."""
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    n = int(vals.size)

    if n < 8:
        return False, f"n={n} (<8): insufficient data, treated as NOT normal", "insufficient_n", n

    if n <= 50:
        stat, pval = stats.shapiro(vals)
        is_normal = bool(pval >= alpha)
        summary = (
            f"Shapiro-Wilk: W={stat:.4f}, p={pval:.4g}"
            + (" -> approximately normal" if is_normal else " -> NOT normal")
        )
        return is_normal, summary, "shapiro", n

    if n <= 300:
        mean_val = np.mean(vals)
        std_val = np.std(vals, ddof=1)
        if (not np.isfinite(std_val)) or (std_val <= 0):
            return False, f"Kolmogorov-Smirnov skipped (std={std_val:.4g}), treated as NOT normal", "ks_invalid_std", n
        stat, pval = stats.kstest(vals, "norm", args=(mean_val, std_val))
        is_normal = bool(pval >= alpha)
        summary = (
            f"Kolmogorov-Smirnov: D={stat:.4f}, p={pval:.4g}"
            + (" -> approximately normal" if is_normal else " -> NOT normal")
        )
        return is_normal, summary, "ks", n

    skew_val = stats.skew(vals)
    kurt_val = stats.kurtosis(vals, fisher=True)
    is_normal = bool(abs(skew_val) < 3 and abs(kurt_val) < 10)
    summary = (
        f"Large sample ({n}): skew={skew_val:.2f}, kurtosis={kurt_val:.2f}, "
        + ("approximately normal" if is_normal else "NOT normal")
    )
    return is_normal, summary, "skew_kurtosis_rule", n


def consensus_significance_from_replicates(
    delta_matrix,
    significance_matrix,
    delta_threshold=0.1,
    require_same_positive_direction=True,
):
    """
    Consensus significance for the normal / CV-binned branch.

    IMPORTANT revision:
    - effect-size gating uses MEAN DeltaRm across replicates;
    - it does NOT require every replicate DeltaRm to exceed the threshold.

    Statistical consistency rule is retained from the original implementation:
    all replicate significance values < 0.10 and at least one < 0.05. The
    significance values are either BH-FDR or raw p-values depending on the
    user-selected multiple-testing option.

    `require_same_positive_direction=True` retains the original normal-branch
    direction-consistency idea, but requires only positive direction, not
    DeltaRm > threshold in every replicate.
    """
    delta_matrix = np.asarray(delta_matrix, dtype=float)
    significance_matrix = np.asarray(significance_matrix, dtype=float)

    finite_delta = np.isfinite(delta_matrix).all(axis=1)
    finite_stat = np.isfinite(significance_matrix).all(axis=1)
    mean_delta = np.nanmean(delta_matrix, axis=1)
    mean_delta_pass = mean_delta > delta_threshold
    same_positive_direction = (delta_matrix > 0).all(axis=1)
    direction_pass = same_positive_direction if require_same_positive_direction else np.ones(len(mean_delta), dtype=bool)
    stat_all_10 = (significance_matrix < 0.10).all(axis=1)
    stat_any_05 = (significance_matrix < 0.05).any(axis=1)

    sig_mask = finite_delta & finite_stat & direction_pass & mean_delta_pass & stat_all_10 & stat_any_05
    return (
        sig_mask,
        mean_delta,
        same_positive_direction,
        mean_delta_pass,
        stat_all_10,
        stat_any_05,
    )

# 9) Function to convert p-values to significance stars
def pval_to_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# ------------------------------
# 1. TMT Protein Quantification QC
# ------------------------------
stage("Stage 1/10 - Protein quantification QC")
file = input("Please provide path of TMT quantification protein CSV file: ").strip()
file = normalize_path(file)
if not os.path.exists(file):
    raise FileNotFoundError(file)

work_dir = os.path.dirname(file)
os.chdir(work_dir)
output_root, tables_dir, figures_dir, metadata_dir, intermediate_dir = build_output_dirs(work_dir)

df = pd.read_csv(file)
protein_col_candidates = ["Protein Group","ProteinGroup","Protein_Group","Accession"]
protein_col = next((c for c in protein_col_candidates if c in df.columns), df.columns[0])

sample_order, sample_cols_map = detect_sample_blocks(df)
if not sample_order:
    tmt_cols = [c for c in df.columns if re.search(r"(126|127|128|129|130|131)", c)]
    num_samples = int(input("Enter number of biological samples: ").strip())
    channels_per_sample = len(tmt_cols)//num_samples
    sample_order = list(range(1,num_samples+1))
    sample_cols_map = {s: tmt_cols[i*channels_per_sample:(i+1)*channels_per_sample] for i,s in enumerate(sample_order)}

sample_names = [input(f"Please enter name for Sample {s} (use 'control' for control group): ").strip() for s in sample_order]

# Rename columns
col_map = {}
sample_block_cols = {}
for i, s in enumerate(sample_order):
    sample_name = sample_names[i]
    cols = sample_cols_map[s]
    sample_block_cols[sample_name] = []
    for old in cols:
        repfrag = extract_reporter_fragment(old)
        newname = f"{sample_name} {repfrag}"
        col_map[old] = newname
        sample_block_cols[sample_name].append(newname)
df.rename(columns=col_map, inplace=True)
all_reporter_cols = [c for cols in sample_block_cols.values() for c in cols]

# Quantified overview
quantified_any = (df[all_reporter_cols] != 0).any(axis=1)
print(f"Total proteins quantified (any non-zero channel): {quantified_any.sum()}")

per_sample_sets = {}
for sample in sample_names:
    cols = sample_block_cols[sample]
    per_sample_sets[sample] = set(df.index[(df[cols] != 0).any(axis=1)])

# Unique peptide filter
unique_col_candidates = ["#Unique","#Unique Peptides","Unique"]
unique_col = next((c for c in unique_col_candidates if c in df.columns), None)
if unique_col:
    choice = input("Do you want to filter proteins by unique peptides ≥2? (YES/NO, default=NO): ").strip().upper()
    if choice=="YES":
        unique_ge2 = (df[unique_col]>=2)
    else:
        unique_ge2 = (df[unique_col]>=1)
else:
    unique_ge2 = pd.Series(True, index=df.index)
    print("No unique peptide column detected, skipping unique filter.")

for sample in sample_names:
    cols = sample_block_cols[sample]
    per_sample_sets[sample] = set(
        df.index[(df[cols] != 0).any(axis=1) & unique_ge2]
    )

# Venn / UpSet plot
num_samples = len(sample_names)
if num_samples>1:
    if num_samples<=3:
        plt.figure(figsize=(6,6))
        if num_samples==2:
            venn2([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
        else:
            venn3([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
        plt.title("Quantified proteins (non-zero & unique filter)")
        plt.savefig(os.path.join(figures_dir, "01_protein_quantification_overlap.tiff"), dpi=300)
        plt.close()
    else:
        upset_data = from_contents(per_sample_sets)
        plt.figure(figsize=(8,6))
        UpSet(upset_data, show_counts=True).plot()
        plt.savefig(os.path.join(figures_dir, "01_protein_quantification_overlap_upset.tiff"), dpi=300)
        plt.close()

# Strict filtering
strict_mask = unique_ge2.copy()
for sample in sample_names:
    cols = sample_block_cols[sample]
    strict_mask &= (df[cols]!=0).all(axis=1)
df_strict = df.loc[strict_mask].copy()

for sample in sample_names:
    cols = sample_block_cols[sample]
    n_complete = (df[cols]!=0).all(axis=1).sum()
    print(f"Sample {sample}: {n_complete} proteins fully quantified in all channels")
n_all_complete = df_strict.shape[0]
print(f"Proteins fully quantified across all samples: {n_all_complete}")

control_sample = next((s for s in sample_names if 'control' in s.lower()), sample_names[0])
control_cols = sample_block_cols[control_sample]
df_control_complete = df.loc[(df[control_cols]!=0).all(axis=1) & unique_ge2]
df_control_complete.to_csv(os.path.join(tables_dir, "01_control_proteins_complete_quantification.csv"), index=False)
df_strict.to_csv(os.path.join(tables_dir, "02_all_samples_proteins_complete_quantification.csv"), index=False)

print("TMT protein QC completed.")

# ------------------------------
# 2. Rm calculation & normalization of control group
# ------------------------------
stage("Stage 2/10 - Control-proteome Rm calculation and replicate-median normalization")
ctrl_cols = sample_block_cols[control_sample]
ch126, ch127, ch128, ch129, ch130, ch131 = ctrl_cols
df_control_complete = df_control_complete.copy()


df_control_complete["Rm1"] = df_control_complete[ch129] / df_control_complete[ch127]
df_control_complete["Rm2"] = df_control_complete[ch130] / df_control_complete[ch128]
df_control_complete["Rm3"] = df_control_complete[ch131] / df_control_complete[ch126]

# Replicate-median multiplicative normalization.
# Each biological replicate has its own control-proteome Rm median. All three
# medians are aligned to a common reference center with a multiplicative factor.
medians, rm_reference_median, CF = compute_replicate_median_normalization(
    df_control_complete, ["Rm1", "Rm2", "Rm3"]
)
for col in ["Rm1", "Rm2", "Rm3"]:
    df_control_complete[f"{col}_normalized"] = df_control_complete[col] * CF[col]

cf_df = pd.DataFrame({
    "Replicate": ["Rm1", "Rm2", "Rm3"],
    "Raw_Rm_Median": medians.values,
    "Reference_Rm_Median": [rm_reference_median] * 3,
    "Correction_Factor": CF.values,
})
cf_df.to_csv(os.path.join(tables_dir, "03_rm_normalization_factors.csv"), index=False)
export_csv(df_control_complete, os.path.join(tables_dir, "04_control_proteins_with_rm.csv"), PROTEIN_EXPORT_RENAME)

print("Rm normalization completed for control group.")

# ------------------------------
# 3. N-glycopeptide analysis
# ------------------------------
stage("Stage 3/10 - N-glycopeptide detection and protein-sequence mapping")
peptide_csv_file = input("Enter path of TMT quantification peptide CSV file: ").strip()
peptide_csv_file = normalize_path(peptide_csv_file)

fasta_file = input("Enter path to UniProt fasta file: ").strip()
fasta_file = normalize_path(fasta_file)

work_dir = os.path.dirname(peptide_csv_file)
os.chdir(work_dir)

fasta_csv_file = os.path.join(intermediate_dir, "uniprot_fasta_parsed.csv")
with open(fasta_file, 'r', encoding='utf-8') as f, open(fasta_csv_file, 'w', newline='', encoding='utf-8') as csvfile:
    import csv
    writer = csv.writer(csvfile)
    writer.writerow(['Accession', 'Protein Name', 'Sequence'])
    accession = ''
    protein_name = ''
    seq_lines = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if accession:
                writer.writerow([accession, protein_name, ''.join(seq_lines)])
            header = line[1:].strip()
            token = header.split()[0] if header else ''
            accession = canonical_protein_id(token) if token else ''

            parts = [p.strip() for p in token.split('|')] if '|' in token else []
            if len(parts) >= 3:
                third_tokens = parts[2].split(' ')
                protein_name = ' '.join(third_tokens[1:]).strip() if len(third_tokens) > 1 else parts[2]
            else:
                tokens = header.split()
                protein_name = ' '.join(tokens[1:]).strip() if len(tokens) > 1 else ''
            seq_lines = []
        else:
            seq_lines.append(line)
    if accession:
        writer.writerow([accession, protein_name, ''.join(seq_lines)])
print(f"Fasta converted to CSV: {fasta_csv_file}")

pep_df = pd.read_csv(peptide_csv_file)
prot_df = pd.read_csv(fasta_csv_file)
pep_df['UniProtID'] = pep_df['Accession'].apply(canonical_protein_id)
prot_dict = dict(zip(prot_df['Accession'], prot_df['Sequence']))

intensity_cols = [c for c in pep_df.columns if re.search(r"Intensity Sample \d+ TMT6-\d+", c)]
if not intensity_cols:
    raise ValueError("No intensity columns detected!")
sample_numbers = sorted(list(set(int(re.search(r"Sample (\d+)", c).group(1)) for c in intensity_cols)))
sample_map = {}
for sn in sample_numbers:
    old_prefix = f"Sample {sn}"
    new_name = input(f"Enter new name for {old_prefix}: ").strip()
    sample_map[old_prefix] = new_name
col_map = {}
for c in intensity_cols:
    sn_match = re.search(r"Sample (\d+)", c)
    if sn_match:
        sn = sn_match.group(1)
        new_prefix = sample_map[f"Sample {sn}"]
        col_map[c] = c.replace(f"Sample {sn}", new_prefix)
pep_df.rename(columns=col_map, inplace=True)
renamed_intensity_cols = list(col_map.values())

# Filter N-glycopeptides
results = []
for idx, row in pep_df.iterrows():
    pep = row['Peptide']
    uni_id = row['UniProtID']
    prot_seq = prot_dict.get(uni_id)
    if not prot_seq:
        continue
    stripped_seq, n_pos_list = parse_peptide(pep)
    if not n_pos_list:
        continue
    start_idx = prot_seq.find(stripped_seq)
    if start_idx == -1:
        continue
    kept_pep_positions = []
    kept_prot_positions = []
    for n_pos in n_pos_list:
        prot_n_pos = start_idx + n_pos
        if prot_n_pos + 2 < len(prot_seq):
            X = prot_seq[prot_n_pos + 1]
            ST = prot_seq[prot_n_pos + 2]
            if X != 'P' and ST in ['S','T']:
                kept_pep_positions.append(n_pos + 1)
                kept_prot_positions.append(prot_n_pos + 1)
    if not kept_prot_positions:
        continue
    out_row = {
        'Protein Accession': uni_id,
        'Peptide': pep,
        'Stripped Sequence': stripped_seq,
        'N-Glycosite in Peptide': ';'.join(map(str, kept_pep_positions)),
        'N-Glycosite in Protein': ';'.join(map(str, kept_prot_positions))
    }
    for c in renamed_intensity_cols:
        out_row[c] = row[c]
    if 'AScore' in pep_df.columns:
        out_row['AScore'] = row['AScore']
    results.append(out_row)

output_csv_file = os.path.join(tables_dir, "05_n_glycopeptides_with_intensity.csv")
result_columns = [
    'Protein Accession',
    'Peptide',
    'Stripped Sequence',
    'N-Glycosite in Peptide',
    'N-Glycosite in Protein',
] + renamed_intensity_cols
if 'AScore' in pep_df.columns:
    result_columns.append('AScore')

results_df = pd.DataFrame(results, columns=result_columns)
results_df.to_csv(output_csv_file, index=False)
print(f"Filtered N-glycopeptides saved: {output_csv_file}")

if results_df.empty:
    raise ValueError(
        "No N-glycopeptides were retrieved after peptide-FASTA mapping. "
        "Please check Accession format consistency between peptide CSV and FASTA, "
        "and verify N(+0.98)-based motif parsing."
    )

# Intensity ratios plot
glyco_intensity = results_df[renamed_intensity_cols].sum().values
total_intensity = pep_df[renamed_intensity_cols].sum().values
non_glyco_intensity = total_intensity - glyco_intensity
ratios = glyco_intensity / total_intensity

fig, ax = plt.subplots(figsize=(10,6))
bar_width = 0.6
x = np.arange(len(renamed_intensity_cols))
ax.bar(x, glyco_intensity, width=bar_width, color='#E63946', label='N-glycopeptides')
ax.bar(x, non_glyco_intensity, bottom=glyco_intensity, width=bar_width, color='#457B9D', label='Other peptides')
for i, ratio in enumerate(ratios):
    ax.text(x[i], glyco_intensity[i]/2, f"{ratio:.2f}", ha='center', va='center', color='white', fontsize=10, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(renamed_intensity_cols, rotation=45, ha='right')
ax.set_ylabel('Total Intensity')
ax.set_title('Intensity Composition of N-glycopeptides per TMT Channel')
ax.legend()
plt.tight_layout()
tiff_file = os.path.join(figures_dir, "02_n_glycopeptide_intensity_composition.tiff")
fig.savefig(tiff_file, dpi=300, format='tiff')
plt.close(fig)
print(f"Stacked bar chart saved: {tiff_file}")

# ----------------------------
# Protein N-glycosite distribution (all filtered N-glycopeptides)
# ----------------------------
stage("Stage 4/10 - Site-centric aggregation and quantification completeness QC")
prot_sites = results_df.groupby('Protein Accession')['N-Glycosite in Protein'].apply(
    lambda x: set(';'.join(x).split(';'))
)
prot_site_counts = prot_sites.apply(len)
total_sites = prot_site_counts.sum()

print(f"Total N-glycopeptides: {len(results_df)}")
print(f"Proteins with N-glycosites: {prot_site_counts.shape[0]}")
print(f"Total N-glycosites: {total_sites}")

# Cap counts >5
prot_site_counts_capped = prot_site_counts.apply(lambda x: x if x <=5 else '>5')
dist = prot_site_counts_capped.value_counts().sort_index(
    key=lambda x: [int(i) if i != '>5' else 6 for i in x]
)

plt.figure(figsize=(8,6))
bars = plt.barh(dist.index.astype(str), dist.values, color='skyblue')

for bar in bars:
    width = bar.get_width()
    plt.text(width + max(dist.values)*0.01, bar.get_y() + bar.get_height()/2,
             str(width), va='center')

plt.xlabel("Protein count")
plt.ylabel("Number of N-glycosites per protein")
plt.title("Distribution of N-glyco site counts per protein (capped at >5)")
plt.tight_layout()

output_file = os.path.join(figures_dir, "03_protein_nglycosite_distribution_all.tiff")
plt.savefig(output_file, dpi=300, format='tiff')
plt.close()
print(f"Protein N-glycosite distribution saved: {output_file}")


# Motif and site-centric intensity
motif_len = int(input("Enter total motif window length (odd number, e.g., 7, 9, 11): ").strip())
if motif_len % 2 == 0:
    raise ValueError("Motif window length must be an odd number.")
half_len = motif_len // 2

site_results = []
for prot, group in results_df.groupby('Protein Accession'):
    prot_seq = prot_dict.get(prot, "")
    prot_len = len(prot_seq)
    all_sites = set(';'.join(group['N-Glycosite in Protein']).split(';'))
    for n_site in all_sites:
        if not n_site:
            continue
        n_site = int(n_site)
        center = n_site - 1
        start = center - half_len
        end = center + half_len
        motif_chars = []
        for i in range(start, end + 1):
            if 0 <= i < prot_len:
                motif_chars.append(prot_seq[i])
            else:
                motif_chars.append('-')
        motif = ''.join(motif_chars)
        pep_mask = group['N-Glycosite in Protein'].apply(lambda x: str(n_site) in x.split(';'))
        pep_subset = group[pep_mask]
        site_row = {'Protein Accession': prot, 'N-Glycosite': n_site, 'Motif': motif}
        for c in renamed_intensity_cols:
            site_row[c] = pep_subset[c].sum()
        site_results.append(site_row)

site_df = pd.DataFrame(site_results)
site_csv_file = os.path.join(tables_dir, "06_nglycosite_site_intensity.csv")
site_df.to_csv(site_csv_file, index=False)
print(f"Site-centric intensity table saved: {site_csv_file}")

# User input for the sample group
chosen_group = input("Enter the sample group to calculate ΔRm (exact name, e.g., 'N-Glyco'): ").strip()

# Map TMT 6-126~131 columns for the selected group
ch126_site = f"Intensity {chosen_group} TMT6-126"
ch127_site = f"Intensity {chosen_group} TMT6-127"
ch128_site = f"Intensity {chosen_group} TMT6-128"
ch129_site = f"Intensity {chosen_group} TMT6-129"
ch130_site = f"Intensity {chosen_group} TMT6-130"
ch131_site = f"Intensity {chosen_group} TMT6-131"

site_cols_for_calc = [ch126_site, ch127_site, ch128_site, ch129_site, ch130_site, ch131_site]

# Check if all columns exist
missing_cols = [c for c in site_cols_for_calc if c not in site_df.columns]
if missing_cols:
    raise ValueError(f"Columns not found for selected group '{chosen_group}': {missing_cols}")

print(f"Selected sample group for ΔRm calculation: '{chosen_group}'")
print("Columns mapped for calculation:", site_cols_for_calc)

print("\n=== Filtering sites with complete quantification in selected-group 6 channels ===")
site_df_clean = site_df.replace(0, np.nan)
complete_mask = site_df_clean[site_cols_for_calc].notna().all(axis=1)
site_df_complete = site_df_clean[complete_mask].copy()
num_sites = site_df_complete.shape[0]
num_proteins = site_df_complete['Protein Accession'].nunique()
print(f"Total quantified sites (complete in selected-group 6 channels): {num_sites}")
print(f"Proteins containing these fully quantified sites: {num_proteins}")
site_complete_csv = os.path.join(tables_dir, "07_nglycosites_complete_quantification.csv")
site_df_complete.to_csv(site_complete_csv, index=False)
print(f"Filtered complete-quantification site table saved: {site_complete_csv}")

prot_site_count = site_df_complete.groupby('Protein Accession')['N-Glycosite'].nunique()
prot_site_counts_capped = prot_site_count.apply(lambda x: x if x <= 5 else '>5')
dist = prot_site_counts_capped.value_counts().sort_index(key=lambda x: [int(i) if i != '>5' else 6 for i in x])

plt.figure(figsize=(8, 6))
bars = plt.barh(dist.index.astype(str), dist.values, color="#A8C9E0")
for bar in bars:
    width = bar.get_width()
    plt.text(width + max(dist.values) * 0.01,
             bar.get_y() + bar.get_height() / 2,
             str(width),
             va='center')
plt.xlabel("Protein count")
plt.ylabel("Number of N-glycosites per protein")
plt.title("Distribution of fully quantified N-glycosites per protein (capped at >5)")
plt.tight_layout()
out_plot = os.path.join(figures_dir, "04_protein_nglycosite_distribution_complete_quant.tiff")
plt.savefig(out_plot, dpi=300, format='tiff')
plt.close()
print(f"Protein N-glycosite distribution saved: {out_plot}")

print("Site-level complete-quantification preprocessing completed.")

# ----------------------------
# Fully quantified N-glycosites with fully quantified proteins (in control group) filtering
# ==============================
stage("Stage 5/10 - Match fully quantified glycosites to fully quantified control proteins")
df_control_complete['UniProtID'] = df_control_complete['Accession'].apply(canonical_protein_id)
proteins_complete = set(df_control_complete['UniProtID'])
site_df_filtered = site_df_complete[site_df_complete['Protein Accession'].isin(proteins_complete)].copy()
df_control_filtered = df_control_complete[df_control_complete['UniProtID'].isin(site_df_filtered['Protein Accession'])].copy()

# Venn Plot (Protein-level)
plt.figure(figsize=(5,5))
v = venn2([proteins_complete, set(site_df_complete['Protein Accession'])],
          set_labels=['Fully quantified proteins in control group', 
                      'Proteins with fully quantified N-Glycosites'],
          set_colors=('#8AB7C2', '#9EAAD4'), alpha=0.5)

# Adjust label font sizes
for text in v.set_labels:
    if text:
        text.set_fontsize(6)
for text in v.subset_labels:
    if text:
        text.set_fontsize(6)
        text.set_fontweight('bold')

plt.title('Protein-level quantification completeness overlap', fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "05_protein_site_quantification_overlap_venn.tiff"), dpi=300, format='tiff')
plt.close()

# Output intersection stats
num_proteins = df_control_filtered.shape[0]
num_sites = site_df_filtered.shape[0]
print(f"Number of proteins with complete control and N-glycosite quantification: {num_proteins}")
print(f"Number of N-glycosites with complete control and N-glycosite quantification: {num_sites}")

# ==============================
# Step 4: Site distribution per protein
# ==============================
prot_site_count = site_df_filtered.groupby('Protein Accession')['N-Glycosite'].nunique()
# Cap counts >5
prot_site_counts_capped = prot_site_count.apply(lambda x: x if x <= 5 else '>5')
dist = prot_site_counts_capped.value_counts().sort_index(key=lambda x: [int(i) if i != '>5' else 6 for i in x])

plt.figure(figsize=(8, 6))
bars = plt.barh(dist.index.astype(str), dist.values, color="#14C0CC")
for bar in bars:
    width = bar.get_width()
    plt.text(width + max(dist.values) * 0.01,
             bar.get_y() + bar.get_height() / 2,
             str(width),
             va='center')
plt.xlabel("Protein count")
plt.ylabel("Number of N-glycosites per protein")
plt.title("Distribution of fully quantified N-glycosites per fully quantified protein (capped at >5)")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "06_protein_nglycosite_distribution_intersection.tiff"), dpi=300, format='tiff')
plt.close()

# ==============================
# Step 5: Calculate Rm, CV, and ΔRm for intersected proteins/sites
# ==============================
stage("Stage 6/10 - Raw Rm CV QC, normalized Rm, and mean DeltaRm calculation")
print(site_df_filtered.columns)
print(df_control_filtered.columns)

# ---------------- Complete sites ----------------
site_df_complete = site_df_complete.copy()
# Calculate raw Rm for each site
site_df_complete['Rm_glyco1'] = site_df_complete[ch129_site] / site_df_complete[ch127_site]
site_df_complete['Rm_glyco2'] = site_df_complete[ch130_site] / site_df_complete[ch128_site]
site_df_complete['Rm_glyco3'] = site_df_complete[ch131_site] / site_df_complete[ch126_site]
# Calculate CV for Rm
site_df_complete['CV_Rm_glyco_complete'] = coefficient_of_variation_from_raw_rm(
    site_df_complete, ['Rm_glyco1', 'Rm_glyco2', 'Rm_glyco3']
)

# Save complete site-level table with Rm and CV
complete_outfile = os.path.join(tables_dir, f"08_nglycosites_rm_raw_cv_{chosen_group}_complete.csv")
export_csv(site_df_complete, complete_outfile, SITE_EXPORT_RENAME)
print(f"Complete site table with Rm and CV saved: {complete_outfile}")

# ---------------- Filtered sites ----------------
site_df_filtered = site_df_filtered.copy()
# Calculate raw Rm for each site
site_df_filtered['Rm_glyco1'] = site_df_filtered[ch129_site] / site_df_filtered[ch127_site]
site_df_filtered['Rm_glyco2'] = site_df_filtered[ch130_site] / site_df_filtered[ch128_site]
site_df_filtered['Rm_glyco3'] = site_df_filtered[ch131_site] / site_df_filtered[ch126_site]
# Calculate CV for Rm
site_df_filtered['CV_Rm_glyco_filtered'] = coefficient_of_variation_from_raw_rm(
    site_df_filtered, ['Rm_glyco1', 'Rm_glyco2', 'Rm_glyco3']
)

# Apply control CF correction
site_df_filtered['Rm_glyco1_norm'] = site_df_filtered['Rm_glyco1'] * CF['Rm1']
site_df_filtered['Rm_glyco2_norm'] = site_df_filtered['Rm_glyco2'] * CF['Rm2']
site_df_filtered['Rm_glyco3_norm'] = site_df_filtered['Rm_glyco3'] * CF['Rm3']

# Protein-level adjusted_Rm_unmod (arithmetic mean of normalized protein Rm) and CV of raw Rm
# Filtered proteins
df_control_filtered = df_control_filtered.copy()
# Arithmetic mean of normalized Rm
df_control_filtered['adjusted_Rm_unmod'] = df_control_filtered[['Rm1_normalized','Rm2_normalized','Rm3_normalized']].mean(axis=1)
# CV of raw Rm (Rm1~3)
df_control_filtered['CV_Rm_filtered'] = coefficient_of_variation_from_raw_rm(
    df_control_filtered, ['Rm1', 'Rm2', 'Rm3']
)

# Complete proteins
df_control_complete = df_control_complete.copy()
df_control_complete['CV_Rm'] = coefficient_of_variation_from_raw_rm(
    df_control_complete, ['Rm1', 'Rm2', 'Rm3']
)

# Save filtered protein-level CSV
new_protein_cols = ['Rm1','Rm2','Rm3','Rm1_normalized','Rm2_normalized','Rm3_normalized','adjusted_Rm_unmod','CV_Rm_filtered']
original_protein_cols = [c for c in df_control_filtered.columns if c not in new_protein_cols]
protein_cols_to_save_filtered = original_protein_cols + new_protein_cols

protein_outfile_filtered = os.path.join(tables_dir, "09_control_proteins_rm_summary_site_matched.csv")
export_csv(df_control_filtered, protein_outfile_filtered, PROTEIN_EXPORT_RENAME, columns=protein_cols_to_save_filtered)
print(f"Filtered protein-level summary saved: {protein_outfile_filtered}")

# Save complete protein-level CSV
new_protein_cols_complete = ['Rm1','Rm2','Rm3','CV_Rm']
original_protein_cols_complete = [c for c in df_control_complete.columns if c not in new_protein_cols_complete]
protein_cols_to_save_complete = original_protein_cols_complete + new_protein_cols_complete

protein_outfile_complete = os.path.join(tables_dir, "10_control_proteins_rm_summary_complete.csv")
export_csv(df_control_complete, protein_outfile_complete, PROTEIN_EXPORT_RENAME, columns=protein_cols_to_save_complete)
print(f"Complete protein-level summary saved: {protein_outfile_complete}")


# ---------------- Site-level ΔRm ----------------
# Arithmetic mean of normalized glycosite Rm
site_df_filtered['adjusted_Rm_glyco'] = site_df_filtered[['Rm_glyco1_norm','Rm_glyco2_norm','Rm_glyco3_norm']].mean(axis=1)
# Map protein-level adjusted Rm_unmod
prot_rm_map = df_control_filtered.set_index('UniProtID')['adjusted_Rm_unmod'].to_dict()
site_df_filtered['adjusted_Rm_unmod'] = site_df_filtered['Protein Accession'].map(prot_rm_map)
# Map protein-level normalized replicate Rm for paired non-parametric testing
prot_rm1_map = df_control_filtered.set_index('UniProtID')['Rm1_normalized'].to_dict()
prot_rm2_map = df_control_filtered.set_index('UniProtID')['Rm2_normalized'].to_dict()
prot_rm3_map = df_control_filtered.set_index('UniProtID')['Rm3_normalized'].to_dict()
site_df_filtered['Rm_unmod1_norm'] = site_df_filtered['Protein Accession'].map(prot_rm1_map)
site_df_filtered['Rm_unmod2_norm'] = site_df_filtered['Protein Accession'].map(prot_rm2_map)
site_df_filtered['Rm_unmod3_norm'] = site_df_filtered['Protein Accession'].map(prot_rm3_map)
# Calculate ΔRm for each biological replicate and their mean
site_df_filtered['ΔRm_rep1'] = site_df_filtered['Rm_glyco1_norm'] - site_df_filtered['Rm_unmod1_norm']
site_df_filtered['ΔRm_rep2'] = site_df_filtered['Rm_glyco2_norm'] - site_df_filtered['Rm_unmod2_norm']
site_df_filtered['ΔRm_rep3'] = site_df_filtered['Rm_glyco3_norm'] - site_df_filtered['Rm_unmod3_norm']
site_df_filtered['ΔRm'] = site_df_filtered[['ΔRm_rep1', 'ΔRm_rep2', 'ΔRm_rep3']].mean(axis=1)
# IMPORTANT: the effect-size threshold is applied to this mean ΔRm later.
# Individual replicate ΔRm values are NOT required to exceed 0.1.

# Save site-level CSV
new_site_cols = ['Rm_glyco1','Rm_glyco2','Rm_glyco3','CV_Rm_glyco_filtered',
                 'Rm_glyco1_norm','Rm_glyco2_norm','Rm_glyco3_norm',
                 'adjusted_Rm_glyco','adjusted_Rm_unmod',
                 'Rm_unmod1_norm','Rm_unmod2_norm','Rm_unmod3_norm',
                 'ΔRm_rep1','ΔRm_rep2','ΔRm_rep3','ΔRm']
original_site_cols = [c for c in site_df_filtered.columns if c not in new_site_cols]
site_cols_to_save = original_site_cols + new_site_cols

site_outfile = os.path.join(tables_dir, "11_nglycosites_delta_rm.csv")
export_csv(site_df_filtered, site_outfile, SITE_EXPORT_RENAME, columns=site_cols_to_save)
print(f"Site-level ΔRm table saved: {site_outfile}")

# ------------------------------
# Step 6: CV QC plots for complete and filtered N-glycosites (TIFF)
# ------------------------------
cv_plot_params = [
    ("All fully quantified sites", site_df_complete['CV_Rm_glyco_complete'], f"CV_distribution_complete_sites_{chosen_group}.tiff"),
    ("Fully quantified sites of fully quantified proteins", site_df_filtered['CV_Rm_glyco_filtered'], f"CV_distribution_filtered_sites_{chosen_group}.tiff")
]

for label, cv_series_raw, fname in cv_plot_params:
    # Remove NaN and sort descending
    cv_series = cv_series_raw.dropna().sort_values(ascending=False).reset_index(drop=True)
    rank = cv_series.index + 1

    plt.figure(figsize=(6,4))
    plt.scatter(rank, cv_series.values, color="#5E556A", alpha=0.7, s=20)

    # Proportion of sites with CV <= 0.3 and <= 0.4
    prop_30 = (cv_series <= 0.3).sum() / len(cv_series) * 100
    prop_40 = (cv_series <= 0.4).sum() / len(cv_series) * 100

    plt.text(0.7*len(cv_series), 0.305, f"≤30%: {prop_30:.1f}%", fontsize=10)
    plt.text(0.7*len(cv_series), 0.405, f"≤40%: {prop_40:.1f}%", fontsize=10)

    plt.xlabel("Sites sorted by CV (high→low)")
    plt.ylabel("CV")
    plt.title(f"CV distribution - {chosen_group} ({label})", fontsize=8)
    plt.tight_layout()

    cv_fig_file = os.path.join(figures_dir, fname)
    plt.savefig(cv_fig_file, dpi=300)
    plt.close()
    print(f"{label} CV distribution TIFF plot saved: {cv_fig_file}")
    
# ------------------------------
# Step 7: CV QC plots for proteins (complete and filtered) (TIFF)
# ------------------------------
protein_cv_plot_params = [
    ("All fully quantified proteins", df_control_complete['CV_Rm'], "CV_distribution_complete_proteins.tiff"),
    ("Fully quantified proteins with fully quantified N-glycosites", df_control_filtered['CV_Rm_filtered'], "CV_distribution_filtered_proteins.tiff")
]

for label, cv_series_raw, fname in protein_cv_plot_params:
    # Remove NaN and sort descending
    cv_series = cv_series_raw.dropna().sort_values(ascending=False).reset_index(drop=True)
    rank = cv_series.index + 1

    plt.figure(figsize=(6,4))
    plt.scatter(rank, cv_series.values, color="#2A9D8F", alpha=0.7, s=20)

    # Proportion of proteins with CV <= 0.3 and <= 0.4
    prop_30 = (cv_series <= 0.3).sum() / len(cv_series) * 100
    prop_40 = (cv_series <= 0.4).sum() / len(cv_series) * 100

    plt.text(0.7*len(cv_series), 0.305, f"≤30%: {prop_30:.1f}%", fontsize=10)
    plt.text(0.7*len(cv_series), 0.405, f"≤40%: {prop_40:.1f}%", fontsize=10)

    plt.xlabel("Proteins sorted by CV (high→low)")
    plt.ylabel("CV")
    plt.title(f"CV distribution - Control Proteins ({label})", fontsize=8)
    plt.tight_layout()

    cv_fig_file = os.path.join(figures_dir, fname)
    plt.savefig(cv_fig_file, dpi=300)
    plt.close()
    print(f"{label} CV distribution TIFF plot saved: {cv_fig_file}")

# ------------------------------
# Step 8: CV-based QC filtering for sites (user threshold)
# ------------------------------
stage("Stage 7/10 - Joint glycosite/protein raw-Rm CV filtering")
# This step should be placed after saving site_outfile, ensuring that
# site_df_filtered and df_control_filtered already contain ΔRm and protein CV columns.

# 1) Read and normalize the CV threshold input (allow "40" or "0.40")
raw_thresh = input("Enter CV threshold for filtering (e.g., '40' for 40% or '0.4' for 0.4): ").strip()
try:
    cv_threshold = float(raw_thresh)
except:
    raise ValueError("CV threshold must be a number like 40 or 0.4")

# Convert percent values (>1) to fractional form
if cv_threshold > 1:
    cv_threshold = cv_threshold / 100.0

print(f"Using CV threshold = {cv_threshold:.3f} (i.e. {cv_threshold*100:.1f}%)")

apply_bh = prompt_yes_no(
    "Apply Benjamini-Hochberg multiple-testing correction? (YES/NO, default=YES): ",
    default=True,
)
print("Multiple-testing mode:", "BH-FDR" if apply_bh else "raw p-values (no BH correction)")

delta_rm_threshold = DELTA_RM_THRESHOLD_DEFAULT
print(f"Mean DeltaRm effect-size threshold = {delta_rm_threshold:.3f}")

# 2) Check required columns (site CV column and protein CV column)
site_cv_col = "CV_Rm_glyco_filtered"   # site-level CV column
protein_cv_col = "CV_Rm_filtered"      # protein-level CV column

missing = []
if site_cv_col not in site_df_filtered.columns:
    missing.append(site_cv_col)
if protein_cv_col not in df_control_filtered.columns:
    missing.append(protein_cv_col)
if "Protein Accession" not in site_df_filtered.columns:
    missing.append("Protein Accession (in site_df_filtered)")
if "UniProtID" not in df_control_filtered.columns:
    missing.append("UniProtID (in df_control_filtered)")

if missing:
    raise KeyError(f"Required columns missing for CV filtering: {missing}")

# 3) Map protein-level CV (UniProtID → CV_Rm_filtered)
protein_cv_series = df_control_filtered.set_index("UniProtID")[protein_cv_col]

# Report match statistics
mapped_proteins = site_df_filtered["Protein Accession"].isin(protein_cv_series.index)
n_mapped = mapped_proteins.sum()
n_total_sites = site_df_filtered.shape[0]
print(f"Site table: {n_total_sites} sites total. {n_mapped} sites have matching protein CV in df_control_filtered.")

# 4) Map protein CV into the site table (keep all original columns including ΔRm)
site_df_filtered = site_df_filtered.copy()
site_df_filtered["Protein_CV"] = site_df_filtered["Protein Accession"].map(protein_cv_series)

# 5) CV filtering (NaN treated as failing)
before_count = site_df_filtered.shape[0]
pass_mask = (
    (site_df_filtered[site_cv_col].notna()) &
    (site_df_filtered["Protein_CV"].notna()) &
    (site_df_filtered[site_cv_col] <= cv_threshold) &
    (site_df_filtered["Protein_CV"] <= cv_threshold)
)
site_df_CV_pass = site_df_filtered.loc[pass_mask].copy()
after_count = site_df_CV_pass.shape[0]

print(f"CV filtering: retained {after_count} / {before_count} sites ({after_count/before_count*100 if before_count>0 else 0:.1f}%).")
print(f"Sites with missing protein CV (excluded): {(~mapped_proteins).sum()}")

# 6) Save CV-filtered site table
cv_filtered_output = os.path.join(tables_dir, f"12_nglycosites_cv_pass_{int(cv_threshold*100)}pct.csv")
export_csv(site_df_CV_pass, cv_filtered_output, SITE_EXPORT_RENAME)
print(f"CV-filtered site table saved: {cv_filtered_output}")

# 7) Optional QC plot: site CV vs protein CV
try:
    plt.figure(figsize=(6,6))
    plt.scatter(site_df_filtered[site_cv_col], site_df_filtered["Protein_CV"], s=10, alpha=0.6)
    plt.axvline(cv_threshold, linestyle='--', linewidth=1, label=f"site CV threshold ({cv_threshold:.2f})")
    plt.axhline(cv_threshold, linestyle='--', linewidth=1, label=f"protein CV threshold ({cv_threshold:.2f})")
    plt.xlabel("Site CV (CV_Rm_glyco_filtered)")
    plt.ylabel("Protein CV (CV_Rm_filtered)")
    plt.title(f"Site CV vs Protein CV (threshold {cv_threshold:.2f})")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    qc_plot_file = os.path.join(figures_dir, f"07_site_vs_protein_raw_rm_cv_{int(cv_threshold*100)}pct.tiff")
    plt.savefig(qc_plot_file, dpi=300, format='tiff')
    plt.close()
    print(f"CV scatter plot saved: {qc_plot_file}")
except Exception as e:
    print("Warning: failed to draw CV scatter plot:", e)

# 8) Final summary: number of proteins with at least one CV-passing site
kept_proteins = site_df_CV_pass["Protein Accession"].nunique()
print(f"Number of proteins with ≥1 site passing CV filter: {kept_proteins}")

# -----------------------------
# Step 9: Visualize Delta_Rm distribution and check normality (robust)
# -----------------------------
stage("Stage 8/10 - Replicate-level DeltaRm distribution and normality assessment")
delta_rep_cols = ["ΔRm_rep1", "ΔRm_rep2", "ΔRm_rep3"]
missing_delta_cols = [c for c in delta_rep_cols if c not in site_df_CV_pass.columns]
if missing_delta_cols:
    raise KeyError(f"Required replicate DeltaRm columns missing: {missing_delta_cols}")

normality_records = []
for rep_idx, rep_col in enumerate(delta_rep_cols, start=1):
    delta_rm_for_plot = site_df_CV_pass[rep_col].dropna().astype(float)

    # --- 1. Histogram + KDE ---
    plt.figure(figsize=(6,4))
    sns.histplot(delta_rm_for_plot, kde=True, color='#A7AED2', bins=30)
    plt.title(f'{rep_col} Distribution (N-glycosites passed CV QC)')
    plt.xlabel(rep_col)
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, f'08_delta_rm_rep{rep_idx}_histogram.tiff'), format='tiff', dpi=300)
    plt.close()

    # --- 2. Q-Q plot ---
    plt.figure(figsize=(6,4))
    sm.qqplot(delta_rm_for_plot, line='s')
    plt.title(f'{rep_col} Q-Q Plot')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, f'09_delta_rm_rep{rep_idx}_qqplot.tiff'), format='tiff', dpi=300)
    plt.close()

    # --- 3. P-P plot ---
    plt.figure(figsize=(6,4))
    sm.ProbPlot(delta_rm_for_plot).ppplot(line='45')
    plt.title(f'{rep_col} P-P Plot')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, f'10_delta_rm_rep{rep_idx}_ppplot.tiff'), format='tiff', dpi=300)
    plt.close()

    # --- 4. Normality test for each biological replicate ---
    is_normal, normality_result, method_used, n_used = assess_distribution_normality(
        delta_rm_for_plot.values,
        alpha=NORMALITY_ALPHA_DEFAULT
    )
    normality_records.append({
        "replicate": rep_col,
        "n": n_used,
        "method": method_used,
        "summary": normality_result,
        "is_normal": is_normal
    })
    print(f"{rep_col}: {normality_result}")

normality_df = pd.DataFrame(normality_records)
normality_outfile = os.path.join(tables_dir, "13_delta_rm_replicate_normality_summary.csv")
normality_df.to_csv(normality_outfile, index=False)
print(f"Replicate-level normality summary saved: {normality_outfile}")

delta_rm_not_normal = not bool(normality_df["is_normal"].all())

# --- 5. Suggest test based on replicate-level normality ---
if delta_rm_not_normal:
    bad_reps = normality_df.loc[~normality_df["is_normal"], "replicate"].tolist()
    print(f"Normality gate failed for replicates: {bad_reps}")
    print("At least one biological replicate is NOT normal; use site-level paired testing (default Wilcoxon, optional paired t-test override)")
else:
    print("All biological replicates are approximately normal, use CV-binned robust sigma testing per replicate")

# --- 6. Debug check: confirm original table unchanged ---
print("Step 9 check: site_df_CV_pass rows =", site_df_CV_pass.shape[0])

# ------------------------------
# Step 10-11: Significant site identification after CV QC (branch by normality)
# ------------------------------
stage("Stage 9/10 - Statistical significance analysis")

# 1) Compute representative CV for each site (max of site CV and protein CV)
site_df_CV_pass = site_df_CV_pass.copy()
site_df_CV_pass["Rep_CV"] = site_df_CV_pass[[site_cv_col, "Protein_CV"]].max(axis=1)

# 2) Sort sites by representative CV (ascending)
site_df_sorted = site_df_CV_pass.sort_values("Rep_CV").reset_index(drop=True)
delta_matrix_sorted = site_df_sorted[delta_rep_cols].to_numpy(dtype=float)

if not delta_rm_not_normal:
    print("Step10 strategy: all replicate ΔRm distributions are approximately normal -> use binning + robust-sigma test per replicate")
    
    # 3) Loop over default bin counts and compute p-values and FDR
    default_bin_counts = [1, 5, 10, 15, 20, 25, 30, 35, 40]
    print(f"Default bin counts for significance testing: {default_bin_counts}")
    
    sig_site_counts = []  # to store number of significant sites for each bin count
    
    for bin_count in default_bin_counts:
        n_sites = len(site_df_sorted)
        bins = get_bins(n_sites, bin_count)

        rep_fdrs = {}
        for rep_col in delta_rep_cols:
            rep_vals = site_df_sorted[rep_col].values
            rep_pvals = np.ones(n_sites)
            for b in bins:
                vals = rep_vals[b]
                # one-tailed test, applied separately in each biological replicate
                rep_pvals[b] = compute_bin_pvals(vals, side="right")
            rep_stat, _ = adjust_pvalues(rep_pvals, apply_bh=apply_bh)
            rep_fdrs[rep_col] = rep_stat

        significance_matrix = np.column_stack([rep_fdrs[c] for c in delta_rep_cols])
        sig_mask, _, _, _, _, _ = consensus_significance_from_replicates(
            delta_matrix=delta_matrix_sorted,
            significance_matrix=significance_matrix,
            delta_threshold=delta_rm_threshold,
            require_same_positive_direction=True,
        )
        sig_count = sig_mask.sum()
        sig_site_counts.append(sig_count)

        print(
            f"Bin count {bin_count}: {sig_count} significant sites "
            f"(same positive direction + mean ΔRm>{delta_rm_threshold:.2f} + statistical consistency rule)"
        )

    # 4) Plot bin count vs significant site number (with cubic spline smoothing)
    try:
        x = np.array(default_bin_counts)
        y = np.array(sig_site_counts)
        
        # Sort x, y for smoothing
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]
        y_sorted = y[sort_idx]
        
        # Generate smooth x
        x_smooth = np.linspace(x_sorted.min(), x_sorted.max(), 300)
        
        # Cubic spline smoothing
        spline = make_interp_spline(x_sorted, y_sorted, k=3)
        y_smooth = spline(x_smooth)
        
        plt.figure(figsize=(6,4))
        
        # Scatter points
        plt.scatter(x, y, label="Observed number", color="#89C9C8")
        
        # Smoothed curve
        plt.plot(x_smooth, y_smooth, color="#F9BEBB", linewidth=2, label="Smoothed spline curve")
        
        plt.xlabel("Bin count")
        plt.ylabel("Number of significant N-glycosites")
        plt.title("Significant N-glycosite number vs bin count")
        plt.legend(frameon=False)
        plt.tight_layout()
        
        plot_file = os.path.join(figures_dir, "11_significant_sites_vs_bin_count.tiff")
        plt.savefig(plot_file, dpi=300, format='tiff')
        plt.close()
        print(f"Scatter plot saved: {plot_file}")
    
    except Exception as e:
        print("Warning: failed to plot significant site scatter:", e)

    # ------------------------------
    # Step 11 (normal branch): User-specified bin size for significance testing
    # ------------------------------
    raw_bin_size_input = input(
        "Enter custom bin size (number of N-glycosites per window) for significance testing, e.g., 50: "
    ).strip()
    try:
        custom_bin_size = int(raw_bin_size_input)
        if custom_bin_size <= 0:
            raise ValueError("Bin size must be positive.")
    except Exception as e:
        raise ValueError(f"Invalid bin size: {e}")
    
    print(f"Using custom bin size = {custom_bin_size} sites per window")
    
    # 2) Compute bins for the user-specified bin size
    n_sites = len(site_df_sorted)
    custom_bins = []
    i = 0
    while i < n_sites:
        j = i + custom_bin_size
        if j >= n_sites:
            if len(custom_bins) > 0 and (n_sites - i) < custom_bin_size:
                custom_bins[-1] = np.concatenate([custom_bins[-1], np.arange(i, n_sites)])
            else:
                custom_bins.append(np.arange(i, n_sites))
            break
        else:
            custom_bins.append(np.arange(i, j))
        i = j
    
    # 3) Compute p-values per bin
    custom_rep_pvals = {}
    custom_rep_fdrs = {}
    for rep_col in delta_rep_cols:
        rep_vals = site_df_sorted[rep_col].values
        rep_pvals = np.ones(n_sites)
        for b in custom_bins:
            vals = rep_vals[b]
            rep_pvals[b] = compute_bin_pvals(vals, side="right")
        rep_stat, stat_label = adjust_pvalues(rep_pvals, apply_bh=apply_bh)
        custom_rep_pvals[rep_col] = rep_pvals
        custom_rep_fdrs[rep_col] = rep_stat

    custom_pval_matrix = np.column_stack([custom_rep_pvals[c] for c in delta_rep_cols])
    custom_significance_matrix = np.column_stack([custom_rep_fdrs[c] for c in delta_rep_cols])
    (
        sig_mask_custom,
        mean_delta_custom,
        same_sign_custom,
        delta_pass_custom,
        stat_all_10_custom,
        stat_any_05_custom,
    ) = consensus_significance_from_replicates(
        delta_matrix=delta_matrix_sorted,
        significance_matrix=custom_significance_matrix,
        delta_threshold=delta_rm_threshold,
        require_same_positive_direction=True,
    )
    
    # 6) Add results to the site table
    site_df_final = site_df_sorted.copy()
    site_df_final["Analysis_method"] = "normal_cv_binned_robust_sigma_right_tailed"
    site_df_final["Multiple_testing"] = "BH" if apply_bh else "none"
    site_df_final["p_value_rep1"] = custom_rep_pvals["ΔRm_rep1"]
    site_df_final["p_value_rep2"] = custom_rep_pvals["ΔRm_rep2"]
    site_df_final["p_value_rep3"] = custom_rep_pvals["ΔRm_rep3"]
    if apply_bh:
        site_df_final["FDR_BH_rep1"] = custom_rep_fdrs["ΔRm_rep1"]
        site_df_final["FDR_BH_rep2"] = custom_rep_fdrs["ΔRm_rep2"]
        site_df_final["FDR_BH_rep3"] = custom_rep_fdrs["ΔRm_rep3"]
        site_df_final["FDR_BH_max_across_reps"] = custom_significance_matrix.max(axis=1)
    site_df_final["p_value_max_across_reps"] = custom_pval_matrix.max(axis=1)
    site_df_final["Significance_metric"] = "FDR_BH" if apply_bh else "p_value"
    site_df_final["Significance_value_rep1"] = custom_rep_fdrs["ΔRm_rep1"]
    site_df_final["Significance_value_rep2"] = custom_rep_fdrs["ΔRm_rep2"]
    site_df_final["Significance_value_rep3"] = custom_rep_fdrs["ΔRm_rep3"]
    site_df_final["Same_positive_direction_across_reps"] = same_sign_custom
    site_df_final["Mean_DeltaRm"] = mean_delta_custom
    site_df_final["Mean_DeltaRm_gt_threshold"] = delta_pass_custom
    site_df_final["Stat_lt0p1_all_reps"] = stat_all_10_custom
    site_df_final["Stat_lt0p05_any_rep"] = stat_any_05_custom
    site_df_final["Significant"] = sig_mask_custom
    
    # 7) Save table
    custom_output_file = os.path.join(tables_dir, f"14_nglycosites_significance_normal_binSize{custom_bin_size}.csv")
    export_csv(site_df_final, custom_output_file, SITE_EXPORT_RENAME)
    n_sig = sig_mask_custom.sum()
    print(f"Consensus rule: same positive direction + mean ΔRm>{delta_rm_threshold:.2f} + all significance values<0.10 + at least one<0.05")
    print(f"Custom bin size {custom_bin_size}: {n_sig} significant sites saved to {custom_output_file}")

    # Diagnostic volcano-like plots for the normal branch.
    # The x-axis is MEAN DeltaRm because the effect-size threshold is defined on
    # mean DeltaRm. The y-axis shows the replicate-specific statistical metric.
    for rep_idx, rep_col in enumerate(delta_rep_cols, start=1):
        stat_col = f"Significance_value_rep{rep_idx}"
        volcano_df_rep = site_df_final[["ΔRm", stat_col, "Significant"]].copy()
        volcano_df_rep = volcano_df_rep.replace([np.inf, -np.inf], np.nan).dropna()
        if volcano_df_rep.shape[0] == 0:
            print(f"Warning: no valid values for normal-branch volcano plot of {rep_col}.")
            continue

        volcano_df_rep["neglog10Stat"] = -np.log10(np.clip(volcano_df_rep[stat_col].values, 1e-300, None))
        significant_mask = volcano_df_rep["Significant"].astype(bool)
        other_mask = ~significant_mask

        plt.figure(figsize=(6,5))
        plt.scatter(volcano_df_rep.loc[other_mask, "ΔRm"], volcano_df_rep.loc[other_mask, "neglog10Stat"], s=14, color="grey", alpha=0.6)
        plt.scatter(volcano_df_rep.loc[significant_mask, "ΔRm"], volcano_df_rep.loc[significant_mask, "neglog10Stat"], s=18, color="red", alpha=0.85, label="Consensus significant")
        plt.axvline(delta_rm_threshold, color="red", linestyle="--", linewidth=1)
        plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
        plt.xlabel("Mean ΔRm across biological replicates")
        plt.ylabel("-log10(FDR_BH)" if apply_bh else "-log10(p-value)")
        plt.title(f"Volcano-like diagnostic (normal branch, statistic from {rep_col})")
        plt.legend(frameon=False, fontsize=8)
        plt.tight_layout()
        volcano_plot_rep = os.path.join(figures_dir, f"12_volcano_normal_rep{rep_idx}.tiff")
        plt.savefig(volcano_plot_rep, dpi=300, format='tiff')
        plt.close()
        print(f"Normal-branch volcano plot saved: {volcano_plot_rep}")

else:
    print("Step10 strategy: replicate-level normality gate failed -> paired one-tailed test (normalized protein Rm vs normalized glycosite Rm) per site")

    site_df_non_normal = site_df_CV_pass.copy().reset_index(drop=True)
    needed_cols = [
        "Rm_unmod1_norm", "Rm_unmod2_norm", "Rm_unmod3_norm",
        "Rm_glyco1_norm", "Rm_glyco2_norm", "Rm_glyco3_norm", "ΔRm"
    ]
    missing_cols = [c for c in needed_cols if c not in site_df_non_normal.columns]
    if missing_cols:
        raise KeyError(f"Required columns missing for non-normal paired branch: {missing_cols}")

    paired_test_input = input(
        "ΔRm not normal. Select paired one-tailed test [wilcoxon/paired_t, default=wilcoxon]: "
    ).strip().lower()
    if paired_test_input == "":
        non_normal_method = "wilcoxon"
    elif paired_test_input in ["wilcoxon", "paired_t"]:
        non_normal_method = paired_test_input
    else:
        raise ValueError("Invalid non-normal paired test. Use 'wilcoxon' or 'paired_t'.")

    print(f"Paired test: {non_normal_method} on normalized Rm (one-tailed H1: glycosite Rm > protein Rm)")

    pvals = []
    n_pairs_used = []
    for _, row in site_df_non_normal.iterrows():
        protein_rm_norm = np.array([row["Rm_unmod1_norm"], row["Rm_unmod2_norm"], row["Rm_unmod3_norm"]], dtype=float)
        glycosite_rm_norm = np.array([row["Rm_glyco1_norm"], row["Rm_glyco2_norm"], row["Rm_glyco3_norm"]], dtype=float)

        p, n_pairs = paired_one_tailed_pvalue(
            protein_vals=protein_rm_norm,
            glycosite_vals=glycosite_rm_norm,
            method=non_normal_method,
            side="right"
        )
        pvals.append(float(p))
        n_pairs_used.append(int(n_pairs))

    pvals = np.array(pvals)
    significance_values, significance_label = adjust_pvalues(pvals, apply_bh=apply_bh)

    site_df_final = site_df_non_normal.copy()
    site_df_final["Analysis_method"] = f"paired_{non_normal_method}_right_tailed_normalized_rm"
    site_df_final["Paired_test"] = non_normal_method
    site_df_final["Multiple_testing"] = "BH" if apply_bh else "none"
    site_df_final["n_pairs_used"] = n_pairs_used
    site_df_final["p_value"] = pvals
    if apply_bh:
        site_df_final["FDR_BH"] = significance_values
    site_df_final["Significance_metric"] = "FDR_BH" if apply_bh else "p_value"
    site_df_final["Significance_value"] = significance_values
    site_df_final["Mean_DeltaRm_gt_threshold"] = site_df_final["ΔRm"] > delta_rm_threshold
    site_df_final["Significant"] = (
        site_df_final["Mean_DeltaRm_gt_threshold"]
        & (site_df_final["Significance_value"] < 0.05)
    )

    if apply_bh and non_normal_method == "wilcoxon" and len(n_pairs_used) > 0:
        max_pairs = int(np.max(n_pairs_used))
        if max_pairs > 0:
            min_theoretical_p = 1.0 / (2 ** max_pairs)
            bh_best_case = min_theoretical_p * len(site_df_final)
            print(f"Wilcoxon one-tailed p-value resolution with {max_pairs} pairs: min theoretical p ≈ {min_theoretical_p:.4g}")
            print(f"Approximate best-case BH floor (p_min * m): {bh_best_case:.4g}")
            if bh_best_case > 0.05:
                print("Warning: with current replicate count and test count, BH<0.05 may be unattainable even for strongest sites.")

    # Save non-normal branch result
    if non_normal_method == "wilcoxon":
        non_normal_output_file = os.path.join(tables_dir, "14_nglycosites_significance_paired_wilcoxon.csv")
    else:
        non_normal_output_file = os.path.join(tables_dir, "14_nglycosites_significance_paired_t.csv")
    export_csv(site_df_final, non_normal_output_file, SITE_EXPORT_RENAME)
    print(f"Non-normal branch results saved: {non_normal_output_file}")
    print(f"Significant up-stabilized sites: {site_df_final['Significant'].sum()}")

    # Top-25 site boxplots
    top_n = min(25, site_df_final.shape[0])
    top_df = site_df_final.sort_values("ΔRm", ascending=False).head(top_n).reset_index(drop=True)
    if top_n > 0:
        ncols = 5
        nrows = int(np.ceil(top_n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.0, nrows * 3.0), squeeze=False)
        axes = axes.flatten()
        rng = np.random.default_rng(2026)

        for i, (_, row) in enumerate(top_df.iterrows()):
            ax = axes[i]
            protein_rm_norm = np.array([row["Rm_unmod1_norm"], row["Rm_unmod2_norm"], row["Rm_unmod3_norm"]], dtype=float)
            glycosite_rm_norm = np.array([row["Rm_glyco1_norm"], row["Rm_glyco2_norm"], row["Rm_glyco3_norm"]], dtype=float)

            ax.boxplot(
                [protein_rm_norm], positions=[1], widths=0.5, showfliers=False, patch_artist=True,
                boxprops=dict(facecolor="#4C78A8", alpha=0.35),
                medianprops=dict(color="#4C78A8"), whiskerprops=dict(color="#4C78A8"), capprops=dict(color="#4C78A8")
            )
            ax.boxplot(
                [glycosite_rm_norm], positions=[2], widths=0.5, showfliers=False, patch_artist=True,
                boxprops=dict(facecolor="#F58518", alpha=0.35),
                medianprops=dict(color="#F58518"), whiskerprops=dict(color="#F58518"), capprops=dict(color="#F58518")
            )

            x1 = 1 + rng.normal(0, 0.04, size=len(protein_rm_norm))
            x2 = 2 + rng.normal(0, 0.04, size=len(glycosite_rm_norm))
            ax.scatter(x1, protein_rm_norm, s=16, color="#4C78A8", zorder=3)
            ax.scatter(x2, glycosite_rm_norm, s=16, color="#F58518", zorder=3)

            combined = np.concatenate([protein_rm_norm, glycosite_rm_norm])
            y_min = np.nanmin(combined)
            y_max = np.nanmax(combined)
            y_span = (y_max - y_min) if (y_max > y_min) else 1.0
            y_text = y_max + 0.18 * y_span
            ax.plot([1, 2], [y_text, y_text], color="black", linewidth=0.8)
            star = pval_to_stars(row["Significance_value"])
            stat_name = "FDR" if apply_bh else "p"
            ax.text(1.5, y_text + 0.03 * y_span, f"{stat_name}={row['Significance_value']:.2e} {star}", ha="center", va="bottom", fontsize=7)

            site_label = f"{row['Protein Accession']}|N{int(row['N-Glycosite'])}"
            ax.set_title(site_label, fontsize=7)
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["Protein", "N-glycosite"], fontsize=7)
            ax.tick_params(axis='y', labelsize=7)
            ax.set_ylabel("Normalized Rm", fontsize=7)

        for j in range(top_n, len(axes)):
            axes[j].axis("off")

        plt.tight_layout()
        top25_plot = os.path.join(figures_dir, f"13_top25_sites_normalized_rm_{non_normal_method}.tiff")
        plt.savefig(top25_plot, dpi=300, format='tiff')
        plt.close()
        print(f"Top-25 site boxplot saved: {top25_plot}")

    # Volcano-like plot for all sites
    volcano_df = site_df_final.copy()
    print("Non-normal volcano uses mean ΔRm across biological replicates (same ΔRm definition as site-level table).")
    volcano_df["neglog10Stat"] = -np.log10(np.clip(volcano_df["Significance_value"].values, 1e-300, None))
    up_mask = (volcano_df["ΔRm"] > delta_rm_threshold) & (volcano_df["Significance_value"] < 0.05)
    down_mask = (volcano_df["ΔRm"] < -delta_rm_threshold) & (volcano_df["Significance_value"] < 0.05)
    other_mask = ~(up_mask | down_mask)

    plt.figure(figsize=(6,5))
    plt.scatter(volcano_df.loc[other_mask, "ΔRm"], volcano_df.loc[other_mask, "neglog10Stat"], s=14, color="grey", alpha=0.6)
    plt.scatter(volcano_df.loc[up_mask, "ΔRm"], volcano_df.loc[up_mask, "neglog10Stat"], s=18, color="red", alpha=0.85, label="Up & significant")
    plt.scatter(volcano_df.loc[down_mask, "ΔRm"], volcano_df.loc[down_mask, "neglog10Stat"], s=18, color="blue", alpha=0.85, label="Down & significant")
    plt.axvline(delta_rm_threshold, color="red", linestyle="--", linewidth=1)
    plt.axvline(-delta_rm_threshold, color="blue", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
    plt.xlabel("ΔRm")
    plt.ylabel("-log10(FDR_BH)" if apply_bh else "-log10(p-value)")
    method_name = "Wilcoxon" if non_normal_method == "wilcoxon" else "paired t-test"
    plt.title(f"Volcano-like plot ({method_name}; {'BH' if apply_bh else 'no BH'})")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    volcano_plot = os.path.join(figures_dir, f"14_volcano_{non_normal_method}.tiff")
    plt.savefig(volcano_plot, dpi=300, format='tiff')
    plt.close()
    print(f"Volcano-like plot saved: {volcano_plot}")

# ------------------------------
# Final run metadata
# ------------------------------
stage("Stage 10/10 - Save run metadata")
run_parameters = {
    "pipeline": PIPELINE_NAME,
    "pipeline_version": PIPELINE_VERSION,
    "run_completed_at": datetime.now().isoformat(timespec="seconds"),
    "protein_csv": file,
    "peptide_csv": peptide_csv_file,
    "fasta_file": fasta_file,
    "control_sample": control_sample,
    "glycosite_sample_group": chosen_group,
    "motif_window_length": motif_len,
    "cv_threshold": float(cv_threshold),
    "cv_definition": "raw_Rm_SD_divided_by_mean_across_biological_replicates",
    "rm_normalization": "replicate_median_multiplicative_alignment",
    "rm_reference_median": float(rm_reference_median),
    "rm_correction_factors": {k: float(v) for k, v in CF.items()},
    "delta_rm_definition": "mean(glycosite_Rm_norm_rep - protein_Rm_norm_rep)",
    "delta_rm_threshold": float(delta_rm_threshold),
    "multiple_testing": "BH" if apply_bh else "none",
    "normality_alpha": NORMALITY_ALPHA_DEFAULT,
    "normal_branch_direction_rule": "all_replicate_DeltaRm_positive",
    "normal_branch_effect_size_rule": "mean_DeltaRm_gt_threshold",
    "n_glycosite_mass_shift_target": 0.98,
    "n_glycosite_mass_shift_tolerance": 0.01,
}
write_json(run_parameters, os.path.join(metadata_dir, "run_parameters.json"))

print(f"\nPipeline completed successfully. Outputs: {output_root}")
