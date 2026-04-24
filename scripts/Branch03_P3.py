#! /usr/env/bin python3

# -*- coding: utf-8 -*-
"""
Refined-TPP plus - Full Preprocessing & Rm Value Analysis Pipeline
---------------------------------------------------------
Created on: Fri Nov  7 11:59:36 2025
Authors: Haoru Song, Chenxin Li, Bin Fu 
Affiliation: Professor Haojie Lu's Lab, Fudan University, Shanghai, China.

Description
-----------
This script implements the complete 6-plex TMT-based quantitative proteomics data 
processing and Rm computation pipeline as designed for the Refined-TPP workflow.
It includes quantification quality control, technical replicate quality control, CV-based binning, robust statistics,
and differential protein identification.

Main Steps
-----------
1. **Automatic channel detection and mapping**
   - Detect TMT reporter ion channels (126-131) automatically.
   - Prompt the user to assign each channel to a biological sample.
   - Ask the user to specify which sample serves as Vehicle (control).

2. **Quantification overview and filtering**
   - Count total quantified proteins (all proteins with any quantitative signal).
   - Count proteins with ≥2 unique peptides.
   - Count proteins quantified in each sample (≥2 unique peptides).
   - Plot a Venn diagram or UpSet plot showing overlap of quantified proteins across samples.
   - Apply strict filtering: retain only proteins with all non-zero TMT intensities 
     in each sample (core 126/127/128 channels must also be non-zero).

3. **Rm computation and normalization**
   - For each protein and sample, compute three replicate Rm:
         Rm1 = 129 / 126, Rm2 = 130 / 127, Rm3 = 131 / 128
   - Compute geometric_mean_Rm as:
         geometric_mean_Rm = exp(mean(log(Rm1, Rm2, Rm3), skipna=True))
     where the log-average protects against outliers and non-positive values.
   - Compute CV for the three raw replicate Rm (CV = std/mean), used for quality control.
   - Median-normalize geometric_mean_Rm values to the Vehicle sample:
         adjusted_Rm = geometric_mean_Rm * (median_geometric_mean_Rm(Vehicle) / median_geometric_mean_Rm(sample))
     to remove global scale differences across TMT channels.
   - **CV-based QC:** 
       • Generate CV distribution plots for each sample (including Vehicle), 
         sorted by descending CV.  
       • Mark cutoffs CV ≤ 20% and CV ≤ 30%, and report percentages of proteins 
         passing these thresholds.  
       • Proteins with CV > 30% in any condition are excluded from subsequent binning and differential analysis.

4. **ΔRm and statistical testing**
   - Compute ΔRm = adjusted_Rm(sample) - adjusted_Rm(Vehicle) for each protein.
   - Perform per-sample normality check on ΔRm distribution.
   - If ΔRm is approximately normal:
       • Use CV-based binning for robust p-value estimation (TPP-TR style):
         proteins are ordered by representative CV (max of Vehicle and experimental sample CVs),
         and per-bin statistics use robust left/right sigma from percentiles (15.87%, 50%, 84.13%).
       • Apply Benjamini-Hochberg (BH) correction for multiple testing (FDR).
   - If ΔRm is NOT approximately normal:
       • Do not report inferential p/FDR due to technical-replicate-only limitation.
       • Provide empirical outlier interpretation using top 5% and bottom 5% ΔRm tails.
       • Generate ΔRm ranked plot (small to large) with dual-tail color coding.

5. **Differential protein output**
   - Normal branch: identify proteins with |ΔRm| > 0.1 and FDR < 0.05.
   - Non-normal branch: identify empirical outliers (top 5% and bottom 5% ΔRm) without p-value/FDR.
   - Output complete results with signed ΔRm values and branch-specific decision basis.

Notes / Remarks
----------------
- CV is computed on raw replicate Rm, not on log-averaged Rm.
- Representative CV for binning is the larger of experimental vs Vehicle.
- Proteins failing CV ≤ 30% in any condition are removed before differential analysis.
- The three Rm values per file are technical replicates; biological replicates are cross-file.
- If ΔRm is non-normal, results are empirical ranking/outlier references rather than inferential statistics.
- All steps are designed to be transparent and allow user customization.

Output
------
- Summary tables of protein counts (total, filtered, per-sample)
- Venn diagram or UpSet plot of quantified proteins (before strict filtering)
- CV distribution plots per sample (with percentages passing 20%/30% cutoffs)
- DeltaRm normality diagnostic plots (Histogram+KDE, Q-Q, P-P) per sample
- Normal branch: CSV with Rm, geometric_mean_Rm, adjusted_Rm, ΔRm, p-value, and FDR
- Normal branch volcano plot (DeltaRm vs -log10(FDR)) per sample
- Non-normal branch: CSV with ΔRm empirical top/bottom 5% outlier labels (p/FDR left as NaN)
- ΔRm ranked empirical plot for non-normal branch
- Differential protein/outlier list with branch-specific criterion labels

Usage
-----
Run the script and follow the on-screen prompts to:
   • Input the path of DBsearch .csv file (please use PEAKS ONLINE)
   • Assign TMT channels to samples
   • Specify which sample is Vehicle (control)

Dependencies
------------
- Python 3.8+
- pandas, numpy, matplotlib, scipy, statsmodels, matplotlib-venn, upsetplot

Installation
------------
    pip install pandas numpy matplotlib scipy statsmodels matplotlib-venn upsetplot

Reference
---------
When using this official script, please cite:
"J. Am. Chem. Soc. 2025, 147, 27, 24127-24139." (continue updating)

---------------------------------------------------------
"""

import os
import re
from collections import defaultdict

import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")

import logging
logging.basicConfig(level = logging.INFO, format = "%(asctime)s - %(levelname)s - %(message)s")
logging.getLogger('numexpr').setLevel(logging.WARNING)
logger = logging.getLogger()

try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    from scipy.stats import norm
    from statsmodels.stats.multitest import multipletests
    from matplotlib_venn import venn2, venn3
    from upsetplot import from_contents, UpSet
except Exception as e:
    raise ImportError(f"{e}\nPlease install required packages: pip install pandas numpy matplotlib scipy statsmodels matplotlib-venn upsetplot")

HAVE_VENN = True
HAVE_UPSET = True

# ------------------------------
# Helper functions
# ------------------------------
def normalize_path(path):
    return os.path.normpath(os.path.abspath(path))

def detect_unique_col(df):
    """Detect #Unique-like column name. Prefer exact '#Unique'."""
    if "#Unique" in df.columns:
        return "#Unique"
    candidates = [c for c in df.columns if re.search(r"#\s*Unique|#Unique|Unique\s*Pept|UniquePept|Unique", c, flags=re.I)]
    return candidates[0] if candidates else None

def detect_sample_blocks(df):
    """
    Detect 'Sample N' grouped columns. Returns ordered sample indices and mapping.
    If headers contain 'Sample <N>' pattern this collects columns for each sample index.
    """
    pattern = re.compile(r"Sample\s*(\d+)\b", flags=re.I)
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
    """Return a short reporter fragment (e.g. 'TMT-126')."""
    m = re.search(r"(126|127|128|129|130|131)", colname)
    if m:
        return f"TMT-{m.group(1)}"
    return colname

def robust_bin_pvals(deltas, cv_sort_order, bin_size=300):
    """
    Compute robust p-values per bin sorted by CV.
    """
    n = len(deltas)
    pvals = np.ones(n)
    idxs = np.array(cv_sort_order)
    bins = []
    i = 0
    while i < n:
        j = i + bin_size
        if j >= n:
            if len(bins) > 0 and (n - i) < bin_size:
                bins[-1] = np.concatenate([bins[-1], idxs[i:n]])
            else:
                bins.append(idxs[i:n])
            break
        else:
            bins.append(idxs[i:j])
        i = j

    for b in bins:
        vals = deltas.loc[b].values
        if len(vals) == 0:
            continue
        p15 = np.percentile(vals, 15.87)
        p50 = np.percentile(vals, 50.0)
        p84 = np.percentile(vals, 84.13)
        left_sigma = p50 - p15
        right_sigma = p84 - p50
        left_sigma = left_sigma if left_sigma > 0 else np.nan
        right_sigma = right_sigma if right_sigma > 0 else np.nan

        for idx in b:
            val = deltas.loc[idx]
            if np.isnan(val):
                pv = 1.0
            else:
                if val >= p50:
                    if np.isnan(right_sigma) or right_sigma == 0:
                        pv = 1.0
                    else:
                        z = (val - p50) / right_sigma
                        pv = 2 * norm.sf(abs(z))
                else:
                    if np.isnan(left_sigma) or left_sigma == 0:
                        pv = 1.0
                    else:
                        z = (p50 - val) / left_sigma
                        pv = 2 * norm.sf(abs(z))
            pvals[np.where(idxs == idx)[0][0]] = pv

    pval_series = pd.Series(index=idxs, data=pvals)
    aligned = pd.Series(1.0, index=deltas.index)
    aligned.loc[pval_series.index] = pval_series.values
    return aligned


def assess_delta_normality(delta_series):
    """Assess normality for one DeltaRm series and return (is_normal, message)."""
    arr = pd.Series(delta_series).replace([np.inf, -np.inf], np.nan).dropna().values
    n = len(arr)

    if n < 3:
        return False, f"n={n}: too few proteins for normality test, treat as non-normal empirical mode"

    if n <= 50:
        stat, p = stats.shapiro(arr)
        is_normal = p >= 0.05
        return is_normal, f"Shapiro-Wilk W={stat:.4f}, p={p:.4g}"

    if n <= 300:
        mean_val = np.mean(arr)
        std_val = np.std(arr, ddof=1)
        if not np.isfinite(std_val) or std_val <= 0:
            return False, "KS skipped because std<=0, treat as non-normal empirical mode"
        stat, p = stats.kstest(arr, 'norm', args=(mean_val, std_val))
        is_normal = p >= 0.05
        return is_normal, f"Kolmogorov-Smirnov D={stat:.4f}, p={p:.4g}"

    skew_val = stats.skew(arr)
    kurt_val = stats.kurtosis(arr, fisher=True)
    is_normal = (abs(skew_val) < 3) and (abs(kurt_val) < 10)
    return is_normal, f"Large-sample rule: skew={skew_val:.2f}, kurtosis={kurt_val:.2f}"


def plot_delta_normality_diagnostics(delta_series, sample_name, out_prefix):
    """Save DeltaRm normality diagnostic plots (Histogram+KDE, Q-Q, P-P)."""
    arr = pd.Series(delta_series).replace([np.inf, -np.inf], np.nan).dropna().values
    if len(arr) < 3:
        return []

    out_files = []

    # Histogram + KDE
    plt.figure(figsize=(6.5, 4.2))
    plt.hist(arr, bins=30, density=True, color="#A7AED2", alpha=0.75, edgecolor="white", linewidth=0.4)
    try:
        pd.Series(arr).plot(kind="kde", color="#374151", linewidth=1.5)
    except Exception:
        pass
    plt.xlabel("DeltaRm")
    plt.ylabel("Density")
    plt.title(f"DeltaRm distribution with KDE - {sample_name}")
    plt.tight_layout()
    hist_kde_file = f"{out_prefix}_DeltaRm_hist_kde_{sample_name}.tiff"
    plt.savefig(hist_kde_file, dpi=300, format='tiff')
    plt.close()
    out_files.append(hist_kde_file)

    # Q-Q
    plt.figure(figsize=(6, 4))
    stats.probplot(arr, dist="norm", plot=plt)
    plt.title(f"DeltaRm Q-Q plot - {sample_name}")
    plt.tight_layout()
    qq_file = f"{out_prefix}_DeltaRm_QQ_{sample_name}.tiff"
    plt.savefig(qq_file, dpi=300, format='tiff')
    plt.close()
    out_files.append(qq_file)

    # P-P
    mean_val = np.mean(arr)
    std_val = np.std(arr, ddof=1)
    if np.isfinite(std_val) and std_val > 0:
        sorted_arr = np.sort(arr)
        n = len(sorted_arr)
        empirical_cdf = (np.arange(1, n + 1) - 0.5) / n
        theoretical_cdf = norm.cdf(sorted_arr, loc=mean_val, scale=std_val)

        plt.figure(figsize=(6, 4))
        plt.scatter(theoretical_cdf, empirical_cdf, s=16, alpha=0.75, color="#4B5563")
        plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1, color="#D7301F")
        plt.xlabel("Theoretical CDF (Normal)")
        plt.ylabel("Empirical CDF")
        plt.title(f"DeltaRm P-P plot - {sample_name}")
        plt.tight_layout()
        pp_file = f"{out_prefix}_DeltaRm_PP_{sample_name}.tiff"
        plt.savefig(pp_file, dpi=300, format='tiff')
        plt.close()
        out_files.append(pp_file)

    return out_files


def plot_delta_volcano(delta_series, fdr_series, sample_name, out_prefix, delta_cutoff=0.1, fdr_cutoff=0.05):
    """Plot volcano for normal-branch DeltaRm analysis."""
    plot_df = pd.DataFrame({
        "DeltaRm": pd.Series(delta_series),
        "FDR": pd.Series(fdr_series)
    }).replace([np.inf, -np.inf], np.nan).dropna()

    if plot_df.empty:
        return None

    plot_df["neglog10FDR"] = -np.log10(np.clip(plot_df["FDR"].values, 1e-300, None))
    up_mask = (plot_df["DeltaRm"] > delta_cutoff) & (plot_df["FDR"] < fdr_cutoff)
    down_mask = (plot_df["DeltaRm"] < -delta_cutoff) & (plot_df["FDR"] < fdr_cutoff)
    other_mask = ~(up_mask | down_mask)

    plt.figure(figsize=(6.2, 5))
    plt.scatter(plot_df.loc[other_mask, "DeltaRm"], plot_df.loc[other_mask, "neglog10FDR"],
                s=14, color="#9E9E9E", alpha=0.65)
    plt.scatter(plot_df.loc[up_mask, "DeltaRm"], plot_df.loc[up_mask, "neglog10FDR"],
                s=18, color="#D7301F", alpha=0.85, label="Up & significant")
    plt.scatter(plot_df.loc[down_mask, "DeltaRm"], plot_df.loc[down_mask, "neglog10FDR"],
                s=18, color="#2C7FB8", alpha=0.85, label="Down & significant")

    plt.axvline(delta_cutoff, color="#D7301F", linestyle="--", linewidth=1)
    plt.axvline(-delta_cutoff, color="#2C7FB8", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(fdr_cutoff), color="black", linestyle="--", linewidth=1)
    plt.xlabel("DeltaRm")
    plt.ylabel("-log10(FDR)")
    plt.title(f"DeltaRm volcano plot - {sample_name}")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()

    out_file = f"{out_prefix}_DeltaRm_volcano_{sample_name}.tiff"
    plt.savefig(out_file, dpi=300, format='tiff')
    plt.close()
    return out_file


def empirical_tail_masks(delta_series, tail_frac=0.05):
    """Return low/high tail masks and tail size for empirical non-normal branch."""
    series = pd.Series(delta_series)
    low_mask = pd.Series(False, index=series.index)
    high_mask = pd.Series(False, index=series.index)

    valid = series.replace([np.inf, -np.inf], np.nan).dropna().sort_values()
    n_valid = len(valid)
    if n_valid == 0:
        return low_mask, high_mask, 0

    n_tail = max(1, int(np.ceil(n_valid * tail_frac)))
    low_mask.loc[valid.index[:n_tail]] = True
    high_mask.loc[valid.index[-n_tail:]] = True
    return low_mask, high_mask, n_tail


def plot_delta_rank_empirical(delta_series, sample_name, out_prefix):
    """Plot DeltaRm ranked (small->large) for empirical non-normal interpretation."""
    ranked = pd.Series(delta_series).replace([np.inf, -np.inf], np.nan).dropna().sort_values()
    if ranked.empty:
        return None

    colors = np.where(
        ranked.values < -0.1,
        "#2C7FB8",
        np.where(ranked.values > 0.1, "#D7301F", "#9E9E9E")
    )

    x = np.arange(1, len(ranked) + 1)
    n_low = int((ranked.values < -0.1).sum())
    n_high = int((ranked.values > 0.1).sum())

    plt.figure(figsize=(8, 4.5))
    plt.scatter(x, ranked.values, c=colors, s=14, alpha=0.9, linewidths=0)
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.axhline(0.1, color="#D7301F", linestyle="--", linewidth=1)
    plt.axhline(-0.1, color="#2C7FB8", linestyle="--", linewidth=1)
    plt.xlabel("Proteins ranked by DeltaRm (small -> large)")
    plt.ylabel("DeltaRm")
    plt.title(f"DeltaRm ranked empirical plot - {sample_name}")
    plt.text(0.02 * len(ranked), np.nanmax(ranked.values), f"< -0.1: {n_low}", color="#2C7FB8", fontsize=8)
    plt.text(0.02 * len(ranked), np.nanmax(ranked.values) - 0.08 * (np.nanmax(ranked.values) - np.nanmin(ranked.values) + 1e-12),
             f"> 0.1: {n_high}", color="#D7301F", fontsize=8)
    plt.tight_layout()

    out_file = f"{out_prefix}_DeltaRm_ranked_empirical_{sample_name}.tiff"
    plt.savefig(out_file, dpi=300, format='tiff')
    plt.close()
    return out_file

# ------------------------------
# 1. Read file & detect sample blocks
# ------------------------------
file = input("Please provide path of TMT quantification Database search CSV file: ").strip()
file = normalize_path(file)
if not os.path.exists(file):
    raise FileNotFoundError(file)

df = pd.read_csv(file)
protein_col = next((c for c in ["Protein Group", "ProteinGroup", "Protein_Group", "Accession", "Accession "] if c in df.columns), df.columns[0])

sample_order, sample_cols_map = detect_sample_blocks(df)
if not sample_order:
    tmt_cols = [c for c in df.columns if re.search(r"(126|127|128|129|130|131)", c)]
    if not tmt_cols:
        raise ValueError("No TMT reporter-like columns (126..131) detected.")
    num_samples = int(input("Enter number of biological samples: ").strip())
    channels_per_sample = len(tmt_cols) // num_samples
    sample_order = list(range(1, num_samples+1))
    sample_cols_map = {s: tmt_cols[i*channels_per_sample:(i+1)*channels_per_sample] for i,s in enumerate(sample_order)}

sample_names = [input(f"Please enter name for Sample {s} (use 'Vehicle' for control): ").strip() for s in sample_order]

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

# ------------------------------
# 2. Quantification QC and Venn / UpSet
# ------------------------------
unique_col = detect_unique_col(df)
if unique_col is None:
    raise ValueError("Could not detect a '#Unique' column.")

quantified_any = (df[all_reporter_cols].fillna(0) != 0).any(axis=1)
unique_ge2 = (df[unique_col].fillna(0) >= 2)

per_sample_sets = {}
for sample in sample_names:
    cols = sample_block_cols[sample]
    has_signal = (df[cols].fillna(0) != 0).any(axis=1)
    ok = has_signal & unique_ge2
    per_sample_sets[sample] = set(df.index[ok])

if len(sample_names) < 4 and HAVE_VENN:
    plt.figure(figsize=(6,6))
    if len(sample_names) == 2:
        venn2([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
    elif len(sample_names) == 3:
        venn3([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
    plt.title("Quantified proteins (unique>=2 & any channel>0) Venn")
    plt.savefig(os.path.splitext(file)[0]+"_quantified_venn.tiff", dpi=300, bbox_inches='tight', format='tiff')
    plt.close()
elif HAVE_UPSET:
    upset_data = from_contents(per_sample_sets)
    plt.figure(figsize=(8,6))
    UpSet(upset_data, show_counts=True).plot()
    plt.savefig(os.path.splitext(file)[0]+"_quantified_upset.tiff", dpi=300, bbox_inches='tight', format='tiff')
    plt.close()

# ------------------------------
# 3. Strict filtering (all channels non-zero & #Unique>=2)
# ------------------------------
final_mask = pd.Series(True, index=df.index)
for sample in sample_names:
    cols = sample_block_cols[sample]
    final_mask &= (df[cols].fillna(0) != 0).all(axis=1)
final_mask &= (df[unique_col].fillna(0) >= 2)
df_final = df.loc[final_mask].copy()
logger.info(f"Protein number after strict filtering: {df_final.shape[0]}")

# ------------------------------
# 4. Compute Rm, CV, adjusted-Rm
# ------------------------------
for sample in sample_names:
    cols = sample_block_cols[sample]
    reporter_map = {int(re.search(r"(126|127|128|129|130|131)", c).group(1)):c for c in cols}
    Rm1 = df_final[reporter_map[129]] / df_final[reporter_map[126]]
    Rm2 = df_final[reporter_map[130]] / df_final[reporter_map[127]]
    Rm3 = df_final[reporter_map[131]] / df_final[reporter_map[128]]
    df_final[f"Rm1_{sample}"] = Rm1.replace([np.inf,-np.inf],np.nan)
    df_final[f"Rm2_{sample}"] = Rm2.replace([np.inf,-np.inf],np.nan)
    df_final[f"Rm3_{sample}"] = Rm3.replace([np.inf,-np.inf],np.nan)
    # log-average
    with np.errstate(divide='ignore', invalid='ignore'):
        log_mean = np.log(df_final[[f"Rm1_{sample}", f"Rm2_{sample}", f"Rm3_{sample}"]]).mean(axis=1, skipna=True)
    df_final[f"geometric_mean_Rm_{sample}"] = np.exp(log_mean)
    # CV from raw Rm
    df_final[f"CV_Rm1-3_{sample}"] = df_final[[f"Rm1_{sample}", f"Rm2_{sample}", f"Rm3_{sample}"]].std(axis=1, ddof=1)/df_final[[f"Rm1_{sample}", f"Rm2_{sample}", f"Rm3_{sample}"]].mean(axis=1)

vehicle = "Vehicle" if "Vehicle" in sample_names else input("Enter Vehicle sample name: ").strip()
vehicle_median = df_final[f"geometric_mean_Rm_{vehicle}"].median()
for sample in sample_names:
    cf = 1.0 if sample==vehicle else float(vehicle_median/df_final[f"geometric_mean_Rm_{sample}"].median())
    df_final[f"adjusted_Rm_{sample}"] = df_final[f"geometric_mean_Rm_{sample}"] * cf

# ------------------------------
# 5. CV QC: plot CV distributions per sample (strict-filtered proteins)
# ------------------------------
cv_cols = [f"CV_Rm1-3_{s}" for s in sample_names]

for s in sample_names:
    plt.figure(figsize=(6,4))
    cv_series = df_final[f"CV_Rm1-3_{s}"].sort_values(ascending=False).reset_index(drop=True)
    rank = cv_series.index + 1
    plt.scatter(rank, cv_series.values, color="#5E556A", alpha=0.7, s=20)
    prop_20 = (cv_series <= 0.2).sum() / len(cv_series) * 100
    prop_30 = (cv_series <= 0.3).sum() / len(cv_series) * 100
    plt.text(0.7*len(cv_series), 0.205, f"≤20%: {prop_20:.1f}%", fontsize=10)
    plt.text(0.7*len(cv_series), 0.305, f"≤30%: {prop_30:.1f}%", fontsize=10)
    plt.xlabel("Proteins sorted by CV (high→low)")
    plt.ylabel("CV")
    plt.title(f"CV distribution - {s} (strict-filtered proteins)")
    plt.tight_layout()
    plt.savefig(os.path.splitext(file)[0]+f"_CV_{s}.tiff", dpi=300, format='tiff')
    plt.close()
# CV QC filtering: remove proteins with CV>30% in any sample (after plotting)
df_final = df_final[(df_final[cv_cols] <= 0.3).all(axis=1)]

# ------------------------------
# 6. ΔRm computation
# ------------------------------
for s in sample_names:
    if s==vehicle: continue
    df_final[f"DeltaRm_{s}"] = df_final[f"adjusted_Rm_{s}"] - df_final[f"adjusted_Rm_{vehicle}"]

# ------------------------------
# 7. ΔRm-Normality-dependent differential analysis
# ------------------------------
delta_cols = [c for c in df_final.columns if c.startswith("DeltaRm_")]
base_prefix = os.path.splitext(file)[0]

for delta_col in delta_cols:
    sample = delta_col.replace("DeltaRm_", "")
    deltas = df_final[delta_col]

    pcol = f"{delta_col}_pval"
    fcol = f"{delta_col}_FDR"
    mode_col = f"{delta_col}_analysis_mode"
    low_col = f"{delta_col}_EmpiricalLow5pct"
    high_col = f"{delta_col}_EmpiricalHigh5pct"
    outlier_col = f"{delta_col}_EmpiricalOutlier"
    tailn_col = f"{delta_col}_EmpiricalTailN"

    normality_plots = plot_delta_normality_diagnostics(deltas, sample, base_prefix)
    if normality_plots:
        logger.info(f"{delta_col}: normality diagnostics saved: {', '.join(normality_plots)}")

    is_normal, normality_msg = assess_delta_normality(deltas)
    branch = "normal_statistical" if is_normal else "non_normal_empirical"
    df_final[mode_col] = branch
    logger.info(f"{delta_col}: {normality_msg}; branch={branch}")

    if is_normal:
        cv_cols_pair = [f"CV_Rm1-3_{sample}", f"CV_Rm1-3_{vehicle}"]
        # Representative CV: max of treatment vs vehicle.
        cv_rep = df_final[cv_cols_pair].max(axis=1).fillna(np.inf)
        order = list(cv_rep.sort_values(ascending=False).index)
        bin_pvals = robust_bin_pvals(deltas, order, bin_size=300)
        pvals_arr = np.where(np.isnan(bin_pvals.loc[df_final.index].values), 1.0, bin_pvals.loc[df_final.index].values)
        _, fdrs, _, _ = multipletests(pvals_arr, method='fdr_bh')

        df_final[pcol] = pvals_arr
        df_final[fcol] = fdrs
        df_final[low_col] = False
        df_final[high_col] = False
        df_final[outlier_col] = False
        df_final[tailn_col] = 0

        volcano_plot = plot_delta_volcano(df_final[delta_col], df_final[fcol], sample, base_prefix)
        if volcano_plot is not None:
            logger.info(f"{delta_col}: volcano plot saved: {volcano_plot}")
    else:
        low_mask, high_mask, n_tail = empirical_tail_masks(deltas, tail_frac=0.05)
        outlier_mask = low_mask | high_mask

        df_final[pcol] = np.nan
        df_final[fcol] = np.nan
        df_final[low_col] = low_mask.reindex(df_final.index, fill_value=False)
        df_final[high_col] = high_mask.reindex(df_final.index, fill_value=False)
        df_final[outlier_col] = outlier_mask.reindex(df_final.index, fill_value=False)
        df_final[tailn_col] = n_tail

        rank_plot = plot_delta_rank_empirical(deltas, sample, base_prefix)
        if rank_plot is not None:
            logger.info(f"{delta_col}: empirical ranked plot saved: {rank_plot}")
        logger.info(
            f"{delta_col}: non-normal branch uses empirical top/bottom 5% only; "
            "p-value/FDR are intentionally omitted due to technical-replicate-only limitation"
        )

# ------------------------------
# 8. Differential proteins/outliers output (full matrix)
# ------------------------------
output_rows = []
for idx,row in df_final.iterrows():
    sig_any = False
    out = {c:row.get(c,np.nan) for c in df.columns}
    for s in sample_names:
        out[f"Rm1_{s}"] = row.get(f"Rm1_{s}",np.nan)
        out[f"Rm2_{s}"] = row.get(f"Rm2_{s}",np.nan)
        out[f"Rm3_{s}"] = row.get(f"Rm3_{s}",np.nan)
        out[f"geometric_mean_Rm_{s}"] = row.get(f"geometric_mean_Rm_{s}",np.nan)
        out[f"adjusted_Rm_{s}"] = row.get(f"adjusted_Rm_{s}",np.nan)
        out[f"CV_Rm1-3_{s}"] = row.get(f"CV_Rm1-3_{s}",np.nan)
        if s==vehicle:
            out[f"DeltaRm_{s}"] = np.nan
            out[f"analysis_mode_{s}"] = np.nan
            out[f"pval_{s}"] = np.nan
            out[f"FDR_{s}"] = np.nan
            out[f"EmpiricalLow5pct_{s}"] = False
            out[f"EmpiricalHigh5pct_{s}"] = False
            out[f"EmpiricalOutlier_{s}"] = False
            out[f"Sig_basis_{s}"] = np.nan
            out[f"Sig_{s}"] = False
        else:
            dcol = f"DeltaRm_{s}"
            pcol = f"{dcol}_pval"
            fcol = f"{dcol}_FDR"
            mode_col = f"{dcol}_analysis_mode"
            low_col = f"{dcol}_EmpiricalLow5pct"
            high_col = f"{dcol}_EmpiricalHigh5pct"
            outlier_col = f"{dcol}_EmpiricalOutlier"

            delta_val = row.get(dcol,np.nan)
            mode_val = row.get(mode_col, "normal_statistical")
            mode_val = mode_val if pd.notna(mode_val) else "normal_statistical"

            emp_low = bool(row.get(low_col, False)) if pd.notna(row.get(low_col, False)) else False
            emp_high = bool(row.get(high_col, False)) if pd.notna(row.get(high_col, False)) else False
            emp_outlier = bool(row.get(outlier_col, False)) if pd.notna(row.get(outlier_col, False)) else False

            out[dcol] = delta_val
            out[f"analysis_mode_{s}"] = mode_val
            out[f"EmpiricalLow5pct_{s}"] = emp_low
            out[f"EmpiricalHigh5pct_{s}"] = emp_high
            out[f"EmpiricalOutlier_{s}"] = emp_outlier

            sig_flag = False
            if mode_val == "normal_statistical":
                fdr_val = row.get(fcol, np.nan)
                out[f"pval_{s}"] = row.get(pcol, np.nan)
                out[f"FDR_{s}"] = fdr_val
                out[f"Sig_basis_{s}"] = "statistical(|DeltaRm|>0.1 & FDR<0.05)"
                if pd.notna(delta_val) and pd.notna(fdr_val):
                    if abs(delta_val) > 0.1 and fdr_val < 0.05:
                        sig_flag = True
                        sig_any = True
            else:
                out[f"pval_{s}"] = np.nan
                out[f"FDR_{s}"] = np.nan
                out[f"Sig_basis_{s}"] = "empirical(top/bottom 5% DeltaRm; no p/FDR)"
                if emp_outlier:
                    sig_flag = True
                    sig_any = True

            out[f"Sig_{s}"] = sig_flag
    if sig_any:
        output_rows.append(out)

output_df = pd.DataFrame(output_rows)
out_name = os.path.splitext(file)[0]+"_Differential_proteins_RefinedTPP_plus_full_matrix.csv"
output_df.to_csv(out_name, index=False)
logger.info(f"Output saved: {out_name}")

# ------------------------------
# 9. Save full results table
# ------------------------------
full_outfile = os.path.splitext(file)[0] + "_RefinedTPP_plus_full_Rm_results.csv"
df_final.to_csv(full_outfile, index=False)
logger.info(f"Full Rm results saved: {full_outfile}")

logger.info("Refined-TPP plus processing completed successfully.")
