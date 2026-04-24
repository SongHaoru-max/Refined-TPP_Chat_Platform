#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Refined-TPP plus - Multi-BioRep + TechRep Analysis Pipeline
-----------------------------------------------------------
Created on: Fri Nov  7 11:59:36 2025
Authors: Haoru Song, Chenxin Li, Bin Fu
Affiliation: Professor Haojie Lu's Lab, Fudan University, Shanghai, China.

Description
-----------
This script extends the Refined-TPP workflow from a single-file (technical-replicate-only)
analysis to a multi-file biological-replicate workflow. Each biological replicate is supplied
as one CSV file, and each file still contains technical replicates (TMT 126-131 channels)
for each experimental condition.

Key updates implemented in this version:
1) Multi-biological-replicate input (at least 3 files)
2) Conservative protein intersection across biological replicates (no missing-value imputation)
3) Batch correction by "median of vehicle medians" across biological replicates
4) Global normality decision:
   - If all (BioRep x Condition) DeltaRm distributions are normal:
       use per-replicate robust-bin p/FDR, then integrate by replicate-level FDR rule
   - If any (BioRep x Condition) is non-normal:
       use paired test across biological replicates on adjusted Rm (default: paired t-test,
       optional: Wilcoxon)
5) Per-BioRep condition naming is supported, with mapping to BR1 canonical conditions

Differential rules
------------------
- Absolute DeltaRm threshold remains fixed at 0.1.
- DeltaRm pass requires abs(DeltaRm) > 0.1 in all biological replicates, and
    DeltaRm must have the same sign in all biological replicates.
- Two-layer CV QC is applied before differential testing using raw Rm values:
    (i) technical-replicate CV QC within each BioRep (default cutoff 0.3), and
    (ii) biological-replicate CV QC computed from pooled technical-replicate
         raw Rm values across all BioReps for each condition (default cutoff 0.4).
    Both cutoffs can be set by user input.
- CV-QC is global for differential analysis: proteins must pass two-layer CV-QC
    across all BioReps and all conditions before any multi-condition testing.
- Global-normal branch significance requires:
    all replicate FDR < 0.1 and at least one replicate FDR < 0.05.
- Global-non-normal branch significance requires:
    paired-test FDR < 0.05.
- In global-non-normal branch, the displayed representative DeltaRm is
    arithmetic mean of per-BioRep DeltaRm values.

Dependencies
------------
- Python 3.8+
- pandas, numpy, matplotlib, scipy, statsmodels, matplotlib-venn, upsetplot

Installation
------------
    pip install pandas numpy matplotlib scipy statsmodels matplotlib-venn upsetplot
"""

import os
import re
from collections import defaultdict

import warnings

warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")

import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logging.getLogger("numexpr").setLevel(logging.WARNING)
logger = logging.getLogger()

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib_venn import venn2, venn3
    from scipy import stats
    from scipy.stats import norm
    from statsmodels.stats.multitest import multipletests
    from upsetplot import UpSet, from_contents
except Exception as e:
    raise ImportError(
        f"{e}\nPlease install required packages: "
        "pip install pandas numpy matplotlib scipy statsmodels matplotlib-venn upsetplot"
    )

HAVE_VENN = True
HAVE_UPSET = True

DEFAULT_TECH_CV_QC_CUTOFF = 0.3
DEFAULT_BIO_CV_QC_CUTOFF = 0.4
TECH_CV_PREVIEW_CANDIDATES = (0.3, 0.35, 0.4)
BIO_CV_PREVIEW_CANDIDATES = (0.3, 0.35, 0.4, 0.45, 0.5)


# ------------------------------
# Helper functions
# ------------------------------
def normalize_path(path):
    return os.path.normpath(os.path.abspath(path))


def parse_csv_paths(raw_text):
    """Parse user input paths split by ';' or ',' and normalize."""
    parts = [p.strip().strip('"').strip("'") for p in re.split(r"[;,]", raw_text) if p.strip()]
    return [normalize_path(p) for p in parts]


def ask_float_with_default(prompt, default, min_value=0.0):
    """Prompt user for a float cutoff with default and basic range validation."""
    while True:
        raw = input(prompt).strip()
        if raw == "":
            return float(default)

        try:
            value = float(raw)
        except ValueError:
            logger.warning(f"Invalid numeric input '{raw}'. Please enter a valid number.")
            continue

        if not np.isfinite(value):
            logger.warning("Input must be a finite number.")
            continue

        if min_value is not None and value <= min_value:
            logger.warning(f"Input must be > {min_value}.")
            continue

        return float(value)


def format_preview_table(df, float_cols=None, digits=3):
    """Render a compact text table for logger output."""
    if df is None or df.empty:
        return "<empty>"

    show_df = df.copy()
    float_cols = float_cols or []
    for col in float_cols:
        if col in show_df.columns:
            show_df[col] = show_df[col].map(lambda x: f"{x:.{digits}f}")
    return show_df.to_string(index=False)


def summarize_global_tech_cv_pass_across_cutoffs(rep_dfs, condition_defs, cutoff_values):
    """Preview global technical CV-QC pass counts across candidate cutoffs."""
    protein_index = next(iter(rep_dfs.values())).index
    rows = []

    for cutoff in cutoff_values:
        global_mask = pd.Series(True, index=protein_index)
        for rep_name, rep_df in rep_dfs.items():
            for cond in condition_defs:
                key = cond["key"]
                cv_col = f"CV_Rm1-3__{key}"
                cond_mask = np.isfinite(rep_df[cv_col]) & (rep_df[cv_col] <= cutoff)
                global_mask &= cond_mask.reindex(protein_index, fill_value=False)

        pass_n = int(global_mask.sum())
        total_n = int(len(global_mask))
        rows.append(
            {
                "Tech_CV_cutoff": float(cutoff),
                "Global_pass_count": pass_n,
                "Total_count": total_n,
                "Global_pass_percent": float(100.0 * pass_n / max(total_n, 1)),
            }
        )

    return pd.DataFrame(rows)


def summarize_global_bio_cv_pass_across_cutoffs(biocv_by_key, condition_defs, cutoff_values):
    """Preview global biological CV-QC pass counts across candidate cutoffs."""
    if not biocv_by_key:
        return pd.DataFrame()

    protein_index = next(iter(biocv_by_key.values())).index
    rows = []

    for cutoff in cutoff_values:
        global_mask = pd.Series(True, index=protein_index)
        for cond in condition_defs:
            key = cond["key"]
            cond_cv = biocv_by_key[key]
            cond_mask = np.isfinite(cond_cv) & (cond_cv <= cutoff)
            global_mask &= cond_mask.reindex(protein_index, fill_value=False)

        pass_n = int(global_mask.sum())
        total_n = int(len(global_mask))
        rows.append(
            {
                "Bio_CV_cutoff": float(cutoff),
                "Global_pass_count": pass_n,
                "Total_count": total_n,
                "Global_pass_percent": float(100.0 * pass_n / max(total_n, 1)),
            }
        )

    return pd.DataFrame(rows)


def make_safe_key(text):
    """Convert display name to an ASCII-safe key for column naming."""
    key = re.sub(r"[^0-9A-Za-z]+", "_", str(text).strip())
    key = key.strip("_")
    return key if key else "Condition"


def detect_unique_col(df):
    """Detect #Unique-like column name. Prefer exact '#Unique'."""
    if "#Unique" in df.columns:
        return "#Unique"
    candidates = [
        c
        for c in df.columns
        if re.search(r"#\s*Unique|#Unique|Unique\s*Pept|UniquePept|Unique", c, flags=re.I)
    ]
    return candidates[0] if candidates else None


def detect_protein_col(df):
    preferred = ["Protein Group", "ProteinGroup", "Protein_Group", "Accession", "Accession "]
    return next((c for c in preferred if c in df.columns), df.columns[0])


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


def draw_quant_overlap(per_sample_sets, sample_names, out_prefix):
    """Draw overlap plot of quantified proteins for one biological replicate file."""
    n = len(sample_names)
    if n < 2:
        return None

    out_file = None
    if n < 4 and HAVE_VENN:
        plt.figure(figsize=(6, 6))
        if n == 2:
            venn2([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
        elif n == 3:
            venn3([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
        plt.title("Quantified proteins (unique>=2 & any channel>0) Venn")
        out_file = f"{out_prefix}_quantified_venn.tiff"
        plt.savefig(out_file, dpi=300, bbox_inches="tight", format="tiff")
        plt.close()
    elif HAVE_UPSET:
        upset_data = from_contents(per_sample_sets)
        plt.figure(figsize=(8, 6))
        UpSet(upset_data, show_counts=True).plot()
        out_file = f"{out_prefix}_quantified_upset.tiff"
        plt.savefig(out_file, dpi=300, bbox_inches="tight", format="tiff")
        plt.close()

    return out_file


def robust_bin_pvals(deltas, cv_sort_order, bin_size=300):
    """Compute robust p-values per bin sorted by representative CV."""
    n = len(deltas)
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
        bins.append(idxs[i:j])
        i = j

    pval_map = {}
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
            elif val >= p50:
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
            pval_map[idx] = pv

    aligned = pd.Series(1.0, index=deltas.index)
    if pval_map:
        pvals = pd.Series(pval_map)
        aligned.loc[pvals.index] = pvals.values
    return aligned


def assess_delta_normality(delta_series):
    """Assess normality for one DeltaRm series and return (is_normal, message)."""
    arr = pd.Series(delta_series).replace([np.inf, -np.inf], np.nan).dropna().values
    n = len(arr)

    if n < 3:
        return False, f"n={n}: too few proteins for normality test, treat as non-normal"

    if n <= 50:
        stat, p = stats.shapiro(arr)
        is_normal = p >= 0.05
        return is_normal, f"Shapiro-Wilk W={stat:.4f}, p={p:.4g}"

    if n <= 300:
        mean_val = np.mean(arr)
        std_val = np.std(arr, ddof=1)
        if not np.isfinite(std_val) or std_val <= 0:
            return False, "KS skipped because std<=0, treat as non-normal"
        stat, p = stats.kstest(arr, "norm", args=(mean_val, std_val))
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

    safe_name = make_safe_key(sample_name)
    out_files = []

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
    hist_kde_file = f"{out_prefix}_DeltaRm_hist_kde_{safe_name}.tiff"
    plt.savefig(hist_kde_file, dpi=300, format="tiff")
    plt.close()
    out_files.append(hist_kde_file)

    plt.figure(figsize=(6, 4))
    stats.probplot(arr, dist="norm", plot=plt)
    plt.title(f"DeltaRm Q-Q plot - {sample_name}")
    plt.tight_layout()
    qq_file = f"{out_prefix}_DeltaRm_QQ_{safe_name}.tiff"
    plt.savefig(qq_file, dpi=300, format="tiff")
    plt.close()
    out_files.append(qq_file)

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
        pp_file = f"{out_prefix}_DeltaRm_PP_{safe_name}.tiff"
        plt.savefig(pp_file, dpi=300, format="tiff")
        plt.close()
        out_files.append(pp_file)

    return out_files


def plot_delta_volcano(delta_series, fdr_series, sample_name, out_prefix, delta_cutoff=0.1, fdr_cutoff=0.05):
    """Plot volcano for robust-bin (global-normal) branch."""
    plot_df = pd.DataFrame({"DeltaRm": pd.Series(delta_series), "FDR": pd.Series(fdr_series)})
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna()
    if plot_df.empty:
        return None

    safe_name = make_safe_key(sample_name)
    plot_df["neglog10FDR"] = -np.log10(np.clip(plot_df["FDR"].values, 1e-300, None))
    up_mask = (plot_df["DeltaRm"] > delta_cutoff) & (plot_df["FDR"] < fdr_cutoff)
    down_mask = (plot_df["DeltaRm"] < -delta_cutoff) & (plot_df["FDR"] < fdr_cutoff)
    other_mask = ~(up_mask | down_mask)

    plt.figure(figsize=(6.2, 5))
    plt.scatter(
        plot_df.loc[other_mask, "DeltaRm"],
        plot_df.loc[other_mask, "neglog10FDR"],
        s=14,
        color="#9E9E9E",
        alpha=0.65,
    )
    plt.scatter(
        plot_df.loc[up_mask, "DeltaRm"],
        plot_df.loc[up_mask, "neglog10FDR"],
        s=18,
        color="#D7301F",
        alpha=0.85,
        label="Up & significant",
    )
    plt.scatter(
        plot_df.loc[down_mask, "DeltaRm"],
        plot_df.loc[down_mask, "neglog10FDR"],
        s=18,
        color="#2C7FB8",
        alpha=0.85,
        label="Down & significant",
    )

    plt.axvline(delta_cutoff, color="#D7301F", linestyle="--", linewidth=1)
    plt.axvline(-delta_cutoff, color="#2C7FB8", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(fdr_cutoff), color="black", linestyle="--", linewidth=1)
    plt.xlabel("DeltaRm")
    plt.ylabel("-log10(FDR)")
    plt.title(f"DeltaRm volcano plot - {sample_name}")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()

    out_file = f"{out_prefix}_DeltaRm_volcano_{safe_name}.tiff"
    plt.savefig(out_file, dpi=300, format="tiff")
    plt.close()
    return out_file


def process_replicate(file_path, rep_name, condition_defs=None, tech_cv_cutoff=0.3):
    """
    Process one biological replicate CSV into a per-protein metric table.
    Returns dict with processed DataFrame indexed by ProteinID and metadata.
    """
    df = pd.read_csv(file_path)
    protein_col = detect_protein_col(df)
    unique_col = detect_unique_col(df)
    if unique_col is None:
        raise ValueError(f"[{rep_name}] Could not detect a '#Unique' column in file: {file_path}")

    sample_order, sample_cols_map = detect_sample_blocks(df)
    if not sample_order:
        tmt_cols = [c for c in df.columns if re.search(r"(126|127|128|129|130|131)", c)]
        if not tmt_cols:
            raise ValueError(f"[{rep_name}] No TMT reporter-like columns (126..131) detected.")

        if condition_defs is None:
            num_samples = int(input(f"[{rep_name}] Enter number of conditions in this file: ").strip())
        else:
            num_samples = len(condition_defs)

        if num_samples <= 0 or (len(tmt_cols) % num_samples != 0):
            raise ValueError(
                f"[{rep_name}] Reporter columns cannot be evenly split by conditions: "
                f"reporter_cols={len(tmt_cols)}, conditions={num_samples}"
            )

        channels_per_sample = len(tmt_cols) // num_samples
        sample_order = list(range(1, num_samples + 1))
        sample_cols_map = {
            s: tmt_cols[i * channels_per_sample : (i + 1) * channels_per_sample]
            for i, s in enumerate(sample_order)
        }

    if condition_defs is None:
        condition_defs = []
        used_keys = set()
        for s in sample_order:
            display_name = input(f"[{rep_name}] Please enter name for Sample {s} (e.g. Vehicle/Treatment): ").strip()
            if not display_name:
                display_name = f"Condition_{s}"
            base_key = make_safe_key(display_name)
            key = base_key
            suffix = 2
            while key in used_keys:
                key = f"{base_key}_{suffix}"
                suffix += 1
            used_keys.add(key)
            condition_defs.append({"name": display_name, "key": key})

        sample_to_cond = {s: condition_defs[i] for i, s in enumerate(sample_order)}
        local_names_by_sample = {s: sample_to_cond[s]["name"] for s in sample_order}
    else:
        if len(condition_defs) != len(sample_order):
            raise ValueError(
                f"[{rep_name}] Condition count mismatch with BR1 mapping: "
                f"expected={len(condition_defs)}, detected={len(sample_order)}"
            )

        canonical_lookup = {}
        for cond in condition_defs:
            canonical_lookup[cond["name"].strip().lower()] = cond
            canonical_lookup[cond["key"].strip().lower()] = cond

        logger.info(
            f"[{rep_name}] Please name each condition in this biological replicate; "
            "names will be mapped to BR1 canonical conditions."
        )

        local_names_by_sample = {}
        for i, s in enumerate(sample_order):
            default_name = condition_defs[i]["name"]
            entered = input(
                f"[{rep_name}] Enter condition name for Sample {s} in this file [default: {default_name}]: "
            ).strip()
            local_names_by_sample[s] = entered if entered else default_name

        sample_to_cond = {}
        used_keys = set()
        unresolved_samples = []

        for s in sample_order:
            token = local_names_by_sample[s].strip().lower()
            cond = canonical_lookup.get(token)
            if cond is not None and cond["key"] not in used_keys:
                sample_to_cond[s] = cond
                used_keys.add(cond["key"])
            else:
                unresolved_samples.append(s)

        if unresolved_samples:
            logger.info(f"[{rep_name}] Auto-mapping incomplete, manual mapping is required for remaining samples.")

        for s in unresolved_samples:
            while True:
                remaining = [c for c in condition_defs if c["key"] not in used_keys]
                if remaining:
                    options = ", ".join([f"{c['name']}[{c['key']}]" for c in remaining])
                else:
                    options = ", ".join([f"{c['name']}[{c['key']}]" for c in condition_defs])

                choice = input(
                    f"[{rep_name}] Map Sample {s} ('{local_names_by_sample[s]}') to canonical condition ({options}): "
                ).strip().lower()

                cond = canonical_lookup.get(choice)
                if cond is None:
                    logger.warning(f"[{rep_name}] Unknown condition '{choice}'. Please input canonical name or key.")
                    continue
                if cond["key"] in used_keys:
                    logger.warning(f"[{rep_name}] Canonical condition '{cond['name']}' already assigned in this BioRep.")
                    continue

                sample_to_cond[s] = cond
                used_keys.add(cond["key"])
                break

        if len(sample_to_cond) != len(sample_order):
            raise ValueError(f"[{rep_name}] Failed to map all samples to canonical conditions.")

    local_label_by_key = {sample_to_cond[s]["key"]: local_names_by_sample[s] for s in sample_order}
    mapping_str = "; ".join(
        [
            f"Sample{s}:{local_names_by_sample[s]}->{sample_to_cond[s]['name']}[{sample_to_cond[s]['key']}]"
            for s in sample_order
        ]
    )
    logger.info(f"[{rep_name}] Condition mapping: {mapping_str}")

    col_map = {}
    cond_block_cols = {}
    for s in sample_order:
        cond = sample_to_cond[s]
        cond_key = cond["key"]
        old_cols = sample_cols_map[s]
        cond_block_cols[cond_key] = []
        for old in old_cols:
            repfrag = extract_reporter_fragment(old)
            newname = f"{cond_key} {repfrag}"
            col_map[old] = newname
            cond_block_cols[cond_key].append(newname)

    df = df.rename(columns=col_map)
    all_reporter_cols = [c for cols in cond_block_cols.values() for c in cols]

    unique_ge2 = df[unique_col].fillna(0) >= 2
    quantified_any = (df[all_reporter_cols].fillna(0) != 0).any(axis=1)

    per_sample_sets = {}
    overlap_labels = [local_label_by_key.get(c["key"], c["name"]) for c in condition_defs]
    for cond in condition_defs:
        key = cond["key"]
        label = local_label_by_key.get(key, cond["name"])
        cols = cond_block_cols[key]
        has_signal = (df[cols].fillna(0) != 0).any(axis=1)
        per_sample_sets[label] = set(df.index[has_signal & unique_ge2])

    overlap_plot = draw_quant_overlap(
        per_sample_sets,
        overlap_labels,
        f"{os.path.splitext(file_path)[0]}_{rep_name}",
    )
    if overlap_plot is not None:
        logger.info(f"[{rep_name}] Quantified overlap plot saved: {overlap_plot}")

    final_mask = pd.Series(True, index=df.index)
    for cond in condition_defs:
        key = cond["key"]
        cols = cond_block_cols[key]
        final_mask &= (df[cols].fillna(0) != 0).all(axis=1)
    final_mask &= unique_ge2

    df_final = df.loc[final_mask].copy()
    df_final["ProteinID"] = df_final[protein_col].astype(str).str.strip()
    df_final = df_final[df_final["ProteinID"] != ""].copy()

    metrics = pd.DataFrame({"ProteinID": df_final["ProteinID"]})

    for cond in condition_defs:
        key = cond["key"]
        cols = cond_block_cols[key]
        reporter_map = {}
        for col in cols:
            m = re.search(r"(126|127|128|129|130|131)", col)
            if m:
                reporter_map[int(m.group(1))] = col

        required = [126, 127, 128, 129, 130, 131]
        missing = [r for r in required if r not in reporter_map]
        if missing:
            raise ValueError(f"[{rep_name}] Condition '{cond['name']}' missing reporter channels: {missing}")

        rm1 = df_final[reporter_map[129]] / df_final[reporter_map[126]]
        rm2 = df_final[reporter_map[130]] / df_final[reporter_map[127]]
        rm3 = df_final[reporter_map[131]] / df_final[reporter_map[128]]

        c_rm1 = f"Rm1__{key}"
        c_rm2 = f"Rm2__{key}"
        c_rm3 = f"Rm3__{key}"
        c_geo = f"geometric_mean_Rm__{key}"
        c_cv = f"CV_Rm1-3__{key}"

        metrics[c_rm1] = rm1.replace([np.inf, -np.inf], np.nan)
        metrics[c_rm2] = rm2.replace([np.inf, -np.inf], np.nan)
        metrics[c_rm3] = rm3.replace([np.inf, -np.inf], np.nan)

        rm_mat = metrics[[c_rm1, c_rm2, c_rm3]]
        with np.errstate(divide="ignore", invalid="ignore"):
            log_mean = np.log(rm_mat).mean(axis=1, skipna=True)
        metrics[c_geo] = np.exp(log_mean)
        metrics[c_cv] = rm_mat.std(axis=1, ddof=1) / rm_mat.mean(axis=1)

    for cond in condition_defs:
        key = cond["key"]
        c_cv = f"CV_Rm1-3__{key}"
        cv_series = (
            metrics[c_cv]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .sort_values(ascending=False)
            .reset_index(drop=True)
        )
        if cv_series.empty:
            continue

        rank = cv_series.index + 1
        plt.figure(figsize=(6, 4))
        plt.scatter(rank, cv_series.values, color="#5E556A", alpha=0.7, s=20)
        prop_20 = (cv_series <= 0.2).sum() / len(cv_series) * 100
        prop_cut = (cv_series <= tech_cv_cutoff).sum() / len(cv_series) * 100
        ymax = max(float(cv_series.max()) * 1.05, max(0.25, tech_cv_cutoff + 0.08))
        plt.ylim(0, ymax)
        plt.text(0.62 * len(cv_series), min(0.205, ymax * 0.95), f"<=20%: {prop_20:.1f}%", fontsize=10)
        plt.text(
            0.62 * len(cv_series),
            min(tech_cv_cutoff + 0.005, ymax * 0.95),
            f"<={tech_cv_cutoff*100:.0f}%: {prop_cut:.1f}%",
            fontsize=10,
        )
        plt.xlabel("Proteins sorted by CV (high->low)")
        plt.ylabel("CV")
        plt.title(f"CV distribution - {rep_name} - {cond['name']} (strict-filtered proteins)")
        plt.tight_layout()
        cv_plot = f"{os.path.splitext(file_path)[0]}_{rep_name}_CV_{key}.tiff"
        plt.savefig(cv_plot, dpi=300, format="tiff")
        plt.close()

    cv_cols = [f"CV_Rm1-3__{cond['key']}" for cond in condition_defs]
    cv_pass_mask = (metrics[cv_cols] <= tech_cv_cutoff).all(axis=1)
    cv_pass_count = int(cv_pass_mask.sum())

    metrics_strict = metrics.copy()
    if metrics_strict["ProteinID"].duplicated().any():
        numeric_cols = [c for c in metrics_strict.columns if c != "ProteinID"]
        metrics_strict = metrics_strict.groupby("ProteinID", as_index=False)[numeric_cols].median()

    metrics_strict = metrics_strict.set_index("ProteinID").sort_index()

    summary = {
        "BioRep": rep_name,
        "File": file_path,
        "Condition_mapping": mapping_str,
        "Protein_column": protein_col,
        "Unique_column": unique_col,
        "Quantified_any_count": int(quantified_any.sum()),
        "Unique_ge2_count": int(unique_ge2.sum()),
        "Strict_filtered_count": int(df_final.shape[0]),
        "TechRep_CV_cutoff": float(tech_cv_cutoff),
        "TechRep_CV_le_cutoff_count": cv_pass_count,
    }

    logger.info(
        f"[{rep_name}] strict-filtered proteins: {summary['Strict_filtered_count']}; "
        f"CV<={tech_cv_cutoff:.3g} proteins (info only): {summary['TechRep_CV_le_cutoff_count']}; "
        "no CV filtering applied before cross-biorep intersection/normalization"
    )

    return {
        "metrics_df": metrics_strict,
        "condition_defs": condition_defs,
        "summary": summary,
    }


def pick_vehicle_key(condition_defs):
    for cond in condition_defs:
        if cond["name"].strip().lower() == "vehicle" or cond["key"].strip().lower() == "vehicle":
            return cond["key"]

    user_text = input("Enter Vehicle condition name (display name or key): ").strip().lower()
    for cond in condition_defs:
        if cond["name"].strip().lower() == user_text or cond["key"].strip().lower() == user_text:
            return cond["key"]
    raise ValueError("Vehicle condition was not found in condition mapping.")


def apply_median_of_medians_normalization(rep_dfs, condition_defs, vehicle_key):
    """
    Batch correction with median-of-medians strategy:
    - Compute vehicle median per BioRep
    - Use median(BioRep vehicle medians) as global reference
    - Normalize each BioRep condition by global_ref / median(condition in that BioRep)

    Note:
    For positive values, median(log(x)) and log(median(x)) are equivalent. We use
    explicit log-domain computation for numerical stability and transparency.
    """
    vehicle_medians = {}
    for rep_name, rep_df in rep_dfs.items():
        vcol = f"geometric_mean_Rm__{vehicle_key}"
        vehicle_medians[rep_name] = float(rep_df[vcol].median())

    vehicle_vals = np.array(list(vehicle_medians.values()), dtype=float)
    if np.any(~np.isfinite(vehicle_vals)) or np.any(vehicle_vals <= 0):
        raise ValueError("Vehicle medians must be finite positive values for log-domain normalization.")
    global_vehicle_median = float(np.exp(np.median(np.log(vehicle_vals))))

    norm_factors = {}
    for rep_name, rep_df in rep_dfs.items():
        norm_factors[rep_name] = {}
        for cond in condition_defs:
            key = cond["key"]
            gcol = f"geometric_mean_Rm__{key}"
            sample_median = float(rep_df[gcol].median())
            if (not np.isfinite(sample_median)) or sample_median == 0:
                cf = np.nan
            else:
                cf = float(np.exp(np.log(global_vehicle_median) - np.log(sample_median)))

            rep_df[f"adjusted_Rm__{key}"] = rep_df[gcol] * cf
            rep_df[f"norm_factor__{key}"] = cf
            norm_factors[rep_name][key] = cf

        for cond in condition_defs:
            key = cond["key"]
            if key == vehicle_key:
                continue
            rep_df[f"DeltaRm__{key}"] = rep_df[f"adjusted_Rm__{key}"] - rep_df[f"adjusted_Rm__{vehicle_key}"]

    return global_vehicle_median, vehicle_medians, norm_factors


def build_cv_qc_masks(rep_dfs, condition_defs, cv_cutoff=0.3):
    """
    Build technical-replicate CV-QC pass masks used before differential analysis.
    Rule per BioRep x condition:
    - CV(condition) <= cutoff

    Global technical CV-QC is later defined as passing this rule for all
    BioReps and all conditions.
    """
    cv_qc_masks = {}
    rows = []

    for rep_name, rep_df in rep_dfs.items():
        for cond in condition_defs:
            key = cond["key"]
            cv_col = f"CV_Rm1-3__{key}"

            finite_mask = np.isfinite(rep_df[cv_col])
            mask = finite_mask & (rep_df[cv_col] <= cv_cutoff)
            mask = mask.fillna(False)

            rep_df[f"CV_QC_pass__{key}"] = mask
            cv_qc_masks[(rep_name, key)] = mask

            total_n = int(mask.shape[0])
            pass_n = int(mask.sum())
            rows.append(
                {
                    "BioRep": rep_name,
                    "Condition_name": cond["name"],
                    "Condition_key": key,
                    "CV_cutoff": cv_cutoff,
                    "Pass_count": pass_n,
                    "Total_count": total_n,
                    "Pass_percent": float(100.0 * pass_n / max(total_n, 1)),
                }
            )

    return cv_qc_masks, pd.DataFrame(rows)


def compute_biorep_pooled_cv_by_condition(rep_dfs, condition_defs):
    """Compute pooled technical raw-Rm CV across biological replicates for each condition."""
    rep_names = list(rep_dfs.keys())
    protein_index = next(iter(rep_dfs.values())).index
    condition_keys = [c["key"] for c in condition_defs]
    required_n = 3 * len(rep_names)

    biocv_by_key = {}
    for key in condition_keys:
        pooled_cols = []
        for rep_name in rep_names:
            rep_df = rep_dfs[rep_name]
            pooled_cols.extend(
                [
                    rep_df[f"Rm1__{key}"].reindex(protein_index).values,
                    rep_df[f"Rm2__{key}"].reindex(protein_index).values,
                    rep_df[f"Rm3__{key}"].reindex(protein_index).values,
                ]
            )

        rm_mat = np.column_stack(pooled_cols)
        valid_n = np.isfinite(rm_mat).sum(axis=1)

        with np.errstate(divide="ignore", invalid="ignore"):
            cv_vals = np.nanstd(rm_mat, axis=1, ddof=1) / np.nanmean(rm_mat, axis=1)

        cv_vals = np.where(valid_n >= required_n, cv_vals, np.nan)
        cv_vals = np.where(np.isfinite(cv_vals), cv_vals, np.nan)
        biocv_by_key[key] = pd.Series(cv_vals, index=protein_index, name=f"BioRep_CV_pooledTech_rawRm__{key}")

    return biocv_by_key


def build_biorep_cv_qc_masks(rep_dfs, condition_defs, out_prefix, cv_cutoff=0.4, biocv_by_key=None):
    """
    Build biological-replicate CV-QC masks using pooled technical-replicate raw Rm across BioReps.

    For each condition, pooled BioRep CV(raw Rm, condition) must be <= cutoff.

    pooled BioRep CV for one protein/condition is calculated from all technical
    raw Rm values across all BioReps, i.e. [Rm1, Rm2, Rm3] x N_BioRep.
    Global biological CV-QC is later defined as passing this rule for all
    conditions.
    """
    key_to_name = {c["key"]: c["name"] for c in condition_defs}
    condition_keys = [c["key"] for c in condition_defs]

    if biocv_by_key is None:
        biocv_by_key = compute_biorep_pooled_cv_by_condition(rep_dfs, condition_defs)

    plot_files = []

    for key in condition_keys:
        cv_series = biocv_by_key[key]

        plot_series = cv_series.dropna().sort_values(ascending=False).reset_index(drop=True)
        if not plot_series.empty:
            rank = plot_series.index + 1
            prop_30 = (plot_series <= 0.3).sum() / len(plot_series) * 100
            prop_40 = (plot_series <= 0.4).sum() / len(plot_series) * 100
            prop_cut = (plot_series <= cv_cutoff).sum() / len(plot_series) * 100
            ymax = max(float(plot_series.max()) * 1.05, 0.5, cv_cutoff + 0.08)

            plt.figure(figsize=(6, 4))
            plt.scatter(rank, plot_series.values, color="#4C6A92", alpha=0.72, s=20)
            if (not np.isclose(cv_cutoff, 0.3)) and (not np.isclose(cv_cutoff, 0.4)):
                plt.axhline(cv_cutoff, linestyle="--", linewidth=1, color="#111827")
            plt.ylim(0, ymax)
            plt.text(0.62 * len(plot_series), min(0.305, ymax * 0.95), f"<=30%: {prop_30:.1f}%", fontsize=10)
            plt.text(0.62 * len(plot_series), min(0.405, ymax * 0.95), f"<=40%: {prop_40:.1f}%", fontsize=10)
            if (not np.isclose(cv_cutoff, 0.3)) and (not np.isclose(cv_cutoff, 0.4)):
                plt.text(
                    0.62 * len(plot_series),
                    min(cv_cutoff + 0.005, ymax * 0.95),
                    f"<={cv_cutoff*100:.0f}%: {prop_cut:.1f}%",
                    fontsize=10,
                )
            plt.xlabel("Proteins sorted by BioRep CV (high->low)")
            plt.ylabel("BioRep CV of pooled technical raw Rm")
            plt.title(f"BioRep pooled-tech CV - {key_to_name.get(key, key)}")
            plt.tight_layout()

            out_plot = f"{out_prefix}_BioRep_CV_{key}.tiff"
            plt.savefig(out_plot, dpi=300, format="tiff")
            plt.close()
            plot_files.append(out_plot)

    bio_qc_masks = {}
    rows = []

    for cond in condition_defs:
        key = cond["key"]
        cond_cv = biocv_by_key[key]
        finite_mask = np.isfinite(cond_cv)
        mask = finite_mask & (cond_cv <= cv_cutoff)
        mask = mask.fillna(False)

        bio_qc_masks[key] = mask
        rows.append(
            {
                "Condition_name": cond["name"],
                "Condition_key": key,
                "CV_cutoff": float(cv_cutoff),
                "Pass_count": int(mask.sum()),
                "Total_count": int(mask.shape[0]),
                "Pass_percent": float(100.0 * mask.sum() / max(mask.shape[0], 1)),
                "Median_BioRep_pooledTech_CV": float(np.nanmedian(cond_cv.values)),
            }
        )

    for rep_name, rep_df in rep_dfs.items():
        for key in condition_keys:
            rep_df[f"BioRep_CV_pooledTech_rawRm__{key}"] = biocv_by_key[key].reindex(rep_df.index)
        for cond in condition_defs:
            key = cond["key"]
            rep_df[f"BioRep_CV_QC_pass__{key}"] = bio_qc_masks[key].reindex(rep_df.index, fill_value=False)

    return bio_qc_masks, pd.DataFrame(rows), plot_files


def build_global_cv_qc_masks(rep_dfs, condition_defs, tech_cv_qc_masks, bio_cv_qc_masks):
    """
    Build global two-layer CV-QC masks.

    Global rule: protein must pass CV-QC in all BioReps and all
    conditions to enter downstream multi-condition differential analysis.
    """
    rep_names = list(rep_dfs.keys())
    protein_index = next(iter(rep_dfs.values())).index

    tech_global = pd.Series(True, index=protein_index)
    for rep_name in rep_names:
        for cond in condition_defs:
            key = cond["key"]
            tech_global &= tech_cv_qc_masks[(rep_name, key)].reindex(protein_index, fill_value=False)

    bio_global = pd.Series(True, index=protein_index)
    for cond in condition_defs:
        key = cond["key"]
        bio_global &= bio_cv_qc_masks[key].reindex(protein_index, fill_value=False)

    two_layer_global = tech_global & bio_global
    total_n = int(len(protein_index))

    summary_df = pd.DataFrame(
        [
            {
                "Mask": "Global_TechRep_CV_QC_all_BioReps_all_conditions",
                "Pass_count": int(tech_global.sum()),
                "Total_count": total_n,
                "Pass_percent": float(100.0 * tech_global.sum() / max(total_n, 1)),
            },
            {
                "Mask": "Global_BioRep_CV_QC_all_conditions",
                "Pass_count": int(bio_global.sum()),
                "Total_count": total_n,
                "Pass_percent": float(100.0 * bio_global.sum() / max(total_n, 1)),
            },
            {
                "Mask": "Global_TwoLayer_CV_QC",
                "Pass_count": int(two_layer_global.sum()),
                "Total_count": total_n,
                "Pass_percent": float(100.0 * two_layer_global.sum() / max(total_n, 1)),
            },
        ]
    )

    for rep_df in rep_dfs.values():
        rep_df["Global_TechRep_CV_QC_pass"] = tech_global.reindex(rep_df.index, fill_value=False)
        rep_df["Global_BioRep_CV_QC_pass"] = bio_global.reindex(rep_df.index, fill_value=False)
        rep_df["Global_TwoLayer_CV_QC_pass"] = two_layer_global.reindex(rep_df.index, fill_value=False)

    return tech_global, bio_global, two_layer_global, summary_df


def evaluate_global_normality(
    rep_dfs,
    exp_conditions,
    out_prefix,
    tech_cv_qc_masks=None,
    bio_cv_qc_masks=None,
    global_cv_qc_mask=None,
):
    """Global decision: all BioRep x Condition DeltaRm distributions must be normal."""
    all_normal = True
    normality_results = {}

    for rep_name, rep_df in rep_dfs.items():
        for cond in exp_conditions:
            key = cond["key"]
            display_name = cond["name"]
            dcol = f"DeltaRm__{key}"
            panel_name = f"{rep_name}_{display_name}"
            if global_cv_qc_mask is not None:
                cv_mask = global_cv_qc_mask.reindex(rep_df.index, fill_value=False)
            else:
                cv_mask = pd.Series(True, index=rep_df.index)
                if tech_cv_qc_masks is not None:
                    cv_mask &= tech_cv_qc_masks[(rep_name, key)].reindex(rep_df.index, fill_value=False)
                if bio_cv_qc_masks is not None:
                    cv_mask &= bio_cv_qc_masks[key].reindex(rep_df.index, fill_value=False)

            deltas_for_test = rep_df.loc[cv_mask, dcol]

            plot_files = plot_delta_normality_diagnostics(deltas_for_test, panel_name, out_prefix)
            if plot_files:
                logger.info(f"{panel_name}: normality diagnostics saved: {', '.join(plot_files)}")

            is_normal, msg = assess_delta_normality(deltas_for_test)
            normality_results[(rep_name, key)] = {"is_normal": is_normal, "message": msg}
            logger.info(f"{panel_name}: {msg}; normal={is_normal}; CV-QC pass n={int(cv_mask.sum())}")
            if not is_normal:
                all_normal = False

    return all_normal, normality_results


def run_global_normal_branch(
    rep_dfs,
    exp_conditions,
    vehicle_key,
    out_prefix,
    tech_cv_qc_masks=None,
    bio_cv_qc_masks=None,
    global_cv_qc_mask=None,
    global_tech_cv_qc_mask=None,
    global_bio_cv_qc_mask=None,
    tech_cv_cutoff=0.3,
    bio_cv_cutoff=0.4,
):
    """
    For global-normal case:
    1) run robust-bin p/FDR within each BioRep
    2) aggregate significance across BioRep using:
       all FDR<0.1 and >=1 BioRep FDR<0.05, plus abs(DeltaRm)>0.1 in all BioRep
    """
    rep_names = list(rep_dfs.keys())
    protein_index = next(iter(rep_dfs.values())).index

    for rep_name, rep_df in rep_dfs.items():
        for cond in exp_conditions:
            key = cond["key"]
            display_name = cond["name"]
            dcol = f"DeltaRm__{key}"
            pcol = f"pval__{key}"
            fcol = f"FDR__{key}"
            if global_cv_qc_mask is not None:
                cv_ok = global_cv_qc_mask.reindex(rep_df.index, fill_value=False)
            else:
                cv_ok = pd.Series(True, index=rep_df.index)
                if tech_cv_qc_masks is not None:
                    cv_ok &= tech_cv_qc_masks[(rep_name, key)].reindex(rep_df.index, fill_value=False)
                if bio_cv_qc_masks is not None:
                    cv_ok &= bio_cv_qc_masks[key].reindex(rep_df.index, fill_value=False)

            cv_pair = [f"CV_Rm1-3__{key}", f"CV_Rm1-3__{vehicle_key}"]
            pvals_series = pd.Series(1.0, index=rep_df.index)
            fdr_series = pd.Series(1.0, index=rep_df.index)

            if int(cv_ok.sum()) >= 3:
                valid_idx = rep_df.index[cv_ok]
                cv_rep = rep_df.loc[valid_idx, cv_pair].max(axis=1).fillna(np.inf)
                order = list(cv_rep.sort_values(ascending=False).index)

                bin_pvals = robust_bin_pvals(rep_df.loc[valid_idx, dcol], order, bin_size=300)
                pvals_valid = np.where(np.isnan(bin_pvals.values), 1.0, bin_pvals.values)
                _, fdrs_valid, _, _ = multipletests(pvals_valid, method="fdr_bh")
                pvals_series.loc[valid_idx] = pvals_valid
                fdr_series.loc[valid_idx] = fdrs_valid

            rep_df[pcol] = pvals_series.values
            rep_df[fcol] = fdr_series.values

            volcano = plot_delta_volcano(
                rep_df.loc[cv_ok, dcol],
                rep_df.loc[cv_ok, fcol],
                f"{rep_name}_{display_name}",
                out_prefix,
            )
            if volcano is not None:
                logger.info(f"{rep_name}_{display_name}: volcano saved: {volcano}")

    aggregate = {}
    for cond in exp_conditions:
        key = cond["key"]
        delta_mat = np.column_stack([rep_dfs[r][f"DeltaRm__{key}"].reindex(protein_index).values for r in rep_names])
        fdr_mat = np.column_stack([rep_dfs[r][f"FDR__{key}"].reindex(protein_index).values for r in rep_names])
        fdr_mat = np.where(np.isfinite(fdr_mat), fdr_mat, 1.0)
        if global_cv_qc_mask is not None:
            pass_cv_all = global_cv_qc_mask.reindex(protein_index, fill_value=False).values
            pass_tech_cv_all = (
                global_tech_cv_qc_mask.reindex(protein_index, fill_value=False).values
                if global_tech_cv_qc_mask is not None
                else pass_cv_all
            )
            pass_bio_cv = (
                global_bio_cv_qc_mask.reindex(protein_index, fill_value=False).values
                if global_bio_cv_qc_mask is not None
                else pass_cv_all
            )
        else:
            tech_cv_pass_mat = np.column_stack(
                [
                    tech_cv_qc_masks[(r, key)].reindex(protein_index, fill_value=False).values
                    if tech_cv_qc_masks is not None
                    else np.ones(len(protein_index), dtype=bool)
                    for r in rep_names
                ]
            )
            pass_tech_cv_all = np.all(tech_cv_pass_mat, axis=1)
            pass_bio_cv = (
                bio_cv_qc_masks[key].reindex(protein_index, fill_value=False).values
                if bio_cv_qc_masks is not None
                else np.ones(len(protein_index), dtype=bool)
            )
            pass_cv_all = pass_tech_cv_all & pass_bio_cv

        pass_delta = np.all(np.isfinite(delta_mat) & (np.abs(delta_mat) > 0.1), axis=1)
        pass_same_sign = np.all(delta_mat > 0, axis=1) | np.all(delta_mat < 0, axis=1)
        pass_fdr_all_01 = np.all(fdr_mat < 0.1, axis=1)
        pass_fdr_one_005 = np.any(fdr_mat < 0.05, axis=1)
        sig = pass_cv_all & pass_delta & pass_same_sign & pass_fdr_all_01 & pass_fdr_one_005

        agg_df = pd.DataFrame(index=protein_index)
        agg_df["DeltaRm_median"] = np.nanmedian(delta_mat, axis=1)
        agg_df["Pass_TechRep_CV_QC_all_reps"] = pass_tech_cv_all
        agg_df["Pass_BioRep_CV_QC"] = pass_bio_cv
        agg_df["Pass_CV_QC_all_reps"] = pass_cv_all
        agg_df["Pass_absDelta_gt_0.1_all_reps"] = pass_delta
        agg_df["Pass_same_sign_all_reps"] = pass_same_sign
        agg_df["Pass_FDR_lt_0.1_all_reps"] = pass_fdr_all_01
        agg_df["Pass_at_least_one_FDR_lt_0.05"] = pass_fdr_one_005
        agg_df["Sig"] = sig
        agg_df["Sig_basis"] = (
            f"global-normal: global TechRep-CV<={tech_cv_cutoff:.3g} across all BioReps x conditions + "
            f"global BioRep-CV<={bio_cv_cutoff:.3g} across all conditions + "
            "abs(DeltaRm)>0.1 all reps + same-sign DeltaRm all reps + all FDR<0.1 + >=1 rep FDR<0.05"
        )
        for i, rep_name in enumerate(rep_names):
            agg_df[f"DeltaRm_{rep_name}"] = delta_mat[:, i]
            agg_df[f"FDR_{rep_name}"] = fdr_mat[:, i]

        aggregate[key] = agg_df

    return aggregate


def run_global_non_normal_branch(
    rep_dfs,
    exp_conditions,
    vehicle_key,
    paired_method="ttest",
    tech_cv_qc_masks=None,
    bio_cv_qc_masks=None,
    global_cv_qc_mask=None,
    global_tech_cv_qc_mask=None,
    global_bio_cv_qc_mask=None,
    tech_cv_cutoff=0.3,
    bio_cv_cutoff=0.4,
):
    """
    For global-non-normal case:
    run paired test across biological replicates on adjusted Rm.
    Each biological replicate contributes one representative adjusted Rm value
    (already summarized from technical replicates via geometric mean Rm).
    Default paired_method='ttest', optional 'wilcoxon'.
    """
    rep_names = list(rep_dfs.keys())
    protein_index = next(iter(rep_dfs.values())).index
    aggregate = {}

    for cond in exp_conditions:
        key = cond["key"]
        delta_mat = np.column_stack([rep_dfs[r][f"DeltaRm__{key}"].reindex(protein_index).values for r in rep_names])
        if global_cv_qc_mask is not None:
            pass_cv_all = global_cv_qc_mask.reindex(protein_index, fill_value=False).values
            pass_tech_cv_all = (
                global_tech_cv_qc_mask.reindex(protein_index, fill_value=False).values
                if global_tech_cv_qc_mask is not None
                else pass_cv_all
            )
            pass_bio_cv = (
                global_bio_cv_qc_mask.reindex(protein_index, fill_value=False).values
                if global_bio_cv_qc_mask is not None
                else pass_cv_all
            )
        else:
            tech_cv_pass_mat = np.column_stack(
                [
                    tech_cv_qc_masks[(r, key)].reindex(protein_index, fill_value=False).values
                    if tech_cv_qc_masks is not None
                    else np.ones(len(protein_index), dtype=bool)
                    for r in rep_names
                ]
            )
            pass_tech_cv_all = np.all(tech_cv_pass_mat, axis=1)
            pass_bio_cv = (
                bio_cv_qc_masks[key].reindex(protein_index, fill_value=False).values
                if bio_cv_qc_masks is not None
                else np.ones(len(protein_index), dtype=bool)
            )
            pass_cv_all = pass_tech_cv_all & pass_bio_cv

        pvals = np.ones(len(protein_index), dtype=float)
        for i, pid in enumerate(protein_index):
            if not pass_cv_all[i]:
                pvals[i] = 1.0
                continue

            x = np.array([rep_dfs[r].at[pid, f"adjusted_Rm__{key}"] for r in rep_names], dtype=float)
            y = np.array([rep_dfs[r].at[pid, f"adjusted_Rm__{vehicle_key}"] for r in rep_names], dtype=float)

            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() < 3:
                pvals[i] = 1.0
                continue

            xv = x[valid]
            yv = y[valid]

            if paired_method == "wilcoxon":
                diff = xv - yv
                if np.allclose(diff, 0):
                    pvals[i] = 1.0
                else:
                    try:
                        pvals[i] = stats.wilcoxon(xv, yv, alternative="two-sided").pvalue
                    except Exception:
                        pvals[i] = 1.0
            else:
                t_res = stats.ttest_rel(xv, yv, nan_policy="omit")
                pv = t_res.pvalue if hasattr(t_res, "pvalue") else t_res[1]
                pvals[i] = pv if np.isfinite(pv) else 1.0

        _, fdrs, _, _ = multipletests(pvals, method="fdr_bh")
        pass_delta = np.all(np.isfinite(delta_mat) & (np.abs(delta_mat) > 0.1), axis=1)
        pass_same_sign = np.all(delta_mat > 0, axis=1) | np.all(delta_mat < 0, axis=1)
        sig = pass_cv_all & pass_delta & pass_same_sign & (fdrs < 0.05)

        delta_mean = np.nanmean(delta_mat, axis=1)
        delta_mean = np.where(pass_cv_all, delta_mean, np.nan)

        agg_df = pd.DataFrame(index=protein_index)
        agg_df["DeltaRm_mean"] = delta_mean
        agg_df["DeltaRm_display"] = agg_df["DeltaRm_mean"]
        agg_df["DeltaRm_median"] = np.nanmedian(delta_mat, axis=1)
        agg_df["Pass_TechRep_CV_QC_all_reps"] = pass_tech_cv_all
        agg_df["Pass_BioRep_CV_QC"] = pass_bio_cv
        agg_df["Pass_CV_QC_all_reps"] = pass_cv_all
        agg_df["Pass_absDelta_gt_0.1_all_reps"] = pass_delta
        agg_df["Pass_same_sign_all_reps"] = pass_same_sign
        agg_df["Paired_test"] = paired_method
        agg_df["pval"] = pvals
        agg_df["FDR"] = fdrs
        agg_df["Sig"] = sig
        agg_df["Sig_basis"] = (
            f"global-non-normal: global TechRep-CV<={tech_cv_cutoff:.3g} across all BioReps x conditions + "
            f"global BioRep-CV<={bio_cv_cutoff:.3g} across all conditions + "
            f"paired-{paired_method} FDR<0.05 + abs(DeltaRm)>0.1 all reps + same-sign DeltaRm all reps; DeltaRm_display=mean"
        )
        for i, rep_name in enumerate(rep_names):
            agg_df[f"DeltaRm_{rep_name}"] = delta_mat[:, i]

        aggregate[key] = agg_df

    return aggregate


def build_full_results(rep_dfs, exp_conditions, aggregate_by_cond, normality_results, analysis_mode):
    protein_index = next(iter(rep_dfs.values())).index
    full_df = pd.DataFrame(index=protein_index)
    full_df["analysis_mode_global"] = analysis_mode

    for (rep_name, cond_key), info in normality_results.items():
        full_df[f"normality__{rep_name}__{cond_key}__is_normal"] = bool(info["is_normal"])
        full_df[f"normality__{rep_name}__{cond_key}__message"] = info["message"]

    for rep_name, rep_df in rep_dfs.items():
        for col in rep_df.columns:
            full_df[f"{rep_name}__{col}"] = rep_df[col].reindex(protein_index)

    for cond in exp_conditions:
        key = cond["key"]
        agg = aggregate_by_cond[key]
        for col in agg.columns:
            full_df[f"AGG__{key}__{col}"] = agg[col].reindex(protein_index)

    return full_df


# ------------------------------
# Main workflow
# ------------------------------
raw_paths = input(
    "Please provide paths of biological replicate CSV files (>=3), separated by ';' or ',': "
).strip()
file_paths = parse_csv_paths(raw_paths)

if len(file_paths) < 3:
    raise ValueError("At least 3 biological replicate CSV files are required.")

for fp in file_paths:
    if not os.path.exists(fp):
        raise FileNotFoundError(fp)

logger.info(
    "Step 1/3: build CV-QC preview first (default cutoffs tech=0.3, bio=0.4), "
    "then ask for final thresholds."
)

rep_dfs = {}
condition_defs = None
summary_rows = []

for i, fp in enumerate(file_paths, start=1):
    rep_name = f"BR{i}"
    rep_result = process_replicate(
        fp,
        rep_name,
        condition_defs=condition_defs,
        tech_cv_cutoff=DEFAULT_TECH_CV_QC_CUTOFF,
    )
    if condition_defs is None:
        condition_defs = rep_result["condition_defs"]
        mapped = ", ".join([f"{c['name']}[{c['key']}]" for c in condition_defs])
        logger.info(f"Condition mapping fixed from BR1: {mapped}")
    rep_dfs[rep_name] = rep_result["metrics_df"]
    summary_rows.append(rep_result["summary"])

vehicle_key = pick_vehicle_key(condition_defs)
exp_conditions = [c for c in condition_defs if c["key"] != vehicle_key]
if not exp_conditions:
    raise ValueError("No experimental condition found after removing Vehicle.")

protein_sets = [set(df.index) for df in rep_dfs.values()]
common_proteins = sorted(set.intersection(*protein_sets))
if len(common_proteins) == 0:
    raise ValueError("No common proteins found across biological replicates after strict quantification filtering.")

for rep_name in rep_dfs:
    rep_dfs[rep_name] = rep_dfs[rep_name].loc[common_proteins].copy()

logger.info(
    "Common proteins retained by conservative intersection before normalization: "
    f"{len(common_proteins)}"
)

global_vehicle_median, vehicle_medians, norm_factors = apply_median_of_medians_normalization(
    rep_dfs, condition_defs, vehicle_key
)
logger.info(f"Global vehicle median (median of BioRep vehicle medians): {global_vehicle_median:.6g}")

for rep_name in rep_dfs:
    fac_str = ", ".join([f"{c['key']}={norm_factors[rep_name][c['key']]:.6g}" for c in condition_defs])
    logger.info(f"{rep_name} normalization factors: {fac_str}")

output_prefix = os.path.splitext(file_paths[0])[0] + "_BioRep"

preview_tech_cv_qc_masks, preview_tech_cv_qc_summary_df = build_cv_qc_masks(
    rep_dfs,
    condition_defs,
    cv_cutoff=DEFAULT_TECH_CV_QC_CUTOFF,
)

preview_biocv_by_key = compute_biorep_pooled_cv_by_condition(rep_dfs, condition_defs)
preview_bio_cv_qc_masks, preview_bio_cv_qc_summary_df, preview_bio_cv_plot_files = build_biorep_cv_qc_masks(
    rep_dfs,
    condition_defs,
    out_prefix=output_prefix,
    cv_cutoff=DEFAULT_BIO_CV_QC_CUTOFF,
    biocv_by_key=preview_biocv_by_key,
)

(
    preview_global_tech_cv_mask,
    preview_global_bio_cv_mask,
    preview_global_two_layer_cv_mask,
    preview_global_cv_qc_summary_df,
) = build_global_cv_qc_masks(
    rep_dfs,
    condition_defs,
    preview_tech_cv_qc_masks,
    preview_bio_cv_qc_masks,
)

tech_cutoff_preview_df = summarize_global_tech_cv_pass_across_cutoffs(
    rep_dfs,
    condition_defs,
    TECH_CV_PREVIEW_CANDIDATES,
)
bio_cutoff_preview_df = summarize_global_bio_cv_pass_across_cutoffs(
    preview_biocv_by_key,
    condition_defs,
    BIO_CV_PREVIEW_CANDIDATES,
)

logger.info("CV-QC preview table (global TechRep CV across all BioReps x conditions):")
logger.info("\n" + format_preview_table(tech_cutoff_preview_df, float_cols=["Tech_CV_cutoff", "Global_pass_percent"]))

logger.info("CV-QC preview table (global BioRep pooled-tech CV across all conditions):")
logger.info("\n" + format_preview_table(bio_cutoff_preview_df, float_cols=["Bio_CV_cutoff", "Global_pass_percent"]))

for _, row in preview_global_cv_qc_summary_df.iterrows():
    logger.info(
        f"[Preview] {row['Mask']}: {int(row['Pass_count'])}/{int(row['Total_count'])} pass "
        f"({float(row['Pass_percent']):.2f}%)"
    )

if preview_bio_cv_plot_files:
    logger.info(f"BioRep CV-QC preview plots saved: {', '.join(preview_bio_cv_plot_files)}")

input(
    "CV-QC preview generated. Please review tables/plots, then press Enter to input final CV-QC cutoffs: "
)

TECH_CV_QC_CUTOFF = ask_float_with_default(
    f"Enter technical-replicate CV-QC cutoff for final analysis (default {DEFAULT_TECH_CV_QC_CUTOFF}): ",
    default=DEFAULT_TECH_CV_QC_CUTOFF,
    min_value=0.0,
)
BIO_CV_QC_CUTOFF = ask_float_with_default(
    f"Enter biological-replicate CV-QC cutoff for final analysis (default {DEFAULT_BIO_CV_QC_CUTOFF}): ",
    default=DEFAULT_BIO_CV_QC_CUTOFF,
    min_value=0.0,
)
logger.info(
    f"Step 2/3 final CV-QC thresholds -> technical: {TECH_CV_QC_CUTOFF:.3g}, biological: {BIO_CV_QC_CUTOFF:.3g}"
)

if np.isclose(TECH_CV_QC_CUTOFF, DEFAULT_TECH_CV_QC_CUTOFF) and np.isclose(
    BIO_CV_QC_CUTOFF, DEFAULT_BIO_CV_QC_CUTOFF
):
    tech_cv_qc_masks = preview_tech_cv_qc_masks
    tech_cv_qc_summary_df = preview_tech_cv_qc_summary_df
    bio_cv_qc_masks = preview_bio_cv_qc_masks
    bio_cv_qc_summary_df = preview_bio_cv_qc_summary_df
    bio_cv_plot_files = preview_bio_cv_plot_files
    global_tech_cv_mask = preview_global_tech_cv_mask
    global_bio_cv_mask = preview_global_bio_cv_mask
    global_two_layer_cv_mask = preview_global_two_layer_cv_mask
    global_cv_qc_summary_df = preview_global_cv_qc_summary_df
else:
    tech_cv_qc_masks, tech_cv_qc_summary_df = build_cv_qc_masks(
        rep_dfs,
        condition_defs,
        cv_cutoff=TECH_CV_QC_CUTOFF,
    )

    bio_cv_qc_masks, bio_cv_qc_summary_df, bio_cv_plot_files = build_biorep_cv_qc_masks(
        rep_dfs,
        condition_defs,
        out_prefix=output_prefix,
        cv_cutoff=BIO_CV_QC_CUTOFF,
        biocv_by_key=preview_biocv_by_key,
    )

    global_tech_cv_mask, global_bio_cv_mask, global_two_layer_cv_mask, global_cv_qc_summary_df = build_global_cv_qc_masks(
        rep_dfs,
        condition_defs,
        tech_cv_qc_masks,
        bio_cv_qc_masks,
    )

for _, row in tech_cv_qc_summary_df.iterrows():
    logger.info(
        f"TechRep CV-QC ({row['BioRep']} | {row['Condition_name']}): "
        f"{int(row['Pass_count'])}/{int(row['Total_count'])} pass (<= {float(row['CV_cutoff']):.2f})"
    )

for _, row in bio_cv_qc_summary_df.iterrows():
    logger.info(
        f"BioRep CV-QC ({row['Condition_name']}): "
        f"{int(row['Pass_count'])}/{int(row['Total_count'])} pass (<= {float(row['CV_cutoff']):.2f})"
    )

if bio_cv_plot_files:
    logger.info(f"BioRep CV-QC plots saved: {', '.join(bio_cv_plot_files)}")

for _, row in global_cv_qc_summary_df.iterrows():
    logger.info(
        f"{row['Mask']}: {int(row['Pass_count'])}/{int(row['Total_count'])} pass "
        f"({float(row['Pass_percent']):.2f}%)"
    )

cv_cols_all_conditions = [f"CV_Rm1-3__{c['key']}" for c in condition_defs]
final_techrep_counts_common = {
    rep_name: int((rep_df[cv_cols_all_conditions] <= TECH_CV_QC_CUTOFF).all(axis=1).sum())
    for rep_name, rep_df in rep_dfs.items()
}

global_normal, normality_results = evaluate_global_normality(
    rep_dfs,
    exp_conditions,
    output_prefix,
    tech_cv_qc_masks=tech_cv_qc_masks,
    bio_cv_qc_masks=bio_cv_qc_masks,
    global_cv_qc_mask=global_two_layer_cv_mask,
)

if global_normal:
    logger.info("Global normality decision: ALL normal -> run robust-bin per BioRep + cross-BioRep FDR integration")
    aggregate_by_cond = run_global_normal_branch(
        rep_dfs,
        exp_conditions,
        vehicle_key,
        output_prefix,
        tech_cv_qc_masks=tech_cv_qc_masks,
        bio_cv_qc_masks=bio_cv_qc_masks,
        global_cv_qc_mask=global_two_layer_cv_mask,
        global_tech_cv_qc_mask=global_tech_cv_mask,
        global_bio_cv_qc_mask=global_bio_cv_mask,
        tech_cv_cutoff=TECH_CV_QC_CUTOFF,
        bio_cv_cutoff=BIO_CV_QC_CUTOFF,
    )
    analysis_mode = "global_normal_robust_per_biorep"
else:
    method_in = input(
        "Global non-normal detected. Choose paired test [ttest/wilcoxon], Enter for default ttest: "
    ).strip().lower()
    paired_method = method_in if method_in in {"ttest", "wilcoxon"} else "ttest"
    logger.info(f"Global non-normal decision -> run paired {paired_method} across biological replicates")
    aggregate_by_cond = run_global_non_normal_branch(
        rep_dfs,
        exp_conditions,
        vehicle_key,
        paired_method=paired_method,
        tech_cv_qc_masks=tech_cv_qc_masks,
        bio_cv_qc_masks=bio_cv_qc_masks,
        global_cv_qc_mask=global_two_layer_cv_mask,
        global_tech_cv_qc_mask=global_tech_cv_mask,
        global_bio_cv_qc_mask=global_bio_cv_mask,
        tech_cv_cutoff=TECH_CV_QC_CUTOFF,
        bio_cv_cutoff=BIO_CV_QC_CUTOFF,
    )
    analysis_mode = f"global_non_normal_paired_{paired_method}"

full_df = build_full_results(rep_dfs, exp_conditions, aggregate_by_cond, normality_results, analysis_mode)

full_outfile = output_prefix + "_full_Rm_results.csv"
full_df.reset_index().rename(columns={"index": "ProteinID"}).to_csv(full_outfile, index=False)
logger.info(f"Full results saved: {full_outfile}")

sig_cols = [f"AGG__{c['key']}__Sig" for c in exp_conditions]
sig_mask = full_df[sig_cols].fillna(False).any(axis=1)
diff_df = full_df.loc[sig_mask].copy()
diff_outfile = output_prefix + "_Differential_proteins.csv"
diff_df.reset_index().rename(columns={"index": "ProteinID"}).to_csv(diff_outfile, index=False)
logger.info(f"Differential proteins saved: {diff_outfile}")

summary_df = pd.DataFrame(summary_rows)
summary_df["Common_intersection_count"] = len(common_proteins)
summary_df["TechRep_CV_cutoff"] = float(TECH_CV_QC_CUTOFF)
summary_df["TechRep_CV_le_cutoff_count"] = (
    summary_df["BioRep"].map(final_techrep_counts_common).fillna(0).astype(int)
)
summary_df["TechRep_CV_cutoff_used"] = float(TECH_CV_QC_CUTOFF)
summary_df["BioRep_CV_cutoff_used"] = float(BIO_CV_QC_CUTOFF)
summary_df["Global_two_layer_CV_QC_pass_count"] = int(global_two_layer_cv_mask.sum())
summary_outfile = output_prefix + "_BioRep_summary.csv"
summary_df.to_csv(summary_outfile, index=False)
logger.info(f"BioRep summary saved: {summary_outfile}")

norm_rows = []
for rep_name, rep_df in rep_dfs.items():
    for cond in condition_defs:
        key = cond["key"]
        gcol = f"geometric_mean_Rm__{key}"
        norm_rows.append(
            {
                "BioRep": rep_name,
                "Condition_name": cond["name"],
                "Condition_key": key,
                "Condition_median": float(rep_df[gcol].median()),
                "Vehicle_median_this_BioRep": float(vehicle_medians[rep_name]),
                "Global_vehicle_median_of_medians": float(global_vehicle_median),
                "Normalization_factor": float(norm_factors[rep_name][key]),
            }
        )

norm_df = pd.DataFrame(norm_rows)
norm_outfile = output_prefix + "_normalization_factors.csv"
norm_df.to_csv(norm_outfile, index=False)
logger.info(f"Normalization factors saved: {norm_outfile}")

cv_qc_outfile = output_prefix + "_cv_qc_summary.csv"
tech_cv_qc_summary_df.to_csv(cv_qc_outfile, index=False)
logger.info(f"TechRep CV-QC summary saved: {cv_qc_outfile}")

bio_cv_qc_outfile = output_prefix + "_biorep_cv_qc_summary.csv"
bio_cv_qc_summary_df.to_csv(bio_cv_qc_outfile, index=False)
logger.info(f"BioRep CV-QC summary saved: {bio_cv_qc_outfile}")

global_cv_qc_outfile = output_prefix + "_global_cv_qc_summary.csv"
global_cv_qc_summary_df.to_csv(global_cv_qc_outfile, index=False)
logger.info(f"Global CV-QC summary saved: {global_cv_qc_outfile}")

sig_summary_rows = []
for cond in exp_conditions:
    key = cond["key"]
    sig_count = int(aggregate_by_cond[key]["Sig"].sum())
    sig_summary_rows.append({"Condition_name": cond["name"], "Condition_key": key, "Significant_protein_count": sig_count})

sig_summary_df = pd.DataFrame(sig_summary_rows)
sig_summary_outfile = output_prefix + "_significance_summary.csv"
sig_summary_df.to_csv(sig_summary_outfile, index=False)
logger.info(f"Significance summary saved: {sig_summary_outfile}")

logger.info("Refined-TPP plus multi-biorep processing completed successfully.")
