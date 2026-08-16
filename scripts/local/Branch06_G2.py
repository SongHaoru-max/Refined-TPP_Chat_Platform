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
4. make later JSON / Agent integration straightforward;
5. standardize PEAKS Protein Group/shared-evidence semantics before statistical
   analysis.

Important definitions in this refined version
-----------------------------------------------
- Protein quantitative unit: one PEAKS Protein Group contributes one reporter
  profile, one Rm profile, and one weight to replicate-median normalization.
  Member accessions are retained for biological identity and FASTA/site mapping.
- Shared peptide/site evidence: identical PEAKS quantitative evidence expanded
  across member accessions is counted once quantitatively, while all valid
  accession/site candidates are retained as annotations. Evidence mapping across
  multiple Protein Groups is retained for QC but excluded from canonical DeltaRm
  because the protein-group reference is not unique.
- N-glycopeptide detection: an N residue is considered deamidated when one of
  its attached modification annotations contains a mass shift around +0.98 Da.
  Other modifications attached to the same N are allowed.
- Rm normalization: replicate-specific multiplicative median normalization.
  Each unique control Protein Group contributes exactly once to M_i. The same
  Protein-Group-derived correction factors are applied to Protein-Group and
  Site-Group Rm.
- Site Group: one site-level quantitative/statistical unit defined by one
  Protein Group plus one candidate-site assignment set. A Site Group may have
  one unique candidate site or multiple ambiguous-within-Protein-Group candidate
  sites; candidate annotations do not multiply quantitative weight.
- CV QC: CV is calculated from RAW Rm values, before normalization, separately
  for Protein Groups and Site Groups.
- DeltaRm: for Site Group s matched to Protein Group g(s), replicate i uses
  DeltaRm_{s,i} = Rm^SG_{s,i,norm} - Rm^PG_{g(s),i,norm}. The effect-size
  threshold is applied to mean_i(DeltaRm_{s,i}), not to every replicate
  separately.
- BH correction: user-selectable. If enabled, multiple-testing correction is
  applied across Site Groups, not candidate accession:site annotation rows.
- Paired tests compare normalized matched Protein-Group Rm and normalized
  Site-Group Rm.
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
try:
    from upsetplot import from_contents, UpSet
except ImportError:
    from_contents = None
    UpSet = None
import statsmodels.api as sm


plt.rcParams.update({'figure.dpi':300})

PIPELINE_NAME = "Branch06_G2"
PIPELINE_VERSION = "refined-3.0"
DELTA_RM_THRESHOLD_DEFAULT = 0.10
NORMALITY_ALPHA_DEFAULT = 0.05
REQUIRED_TMT_CHANNELS = (126, 127, 128, 129, 130, 131)
PROTEIN_GROUP_CANDIDATES = ("Protein Group", "ProteinGroup", "Protein_Group")
ACCESSION_CANDIDATES = ("Accession", "Protein Accession")
PEPTIDE_CANDIDATES = ("Peptide", "Peptide Sequence", "Sequence")

# Formal PG/SG notation recorded in run metadata and mirrored by exported column names.
# PG = Protein Group quantitative unit; SG = Site Group quantitative/statistical unit.
METHOD_FORMULAS = {
    "protein_group_rm_rep1": "Rm^PG_{g,1} = I^PG_{g,129} / I^PG_{g,127}",
    "protein_group_rm_rep2": "Rm^PG_{g,2} = I^PG_{g,130} / I^PG_{g,128}",
    "protein_group_rm_rep3": "Rm^PG_{g,3} = I^PG_{g,131} / I^PG_{g,126}",
    "normalization_median": "M_i = median_g(Rm^PG_{g,i}) over unique control Protein Groups",
    "normalization_reference": "T = exp(median_i(log(M_i)))",
    "normalization_factor": "CF_i = T / M_i",
    "protein_group_rm_normalized": "Rm^PG_{g,i,norm} = Rm^PG_{g,i} * CF_i",
    "site_group_intensity": "I^SG_{s,c} = sum_{e in E_s} I_{e,c}; each peptide_evidence_id contributes once",
    "site_group_rm_rep1": "Rm^SG_{s,1} = I^SG_{s,129} / I^SG_{s,127}",
    "site_group_rm_rep2": "Rm^SG_{s,2} = I^SG_{s,130} / I^SG_{s,128}",
    "site_group_rm_rep3": "Rm^SG_{s,3} = I^SG_{s,131} / I^SG_{s,126}",
    "site_group_rm_normalized": "Rm^SG_{s,i,norm} = Rm^SG_{s,i} * CF_i",
    "protein_group_cv": "CV^PG_g = SD_i(Rm^PG_{g,i}) / mean_i(Rm^PG_{g,i}) using raw Rm",
    "site_group_cv": "CV^SG_s = SD_i(Rm^SG_{s,i}) / mean_i(Rm^SG_{s,i}) using raw Rm",
    "delta_rm_replicate": "DeltaRm_{s,i} = Rm^SG_{s,i,norm} - Rm^PG_{g(s),i,norm}",
    "delta_rm_mean": "mean_DeltaRm_s = mean_i(DeltaRm_{s,i})",
    "multiple_testing_unit": "one hypothesis per SiteGroupID",
}


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
        M_i = median_g(Rm^PG_{g,i} across unique control Protein Groups)
        reference = exp(median(log(M_i)))
        CF_i = reference / M_i
        Rm^PG_{g,i,norm} = Rm^PG_{g,i} * CF_i

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
    "Rm1": "protein_group_rm_rep1_raw",
    "Rm2": "protein_group_rm_rep2_raw",
    "Rm3": "protein_group_rm_rep3_raw",
    "Rm1_normalized": "protein_group_rm_rep1_norm",
    "Rm2_normalized": "protein_group_rm_rep2_norm",
    "Rm3_normalized": "protein_group_rm_rep3_norm",
    "mean_Rm_pg_norm": "protein_group_rm_norm_mean",
    "CV_Rm_filtered": "protein_group_rm_raw_cv",
    "CV_Rm": "protein_group_rm_raw_cv",
}

SITE_EXPORT_RENAME = {
    "Rm_sg1": "site_group_rm_rep1_raw",
    "Rm_sg2": "site_group_rm_rep2_raw",
    "Rm_sg3": "site_group_rm_rep3_raw",
    "CV_Rm_sg_complete": "site_group_rm_raw_cv",
    "CV_Rm_sg_filtered": "site_group_rm_raw_cv",
    "Rm_sg1_norm": "site_group_rm_rep1_norm",
    "Rm_sg2_norm": "site_group_rm_rep2_norm",
    "Rm_sg3_norm": "site_group_rm_rep3_norm",
    "mean_Rm_sg_norm": "site_group_rm_norm_mean",
    "mean_Rm_pg_norm": "matched_protein_group_rm_norm_mean",
    "Rm_pg1_norm_ref": "matched_protein_group_rm_rep1_norm",
    "Rm_pg2_norm_ref": "matched_protein_group_rm_rep2_norm",
    "Rm_pg3_norm_ref": "matched_protein_group_rm_rep3_norm",
    "ΔRm_rep1": "delta_rm_rep1",
    "ΔRm_rep2": "delta_rm_rep2",
    "ΔRm_rep3": "delta_rm_rep3",
    "ΔRm": "delta_rm_mean",
    "Matched_PG_CV": "matched_protein_group_rm_raw_cv",
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

def canonical_group_id(value):
    """Normalize PEAKS Protein Group identifiers to stable strings."""
    if pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)) and np.isfinite(value) and float(value).is_integer():
        return str(int(value))
    text = str(value).strip()
    if re.fullmatch(r"\d+\.0", text):
        text = text[:-2]
    return text


def resolve_column(df, candidates, label, required=True):
    """Resolve one semantic column without mixing Protein Group and Accession roles."""
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    if required:
        raise KeyError(
            f"Required {label} column not found. Expected one of: {list(candidates)}"
        )
    return None


def reporter_channel_from_col(colname):
    """Return reporter channel number 126..131 when uniquely encoded in a column name."""
    matches = re.findall(r"(?<!\d)(126|127|128|129|130|131)(?!\d)", str(colname))
    unique = sorted(set(int(x) for x in matches))
    if len(unique) == 1:
        return unique[0]
    return None


def detect_sample_channel_map(df):
    """
    Detect Sample N reporter columns and map them explicitly by reporter channel.

    This deliberately does not infer channel identity from column position.
    """
    sample_pattern = re.compile(r"Sample\s*(\d+)", flags=re.I)
    sample_map = defaultdict(dict)
    for col in df.columns:
        sample_match = sample_pattern.search(str(col))
        channel = reporter_channel_from_col(col)
        if not sample_match or channel is None:
            continue
        sample_idx = int(sample_match.group(1))
        if channel in sample_map[sample_idx]:
            raise ValueError(
                f"Duplicate reporter channel {channel} detected for Sample {sample_idx}: "
                f"'{sample_map[sample_idx][channel]}' and '{col}'"
            )
        sample_map[sample_idx][channel] = col

    if not sample_map:
        return [], {}

    expected = set(REQUIRED_TMT_CHANNELS)
    for sample_idx, channel_map in sample_map.items():
        observed = set(channel_map)
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        if missing or extra:
            raise ValueError(
                f"Sample {sample_idx} reporter-channel mapping is incomplete. "
                f"Missing={missing}, unexpected={extra}. "
                "Expected exactly TMT 126, 127, 128, 129, 130, 131."
            )

    return sorted(sample_map), sample_map


def rename_sample_channels(df, sample_order, sample_channel_map, sample_names):
    """Rename sample reporter columns and return both ordered lists and channel maps."""
    if len(sample_order) != len(sample_names):
        raise ValueError("Sample-order and sample-name lengths do not match.")
    if len(set(sample_names)) != len(sample_names):
        raise ValueError("Sample names must be unique.")

    rename_map = {}
    sample_block_cols = {}
    sample_channel_cols = {}
    for sample_idx, sample_name in zip(sample_order, sample_names):
        sample_block_cols[sample_name] = []
        sample_channel_cols[sample_name] = {}
        for channel in REQUIRED_TMT_CHANNELS:
            old = sample_channel_map[sample_idx][channel]
            new = f"{sample_name} TMT-{channel}"
            rename_map[old] = new
            sample_block_cols[sample_name].append(new)
            sample_channel_cols[sample_name][channel] = new

    df = df.rename(columns=rename_map).copy()
    return df, sample_block_cols, sample_channel_cols, rename_map


def ensure_numeric_columns(df, columns, label):
    """Fail early when reporter-intensity columns contain non-numeric values."""
    out = df.copy()
    for col in columns:
        try:
            out[col] = pd.to_numeric(out[col], errors="raise")
        except Exception as exc:
            raise ValueError(f"Non-numeric values detected in {label} column '{col}': {exc}") from exc
    return out


def validate_group_quantitative_profiles(df, group_col, reporter_cols):
    """
    Verify that all accession rows within a PEAKS Protein Group share one reporter profile.

    The refined workflow treats Protein Group as the protein-level quantitative unit.
    """
    inconsistent = []
    for group_id, group in df.groupby(group_col, dropna=False, sort=False):
        profiles = group[reporter_cols].drop_duplicates()
        if len(profiles) > 1:
            inconsistent.append(str(group_id))
            if len(inconsistent) >= 10:
                break
    if inconsistent:
        raise ValueError(
            "Inconsistent reporter-ion profiles were found within PEAKS Protein Group(s): "
            + ", ".join(inconsistent)
            + ". Protein Group can only be used as one quantitative unit when member rows "
              "share the same reporter profile."
        )


def find_all_occurrences(sequence, peptide):
    """Return every start index where peptide occurs in sequence."""
    starts = []
    start = 0
    while True:
        idx = sequence.find(peptide, start)
        if idx < 0:
            break
        starts.append(idx)
        start = idx + 1
    return starts


def parse_site_label(label):
    """Parse a canonical site label such as 'P36507:N123'."""
    accession, site_text = str(label).rsplit(":N", 1)
    return accession, int(site_text)


def build_site_motif(prot_dict, site_label, half_len):
    accession, site_pos = parse_site_label(site_label)
    sequence = prot_dict.get(accession, "")
    center = site_pos - 1
    chars = []
    for idx in range(center - half_len, center + half_len + 1):
        chars.append(sequence[idx] if 0 <= idx < len(sequence) else "-")
    return "".join(chars)


def parse_peptide(peptide):
    """
    Strip modification annotations and return zero-based N positions carrying +0.98-like mass shifts.

    Malformed unmatched parentheses raise a clear error instead of risking a non-advancing loop.
    """
    peptide = str(peptide)
    n_positions = []
    bare_seq = ""
    i = 0
    while i < len(peptide):
        if peptide[i] == "N":
            mods = []
            j = i + 1
            while j < len(peptide) and peptide[j] == "(":
                k = peptide.find(")", j)
                if k == -1:
                    raise ValueError(f"Malformed peptide modification annotation: {peptide}")
                mods.append(peptide[j + 1:k])
                j = k + 1
            if contains_target_mass_shift(mods, target=0.98, tolerance=0.01):
                n_positions.append(len(bare_seq))
            bare_seq += "N"
            i = j
        elif peptide[i] == "(":
            k = peptide.find(")", i)
            if k == -1:
                raise ValueError(f"Malformed peptide modification annotation: {peptide}")
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
    """Return non-empty index bins; safe when requested bin_count exceeds n_sites."""
    n_sites = int(n_sites)
    bin_count = int(bin_count)
    if n_sites <= 0:
        return []
    if bin_count <= 0:
        raise ValueError("bin_count must be positive")
    effective_bins = min(bin_count, n_sites)
    return [arr for arr in np.array_split(np.arange(n_sites), effective_bins) if len(arr) > 0]

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


def paired_one_tailed_pvalue(protein_group_vals, site_group_vals, method="wilcoxon", side="right"):
    """Compute a one-tailed paired p-value from matched Protein-Group/Site-Group normalized Rm pairs."""
    protein_group_vals = np.asarray(protein_group_vals, dtype=float)
    site_group_vals = np.asarray(site_group_vals, dtype=float)

    valid_mask = (~np.isnan(protein_group_vals)) & (~np.isnan(site_group_vals))
    protein_valid = protein_group_vals[valid_mask]
    site_group_valid = site_group_vals[valid_mask]
    n_pairs = len(protein_valid)

    if n_pairs < 2:
        return 1.0, n_pairs

    diffs = site_group_valid - protein_valid
    if np.allclose(diffs, 0):
        return 1.0, n_pairs

    if method == "paired_t":
        t_stat, p_two = stats.ttest_rel(site_group_valid, protein_valid, nan_policy="omit")
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
                site_group_valid,
                protein_valid,
                alternative=alt,
                zero_method="wilcox",
                method="auto"
            )
        except TypeError:
            _, p_one = stats.wilcoxon(
                site_group_valid,
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
stage("Stage 1/10 - Protein Group quantification QC")
file = normalize_path(input("Please provide path of TMT quantification protein CSV file: ").strip())
if not os.path.exists(file):
    raise FileNotFoundError(file)

work_dir = os.path.dirname(file)
os.chdir(work_dir)
output_root, tables_dir, figures_dir, metadata_dir, intermediate_dir = build_output_dirs(work_dir)

df_raw = pd.read_csv(file)
protein_group_col = resolve_column(df_raw, PROTEIN_GROUP_CANDIDATES, "Protein Group")
accession_col = resolve_column(df_raw, ACCESSION_CANDIDATES, "protein Accession")

df_raw = df_raw.copy()
df_raw["protein_group_id"] = df_raw[protein_group_col].apply(canonical_group_id)
df_raw["uniprot_id"] = df_raw[accession_col].apply(canonical_protein_id)
if (df_raw["protein_group_id"] == "").any():
    raise ValueError("Empty Protein Group identifiers detected in protein CSV.")
if (df_raw["uniprot_id"] == "").any():
    raise ValueError("Empty/invalid Accession identifiers detected in protein CSV.")

sample_order, protein_sample_channel_map = detect_sample_channel_map(df_raw)
if not sample_order:
    raise ValueError(
        "Could not detect PEAKS Sample N reporter columns in the protein CSV. "
        "For deterministic G2 analysis, columns must explicitly encode Sample N and TMT 126-131; "
        "positional fallback grouping has been removed."
    )

sample_names = [
    input(f"Please enter name for Sample {s} (use 'control' for control group): ").strip()
    for s in sample_order
]
if any(not name for name in sample_names):
    raise ValueError("Sample names cannot be empty.")

df_raw, sample_block_cols, sample_channel_cols, protein_col_map = rename_sample_channels(
    df_raw, sample_order, protein_sample_channel_map, sample_names
)
all_reporter_cols = [c for cols in sample_block_cols.values() for c in cols]
df_raw = ensure_numeric_columns(df_raw, all_reporter_cols, "protein reporter intensity")

# Protein Group is the protein-level quantitative unit. Member accessions are retained
# as biological/sequence annotations but do not independently weight normalization.
validate_group_quantitative_profiles(df_raw, "protein_group_id", all_reporter_cols)
protein_member_map = (
    df_raw[["protein_group_id", "uniprot_id", accession_col]]
    .drop_duplicates()
    .rename(columns={accession_col: "accession_raw"})
    .sort_values(["protein_group_id", "uniprot_id"])
    .reset_index(drop=True)
)
protein_member_map.to_csv(
    os.path.join(tables_dir, "00_protein_group_member_map.csv"), index=False
)
member_accessions = protein_member_map.groupby("protein_group_id")["uniprot_id"].apply(
    lambda s: ";".join(sorted(set(str(x) for x in s if str(x))))
)
member_counts = protein_member_map.groupby("protein_group_id")["uniprot_id"].nunique()

# Unique-peptide threshold is applied to PEAKS accession rows first; a Protein Group is
# retained when at least one member row satisfies the threshold. Quantification is then
# collapsed to one row per Protein Group.
unique_col_candidates = ["#Unique", "#Unique Peptides", "Unique"]
unique_col = next((c for c in unique_col_candidates if c in df_raw.columns), None)
min_unique_peptides = None
if unique_col:
    df_raw[unique_col] = pd.to_numeric(df_raw[unique_col], errors="raise")
    raw_min_unique = input(
        "Minimum number of unique peptides required per PEAKS accession row "
        "before Protein Group collapse (default=1; use 0 to disable): "
    ).strip()
    min_unique_peptides = 1 if raw_min_unique == "" else int(raw_min_unique)
    if min_unique_peptides < 0:
        raise ValueError("Minimum unique-peptide count cannot be negative.")
    unique_pass = df_raw[unique_col] >= min_unique_peptides
else:
    unique_pass = pd.Series(True, index=df_raw.index)
    print("No unique peptide column detected; unique-peptide filtering skipped.")

eligible_rows = df_raw.loc[unique_pass].copy()
if eligible_rows.empty:
    raise ValueError("No protein rows remain after the unique-peptide filter.")

df = eligible_rows.drop_duplicates(subset=["protein_group_id"], keep="first").copy()
df["member_accessions"] = df["protein_group_id"].map(member_accessions)
df["member_accession_count"] = df["protein_group_id"].map(member_counts).astype(int)

print(f"Raw PEAKS protein rows: {len(df_raw)}")
print(f"Unique PEAKS Protein Groups represented: {df_raw['protein_group_id'].nunique()}")
print(f"Protein Groups after unique-peptide filtering: {len(df)}")
print(f"Unique member UniProt accessions: {protein_member_map['uniprot_id'].nunique()}")

# Quantified overview using Protein Group counts, not accession-row counts.
quantified_any = (df[all_reporter_cols] != 0).any(axis=1)
print(f"Protein Groups quantified in any reporter channel: {int(quantified_any.sum())}")

per_sample_sets = {}
for sample in sample_names:
    cols = sample_block_cols[sample]
    per_sample_sets[sample] = set(df.loc[(df[cols] != 0).any(axis=1), "protein_group_id"])

num_samples = len(sample_names)
if num_samples > 1:
    if num_samples <= 3:
        plt.figure(figsize=(6, 6))
        if num_samples == 2:
            venn2([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
        else:
            venn3([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
        plt.title("Quantified PEAKS Protein Groups")
        plt.savefig(os.path.join(figures_dir, "01_protein_group_quantification_overlap.tiff"), dpi=300)
        plt.close()
    else:
        if from_contents is None or UpSet is None:
            raise ImportError(
                "upsetplot is required when more than three samples are present. "
                "Install it with: pip install upsetplot"
            )
        upset_data = from_contents(per_sample_sets)
        plt.figure(figsize=(8, 6))
        UpSet(upset_data, show_counts=True).plot()
        plt.savefig(os.path.join(figures_dir, "01_protein_group_quantification_overlap_upset.tiff"), dpi=300)
        plt.close()

strict_mask = pd.Series(True, index=df.index)
for sample in sample_names:
    strict_mask &= (df[sample_block_cols[sample]] != 0).all(axis=1)
df_strict = df.loc[strict_mask].copy()

for sample in sample_names:
    n_complete = int((df[sample_block_cols[sample]] != 0).all(axis=1).sum())
    print(f"Sample {sample}: {n_complete} Protein Groups fully quantified in all six channels")
print(f"Protein Groups fully quantified across all samples: {len(df_strict)}")

control_sample = next((s for s in sample_names if "control" in s.lower()), sample_names[0])
control_cols = sample_block_cols[control_sample]
df_control_complete = df.loc[(df[control_cols] != 0).all(axis=1)].copy()
df_control_complete.to_csv(
    os.path.join(tables_dir, "01_control_protein_groups_complete_quantification.csv"), index=False
)
df_strict.to_csv(
    os.path.join(tables_dir, "02_all_samples_protein_groups_complete_quantification.csv"), index=False
)
print("TMT Protein Group QC completed.")

# ------------------------------
# 2. Rm calculation & normalization of control group
# ------------------------------
stage("Stage 2/10 - Protein Group Rm calculation and replicate-median normalization")
ctrl_channel_cols = sample_channel_cols[control_sample]
ch126 = ctrl_channel_cols[126]
ch127 = ctrl_channel_cols[127]
ch128 = ctrl_channel_cols[128]
ch129 = ctrl_channel_cols[129]
ch130 = ctrl_channel_cols[130]
ch131 = ctrl_channel_cols[131]

df_control_complete["Rm1"] = df_control_complete[ch129] / df_control_complete[ch127]
df_control_complete["Rm2"] = df_control_complete[ch130] / df_control_complete[ch128]
df_control_complete["Rm3"] = df_control_complete[ch131] / df_control_complete[ch126]

# Each unique PEAKS Protein Group contributes once to the replicate median.
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
    "Protein_Quantitative_Unit": ["PEAKS Protein Group"] * 3,
})
cf_df.to_csv(os.path.join(tables_dir, "03_rm_normalization_factors.csv"), index=False)
export_csv(
    df_control_complete,
    os.path.join(tables_dir, "04_control_protein_groups_with_rm.csv"),
    PROTEIN_EXPORT_RENAME,
)
print("Rm normalization completed using unique PEAKS Protein Groups.")

# ------------------------------
# 3. N-glycopeptide analysis
# ------------------------------
stage("Stage 3/10 - N-glycopeptide evidence collapse and sequence/site mapping")
peptide_csv_file = normalize_path(input("Enter path of TMT quantification peptide CSV file: ").strip())
fasta_file = normalize_path(input("Enter path to UniProt fasta file: ").strip())
if not os.path.exists(peptide_csv_file):
    raise FileNotFoundError(peptide_csv_file)
if not os.path.exists(fasta_file):
    raise FileNotFoundError(fasta_file)

fasta_csv_file = os.path.join(intermediate_dir, "uniprot_fasta_parsed.csv")
with open(fasta_file, "r", encoding="utf-8") as f, open(fasta_csv_file, "w", newline="", encoding="utf-8") as csvfile:
    import csv
    writer = csv.writer(csvfile)
    writer.writerow(["Accession", "Protein Name", "Sequence"])
    accession = ""
    protein_name = ""
    seq_lines = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if accession:
                writer.writerow([accession, protein_name, "".join(seq_lines)])
            header = line[1:].strip()
            token = header.split()[0] if header else ""
            accession = canonical_protein_id(token) if token else ""
            parts = [p.strip() for p in token.split("|")] if "|" in token else []
            if len(parts) >= 3:
                protein_name = parts[2]
            else:
                tokens = header.split()
                protein_name = " ".join(tokens[1:]).strip() if len(tokens) > 1 else ""
            seq_lines = []
        else:
            seq_lines.append(line)
    if accession:
        writer.writerow([accession, protein_name, "".join(seq_lines)])
print(f"Fasta converted to CSV: {fasta_csv_file}")

pep_df = pd.read_csv(peptide_csv_file)
prot_df = pd.read_csv(fasta_csv_file)
pep_accession_col = resolve_column(pep_df, ACCESSION_CANDIDATES, "peptide Accession")
pep_peptide_col = resolve_column(pep_df, PEPTIDE_CANDIDATES, "peptide sequence")
pep_group_col = resolve_column(pep_df, PROTEIN_GROUP_CANDIDATES, "peptide Protein Group", required=False)
pep_df = pep_df.copy()
pep_df["uniprot_id"] = pep_df[pep_accession_col].apply(canonical_protein_id)
if (pep_df["uniprot_id"] == "").any():
    raise ValueError("Empty/invalid Accession identifiers detected in peptide CSV.")

if pep_group_col is not None:
    pep_df["protein_group_id"] = pep_df[pep_group_col].apply(canonical_group_id)
else:
    accession_group_counts = protein_member_map.groupby("uniprot_id")["protein_group_id"].nunique()
    ambiguous_accessions = accession_group_counts[accession_group_counts > 1]
    if not ambiguous_accessions.empty:
        raise ValueError(
            "Peptide CSV has no Protein Group column and some accessions map to multiple Protein Groups; "
            "cannot infer a unique protein-group reference."
        )
    accession_to_group = (
        protein_member_map.drop_duplicates("uniprot_id").set_index("uniprot_id")["protein_group_id"].to_dict()
    )
    pep_df["protein_group_id"] = pep_df["uniprot_id"].map(accession_to_group)
    if pep_df["protein_group_id"].isna().any():
        missing = pep_df.loc[pep_df["protein_group_id"].isna(), "uniprot_id"].drop_duplicates().head(10).tolist()
        raise ValueError(f"Unable to map peptide accessions to Protein Groups: {missing}")

pep_sample_order, peptide_sample_channel_map = detect_sample_channel_map(pep_df)
if not pep_sample_order:
    raise ValueError(
        "Could not detect Sample N / TMT126-131 reporter columns in peptide CSV. "
        "Expected explicit PEAKS-style sample/channel names."
    )
peptide_sample_names = [
    input(f"Enter new name for Sample {sn} in peptide table: ").strip()
    for sn in pep_sample_order
]
if any(not name for name in peptide_sample_names):
    raise ValueError("Peptide sample names cannot be empty.")
pep_df, peptide_sample_block_cols, peptide_sample_channel_cols, peptide_col_map = rename_sample_channels(
    pep_df, pep_sample_order, peptide_sample_channel_map, peptide_sample_names
)
renamed_intensity_cols = [c for cols in peptide_sample_block_cols.values() for c in cols]
pep_df = ensure_numeric_columns(pep_df, renamed_intensity_cols, "peptide reporter intensity")
prot_dict = dict(zip(prot_df["Accession"].astype(str), prot_df["Sequence"].astype(str)))

# Collapse duplicated accession assignments when they represent the same PEAKS
# quantitative evidence: same modified peptide and same complete reporter profile.
# If that one evidence record is assigned across multiple Protein Groups, it is retained
# for provenance/QC but excluded from canonical DeltaRm because the protein reference is not unique.
peptide_evidence_rows = []
evidence_counter = 1
for modified_peptide, peptide_group in pep_df.groupby(pep_peptide_col, dropna=False, sort=False):
    profile_groups = peptide_group.groupby(renamed_intensity_cols, dropna=False, sort=False)
    for profile_values, profile_group in profile_groups:
        if not isinstance(profile_values, tuple):
            profile_values = (profile_values,)
        assignment_pairs = sorted(set(
            (canonical_group_id(pg), str(acc))
            for pg, acc in zip(profile_group["protein_group_id"], profile_group["uniprot_id"])
            if canonical_group_id(pg) and str(acc)
        ))
        candidate_groups = sorted(set(pg for pg, _ in assignment_pairs))
        candidate_accessions = sorted(set(acc for _, acc in assignment_pairs))
        evidence_row = {
            "peptide_evidence_id": f"PE{evidence_counter:06d}",
            "protein_group_id": candidate_groups[0] if len(candidate_groups) == 1 else "",
            "candidate_protein_groups": ";".join(candidate_groups),
            "candidate_protein_group_count": len(candidate_groups),
            "modified_peptide": str(modified_peptide),
            "candidate_accessions": ";".join(candidate_accessions),
            "candidate_accession_count": len(candidate_accessions),
            "candidate_group_accession_assignments": ";".join(
                f"{pg}|{acc}" for pg, acc in assignment_pairs
            ),
        }
        for col, value in zip(renamed_intensity_cols, profile_values):
            evidence_row[col] = value
        if "AScore" in profile_group.columns:
            ascore_values = pd.to_numeric(profile_group["AScore"], errors="coerce")
            evidence_row["AScore"] = float(ascore_values.max()) if ascore_values.notna().any() else np.nan
        peptide_evidence_rows.append(evidence_row)
        evidence_counter += 1

peptide_evidence_df = pd.DataFrame(peptide_evidence_rows)
if peptide_evidence_df.empty:
    raise ValueError("No peptide quantitative evidence could be constructed from peptide CSV.")

results = []
for _, row in peptide_evidence_df.iterrows():
    peptide = row["modified_peptide"]
    stripped_seq, n_pos_list = parse_peptide(peptide)
    if not n_pos_list:
        continue

    candidate_site_assignments = set()
    valid_accessions = set()
    valid_groups = set()
    assignment_tokens = [
        x for x in str(row["candidate_group_accession_assignments"]).split(";") if x
    ]
    for token in assignment_tokens:
        protein_group_id, uni_id = token.split("|", 1)
        prot_seq = prot_dict.get(uni_id)
        if not prot_seq:
            continue
        starts = find_all_occurrences(prot_seq, stripped_seq)
        for start_idx in starts:
            for n_pos in n_pos_list:
                prot_n_pos = start_idx + n_pos
                if prot_n_pos + 2 >= len(prot_seq):
                    continue
                x_residue = prot_seq[prot_n_pos + 1]
                st_residue = prot_seq[prot_n_pos + 2]
                if x_residue != "P" and st_residue in {"S", "T"}:
                    candidate_site_assignments.add(
                        f"{protein_group_id}|{uni_id}:N{prot_n_pos + 1}"
                    )
                    valid_accessions.add(uni_id)
                    valid_groups.add(protein_group_id)

    if not candidate_site_assignments:
        continue

    site_assignment_sorted = sorted(candidate_site_assignments)
    candidate_sites_sorted = sorted(set(x.split("|", 1)[1] for x in site_assignment_sorted))
    valid_groups_sorted = sorted(valid_groups)
    if len(valid_groups_sorted) > 1:
        mapping_status = "ambiguous_across_protein_groups"
        canonical_group = ""
    elif len(candidate_sites_sorted) == 1:
        mapping_status = "unique"
        canonical_group = valid_groups_sorted[0]
    else:
        mapping_status = "ambiguous_within_protein_group"
        canonical_group = valid_groups_sorted[0]

    out_row = {
        "peptide_evidence_id": row["peptide_evidence_id"],
        "protein_group_id": canonical_group,
        "candidate_protein_groups": ";".join(valid_groups_sorted),
        "candidate_protein_group_count": len(valid_groups_sorted),
        "candidate_accessions": ";".join(sorted(valid_accessions)),
        "candidate_accession_count": len(valid_accessions),
        "candidate_sites": ";".join(candidate_sites_sorted),
        "candidate_site_assignments": ";".join(site_assignment_sorted),
        "candidate_site_count": len(candidate_sites_sorted),
        "mapping_status": mapping_status,
        "modified_peptide": peptide,
        "stripped_sequence": stripped_seq,
        "modified_n_positions_in_peptide": ";".join(str(x + 1) for x in n_pos_list),
    }
    for col in renamed_intensity_cols:
        out_row[col] = row[col]
    if "AScore" in peptide_evidence_df.columns:
        out_row["AScore"] = row.get("AScore", np.nan)
    results.append(out_row)

results_df = pd.DataFrame(results)
if results_df.empty:
    raise ValueError(
        "No N-glycopeptide evidence remained after +0.98 parsing, FASTA mapping, and N-X-S/T validation."
    )

output_csv_file = os.path.join(tables_dir, "05_n_glycopeptide_evidence_with_intensity.csv")
results_df.to_csv(output_csv_file, index=False)
print(f"N-glycopeptide evidence table saved: {output_csv_file}")
print(f"Raw peptide rows: {len(pep_df)}")
print(f"Unique peptide quantitative evidence units: {len(peptide_evidence_df)}")
print(f"Valid N-glycopeptide evidence units: {len(results_df)}")
print(
    "Ambiguous within-Protein-Group evidence units: "
    f"{int((results_df['mapping_status'] == 'ambiguous_within_protein_group').sum())}"
)
print(
    "Ambiguous across-Protein-Group evidence units (retained for QC, excluded from DeltaRm): "
    f"{int((results_df['mapping_status'] == 'ambiguous_across_protein_groups').sum())}"
)

# Intensity composition uses de-duplicated quantitative evidence, not accession-expanded rows.
glyco_intensity = results_df[renamed_intensity_cols].sum().values
total_intensity = peptide_evidence_df[renamed_intensity_cols].sum().values
non_glyco_intensity = np.maximum(total_intensity - glyco_intensity, 0)
ratios = np.divide(glyco_intensity, total_intensity, out=np.zeros_like(glyco_intensity, dtype=float), where=total_intensity != 0)
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(renamed_intensity_cols))
ax.bar(x, glyco_intensity, width=0.6, color="#E63946", label="N-glycopeptide evidence")
ax.bar(x, non_glyco_intensity, bottom=glyco_intensity, width=0.6, color="#457B9D", label="Other peptide evidence")
for i, ratio in enumerate(ratios):
    ax.text(x[i], glyco_intensity[i] / 2 if glyco_intensity[i] else 0, f"{ratio:.2f}", ha="center", va="center", color="white", fontsize=10, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(renamed_intensity_cols, rotation=45, ha="right")
ax.set_ylabel("Total Intensity")
ax.set_title("Intensity Composition of N-glycopeptide Quantitative Evidence")
ax.legend()
plt.tight_layout()
fig.savefig(os.path.join(figures_dir, "02_n_glycopeptide_intensity_composition.tiff"), dpi=300, format="tiff")
plt.close(fig)

# ----------------------------
# 4. Site Groups
# ----------------------------
stage("Stage 4/10 - Site Group aggregation and completeness QC")
motif_len = int(input("Enter total motif window length (odd number, e.g., 7, 9, 11): ").strip())
if motif_len <= 0 or motif_len % 2 == 0:
    raise ValueError("Motif window length must be a positive odd number.")
half_len = motif_len // 2

site_results = []
site_counter = 1
site_mappable_df = results_df[
    results_df["mapping_status"] != "ambiguous_across_protein_groups"
].copy()
for (protein_group_id, candidate_sites), group in site_mappable_df.groupby(
    ["protein_group_id", "candidate_sites"], dropna=False, sort=True
):
    site_labels = [x for x in str(candidate_sites).split(";") if x]
    if not site_labels:
        continue
    candidate_accessions = sorted(set(parse_site_label(x)[0] for x in site_labels))
    motif_annotations = [
        f"{label}:{build_site_motif(prot_dict, label, half_len)}" for label in site_labels
    ]
    site_row = {
        "site_group_id": f"SG{site_counter:06d}",
        "protein_group_id": canonical_group_id(protein_group_id),
        "candidate_sites": ";".join(site_labels),
        "candidate_site_count": len(site_labels),
        "candidate_accessions": ";".join(candidate_accessions),
        "candidate_accession_count": len(candidate_accessions),
        "mapping_status": "unique" if len(site_labels) == 1 else "ambiguous_within_protein_group",
        "candidate_motifs": ";".join(motif_annotations),
        "supporting_evidence_ids": ";".join(sorted(set(group["peptide_evidence_id"].astype(str)))),
        "supporting_evidence_count": int(group["peptide_evidence_id"].nunique()),
        "supporting_peptides": ";".join(sorted(set(group["modified_peptide"].astype(str)))),
    }
    for col in renamed_intensity_cols:
        site_row[col] = group[col].sum()
    site_results.append(site_row)
    site_counter += 1

site_df = pd.DataFrame(site_results)
if site_df.empty:
    raise ValueError("No Site Groups were constructed.")
site_csv_file = os.path.join(tables_dir, "06_nglycosite_site_groups.csv")
site_df.to_csv(site_csv_file, index=False)
print(f"Site Group table saved: {site_csv_file}")
print(f"Protein Groups with mapped N-glycosite evidence: {site_df['protein_group_id'].nunique()}")
print(f"Site Groups: {len(site_df)}")
print(f"Candidate site annotations represented: {int(site_df['candidate_site_count'].sum())}")
print(f"Unique Site Groups: {int((site_df['mapping_status'] == 'unique').sum())}")
print(f"Ambiguous within-group Site Groups: {int((site_df['mapping_status'] != 'unique').sum())}")

# Distribution of Site Groups per Protein Group.
pg_site_counts = site_df.groupby("protein_group_id")["site_group_id"].nunique()
pg_site_counts_capped = pg_site_counts.apply(lambda x: x if x <= 5 else ">5")
dist = pg_site_counts_capped.value_counts().sort_index(
    key=lambda x: [int(i) if i != ">5" else 6 for i in x]
)
plt.figure(figsize=(8, 6))
bars = plt.barh(dist.index.astype(str), dist.values, color="skyblue")
for bar in bars:
    width = bar.get_width()
    plt.text(width + max(dist.values) * 0.01, bar.get_y() + bar.get_height() / 2, str(width), va="center")
plt.xlabel("Protein Group count")
plt.ylabel("Number of Site Groups per Protein Group")
plt.title("Distribution of N-glycosite Site Groups per Protein Group")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "03_protein_group_nglycosite_group_distribution_all.tiff"), dpi=300, format="tiff")
plt.close()

chosen_group = input("Enter the peptide sample group to calculate DeltaRm (exact renamed sample name, e.g., 'N-Glyco'): ").strip()
if chosen_group not in peptide_sample_channel_cols:
    raise ValueError(
        f"Unknown peptide sample group '{chosen_group}'. Available groups: {list(peptide_sample_channel_cols)}"
    )
site_channel_map = peptide_sample_channel_cols[chosen_group]
ch126_site = site_channel_map[126]
ch127_site = site_channel_map[127]
ch128_site = site_channel_map[128]
ch129_site = site_channel_map[129]
ch130_site = site_channel_map[130]
ch131_site = site_channel_map[131]
site_cols_for_calc = [ch126_site, ch127_site, ch128_site, ch129_site, ch130_site, ch131_site]

site_df_clean = site_df.copy()
site_df_clean[site_cols_for_calc] = site_df_clean[site_cols_for_calc].replace(0, np.nan)
complete_mask = site_df_clean[site_cols_for_calc].notna().all(axis=1)
site_df_complete = site_df_clean.loc[complete_mask].copy()
print(f"Fully quantified Site Groups: {len(site_df_complete)}")
print(f"Protein Groups containing fully quantified Site Groups: {site_df_complete['protein_group_id'].nunique()}")
site_complete_csv = os.path.join(tables_dir, "07_nglycosite_site_groups_complete_quantification.csv")
site_df_complete.to_csv(site_complete_csv, index=False)

pg_site_count_complete = site_df_complete.groupby("protein_group_id")["site_group_id"].nunique()
pg_site_counts_capped = pg_site_count_complete.apply(lambda x: x if x <= 5 else ">5")
dist = pg_site_counts_capped.value_counts().sort_index(key=lambda x: [int(i) if i != ">5" else 6 for i in x])
plt.figure(figsize=(8, 6))
bars = plt.barh(dist.index.astype(str), dist.values, color="#A8C9E0")
for bar in bars:
    width = bar.get_width()
    plt.text(width + max(dist.values) * 0.01, bar.get_y() + bar.get_height() / 2, str(width), va="center")
plt.xlabel("Protein Group count")
plt.ylabel("Fully quantified Site Groups per Protein Group")
plt.title("Distribution of fully quantified N-glycosite Site Groups per Protein Group")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "04_protein_group_nglycosite_group_distribution_complete.tiff"), dpi=300, format="tiff")
plt.close()

# ----------------------------
# 5. Match Site Groups to fully quantified control Protein Groups
# ----------------------------
stage("Stage 5/10 - Match fully quantified Site Groups to fully quantified control Protein Groups")
protein_groups_complete = set(df_control_complete["protein_group_id"].astype(str))
site_df_filtered = site_df_complete[
    site_df_complete["protein_group_id"].astype(str).isin(protein_groups_complete)
].copy()
matched_groups = set(site_df_filtered["protein_group_id"].astype(str))
df_control_filtered = df_control_complete[
    df_control_complete["protein_group_id"].astype(str).isin(matched_groups)
].copy()

plt.figure(figsize=(5, 5))
v = venn2(
    [protein_groups_complete, set(site_df_complete["protein_group_id"].astype(str))],
    set_labels=["Fully quantified control Protein Groups", "Protein Groups with fully quantified Site Groups"],
    set_colors=("#8AB7C2", "#9EAAD4"),
    alpha=0.5,
)
for text_obj in list(v.set_labels) + list(v.subset_labels):
    if text_obj:
        text_obj.set_fontsize(6)
plt.title("Protein Group / Site Group quantification completeness overlap", fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "05_protein_group_site_quantification_overlap_venn.tiff"), dpi=300, format="tiff")
plt.close()

print(f"Matched Protein Groups: {len(df_control_filtered)}")
print(f"Matched Site Groups: {len(site_df_filtered)}")
print(f"Candidate site annotations in matched units: {int(site_df_filtered['candidate_site_count'].sum())}")

pg_site_count_intersection = site_df_filtered.groupby("protein_group_id")["site_group_id"].nunique()
pg_site_counts_capped = pg_site_count_intersection.apply(lambda x: x if x <= 5 else ">5")
dist = pg_site_counts_capped.value_counts().sort_index(key=lambda x: [int(i) if i != ">5" else 6 for i in x])
plt.figure(figsize=(8, 6))
bars = plt.barh(dist.index.astype(str), dist.values, color="#14C0CC")
for bar in bars:
    width = bar.get_width()
    plt.text(width + max(dist.values) * 0.01, bar.get_y() + bar.get_height() / 2, str(width), va="center")
plt.xlabel("Protein Group count")
plt.ylabel("Matched Site Groups per Protein Group")
plt.title("Distribution of matched N-glycosite Site Groups per Protein Group")
plt.tight_layout()
plt.savefig(os.path.join(figures_dir, "06_protein_group_nglycosite_group_distribution_intersection.tiff"), dpi=300, format="tiff")
plt.close()

# ----------------------------
# 6. Rm, CV, normalization and DeltaRm
# ----------------------------
stage("Stage 6/10 - Raw Rm CV QC, normalized Rm, and mean DeltaRm calculation")
site_df_complete = site_df_complete.copy()
site_df_complete["Rm_sg1"] = site_df_complete[ch129_site] / site_df_complete[ch127_site]
site_df_complete["Rm_sg2"] = site_df_complete[ch130_site] / site_df_complete[ch128_site]
site_df_complete["Rm_sg3"] = site_df_complete[ch131_site] / site_df_complete[ch126_site]
site_df_complete["CV_Rm_sg_complete"] = coefficient_of_variation_from_raw_rm(
    site_df_complete, ["Rm_sg1", "Rm_sg2", "Rm_sg3"]
)
complete_outfile = os.path.join(tables_dir, f"08_nglycosite_site_groups_rm_raw_cv_{chosen_group}_complete.csv")
export_csv(site_df_complete, complete_outfile, SITE_EXPORT_RENAME)

site_df_filtered = site_df_filtered.copy()
site_df_filtered["Rm_sg1"] = site_df_filtered[ch129_site] / site_df_filtered[ch127_site]
site_df_filtered["Rm_sg2"] = site_df_filtered[ch130_site] / site_df_filtered[ch128_site]
site_df_filtered["Rm_sg3"] = site_df_filtered[ch131_site] / site_df_filtered[ch126_site]
site_df_filtered["CV_Rm_sg_filtered"] = coefficient_of_variation_from_raw_rm(
    site_df_filtered, ["Rm_sg1", "Rm_sg2", "Rm_sg3"]
)
site_df_filtered["Rm_sg1_norm"] = site_df_filtered["Rm_sg1"] * CF["Rm1"]
site_df_filtered["Rm_sg2_norm"] = site_df_filtered["Rm_sg2"] * CF["Rm2"]
site_df_filtered["Rm_sg3_norm"] = site_df_filtered["Rm_sg3"] * CF["Rm3"]

# Protein Group Rm/CV summary.
df_control_filtered = df_control_filtered.copy()
df_control_filtered["mean_Rm_pg_norm"] = df_control_filtered[
    ["Rm1_normalized", "Rm2_normalized", "Rm3_normalized"]
].mean(axis=1)
df_control_filtered["CV_Rm_filtered"] = coefficient_of_variation_from_raw_rm(
    df_control_filtered, ["Rm1", "Rm2", "Rm3"]
)
df_control_complete = df_control_complete.copy()
df_control_complete["CV_Rm"] = coefficient_of_variation_from_raw_rm(
    df_control_complete, ["Rm1", "Rm2", "Rm3"]
)

new_protein_cols = ["Rm1", "Rm2", "Rm3", "Rm1_normalized", "Rm2_normalized", "Rm3_normalized", "mean_Rm_pg_norm", "CV_Rm_filtered"]
original_protein_cols = [c for c in df_control_filtered.columns if c not in new_protein_cols]
export_csv(
    df_control_filtered,
    os.path.join(tables_dir, "09_control_protein_groups_rm_summary_site_matched.csv"),
    PROTEIN_EXPORT_RENAME,
    columns=original_protein_cols + new_protein_cols,
)
new_protein_cols_complete = ["Rm1", "Rm2", "Rm3", "CV_Rm"]
original_protein_cols_complete = [c for c in df_control_complete.columns if c not in new_protein_cols_complete]
export_csv(
    df_control_complete,
    os.path.join(tables_dir, "10_control_protein_groups_rm_summary_complete.csv"),
    PROTEIN_EXPORT_RENAME,
    columns=original_protein_cols_complete + new_protein_cols_complete,
)

# Site Groups use their matched Protein Group as the protein reference.
site_df_filtered["mean_Rm_sg_norm"] = site_df_filtered[
    ["Rm_sg1_norm", "Rm_sg2_norm", "Rm_sg3_norm"]
].mean(axis=1)
pg_indexed = df_control_filtered.set_index("protein_group_id")
site_df_filtered["mean_Rm_pg_norm"] = site_df_filtered["protein_group_id"].map(pg_indexed["mean_Rm_pg_norm"])
site_df_filtered["Rm_pg1_norm_ref"] = site_df_filtered["protein_group_id"].map(pg_indexed["Rm1_normalized"])
site_df_filtered["Rm_pg2_norm_ref"] = site_df_filtered["protein_group_id"].map(pg_indexed["Rm2_normalized"])
site_df_filtered["Rm_pg3_norm_ref"] = site_df_filtered["protein_group_id"].map(pg_indexed["Rm3_normalized"])
site_df_filtered["ΔRm_rep1"] = site_df_filtered["Rm_sg1_norm"] - site_df_filtered["Rm_pg1_norm_ref"]
site_df_filtered["ΔRm_rep2"] = site_df_filtered["Rm_sg2_norm"] - site_df_filtered["Rm_pg2_norm_ref"]
site_df_filtered["ΔRm_rep3"] = site_df_filtered["Rm_sg3_norm"] - site_df_filtered["Rm_pg3_norm_ref"]
site_df_filtered["ΔRm"] = site_df_filtered[["ΔRm_rep1", "ΔRm_rep2", "ΔRm_rep3"]].mean(axis=1)

new_site_cols = [
    "Rm_sg1", "Rm_sg2", "Rm_sg3", "CV_Rm_sg_filtered",
    "Rm_sg1_norm", "Rm_sg2_norm", "Rm_sg3_norm",
    "mean_Rm_sg_norm", "mean_Rm_pg_norm",
    "Rm_pg1_norm_ref", "Rm_pg2_norm_ref", "Rm_pg3_norm_ref",
    "ΔRm_rep1", "ΔRm_rep2", "ΔRm_rep3", "ΔRm",
]
original_site_cols = [c for c in site_df_filtered.columns if c not in new_site_cols]
site_outfile = os.path.join(tables_dir, "11_nglycosite_site_groups_delta_rm.csv")
export_csv(site_df_filtered, site_outfile, SITE_EXPORT_RENAME, columns=original_site_cols + new_site_cols)
print(f"Site Group DeltaRm table saved: {site_outfile}")

# CV QC plots.
cv_plot_params = [
    ("All fully quantified Site Groups", site_df_complete["CV_Rm_sg_complete"], f"CV_distribution_complete_site_groups_{chosen_group}.tiff"),
    ("Matched fully quantified Site Groups", site_df_filtered["CV_Rm_sg_filtered"], f"CV_distribution_matched_site_groups_{chosen_group}.tiff"),
]
for label, cv_series_raw, fname in cv_plot_params:
    cv_series = cv_series_raw.dropna().sort_values(ascending=False).reset_index(drop=True)
    if cv_series.empty:
        continue
    rank = cv_series.index + 1
    plt.figure(figsize=(6, 4))
    plt.scatter(rank, cv_series.values, color="#5E556A", alpha=0.7, s=20)
    prop_30 = (cv_series <= 0.3).mean() * 100
    prop_40 = (cv_series <= 0.4).mean() * 100
    plt.text(0.7 * len(cv_series), 0.305, f"≤30%: {prop_30:.1f}%", fontsize=10)
    plt.text(0.7 * len(cv_series), 0.405, f"≤40%: {prop_40:.1f}%", fontsize=10)
    plt.xlabel("Site Groups sorted by CV (high→low)")
    plt.ylabel("CV")
    plt.title(f"CV distribution - {chosen_group} ({label})", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, fname), dpi=300)
    plt.close()

protein_cv_plot_params = [
    ("All fully quantified Protein Groups", df_control_complete["CV_Rm"], "CV_distribution_complete_protein_groups.tiff"),
    ("Protein Groups with matched Site Groups", df_control_filtered["CV_Rm_filtered"], "CV_distribution_matched_protein_groups.tiff"),
]
for label, cv_series_raw, fname in protein_cv_plot_params:
    cv_series = cv_series_raw.dropna().sort_values(ascending=False).reset_index(drop=True)
    if cv_series.empty:
        continue
    rank = cv_series.index + 1
    plt.figure(figsize=(6, 4))
    plt.scatter(rank, cv_series.values, color="#2A9D8F", alpha=0.7, s=20)
    prop_30 = (cv_series <= 0.3).mean() * 100
    prop_40 = (cv_series <= 0.4).mean() * 100
    plt.text(0.7 * len(cv_series), 0.305, f"≤30%: {prop_30:.1f}%", fontsize=10)
    plt.text(0.7 * len(cv_series), 0.405, f"≤40%: {prop_40:.1f}%", fontsize=10)
    plt.xlabel("Protein Groups sorted by CV (high→low)")
    plt.ylabel("CV")
    plt.title(f"CV distribution - Control Protein Groups ({label})", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, fname), dpi=300)
    plt.close()

# ------------------------------
# Step 8: CV-based QC filtering for sites (user threshold)
# ------------------------------
stage("Stage 7/10 - Joint Site Group / Protein-Group raw-Rm CV filtering")
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

# 2) Check required columns (Site-Group CV column and matched Protein-Group CV column)
site_cv_col = "CV_Rm_sg_filtered"   # Site-Group-level raw-Rm CV column
protein_cv_col = "CV_Rm_filtered"      # Protein-Group-level raw-Rm CV column

missing = []
if site_cv_col not in site_df_filtered.columns:
    missing.append(site_cv_col)
if protein_cv_col not in df_control_filtered.columns:
    missing.append(protein_cv_col)
if "protein_group_id" not in site_df_filtered.columns:
    missing.append("protein_group_id (in site_df_filtered)")
if "protein_group_id" not in df_control_filtered.columns:
    missing.append("protein_group_id (in df_control_filtered)")

if missing:
    raise KeyError(f"Required columns missing for CV filtering: {missing}")

# 3) Map Protein-Group-level CV
protein_cv_series = df_control_filtered.set_index("protein_group_id")[protein_cv_col]

# Report match statistics
mapped_proteins = site_df_filtered["protein_group_id"].isin(protein_cv_series.index)
n_mapped = mapped_proteins.sum()
n_total_sites = site_df_filtered.shape[0]
print(f"Site table: {n_total_sites} Site Groups total. {n_mapped} have matching Protein Group CV.")

# 4) Map matched Protein-Group CV into the Site Group table (keep all original columns including ΔRm)
site_df_filtered = site_df_filtered.copy()
site_df_filtered["Matched_PG_CV"] = site_df_filtered["protein_group_id"].map(protein_cv_series)

# 5) CV filtering (NaN treated as failing)
before_count = site_df_filtered.shape[0]
pass_mask = (
    (site_df_filtered[site_cv_col].notna()) &
    (site_df_filtered["Matched_PG_CV"].notna()) &
    (site_df_filtered[site_cv_col] <= cv_threshold) &
    (site_df_filtered["Matched_PG_CV"] <= cv_threshold)
)
site_df_CV_pass = site_df_filtered.loc[pass_mask].copy()
after_count = site_df_CV_pass.shape[0]

print(f"CV filtering: retained {after_count} / {before_count} Site Groups ({after_count/before_count*100 if before_count>0 else 0:.1f}%).")
print(f"Site units with missing Protein Group CV (excluded): {(~mapped_proteins).sum()}")

# 6) Save CV-filtered site table
cv_filtered_output = os.path.join(tables_dir, f"12_nglycosite_site_groups_cv_pass_{int(cv_threshold*100)}pct.csv")
export_csv(site_df_CV_pass, cv_filtered_output, SITE_EXPORT_RENAME)
print(f"CV-filtered site table saved: {cv_filtered_output}")

# 7) Optional QC plot: Site Group CV vs matched Protein Group CV
try:
    plt.figure(figsize=(6,6))
    plt.scatter(site_df_filtered[site_cv_col], site_df_filtered["Matched_PG_CV"], s=10, alpha=0.6)
    plt.axvline(cv_threshold, linestyle='--', linewidth=1, label=f"Site Group CV threshold ({cv_threshold:.2f})")
    plt.axhline(cv_threshold, linestyle='--', linewidth=1, label=f"Protein Group CV threshold ({cv_threshold:.2f})")
    plt.xlabel("Site Group CV (raw Rm)")
    plt.ylabel("Matched Protein Group CV (raw Rm)")
    plt.title(f"Site Group CV vs matched Protein Group CV (threshold {cv_threshold:.2f})")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    qc_plot_file = os.path.join(figures_dir, f"07_site_vs_protein_raw_rm_cv_{int(cv_threshold*100)}pct.tiff")
    plt.savefig(qc_plot_file, dpi=300, format='tiff')
    plt.close()
    print(f"CV scatter plot saved: {qc_plot_file}")
except Exception as e:
    print("Warning: failed to draw CV scatter plot:", e)

# 8) Final summary: number of Protein Groups with at least one CV-passing Site Group
kept_proteins = site_df_CV_pass["protein_group_id"].nunique()
print(f"Number of Protein Groups with ≥1 Site Group passing CV filter: {kept_proteins}")

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
    plt.title(f'{rep_col} Distribution (Site Groups passed CV QC)')
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
    print("At least one biological replicate is NOT normal; use Site-Group-level paired testing (default Wilcoxon, optional paired t-test override)")
else:
    print("All biological replicates are approximately normal, use CV-binned robust sigma testing per replicate")

# --- 6. Debug check: confirm original table unchanged ---
print("Step 9 check: site_df_CV_pass rows =", site_df_CV_pass.shape[0])

# ------------------------------
# Step 10-11: Significant site identification after CV QC (branch by normality)
# ------------------------------
stage("Stage 9/10 - Statistical significance analysis")

# 1) Compute representative CV for each Site Group (max of site CV and protein CV)
site_df_CV_pass = site_df_CV_pass.copy()
site_df_CV_pass["Rep_CV"] = site_df_CV_pass[[site_cv_col, "Matched_PG_CV"]].max(axis=1)

# 2) Sort Site Groups by representative CV (ascending)
site_df_sorted = site_df_CV_pass.sort_values("Rep_CV").reset_index(drop=True)
delta_matrix_sorted = site_df_sorted[delta_rep_cols].to_numpy(dtype=float)

if not delta_rm_not_normal:
    print("Step10 strategy: all replicate ΔRm distributions are approximately normal -> use binning + robust-sigma test per replicate")
    
    # 3) Loop over default bin counts and compute p-values and FDR
    default_bin_counts = [1, 5, 10, 15, 20, 25, 30, 35, 40]
    print(f"Default bin counts for significance testing: {default_bin_counts}")
    
    sig_site_counts = []  # to store number of significant Site Groups for each bin count
    
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
            f"Bin count {bin_count}: {sig_count} significant Site Groups "
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
        plt.ylabel("Number of significant N-glycosite Site Groups")
        plt.title("Significant Site Group count vs bin count")
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
        "Enter custom bin size (number of N-glycosite Site Groups per window) for significance testing, e.g., 50: "
    ).strip()
    try:
        custom_bin_size = int(raw_bin_size_input)
        if custom_bin_size <= 0:
            raise ValueError("Bin size must be positive.")
    except Exception as e:
        raise ValueError(f"Invalid bin size: {e}")
    
    print(f"Using custom bin size = {custom_bin_size} Site Groups per window")
    
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
    custom_output_file = os.path.join(tables_dir, f"14_nglycosite_site_groups_significance_normal_binSize{custom_bin_size}.csv")
    export_csv(site_df_final, custom_output_file, SITE_EXPORT_RENAME)
    n_sig = sig_mask_custom.sum()
    print(f"Consensus rule: same positive direction + mean ΔRm>{delta_rm_threshold:.2f} + all significance values<0.10 + at least one<0.05")
    print(f"Custom bin size {custom_bin_size}: {n_sig} significant Site Groups saved to {custom_output_file}")

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
    print("Step10 strategy: replicate-level normality gate failed -> paired one-tailed test (normalized matched Protein-Group Rm vs normalized Site-Group Rm) per Site Group")

    site_df_non_normal = site_df_CV_pass.copy().reset_index(drop=True)
    needed_cols = [
        "Rm_pg1_norm_ref", "Rm_pg2_norm_ref", "Rm_pg3_norm_ref",
        "Rm_sg1_norm", "Rm_sg2_norm", "Rm_sg3_norm", "ΔRm"
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

    print(f"Paired test: {non_normal_method} on normalized Rm (one-tailed H1: Site-Group Rm > matched Protein-Group Rm)")

    pvals = []
    n_pairs_used = []
    for _, row in site_df_non_normal.iterrows():
        protein_rm_norm = np.array([row["Rm_pg1_norm_ref"], row["Rm_pg2_norm_ref"], row["Rm_pg3_norm_ref"]], dtype=float)
        glycosite_rm_norm = np.array([row["Rm_sg1_norm"], row["Rm_sg2_norm"], row["Rm_sg3_norm"]], dtype=float)

        p, n_pairs = paired_one_tailed_pvalue(
            protein_group_vals=protein_rm_norm,
            site_group_vals=glycosite_rm_norm,
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
        non_normal_output_file = os.path.join(tables_dir, "14_nglycosite_site_groups_significance_paired_wilcoxon.csv")
    else:
        non_normal_output_file = os.path.join(tables_dir, "14_nglycosite_site_groups_significance_paired_t.csv")
    export_csv(site_df_final, non_normal_output_file, SITE_EXPORT_RENAME)
    print(f"Non-normal branch results saved: {non_normal_output_file}")
    print(f"Significant up-stabilized Site Groups: {site_df_final['Significant'].sum()}")

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
            protein_rm_norm = np.array([row["Rm_pg1_norm_ref"], row["Rm_pg2_norm_ref"], row["Rm_pg3_norm_ref"]], dtype=float)
            glycosite_rm_norm = np.array([row["Rm_sg1_norm"], row["Rm_sg2_norm"], row["Rm_sg3_norm"]], dtype=float)

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

            site_label = str(row.get("candidate_sites", row.get("site_group_id", "site_group")))
            if len(site_label) > 42:
                site_label = site_label[:39] + "..."
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
    print("Non-normal volcano uses mean ΔRm across biological replicates (same DeltaRm definition as the Site Group table).")
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
    "protein_quantitative_unit": "PEAKS_Protein_Group",
    "protein_member_identity_unit": "UniProt_accession",
    "site_group_definition": "one_protein_group_plus_one_candidate_site_assignment_set",
    "shared_peptide_policy": "collapse_identical_quantitative_evidence_within_protein_group",
    "candidate_annotation_policy": "retain_all_valid_accession_site_candidates",
    "normalization_weighting": "one_vote_per_unique_protein_group",
    "normalization_population": "unique_control_protein_groups_after_qc",
    "multiple_testing_unit": "SiteGroupID",
    "candidate_annotations_are_independent_tests": False,
    "method_formulas": METHOD_FORMULAS,
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
    "delta_rm_definition": "mean_i(site_group_Rm_norm_i - matched_protein_group_Rm_norm_i)",
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
