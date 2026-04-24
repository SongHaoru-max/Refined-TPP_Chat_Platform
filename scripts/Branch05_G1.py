# -*- coding: utf-8 -*-
"""
Branch05_G1 - Multi-BioRep paired Site/Protein pipeline
--------------------------------------------------------
Target scenario:
- Biological replicates >= 3
- Each biological replicate contains 3 technical repeats in TMT-6plex
- Input is paired CSV files per BioRep: one site-group CSV + one protein-group CSV

Core logic:
0) FASTA-based site retrieval for each BioRep site CSV (peptide-level -> site-level)
1) Protein quantification overlap plots (per BioRep, all conditions, quantified-any rule)
2) Strict completeness gate across BioReps for downstream analysis:
   - site group: selected condition must be fully quantified in all 6 channels
   - protein group: control condition must be fully quantified in all 6 channels
3) Rm workflow:
   - raw Rm1-3 from TMT channel pairs (129/127, 130/128, 131/126)
   - geometric mean within each BioRep
   - median-of-medians normalization in log space across BioReps
   - DeltaRm = adjusted_site_Rm - adjusted_protein_control_Rm
4) Two-layer CV-QC:
   - technical CV gate (default 0.30): site CV and protein CV both pass in every BioRep
   - biological CV gate (default 0.40): pooled raw-Rm CV across BioReps for site and protein both pass
   - CV distribution plots generated before final threshold input
5) Branching differential analysis:
   - If DeltaRm distributions are normal in all BioReps:
     robust sigma by CV bin per BioRep; consensus significance rule:
     same DeltaRm sign across BioReps + abs(DeltaRm)>0.1 in all BioReps +
     all replicate FDR<0.1 + at least one replicate FDR<0.05
     output one volcano per BioRep
   - Otherwise:
     paired one-tailed test across BioReps on adjusted Rm (paired t-test or Wilcoxon),
     BH correction, significance requires FDR<0.05 + same sign + abs(DeltaRm)>0.1 in all BioReps
     output one volcano using mean DeltaRm and TOP25 boxplots
"""

import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
from scipy import stats
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from upsetplot import UpSet, from_contents


plt.rcParams.update({"figure.dpi": 300})

DEFAULT_TECH_CV_CUTOFF = 0.30
DEFAULT_BIO_CV_CUTOFF = 0.40
TECH_CV_PREVIEW_CANDIDATES = (0.25, 0.30, 0.35, 0.40)
BIO_CV_PREVIEW_CANDIDATES = (0.30, 0.35, 0.40, 0.45, 0.50)


def normalize_path(path):
    return os.path.normpath(os.path.abspath(path))


def canonical_protein_id(value):
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""

    text = text.split(";")[0].strip()
    if "|" in text:
        parts = [p.strip() for p in text.split("|") if p.strip()]
        if len(parts) >= 2 and parts[0].lower() in {"sp", "tr", "up", "ref", "gb", "emb", "dbj", "gi"}:
            text = parts[1]
        elif parts:
            text = parts[0]

    text = text.split()[0].strip()
    return text


def make_safe_key(text):
    key = re.sub(r"[^0-9A-Za-z]+", "_", str(text).strip())
    key = key.strip("_")
    return key if key else "Condition"


def detect_sample_blocks(df):
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


def select_tmt_quant_cols(df):
    # Branch06/07/08-aligned behavior: prefer strict Intensity-Sample-TMT6 columns.
    strict_intensity_pat = re.compile(r"Intensity\s+Sample\s*\d+\s+TMT6-(126|127|128|129|130|131)", flags=re.I)
    strict_intensity_cols = [c for c in df.columns if strict_intensity_pat.search(str(c))]
    if strict_intensity_cols:
        return strict_intensity_cols

    # Fallback for variant headers: still exclude ratio-like columns.
    reporter_pat = re.compile(r"(126|127|128|129|130|131)")
    sample_pat = re.compile(r"Sample\s*\d+", flags=re.I)

    def is_ratio_like(colname):
        low = str(colname).lower()
        return ("ratio" in low) or ("fold" in low and "change" in low)

    base = [
        c
        for c in df.columns
        if sample_pat.search(str(c)) and reporter_pat.search(str(c)) and (not is_ratio_like(c))
    ]
    intensity_cols = [c for c in base if "intensity" in str(c).lower()]
    return intensity_cols if intensity_cols else base


def extract_reporter(colname):
    m = re.search(r"(126|127|128|129|130|131)", str(colname))
    if m:
        return int(m.group(1))
    return None


def detect_channel_map(sample_cols):
    channel_map = {}
    for col in sample_cols:
        ch = extract_reporter(col)
        if ch is not None and ch not in channel_map:
            channel_map[ch] = col
    required = [126, 127, 128, 129, 130, 131]
    missing = [c for c in required if c not in channel_map]
    if missing:
        raise ValueError(f"Missing TMT channels {missing} from sample block columns: {sample_cols}")
    return channel_map


def detect_unique_col(df):
    if "#Unique" in df.columns:
        return "#Unique"
    candidates = [
        c
        for c in df.columns
        if re.search(r"#\s*Unique|#Unique|Unique\s*Pept|UniquePept|Unique", c, flags=re.I)
    ]
    return candidates[0] if candidates else None


def detect_protein_col(df):
    # Keep consistent with Branch06 matching: use accession-like identifiers first.
    preferred = ["Accession", "Accession ", "Protein Group", "ProteinGroup", "Protein_Group"]
    for c in preferred:
        if c in df.columns:
            return c
    return df.columns[0]


def detect_site_columns(df):
    prot_candidates = ["Protein Accession", "Accession", "Protein Group", "ProteinGroup", "Protein_Group"]
    nsite_candidates = ["N-Glycosite", "N-Glycosite in Protein", "N-glycosite", "N_site", "Site"]

    prot_col = next((c for c in prot_candidates if c in df.columns), None)
    nsite_col = next((c for c in nsite_candidates if c in df.columns), None)

    if prot_col is None:
        raise ValueError("Cannot detect protein column in site CSV.")

    return prot_col, nsite_col


def load_fasta_sequences(fasta_file):
    seq_dict = {}
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(fasta_file)

    accession = ""
    seq_lines = []
    with open(fasta_file, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession and seq_lines:
                    seq_dict[accession] = "".join(seq_lines)
                header = line[1:].strip()
                first = header.split()[0] if header else ""
                accession = canonical_protein_id(first)
                seq_lines = []
            else:
                seq_lines.append(line)

    if accession and seq_lines:
        seq_dict[accession] = "".join(seq_lines)

    if not seq_dict:
        raise ValueError("No sequences parsed from FASTA file.")

    return seq_dict


def parse_peptide_for_nglyco(peptide):
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
                    break
                mods.append(peptide[j + 1 : k])
                j = k + 1
            if "+0.98" in mods:
                n_positions.append(len(bare_seq))
            bare_seq += "N"
            i = j
        elif peptide[i] == "(":
            k = peptide.find(")", i)
            i = k + 1 if k != -1 else i + 1
        else:
            bare_seq += peptide[i]
            i += 1
    return bare_seq, n_positions


def build_sites_from_peptide_table(peptide_df, fasta_dict, motif_len=9):
    if "Peptide" not in peptide_df.columns:
        raise ValueError("Site CSV lacks 'Peptide' column for FASTA-based site retrieval.")

    accession_col = None
    for c in ["Protein Accession", "Accession", "Protein Group", "ProteinGroup", "Protein_Group"]:
        if c in peptide_df.columns:
            accession_col = c
            break
    if accession_col is None:
        raise ValueError("Site CSV lacks accession column for FASTA mapping.")

    quant_cols = select_tmt_quant_cols(peptide_df)
    if not quant_cols:
        raise ValueError("No quant columns (Sample N + reporter channel) found in site CSV.")

    half_len = int(motif_len // 2)
    records = []

    for _, row in peptide_df.iterrows():
        protein_id = canonical_protein_id(row.get(accession_col, np.nan))
        if not protein_id:
            continue

        prot_seq = fasta_dict.get(protein_id, "")
        if not prot_seq:
            continue

        pep = str(row.get("Peptide", "")).strip()
        if not pep:
            continue

        stripped_seq, n_pos_list = parse_peptide_for_nglyco(pep)
        if not n_pos_list or not stripped_seq:
            continue

        start_idx = prot_seq.find(stripped_seq)
        if start_idx == -1:
            continue

        for n_pos in n_pos_list:
            prot_n_pos = start_idx + n_pos
            if prot_n_pos + 2 >= len(prot_seq):
                continue

            x = prot_seq[prot_n_pos + 1]
            st = prot_seq[prot_n_pos + 2]
            if x == "P" or st not in {"S", "T"}:
                continue

            center = prot_n_pos
            motif_chars = []
            for k in range(center - half_len, center + half_len + 1):
                motif_chars.append(prot_seq[k] if 0 <= k < len(prot_seq) else "-")
            motif = "".join(motif_chars)

            rec = {
                "ProteinID": protein_id,
                "Protein Accession": protein_id,
                "N-Glycosite": int(prot_n_pos + 1),
                "Motif": motif,
                "SiteID": f"{protein_id}|N{int(prot_n_pos + 1)}",
                "Peptide": pep,
            }
            for c in quant_cols:
                rec[c] = row.get(c, np.nan)
            records.append(rec)

    if not records:
        raise ValueError("No N-glycosites retrieved from peptide table with current FASTA mapping.")

    retrieved_peptide_df = pd.DataFrame(records)
    site_level_df = (
        retrieved_peptide_df.groupby(["ProteinID", "Protein Accession", "N-Glycosite", "Motif", "SiteID"], as_index=False)[quant_cols]
        .sum(min_count=1)
    )
    return retrieved_peptide_df, site_level_df


def build_site_id(df, prot_col, nsite_col):
    prot = df[prot_col].map(canonical_protein_id)
    if nsite_col is not None:
        nsite_text = df[nsite_col].astype(str).str.strip()
        return prot + "|N" + nsite_text

    if "Peptide" in df.columns:
        pep = df["Peptide"].astype(str).str.strip()
        return prot + "|" + pep

    return prot + "|ROW" + df.index.astype(str)


def geometric_mean_safe(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0 or np.any(arr <= 0):
        return np.nan
    return float(np.exp(np.mean(np.log(arr))))


def compute_sigma(vals):
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


def compute_bin_pvals(vals, side="right"):
    vals = np.asarray(vals, dtype=float)
    pvals = np.ones(vals.shape[0], dtype=float)

    valid_vals = vals[~np.isnan(vals)]
    if valid_vals.size == 0:
        return pvals

    p50, left_sigma, right_sigma = compute_sigma(valid_vals)

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


def robust_bin_pvals(deltas, cv_sort_order, bin_size=300):
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
                    pv = 1 - norm.cdf((val - p50) / right_sigma)
            else:
                if np.isnan(left_sigma) or left_sigma == 0:
                    pv = 1.0
                else:
                    pv = 1 - norm.cdf((val - p50) / left_sigma)
            pval_map[idx] = pv

    aligned = pd.Series(1.0, index=deltas.index)
    for idx, pv in pval_map.items():
        aligned.loc[idx] = pv
    return aligned


def paired_one_tailed_pvalue(bulk_vals, glyco_vals, method="paired_t", side="right"):
    bulk_vals = np.asarray(bulk_vals, dtype=float)
    glyco_vals = np.asarray(glyco_vals, dtype=float)

    valid_mask = (~np.isnan(bulk_vals)) & (~np.isnan(glyco_vals))
    bulk_valid = bulk_vals[valid_mask]
    glyco_valid = glyco_vals[valid_mask]
    n_pairs = len(bulk_valid)

    if n_pairs < 2:
        return 1.0, n_pairs

    diffs = glyco_valid - bulk_valid
    if np.allclose(diffs, 0):
        return 1.0, n_pairs

    if method == "paired_t":
        t_stat, p_two = stats.ttest_rel(glyco_valid, bulk_valid, nan_policy="omit")
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
                glyco_valid,
                bulk_valid,
                alternative=alt,
                zero_method="wilcox",
                method="auto",
            )
        except TypeError:
            _, p_one = stats.wilcoxon(
                glyco_valid,
                bulk_valid,
                alternative=alt,
                zero_method="wilcox",
                mode="auto",
            )
        except Exception:
            p_one = 1.0
        return float(np.clip(p_one, 0.0, 1.0)), n_pairs

    raise ValueError("method must be 'paired_t' or 'wilcoxon'")


def assess_distribution_normality(values, alpha=0.05):
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    n = int(vals.size)

    if n < 8:
        return False, f"n={n} (<8): insufficient data, treated as NOT normal", "insufficient_n", n

    if n <= 50:
        stat, pval = stats.shapiro(vals)
        is_normal = bool(pval >= alpha)
        summary = f"Shapiro-Wilk: W={stat:.4f}, p={pval:.4g}" + (
            " -> approximately normal" if is_normal else " -> NOT normal"
        )
        return is_normal, summary, "shapiro", n

    if n <= 300:
        mean_val = np.mean(vals)
        std_val = np.std(vals, ddof=1)
        if (not np.isfinite(std_val)) or (std_val <= 0):
            return False, f"Kolmogorov-Smirnov skipped (std={std_val:.4g}), treated as NOT normal", "ks_invalid_std", n
        stat, pval = stats.kstest(vals, "norm", args=(mean_val, std_val))
        is_normal = bool(pval >= alpha)
        summary = f"Kolmogorov-Smirnov: D={stat:.4f}, p={pval:.4g}" + (
            " -> approximately normal" if is_normal else " -> NOT normal"
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


def pval_to_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def draw_quant_overlap(per_sample_sets, sample_names, out_prefix):
    n = len(sample_names)
    if n < 2:
        return None

    out_file = None
    if n < 4:
        plt.figure(figsize=(6, 6))
        if n == 2:
            venn2([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
        elif n == 3:
            venn3([per_sample_sets[s] for s in sample_names], set_labels=sample_names)
        plt.title("Quantified proteins overlap")
        out_file = f"{out_prefix}_quantified_venn.tiff"
        plt.savefig(out_file, dpi=300, bbox_inches="tight", format="tiff")
        plt.close()
    else:
        upset_data = from_contents(per_sample_sets)
        plt.figure(figsize=(8, 6))
        UpSet(upset_data, show_counts=True).plot()
        out_file = f"{out_prefix}_quantified_upset.tiff"
        plt.savefig(out_file, dpi=300, bbox_inches="tight", format="tiff")
        plt.close()

    return out_file


def draw_site_retrieval_reports(rep_name, raw_site_df, retrieved_peptide_df, retrieved_site_df, quant_cols, out_prefix):
    if quant_cols:
        # Coerce to numeric and align by quant_cols to avoid length mismatch in plotting.
        raw_cols = [c for c in quant_cols if c in raw_site_df.columns]
        glyco_cols = [c for c in quant_cols if c in retrieved_peptide_df.columns]

        raw_num = raw_site_df[raw_cols].apply(pd.to_numeric, errors="coerce") if raw_cols else pd.DataFrame(index=raw_site_df.index)
        glyco_num = (
            retrieved_peptide_df[glyco_cols].apply(pd.to_numeric, errors="coerce")
            if glyco_cols
            else pd.DataFrame(index=retrieved_peptide_df.index)
        )

        total_intensity = raw_num.sum(axis=0).reindex(quant_cols, fill_value=0.0).values.astype(float)
        glyco_intensity = glyco_num.sum(axis=0).reindex(quant_cols, fill_value=0.0).values.astype(float)
        non_glyco_intensity = np.maximum(total_intensity - glyco_intensity, 0.0)
        ratio = np.divide(glyco_intensity, total_intensity, out=np.zeros_like(glyco_intensity), where=total_intensity > 0)

        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(quant_cols))
        ax.bar(x, glyco_intensity, width=0.6, color="#E63946", label="N-glycopeptides")
        ax.bar(x, non_glyco_intensity, bottom=glyco_intensity, width=0.6, color="#457B9D", label="Other peptides")
        for i, r in enumerate(ratio):
            ax.text(x[i], glyco_intensity[i] / 2.0 if glyco_intensity[i] > 0 else 0.0, f"{r:.2f}",
                    ha="center", va="center", color="white", fontsize=8, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(quant_cols, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("Total Intensity")
        ax.set_title(f"{rep_name} - Intensity composition of N-glycopeptides")
        ax.legend(frameon=False)
        plt.tight_layout()
        intensity_plot = f"{out_prefix}_intensity_composition.tiff"
        plt.savefig(intensity_plot, dpi=300, format="tiff")
        plt.close(fig)

    site_count_series = retrieved_site_df.groupby("ProteinID")["N-Glycosite"].nunique()
    capped = site_count_series.apply(lambda x: x if x <= 5 else ">5")
    dist = capped.value_counts().sort_index(key=lambda v: [int(i) if i != ">5" else 6 for i in v])

    if not dist.empty:
        plt.figure(figsize=(8, 6))
        bars = plt.barh(dist.index.astype(str), dist.values, color="#8CB9D0")
        for bar in bars:
            width = bar.get_width()
            plt.text(width + max(dist.values) * 0.01, bar.get_y() + bar.get_height() / 2.0, str(int(width)), va="center")
        plt.xlabel("Protein count")
        plt.ylabel("Number of N-glycosites per protein")
        plt.title(f"{rep_name} - N-glycosite count distribution per protein")
        plt.tight_layout()
        dist_plot = f"{out_prefix}_protein_nsite_distribution.tiff"
        plt.savefig(dist_plot, dpi=300, format="tiff")
        plt.close()

    summary = pd.DataFrame([
        {
            "BioRep": rep_name,
            "Retrieved_glycopeptide_rows": int(retrieved_peptide_df.shape[0]),
            "Retrieved_site_rows": int(retrieved_site_df.shape[0]),
            "Retrieved_protein_count": int(retrieved_site_df["ProteinID"].nunique()),
            "Total_unique_sites": int(retrieved_site_df["N-Glycosite"].count()),
        }
    ])
    summary_file = f"{out_prefix}_retrieval_summary.csv"
    summary.to_csv(summary_file, index=False)


def draw_biorep_cv_reports(
    rep_tag,
    out_dir,
    site_cv_complete_target,
    site_cv_intersect,
    protein_cv_complete_control,
    protein_cv_intersect,
    scatter_site_cv=None,
    scatter_protein_cv=None,
    target_label="N-Glyco",
    control_label="Control",
    tech_cutoff=0.30,
):
    cv_panels = [
        {
            "series": site_cv_complete_target,
            "color": "#5E556A",
            "xlabel": "Sites sorted by CV (high->low)",
            "title": (
                f"{rep_tag} - Site CV distribution\n"
                f"{target_label} (All fully quantified sites)"
            ),
            "out_name": f"Branch05_G1_{rep_tag}_CV_sites_complete_target.tiff",
        },
        {
            "series": site_cv_intersect,
            "color": "#5E556A",
            "xlabel": "Sites sorted by CV (high->low)",
            "title": (
                f"{rep_tag} - Site CV distribution\n"
                f"{target_label} (Fully quantified sites of fully quantified {control_label} proteins)"
            ),
            "out_name": f"Branch05_G1_{rep_tag}_CV_sites_target_on_complete_protein.tiff",
        },
        {
            "series": protein_cv_complete_control,
            "color": "#2A9D8F",
            "xlabel": f"Fully quantified {control_label} proteins sorted by CV (high->low)",
            "title": (
                f"{rep_tag} - Protein CV distribution\n"
                f"{control_label} Proteins (All fully quantified proteins)"
            ),
            "out_name": f"Branch05_G1_{rep_tag}_CV_proteins_complete_control.tiff",
        },
        {
            "series": protein_cv_intersect,
            "color": "#2A9D8F",
            "xlabel": f"Fully quantified {control_label} proteins sorted by CV (high->low)",
            "title": (
                f"{rep_tag} - Protein CV distribution\n"
                f"{control_label} Proteins (Fully quantified proteins with fully quantified {target_label} sites)"
            ),
            "out_name": f"Branch05_G1_{rep_tag}_CV_proteins_complete_control_with_complete_sites.tiff",
        },
    ]

    for panel in cv_panels:
        series = pd.Series(panel["series"]).replace([np.inf, -np.inf], np.nan).dropna()
        if series.empty:
            continue

        plot_series = series.sort_values(ascending=False).reset_index(drop=True)
        rank = plot_series.index + 1
        prop_30 = float((plot_series <= 0.30).sum() / len(plot_series) * 100.0)
        prop_40 = float((plot_series <= 0.40).sum() / len(plot_series) * 100.0)
        prop_cut = float((plot_series <= tech_cutoff).sum() / len(plot_series) * 100.0)
        ymax = max(float(plot_series.max()) * 1.05, max(0.25, tech_cutoff + 0.08))

        plt.figure(figsize=(6, 4))
        plt.scatter(rank, plot_series.values, color=panel["color"], alpha=0.7, s=18)
        plt.text(0.62 * len(plot_series), min(0.305, ymax * 0.95), f"<=30%: {prop_30:.1f}%", fontsize=9)
        plt.text(0.62 * len(plot_series), min(0.405, ymax * 0.95), f"<=40%: {prop_40:.1f}%", fontsize=9)
        if (not np.isclose(tech_cutoff, 0.30)) and (not np.isclose(tech_cutoff, 0.40)):
            plt.text(0.62 * len(plot_series), min(tech_cutoff + 0.005, ymax * 0.95), f"<={tech_cutoff*100:.0f}%: {prop_cut:.1f}%", fontsize=9)
        plt.axhline(tech_cutoff, linestyle="--", linewidth=1, color="#111827")
        plt.ylim(0, ymax)
        plt.xlabel(panel["xlabel"])
        plt.ylabel("CV")
        plt.title(panel["title"])
        plt.tight_layout()
        out_file = os.path.join(out_dir, panel["out_name"])
        plt.savefig(out_file, dpi=300, format="tiff")
        plt.close()
        print(f"{rep_tag} CV distribution saved: {out_file}")

    if scatter_site_cv is not None and scatter_protein_cv is not None:
        scatter_df = pd.DataFrame(
            {
                "CV_site_tech": pd.Series(scatter_site_cv),
                "CV_ctrl_tech": pd.Series(scatter_protein_cv),
            }
        ).replace([np.inf, -np.inf], np.nan).dropna()
        if not scatter_df.empty:
            plt.figure(figsize=(6, 6))
            plt.scatter(scatter_df["CV_site_tech"], scatter_df["CV_ctrl_tech"], s=10, alpha=0.6)
            plt.axvline(tech_cutoff, linestyle="--", linewidth=1, label=f"site CV cutoff ({tech_cutoff:.2f})")
            plt.axhline(tech_cutoff, linestyle="--", linewidth=1, label=f"protein CV cutoff ({tech_cutoff:.2f})")
            plt.xlabel("Site CV (raw Rm technical repeats)")
            plt.ylabel("Protein CV (raw Rm technical repeats)")
            plt.title(f"{rep_tag} - Site CV vs Protein CV")
            plt.legend(frameon=False, fontsize=8)
            plt.tight_layout()
            out_file = os.path.join(out_dir, f"Branch05_G1_{rep_tag}_Site_vs_Protein_CV_scatter.tiff")
            plt.savefig(out_file, dpi=300, format="tiff")
            plt.close()
            print(f"{rep_tag} site-vs-protein CV scatter saved: {out_file}")


def plot_protein_site_distribution(site_df, title, out_file, color="#14C0CC"):
    if site_df is None or site_df.empty:
        return None

    prot_site_count = site_df.groupby("ProteinID")["SiteID"].nunique()
    if prot_site_count.empty:
        return None

    capped = prot_site_count.apply(lambda x: x if x <= 5 else ">5")
    dist = capped.value_counts().sort_index(key=lambda v: [int(i) if i != ">5" else 6 for i in v])
    if dist.empty:
        return None

    plt.figure(figsize=(8, 6))
    bars = plt.barh(dist.index.astype(str), dist.values, color=color)
    for bar in bars:
        width = bar.get_width()
        plt.text(width + max(dist.values) * 0.01, bar.get_y() + bar.get_height() / 2.0, str(int(width)), va="center")
    plt.xlabel("Protein count")
    plt.ylabel("Number of N-glycosites per protein")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, format="tiff")
    plt.close()
    return out_file


def plot_cv_distribution(series, title, out_file, threshold_lines=None):
    vals = pd.Series(series).replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False).reset_index(drop=True)
    if vals.empty:
        return None

    rank = vals.index + 1
    plt.figure(figsize=(6, 4))
    plt.scatter(rank, vals.values, color="#4C6A92", alpha=0.72, s=18)

    ymax = max(float(vals.max()) * 1.05, 0.5)
    if threshold_lines:
        for t in threshold_lines:
            plt.axhline(float(t), linestyle="--", linewidth=1)
            prop = float((vals <= float(t)).sum() / len(vals) * 100.0)
            plt.text(0.62 * len(vals), min(float(t) + 0.005, ymax * 0.95), f"<={float(t)*100:.0f}%: {prop:.1f}%", fontsize=9)

    plt.ylim(0, ymax)
    plt.xlabel("Features sorted by CV (high->low)")
    plt.ylabel("CV")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, format="tiff")
    plt.close()
    return out_file


def plot_delta_normality_diagnostics(delta_series, sample_name, out_prefix):
    arr = pd.Series(delta_series).replace([np.inf, -np.inf], np.nan).dropna().values
    if arr.size < 3:
        return []

    safe_name = re.sub(r"[^0-9A-Za-z_-]+", "_", str(sample_name))
    out_files = []

    plt.figure(figsize=(6, 4))
    sns.histplot(arr, kde=True, bins=30, color="#A7AED2")
    plt.xlabel("DeltaRm")
    plt.ylabel("Count")
    plt.title(f"DeltaRm distribution with KDE - {sample_name}")
    plt.tight_layout()
    hist_file = f"{out_prefix}_DeltaRm_hist_kde_{safe_name}.tiff"
    plt.savefig(hist_file, dpi=300, format="tiff")
    plt.close()
    out_files.append(hist_file)

    plt.figure(figsize=(6, 4))
    sm.qqplot(arr, line="s")
    plt.title(f"DeltaRm Q-Q plot - {sample_name}")
    plt.tight_layout()
    qq_file = f"{out_prefix}_DeltaRm_QQ_{safe_name}.tiff"
    plt.savefig(qq_file, dpi=300, format="tiff")
    plt.close()
    out_files.append(qq_file)

    plt.figure(figsize=(6, 4))
    sm.ProbPlot(arr).ppplot(line="45")
    plt.title(f"DeltaRm P-P plot - {sample_name}")
    plt.tight_layout()
    pp_file = f"{out_prefix}_DeltaRm_PP_{safe_name}.tiff"
    plt.savefig(pp_file, dpi=300, format="tiff")
    plt.close()
    out_files.append(pp_file)

    return out_files


def plot_volcano(delta_series, fdr_series, title, out_file, delta_cut=0.1, fdr_cut=0.05):
    plot_df = pd.DataFrame({"DeltaRm": pd.Series(delta_series), "FDR": pd.Series(fdr_series)})
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna()
    if plot_df.empty:
        return None

    plot_df["neglog10FDR"] = -np.log10(np.clip(plot_df["FDR"].values, 1e-300, None))
    up_mask = (plot_df["DeltaRm"] > delta_cut) & (plot_df["FDR"] < fdr_cut)
    down_mask = (plot_df["DeltaRm"] < -delta_cut) & (plot_df["FDR"] < fdr_cut)
    other_mask = ~(up_mask | down_mask)

    plt.figure(figsize=(6, 5))
    plt.scatter(plot_df.loc[other_mask, "DeltaRm"], plot_df.loc[other_mask, "neglog10FDR"], s=14, color="grey", alpha=0.6)
    plt.scatter(plot_df.loc[up_mask, "DeltaRm"], plot_df.loc[up_mask, "neglog10FDR"], s=18, color="red", alpha=0.85, label="Up & significant")
    plt.scatter(plot_df.loc[down_mask, "DeltaRm"], plot_df.loc[down_mask, "neglog10FDR"], s=18, color="blue", alpha=0.85, label="Down & significant")

    plt.axvline(delta_cut, color="red", linestyle="--", linewidth=1)
    plt.axvline(-delta_cut, color="blue", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(fdr_cut), color="black", linestyle="--", linewidth=1)

    plt.xlabel("DeltaRm")
    plt.ylabel("-log10(FDR)")
    plt.title(title)
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, format="tiff")
    plt.close()
    return out_file


def ask_float_with_default(prompt, default, min_value=0.0):
    while True:
        raw = input(prompt).strip()
        if raw == "":
            return float(default)
        try:
            val = float(raw)
        except ValueError:
            print(f"Invalid numeric input: {raw}")
            continue
        if not np.isfinite(val) or val <= min_value:
            print(f"Input must be finite and > {min_value}")
            continue
        return float(val)


def prepare_replicate(
    rep_name,
    protein_path,
    site_path,
    protein_condition_defs,
    site_condition_defs,
    use_unique_ge2,
    fasta_dict,
    motif_len,
):
    protein_df = pd.read_csv(protein_path)
    site_df_raw = pd.read_csv(site_path)

    # Mandatory site retrieval and analysis from peptide-level table using FASTA.
    retrieved_peptide_df, site_df = build_sites_from_peptide_table(site_df_raw, fasta_dict, motif_len=motif_len)

    site_base = os.path.splitext(os.path.basename(site_path))[0]
    site_dir = os.path.dirname(site_path)
    retrieved_pep_file = os.path.join(site_dir, f"{site_base}_{rep_name}_retrieved_nglycopeptides.csv")
    retrieved_site_file = os.path.join(site_dir, f"{site_base}_{rep_name}_retrieved_sites.csv")
    retrieved_peptide_df.to_csv(retrieved_pep_file, index=False)
    site_df.to_csv(retrieved_site_file, index=False)
    print(f"{rep_name} FASTA retrieval peptide table saved: {retrieved_pep_file}")
    print(f"{rep_name} FASTA retrieval site table saved: {retrieved_site_file}")

    quant_cols = select_tmt_quant_cols(site_df_raw)
    draw_site_retrieval_reports(
        rep_name=rep_name,
        raw_site_df=site_df_raw,
        retrieved_peptide_df=retrieved_peptide_df,
        retrieved_site_df=site_df,
        quant_cols=quant_cols,
        out_prefix=os.path.join(site_dir, f"{site_base}_{rep_name}"),
    )

    p_prot_col = detect_protein_col(protein_df)
    p_sample_order, p_sample_cols = detect_sample_blocks(protein_df)
    if not p_sample_order:
        raise ValueError(f"No 'Sample N' block detected in protein CSV: {protein_path}")

    s_sample_order, s_sample_cols = detect_sample_blocks(site_df)
    if not s_sample_order:
        raise ValueError(f"No 'Sample N' block detected in site CSV: {site_path}")

    if protein_condition_defs is None:
        protein_condition_defs = []
        if len(p_sample_order) == 1:
            raw_name = input(
                f"{rep_name} - Enter condition name for protein Sample {p_sample_order[0]} (default=control): "
            ).strip()
            cond_name = raw_name if raw_name else "control"
            protein_condition_defs.append({"name": cond_name, "key": make_safe_key(cond_name)})
        else:
            for s in p_sample_order:
                cond_name = input(f"{rep_name} - Enter condition name for protein Sample {s}: ").strip()
                if not cond_name:
                    raise ValueError("Condition name cannot be empty.")
                protein_condition_defs.append({"name": cond_name, "key": make_safe_key(cond_name)})

    if site_condition_defs is None:
        site_condition_defs = []
        for s in s_sample_order:
            cond_name = input(f"{rep_name} - Enter condition name for site Sample {s}: ").strip()
            if not cond_name:
                raise ValueError("Condition name cannot be empty.")
            site_condition_defs.append({"name": cond_name, "key": make_safe_key(cond_name)})

    if len(protein_condition_defs) != len(p_sample_order):
        raise ValueError(f"{rep_name}: protein sample count mismatch against BR1 protein condition definitions.")
    if len(site_condition_defs) != len(s_sample_order):
        raise ValueError(f"{rep_name}: site sample count mismatch against BR1 site condition definitions.")

    p_cond_to_cols = {}
    s_cond_to_cols = {}
    for i, cond in enumerate(protein_condition_defs):
        p_cond_to_cols[cond["key"]] = p_sample_cols[p_sample_order[i]]
    for i, cond in enumerate(site_condition_defs):
        s_cond_to_cols[cond["key"]] = s_sample_cols[s_sample_order[i]]

    unique_col = detect_unique_col(protein_df)
    if unique_col is not None:
        unique_mask = protein_df[unique_col] >= (2 if use_unique_ge2 else 1)
    else:
        unique_mask = pd.Series(True, index=protein_df.index)

    quantified_sets = {}
    for cond in protein_condition_defs:
        key = cond["key"]
        ccols = p_cond_to_cols[key]
        quantified = (protein_df[ccols] != 0).any(axis=1) & unique_mask
        protein_ids = protein_df.loc[quantified, p_prot_col].map(canonical_protein_id)
        quantified_sets[cond["name"]] = set(protein_ids[protein_ids != ""])

    p_id = protein_df[p_prot_col].map(canonical_protein_id)
    protein_df = protein_df.copy()
    protein_df["ProteinID"] = p_id
    protein_df["Unique_pass"] = unique_mask.values
    # Align with Branch06 behavior: unique filter is applied before downstream completeness.
    protein_df = protein_df[(protein_df["ProteinID"] != "") & (protein_df["Unique_pass"])].copy()

    site_prot_col, site_nsite_col = detect_site_columns(site_df)
    site_df = site_df.copy()
    site_df["ProteinID"] = site_df[site_prot_col].map(canonical_protein_id)
    site_df["SiteID"] = build_site_id(site_df, site_prot_col, site_nsite_col)

    return {
        "protein_df": protein_df,
        "site_df": site_df,
        "protein_condition_defs": protein_condition_defs,
        "site_condition_defs": site_condition_defs,
        "protein_cond_cols": p_cond_to_cols,
        "site_cond_cols": s_cond_to_cols,
        "quantified_sets": quantified_sets,
        "protein_path": protein_path,
        "site_path": site_path,
        "retrieved_peptide_file": retrieved_pep_file,
        "retrieved_site_file": retrieved_site_file,
    }


def calc_rm_triplet_from_cols(df, cols):
    cmap = detect_channel_map(cols)
    rm1 = df[cmap[129]] / df[cmap[127]]
    rm2 = df[cmap[130]] / df[cmap[128]]
    rm3 = df[cmap[131]] / df[cmap[126]]
    return rm1, rm2, rm3


def compute_replicate_site_metrics(
    rep_obj,
    control_key,
    target_key,
    rep_name=None,
    control_label="Control",
    target_label="N-Glyco",
):
    p_df = rep_obj["protein_df"].copy()
    s_df = rep_obj["site_df"].copy()

    p_rm1, p_rm2, p_rm3 = calc_rm_triplet_from_cols(p_df, rep_obj["protein_cond_cols"][control_key])
    p_df["Rm_ctrl1"] = p_rm1
    p_df["Rm_ctrl2"] = p_rm2
    p_df["Rm_ctrl3"] = p_rm3
    p_df["CV_ctrl_tech"] = p_df[["Rm_ctrl1", "Rm_ctrl2", "Rm_ctrl3"]].std(axis=1, ddof=1) / p_df[["Rm_ctrl1", "Rm_ctrl2", "Rm_ctrl3"]].mean(axis=1)
    p_df["Rm_ctrl_geom"] = p_df[["Rm_ctrl1", "Rm_ctrl2", "Rm_ctrl3"]].apply(geometric_mean_safe, axis=1)

    s_rm1, s_rm2, s_rm3 = calc_rm_triplet_from_cols(s_df, rep_obj["site_cond_cols"][target_key])
    s_df["Rm_site1"] = s_rm1
    s_df["Rm_site2"] = s_rm2
    s_df["Rm_site3"] = s_rm3
    s_df["CV_site_tech"] = s_df[["Rm_site1", "Rm_site2", "Rm_site3"]].std(axis=1, ddof=1) / s_df[["Rm_site1", "Rm_site2", "Rm_site3"]].mean(axis=1)
    s_df["Rm_site_geom"] = s_df[["Rm_site1", "Rm_site2", "Rm_site3"]].apply(geometric_mean_safe, axis=1)

    # Completeness gate at replicate level:
    # site target complete (all 6 channels) + mapped protein control complete (all 6 channels)
    site_target_cols = rep_obj["site_cond_cols"][target_key]
    s_complete = (s_df[site_target_cols] != 0).all(axis=1)

    ctrl_cols = rep_obj["protein_cond_cols"][control_key]
    p_complete = (p_df[ctrl_cols] != 0).all(axis=1)

    p_df["Control_complete"] = p_complete.values
    p_pass = p_df.loc[p_df["Control_complete"], ["ProteinID", "Rm_ctrl1", "Rm_ctrl2", "Rm_ctrl3", "CV_ctrl_tech", "Rm_ctrl_geom"]]
    p_pass = p_pass.drop_duplicates(subset=["ProteinID"]).set_index("ProteinID")
    p_cv_complete_control = p_pass["CV_ctrl_tech"].copy()

    proteins_ctrl_complete = set(p_df.loc[p_df["Control_complete"], "ProteinID"])
    proteins_site_complete = set(s_df.loc[s_complete, "ProteinID"])

    s_df_all_retrieved = s_df.copy()
    s_df_complete_target = s_df[s_complete].copy()
    s_df = s_df_complete_target[s_df_complete_target["ProteinID"].isin(p_pass.index)].copy()

    s_df["Rm_ctrl1"] = s_df["ProteinID"].map(p_pass["Rm_ctrl1"])
    s_df["Rm_ctrl2"] = s_df["ProteinID"].map(p_pass["Rm_ctrl2"])
    s_df["Rm_ctrl3"] = s_df["ProteinID"].map(p_pass["Rm_ctrl3"])
    s_df["CV_ctrl_tech"] = s_df["ProteinID"].map(p_pass["CV_ctrl_tech"])
    s_df["Rm_ctrl_geom"] = s_df["ProteinID"].map(p_pass["Rm_ctrl_geom"])

    intersect_protein_ids = set(s_df["ProteinID"].dropna().astype(str))
    p_cv_intersect = p_pass.loc[p_pass.index.isin(intersect_protein_ids), "CV_ctrl_tech"].copy()

    rep_tag = rep_name if rep_name else "BR"
    out_dir = os.path.dirname(rep_obj["site_path"])

    plot_protein_site_distribution(
        s_df_all_retrieved,
        f"{rep_tag} - Distribution of {target_label} sites per protein (all retrieved sites)",
        os.path.join(out_dir, f"Branch05_G1_{rep_tag}_site_dist_all_retrieved.tiff"),
        color="#8CB9D0",
    )
    plot_protein_site_distribution(
        s_df_complete_target,
        f"{rep_tag} - Distribution of fully quantified {target_label} sites per protein",
        os.path.join(out_dir, f"Branch05_G1_{rep_tag}_site_dist_complete_target.tiff"),
        color="#A8C9E0",
    )
    plot_protein_site_distribution(
        s_df,
        (
            f"{rep_tag} - Distribution of fully quantified {target_label} sites\n"
            f"on fully quantified {control_label} proteins"
        ),
        os.path.join(out_dir, f"Branch05_G1_{rep_tag}_site_dist_target_on_complete_protein.tiff"),
        color="#14C0CC",
    )

    draw_biorep_cv_reports(
        rep_tag=rep_tag,
        out_dir=out_dir,
        site_cv_complete_target=s_df_complete_target["CV_site_tech"],
        site_cv_intersect=s_df["CV_site_tech"],
        protein_cv_complete_control=p_cv_complete_control,
        protein_cv_intersect=p_cv_intersect,
        scatter_site_cv=s_df["CV_site_tech"],
        scatter_protein_cv=s_df["CV_ctrl_tech"],
        target_label=target_label,
        control_label=control_label,
        tech_cutoff=DEFAULT_TECH_CV_CUTOFF,
    )

    venn_file = os.path.join(out_dir, f"Branch05_G1_{rep_tag}_protein_site_completeness_venn.tiff")
    if len(proteins_ctrl_complete) > 0 and len(proteins_site_complete) > 0:
        plt.figure(figsize=(5, 5))
        v = venn2(
            [proteins_ctrl_complete, proteins_site_complete],
            set_labels=["Control-complete proteins", "Target-site-complete proteins"],
            set_colors=("#8AB7C2", "#9EAAD4"),
            alpha=0.5,
        )
        for text in v.set_labels:
            if text:
                text.set_fontsize(6)
        for text in v.subset_labels:
            if text:
                text.set_fontsize(6)
                text.set_fontweight("bold")

        plt.title(f"{rep_tag} - Protein/site completeness overlap", fontsize=8)
        plt.tight_layout()
        plt.savefig(venn_file, dpi=300, format="tiff")
        plt.close()

    metric_cols = [
        "SiteID",
        "ProteinID",
        "Rm_site1",
        "Rm_site2",
        "Rm_site3",
        "CV_site_tech",
        "Rm_site_geom",
        "Rm_ctrl1",
        "Rm_ctrl2",
        "Rm_ctrl3",
        "CV_ctrl_tech",
        "Rm_ctrl_geom",
    ]
    metrics = s_df[metric_cols].copy().drop_duplicates(subset=["SiteID"]).set_index("SiteID")

    return metrics


def apply_median_of_medians_normalization(site_metrics_by_rep):
    rep_names = list(site_metrics_by_rep.keys())

    site_medians = {}
    ctrl_medians = {}
    for rep in rep_names:
        rdf = site_metrics_by_rep[rep]
        site_medians[rep] = float(rdf["Rm_site_geom"].median())
        ctrl_medians[rep] = float(rdf["Rm_ctrl_geom"].median())

    site_vals = np.array(list(site_medians.values()), dtype=float)
    ctrl_vals = np.array(list(ctrl_medians.values()), dtype=float)
    if np.any(site_vals <= 0) or np.any(ctrl_vals <= 0) or np.any(~np.isfinite(site_vals)) or np.any(~np.isfinite(ctrl_vals)):
        raise ValueError("Median-of-medians normalization requires finite positive medians.")

    global_site_ref = float(np.exp(np.median(np.log(site_vals))))
    global_ctrl_ref = float(np.exp(np.median(np.log(ctrl_vals))))

    norm_rows = []
    for rep in rep_names:
        site_cf = float(np.exp(np.log(global_site_ref) - np.log(site_medians[rep])))
        ctrl_cf = float(np.exp(np.log(global_ctrl_ref) - np.log(ctrl_medians[rep])))

        rdf = site_metrics_by_rep[rep]
        rdf["site_norm_factor"] = site_cf
        rdf["ctrl_norm_factor"] = ctrl_cf
        rdf["Rm_site_adjusted"] = rdf["Rm_site_geom"] * site_cf
        rdf["Rm_ctrl_adjusted"] = rdf["Rm_ctrl_geom"] * ctrl_cf
        rdf["DeltaRm"] = rdf["Rm_site_adjusted"] - rdf["Rm_ctrl_adjusted"]

        norm_rows.append(
            {
                "BioRep": rep,
                "Site_median_geom_Rm": site_medians[rep],
                "Ctrl_median_geom_Rm": ctrl_medians[rep],
                "Global_site_median_of_medians": global_site_ref,
                "Global_ctrl_median_of_medians": global_ctrl_ref,
                "Site_norm_factor": site_cf,
                "Ctrl_norm_factor": ctrl_cf,
            }
        )

    return pd.DataFrame(norm_rows)


def compute_bio_cv_pooled(site_metrics_by_rep):
    rep_names = list(site_metrics_by_rep.keys())
    site_index = next(iter(site_metrics_by_rep.values())).index

    site_raw_cols = []
    ctrl_raw_cols = []
    for rep in rep_names:
        rdf = site_metrics_by_rep[rep]
        site_raw_cols.extend([rdf["Rm_site1"].reindex(site_index).values, rdf["Rm_site2"].reindex(site_index).values, rdf["Rm_site3"].reindex(site_index).values])
        ctrl_raw_cols.extend([rdf["Rm_ctrl1"].reindex(site_index).values, rdf["Rm_ctrl2"].reindex(site_index).values, rdf["Rm_ctrl3"].reindex(site_index).values])

    site_mat = np.column_stack(site_raw_cols)
    ctrl_mat = np.column_stack(ctrl_raw_cols)

    required_n = 3 * len(rep_names)

    with np.errstate(divide="ignore", invalid="ignore"):
        site_cv = np.nanstd(site_mat, axis=1, ddof=1) / np.nanmean(site_mat, axis=1)
        ctrl_cv = np.nanstd(ctrl_mat, axis=1, ddof=1) / np.nanmean(ctrl_mat, axis=1)

    site_valid_n = np.isfinite(site_mat).sum(axis=1)
    ctrl_valid_n = np.isfinite(ctrl_mat).sum(axis=1)

    site_cv = np.where(site_valid_n >= required_n, site_cv, np.nan)
    ctrl_cv = np.where(ctrl_valid_n >= required_n, ctrl_cv, np.nan)

    return pd.Series(site_cv, index=site_index, name="CV_site_biorep_pooled"), pd.Series(ctrl_cv, index=site_index, name="CV_ctrl_biorep_pooled")


def build_global_cv_masks(site_metrics_by_rep, tech_cutoff, bio_cutoff, bio_site_cv, bio_ctrl_cv):
    rep_names = list(site_metrics_by_rep.keys())
    site_index = next(iter(site_metrics_by_rep.values())).index

    tech_mask = pd.Series(True, index=site_index)
    for rep in rep_names:
        rdf = site_metrics_by_rep[rep]
        rep_mask = (
            np.isfinite(rdf["CV_site_tech"]) &
            np.isfinite(rdf["CV_ctrl_tech"]) &
            (rdf["CV_site_tech"] <= tech_cutoff) &
            (rdf["CV_ctrl_tech"] <= tech_cutoff)
        )
        tech_mask &= rep_mask.reindex(site_index, fill_value=False)

    bio_mask = (
        np.isfinite(bio_site_cv) &
        np.isfinite(bio_ctrl_cv) &
        (bio_site_cv <= bio_cutoff) &
        (bio_ctrl_cv <= bio_cutoff)
    )

    global_mask = tech_mask & bio_mask
    return tech_mask, bio_mask, global_mask


def main():
    print("=== Branch05_G1: Multi-BioRep paired Site/Protein pipeline ===")

    fasta_file = normalize_path(input("Enter UniProt FASTA path (used for site retrieval in each BioRep): ").strip())
    fasta_dict = load_fasta_sequences(fasta_file)
    motif_len_in = input("Enter motif window length for site retrieval (odd, default=9): ").strip()
    motif_len = 9 if motif_len_in == "" else int(motif_len_in)
    if motif_len <= 0 or motif_len % 2 == 0:
        raise ValueError("Motif window length must be a positive odd number.")

    n_reps = int(input("Enter number of biological replicates (>=3): ").strip())
    if n_reps < 3:
        raise ValueError("At least 3 biological replicates are required.")

    use_unique_in = input("Apply protein unique-peptide >=2 filter for quantified overlap plot? (YES/NO, default=NO): ").strip().upper()
    use_unique_ge2 = (use_unique_in == "YES")

    reps_raw = {}
    protein_condition_defs = None
    site_condition_defs = None

    for i in range(1, n_reps + 1):
        rep_name = f"BR{i}"
        p_path = normalize_path(input(f"{rep_name} protein CSV path: ").strip())
        s_path = normalize_path(input(f"{rep_name} site CSV path: ").strip())
        if not os.path.exists(p_path):
            raise FileNotFoundError(p_path)
        if not os.path.exists(s_path):
            raise FileNotFoundError(s_path)

        rep_obj = prepare_replicate(
            rep_name=rep_name,
            protein_path=p_path,
            site_path=s_path,
            protein_condition_defs=protein_condition_defs,
            site_condition_defs=site_condition_defs,
            use_unique_ge2=use_unique_ge2,
            fasta_dict=fasta_dict,
            motif_len=motif_len,
        )
        protein_condition_defs = rep_obj["protein_condition_defs"]
        site_condition_defs = rep_obj["site_condition_defs"]
        reps_raw[rep_name] = rep_obj

    print("Detected canonical protein conditions from BR1 (protein Sample blocks):")
    for cond in protein_condition_defs:
        print(f"- {cond['name']} [{cond['key']}]")

    print("Detected canonical site conditions from BR1 (site Sample blocks):")
    for cond in site_condition_defs:
        print(f"- {cond['name']} [{cond['key']}]")

    if len(protein_condition_defs) < 1:
        raise ValueError("At least 1 protein condition is required to choose control.")
    if len(site_condition_defs) < 1:
        raise ValueError("At least 1 site condition is required to choose target.")

    available_protein_names = [c["name"] for c in protein_condition_defs]
    available_site_names = [c["name"] for c in site_condition_defs]

    print("Available protein condition names for control input:")
    for name in available_protein_names:
        print(f"- {name}")

    print("Available site condition names for target input:")
    for name in available_site_names:
        print(f"- {name}")

    control_name = input("Enter control condition name (must match canonical name): ").strip().lower()
    protein_key_map = {c["name"].strip().lower(): c["key"] for c in protein_condition_defs}
    protein_name_map = {c["name"].strip().lower(): c["name"] for c in protein_condition_defs}
    if control_name not in protein_key_map:
        raise ValueError(
            "Control condition not found. "
            f"Input='{control_name}'. Available protein conditions={available_protein_names}"
        )
    control_key = protein_key_map[control_name]
    control_display_name = protein_name_map[control_name]

    target_name = input("Enter target site condition name for completeness and downstream analysis: ").strip().lower()
    site_key_map = {c["name"].strip().lower(): c["key"] for c in site_condition_defs}
    site_name_map = {c["name"].strip().lower(): c["name"] for c in site_condition_defs}
    if target_name not in site_key_map:
        raise ValueError(
            "Target site condition not found. "
            f"Input='{target_name}'. Available site conditions={available_site_names}."
        )
    target_key = site_key_map[target_name]
    target_display_name = site_name_map[target_name]

    print(f"Condition scope fixed: protein-group uses control '{control_key}' only; site-group uses target '{target_key}' only.")

    out_dir = os.path.dirname(next(iter(reps_raw.values()))["site_path"])
    out_prefix = os.path.join(out_dir, "Branch05_G1")

    # Protein overlap plots: all conditions, quantified-any rule, per BioRep
    for rep_name, rep_obj in reps_raw.items():
        overlap_file = draw_quant_overlap(
            rep_obj["quantified_sets"],
            list(rep_obj["quantified_sets"].keys()),
            f"{out_prefix}_{rep_name}",
        )
        if overlap_file is not None:
            print(f"{rep_name} overlap plot saved: {overlap_file}")

    # Build per-rep site metrics then intersect sites across all BioReps
    print("Generating per-BioRep CV distributions and site/protein CV scatter BEFORE asking final CV thresholds...")
    site_metrics_by_rep = {}
    for rep_name, rep_obj in reps_raw.items():
        site_metrics_by_rep[rep_name] = compute_replicate_site_metrics(
            rep_obj,
            control_key,
            target_key,
            rep_name=rep_name,
            control_label=control_display_name,
            target_label=target_display_name,
        )

    common_sites = sorted(set.intersection(*[set(df.index) for df in site_metrics_by_rep.values()]))
    if len(common_sites) == 0:
        raise ValueError("No common sites remain across BioReps after strict completeness gate.")

    for rep_name in site_metrics_by_rep:
        site_metrics_by_rep[rep_name] = site_metrics_by_rep[rep_name].loc[common_sites].copy()

    rep_names = list(site_metrics_by_rep.keys())

    print(f"Common strict-complete sites across all BioReps: {len(common_sites)}")

    norm_df = apply_median_of_medians_normalization(site_metrics_by_rep)
    norm_file = f"{out_prefix}_normalization_factors.csv"
    norm_df.to_csv(norm_file, index=False)
    print(f"Normalization factors saved: {norm_file}")

    bio_site_cv, bio_ctrl_cv = compute_bio_cv_pooled(site_metrics_by_rep)

    cv_method_note = pd.DataFrame(
        [
            {
                "Layer": "Technical_CV",
                "Definition": "CV within each BioRep using 3 raw Rm values (Rm1,Rm2,Rm3)",
                "Raw_value_count_per_feature": 3,
            },
            {
                "Layer": "Biological_pooled_CV",
                "Definition": "CV across BioReps using pooled raw Rm values",
                "Raw_value_count_per_feature": int(3 * len(rep_names)),
            },
        ]
    )
    cv_method_file = f"{out_prefix}_cv_method_note.csv"
    cv_method_note.to_csv(cv_method_file, index=False)
    print(
        f"CV method note saved: {cv_method_file} "
        f"(Tech uses 3 values/BioRep; Bio pooled uses 3*{len(rep_names)}={3*len(rep_names)} values)"
    )

    # CV preview plots and pass preview
    tech_site_max = pd.concat([site_metrics_by_rep[r]["CV_site_tech"] for r in site_metrics_by_rep], axis=1).max(axis=1)
    tech_ctrl_max = pd.concat([site_metrics_by_rep[r]["CV_ctrl_tech"] for r in site_metrics_by_rep], axis=1).max(axis=1)

    preview_files = [
        plot_cv_distribution(tech_site_max, "Tech CV (site, max across BioReps)", f"{out_prefix}_CV_preview_tech_site.tiff", [0.30, 0.40]),
        plot_cv_distribution(tech_ctrl_max, "Tech CV (protein control, max across BioReps)", f"{out_prefix}_CV_preview_tech_ctrl.tiff", [0.30, 0.40]),
        plot_cv_distribution(bio_site_cv, "BioRep pooled CV (site raw Rm)", f"{out_prefix}_CV_preview_bio_site.tiff", [0.30, 0.40]),
        plot_cv_distribution(bio_ctrl_cv, "BioRep pooled CV (protein control raw Rm)", f"{out_prefix}_CV_preview_bio_ctrl.tiff", [0.30, 0.40]),
    ]
    for rep in rep_names:
        rdf = site_metrics_by_rep[rep]
        preview_files.extend(
            [
                plot_cv_distribution(
                    rdf["CV_site_tech"],
                    (
                        f"Tech CV (site, {rep}, 3 raw Rm within BioRep)\n"
                        f"{target_display_name} (Fully quantified sites of fully quantified proteins)"
                    ),
                    f"{out_prefix}_CV_preview_tech_site_{rep}.tiff",
                    [0.30, 0.40],
                ),
                plot_cv_distribution(
                    rdf["CV_ctrl_tech"],
                    (
                        f"Tech CV (protein control, {rep}, 3 raw Rm within BioRep)\n"
                        f"{control_display_name} Proteins (Fully quantified proteins with fully quantified N-glycosites)"
                    ),
                    f"{out_prefix}_CV_preview_tech_ctrl_{rep}.tiff",
                    [0.30, 0.40],
                ),
            ]
        )
    preview_files = [f for f in preview_files if f is not None]
    if preview_files:
        print("CV preview plots saved:")
        for f in preview_files:
            print(f"- {f}")

    preview_rows = []
    for tc in TECH_CV_PREVIEW_CANDIDATES:
        for bc in BIO_CV_PREVIEW_CANDIDATES:
            _, _, gmask = build_global_cv_masks(site_metrics_by_rep, tc, bc, bio_site_cv, bio_ctrl_cv)
            preview_rows.append(
                {
                    "Tech_CV_cutoff": float(tc),
                    "Bio_CV_cutoff": float(bc),
                    "Global_pass_count": int(gmask.sum()),
                    "Total_count": int(len(gmask)),
                    "Global_pass_percent": float(100.0 * gmask.sum() / max(len(gmask), 1)),
                }
            )
    cv_preview_df = pd.DataFrame(preview_rows)
    cv_preview_file = f"{out_prefix}_cv_preview_grid.csv"
    cv_preview_df.to_csv(cv_preview_file, index=False)
    print(f"CV threshold preview table saved: {cv_preview_file}")

    input("Review the generated CV plots/tables above, then press Enter to set final CV thresholds... ")

    tech_cutoff = ask_float_with_default(
        f"Enter final technical CV cutoff (default {DEFAULT_TECH_CV_CUTOFF}): ",
        DEFAULT_TECH_CV_CUTOFF,
        min_value=0.0,
    )
    bio_cutoff = ask_float_with_default(
        f"Enter final biological CV cutoff (default {DEFAULT_BIO_CV_CUTOFF}): ",
        DEFAULT_BIO_CV_CUTOFF,
        min_value=0.0,
    )

    tech_mask, bio_mask, global_cv_mask = build_global_cv_masks(site_metrics_by_rep, tech_cutoff, bio_cutoff, bio_site_cv, bio_ctrl_cv)

    cv_summary = pd.DataFrame(
        [
            {
                "Mask": "Global_Tech_CV_QC",
                "Pass_count": int(tech_mask.sum()),
                "Total_count": int(len(tech_mask)),
                "Pass_percent": float(100.0 * tech_mask.sum() / max(len(tech_mask), 1)),
                "Tech_CV_cutoff": float(tech_cutoff),
                "Bio_CV_cutoff": float(bio_cutoff),
            },
            {
                "Mask": "Global_Bio_CV_QC",
                "Pass_count": int(bio_mask.sum()),
                "Total_count": int(len(bio_mask)),
                "Pass_percent": float(100.0 * bio_mask.sum() / max(len(bio_mask), 1)),
                "Tech_CV_cutoff": float(tech_cutoff),
                "Bio_CV_cutoff": float(bio_cutoff),
            },
            {
                "Mask": "Global_TwoLayer_CV_QC",
                "Pass_count": int(global_cv_mask.sum()),
                "Total_count": int(len(global_cv_mask)),
                "Pass_percent": float(100.0 * global_cv_mask.sum() / max(len(global_cv_mask), 1)),
                "Tech_CV_cutoff": float(tech_cutoff),
                "Bio_CV_cutoff": float(bio_cutoff),
            },
        ]
    )
    cv_summary_file = f"{out_prefix}_cv_qc_summary.csv"
    cv_summary.to_csv(cv_summary_file, index=False)
    print(f"CV QC summary saved: {cv_summary_file}")

    # Normality by BioRep on CV-passed sites
    normality_rows = []
    all_normal = True
    for rep in rep_names:
        rdf = site_metrics_by_rep[rep]
        deltas = rdf.loc[global_cv_mask, "DeltaRm"]
        plot_delta_normality_diagnostics(deltas, rep, out_prefix)
        is_normal, msg, method_used, n_used = assess_distribution_normality(deltas.values, alpha=0.05)
        normality_rows.append({"BioRep": rep, "n": n_used, "method": method_used, "summary": msg, "is_normal": is_normal})
        print(f"{rep}: {msg}")
        if not is_normal:
            all_normal = False

    normality_df = pd.DataFrame(normality_rows)
    normality_file = f"{out_prefix}_normality_summary.csv"
    normality_df.to_csv(normality_file, index=False)
    print(f"Normality summary saved: {normality_file}")

    # Build consolidated table
    result = pd.DataFrame(index=common_sites)
    result["Pass_Tech_CV_QC"] = tech_mask
    result["Pass_Bio_CV_QC"] = bio_mask
    result["Pass_TwoLayer_CV_QC"] = global_cv_mask
    result["CV_site_biorep_pooled"] = bio_site_cv
    result["CV_ctrl_biorep_pooled"] = bio_ctrl_cv

    delta_mat = []
    for rep in rep_names:
        rdf = site_metrics_by_rep[rep]
        for c in [
            "ProteinID",
            "Rm_site1",
            "Rm_site2",
            "Rm_site3",
            "CV_site_tech",
            "Rm_site_geom",
            "Rm_site_adjusted",
            "Rm_ctrl1",
            "Rm_ctrl2",
            "Rm_ctrl3",
            "CV_ctrl_tech",
            "Rm_ctrl_geom",
            "Rm_ctrl_adjusted",
            "DeltaRm",
            "site_norm_factor",
            "ctrl_norm_factor",
        ]:
            result[f"{rep}__{c}"] = rdf[c].reindex(common_sites)
        delta_mat.append(rdf["DeltaRm"].reindex(common_sites).values)

    delta_mat = np.column_stack(delta_mat)
    pass_same_sign = (np.all(delta_mat > 0, axis=1) | np.all(delta_mat < 0, axis=1))
    pass_abs_delta = np.all(np.isfinite(delta_mat) & (np.abs(delta_mat) > 0.1), axis=1)

    if all_normal:
        print("Global normality decision: all BioRep DeltaRm distributions are approximately normal.")
        print("Running robust sigma by CV-bin in each BioRep, then consensus integration.")

        fdr_mat = []
        pval_mat = []
        for rep in rep_names:
            rdf = site_metrics_by_rep[rep]
            cv_ok_idx = list(result.index[result["Pass_TwoLayer_CV_QC"]])
            pvals = pd.Series(1.0, index=result.index)
            fdrs = pd.Series(1.0, index=result.index)
            if len(cv_ok_idx) >= 3:
                rep_cv = rdf.loc[cv_ok_idx, ["CV_site_tech", "CV_ctrl_tech"]].max(axis=1)
                order = list(rep_cv.sort_values(ascending=False).index)
                bin_pvals = robust_bin_pvals(rdf.loc[cv_ok_idx, "DeltaRm"], order, bin_size=300)
                pv = np.where(np.isfinite(bin_pvals.values), bin_pvals.values, 1.0)
                _, fd, _, _ = multipletests(pv, method="fdr_bh")
                pvals.loc[cv_ok_idx] = pv
                fdrs.loc[cv_ok_idx] = fd

            result[f"{rep}__p_value"] = pvals.values
            result[f"{rep}__FDR"] = fdrs.values
            pval_mat.append(pvals.values)
            fdr_mat.append(fdrs.values)

            volcano_file = f"{out_prefix}_volcano_normal_{rep}.tiff"
            plot_volcano(
                delta_series=rdf.loc[result["Pass_TwoLayer_CV_QC"], "DeltaRm"],
                fdr_series=fdrs.loc[result["Pass_TwoLayer_CV_QC"]],
                title=f"Volcano (normal branch, {rep})",
                out_file=volcano_file,
            )

        pval_mat = np.column_stack(pval_mat)
        fdr_mat = np.column_stack(fdr_mat)

        pass_fdr_all_01 = np.all(fdr_mat < 0.1, axis=1)
        pass_fdr_one_005 = np.any(fdr_mat < 0.05, axis=1)
        sig = result["Pass_TwoLayer_CV_QC"].values & pass_same_sign & pass_abs_delta & pass_fdr_all_01 & pass_fdr_one_005

        result["Analysis_mode"] = "global_normal_robust_sigma"
        result["Pass_same_sign_all_BioReps"] = pass_same_sign
        result["Pass_absDelta_gt_0p1_all_BioReps"] = pass_abs_delta
        result["Pass_all_FDR_lt_0p1"] = pass_fdr_all_01
        result["Pass_any_FDR_lt_0p05"] = pass_fdr_one_005
        result["Significant"] = sig
        result["DeltaRm_display"] = np.nanmedian(delta_mat, axis=1)

    else:
        print("Global normality decision: at least one BioRep DeltaRm distribution is non-normal.")
        method_in = input("Choose paired one-tailed test across BioReps [paired_t/wilcoxon, default=paired_t]: ").strip().lower()
        paired_method = method_in if method_in in {"paired_t", "wilcoxon"} else "paired_t"
        print(f"Using non-normal branch paired method: {paired_method}, one-tailed (site > control)")

        pvals = np.ones(len(result), dtype=float)
        n_pairs = np.zeros(len(result), dtype=int)
        for i, sid in enumerate(result.index):
            if not bool(result.at[sid, "Pass_TwoLayer_CV_QC"]):
                pvals[i] = 1.0
                n_pairs[i] = 0
                continue

            x = np.array([site_metrics_by_rep[r].at[sid, "Rm_site_adjusted"] for r in rep_names], dtype=float)
            y = np.array([site_metrics_by_rep[r].at[sid, "Rm_ctrl_adjusted"] for r in rep_names], dtype=float)
            p_one, n_valid = paired_one_tailed_pvalue(y, x, method=paired_method, side="right")
            pvals[i] = p_one
            n_pairs[i] = n_valid

        _, fdrs, _, _ = multipletests(pvals, method="fdr_bh")
        sig = result["Pass_TwoLayer_CV_QC"].values & pass_same_sign & pass_abs_delta & (fdrs < 0.05)

        delta_mean = np.nanmean(delta_mat, axis=1)

        result["Analysis_mode"] = f"global_non_normal_paired_{paired_method}_right_tailed"
        result["Paired_test"] = paired_method
        result["n_pairs_used"] = n_pairs
        result["p_value"] = pvals
        result["FDR"] = fdrs
        result["Pass_same_sign_all_BioReps"] = pass_same_sign
        result["Pass_absDelta_gt_0p1_all_BioReps"] = pass_abs_delta
        result["Significant"] = sig
        result["DeltaRm_display"] = delta_mean

        volcano_file = f"{out_prefix}_volcano_non_normal_{paired_method}.tiff"
        plot_volcano(
            delta_series=result.loc[result["Pass_TwoLayer_CV_QC"], "DeltaRm_display"],
            fdr_series=result.loc[result["Pass_TwoLayer_CV_QC"], "FDR"],
            title=f"Volcano ({paired_method}, one-tailed BH)",
            out_file=volcano_file,
        )

        # TOP25 boxplots
        sig_df = result[result["Significant"]].copy()
        top_pool = sig_df if not sig_df.empty else result[result["Pass_TwoLayer_CV_QC"]].copy()
        top_df = top_pool.sort_values("DeltaRm_display", ascending=False).head(25)

        if not top_df.empty:
            n_panels = top_df.shape[0]
            ncols = 5
            nrows = int(np.ceil(n_panels / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.0, nrows * 3.0), squeeze=False)
            axes = axes.flatten()
            rng = np.random.default_rng(2026)

            for i, sid in enumerate(top_df.index):
                ax = axes[i]
                site_vals = np.array([site_metrics_by_rep[r].at[sid, "Rm_site_adjusted"] for r in rep_names], dtype=float)
                ctrl_vals = np.array([site_metrics_by_rep[r].at[sid, "Rm_ctrl_adjusted"] for r in rep_names], dtype=float)

                ax.boxplot([ctrl_vals], positions=[1], widths=0.5, showfliers=False, patch_artist=True,
                           boxprops=dict(facecolor="#4C78A8", alpha=0.35), medianprops=dict(color="#4C78A8"))
                ax.boxplot([site_vals], positions=[2], widths=0.5, showfliers=False, patch_artist=True,
                           boxprops=dict(facecolor="#F58518", alpha=0.35), medianprops=dict(color="#F58518"))

                x1 = 1 + rng.normal(0, 0.04, size=len(ctrl_vals))
                x2 = 2 + rng.normal(0, 0.04, size=len(site_vals))
                ax.scatter(x1, ctrl_vals, s=16, color="#4C78A8", zorder=3)
                ax.scatter(x2, site_vals, s=16, color="#F58518", zorder=3)

                sid_short = sid if len(str(sid)) <= 28 else str(sid)[:25] + "..."
                fdr_val = float(top_df.at[sid, "FDR"]) if "FDR" in top_df.columns else np.nan
                star = pval_to_stars(fdr_val) if np.isfinite(fdr_val) else ""
                ax.set_title(f"{sid_short}\nFDR={fdr_val:.2e} {star}" if np.isfinite(fdr_val) else sid_short, fontsize=7)
                ax.set_xticks([1, 2])
                ax.set_xticklabels(["Ctrl", "Site"], fontsize=7)
                ax.tick_params(axis="y", labelsize=7)

            for j in range(n_panels, len(axes)):
                axes[j].axis("off")

            plt.tight_layout()
            top25_file = f"{out_prefix}_Top25_boxplot_{paired_method}.tiff"
            plt.savefig(top25_file, dpi=300, format="tiff")
            plt.close()
            print(f"Top25 boxplot saved: {top25_file}")

    # Final outputs
    full_file = f"{out_prefix}_full_results.csv"
    result.reset_index().rename(columns={"index": "SiteID"}).to_csv(full_file, index=False)
    print(f"Full results saved: {full_file}")

    diff = result[result["Significant"]].copy()
    diff_file = f"{out_prefix}_Differential_sites.csv"
    diff.reset_index().rename(columns={"index": "SiteID"}).to_csv(diff_file, index=False)
    print(f"Differential sites saved: {diff_file}")

    sig_summary = pd.DataFrame(
        [
            {
                "Target_condition": target_key,
                "Control_condition": control_key,
                "Common_sites_after_strict_gate": int(len(common_sites)),
                "Global_two_layer_cv_pass": int(global_cv_mask.sum()),
                "Significant_sites": int(result["Significant"].sum()),
                "Tech_CV_cutoff": float(tech_cutoff),
                "Bio_CV_cutoff": float(bio_cutoff),
            }
        ]
    )
    sig_summary_file = f"{out_prefix}_significance_summary.csv"
    sig_summary.to_csv(sig_summary_file, index=False)
    print(f"Significance summary saved: {sig_summary_file}")

    print("Branch05_G1 pipeline finished.")


if __name__ == "__main__":
    main()
