#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Refined-TPP plus - BioRep-only Analysis Pipeline (No Technical Replicates)
-------------------------------------------------------------------------
Created on: Thu Mar 19 2026
Authors: Haoru Song, Chenxin Li, Bin Fu
Affiliation: Professor Haojie Lu's Lab, Fudan University, Shanghai, China.

Description
-----------
This script is designed for experiments with biological replicates only.
Unlike Branch01_P1_update.py, no within-file technical replicates are used.

Data design in this script:
1) One CSV file contains all conditions.
2) For each condition, TMT-6plex channels encode 3 biological replicates using
   fixed reporter pairs:
   - BR1: 129 / 126
   - BR2: 130 / 127
   - BR3: 131 / 128
3) Quant completeness for downstream analysis requires proteins to have complete
   positive quantification in all channels of all conditions.
4) Pre-analysis overlap plot (Venn/UpSet) is based on quantified proteins per
   biological replicate and does not require complete quantification. A
   user switch controls whether "unique peptide >= 2" is applied for this plot.

Normalization and QC:
---------------------
- No geometric averaging across technical replicates is performed.
- Raw Rm values are normalized by median-of-medians in log space using Vehicle
  medians across BR1/BR2/BR3.
- Biological-replicate CV-QC is computed per condition from raw Rm across BR1-3.
- Global BioRep CV-QC requires passing CV cutoff in all conditions before
  differential analysis.

Differential strategy:
----------------------
- Global DeltaRm normality decision over all (BioRep x condition) panels.

- If all normal:
	global robust sigma branch (no CV binning).
	Replicate-level FDR integration rule:
	all FDR < 0.1 and at least one FDR < 0.05,
	plus abs(DeltaRm) > 0.1 and same-sign DeltaRm in all BioReps.
- If any non-normal:
	paired test across adjusted Rm of BR1/BR2/BR3 (ttest default; optional wilcoxon),
	with FDR < 0.05 plus DeltaRm consistency filters.

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

REQUIRED_REPORTERS = (126, 127, 128, 129, 130, 131)
BR_LABELS = ("BR1", "BR2", "BR3")
BR_REPORTER_PAIRS = {
	"BR1": (129, 126),
	"BR2": (130, 127),
	"BR3": (131, 128),
}

DEFAULT_BIO_CV_QC_CUTOFF = 0.4
BIO_CV_PREVIEW_CANDIDATES = (0.3, 0.35, 0.4, 0.45, 0.5)
DEFAULT_DELTA_CUTOFF = 0.1


# ------------------------------
# Helper functions
# ------------------------------
def normalize_path(path):
	return os.path.normpath(os.path.abspath(path))


def parse_single_csv_path(raw_text):
	"""Parse a single CSV path, tolerant to quotes/spaces."""
	cleaned = raw_text.strip().strip('"').strip("'")
	if not cleaned:
		raise ValueError("Input path is empty.")
	return normalize_path(cleaned)


def ask_yes_no(prompt, default=False):
	"""Prompt user for yes/no with default value."""
	default_token = "y" if default else "n"
	while True:
		raw = input(prompt).strip().lower()
		if raw == "":
			return default
		if raw in {"y", "yes", "1", "true", "t"}:
			return True
		if raw in {"n", "no", "0", "false", "f"}:
			return False
		logger.warning(f"Invalid input '{raw}'. Please enter y or n.")
		logger.info(f"Press Enter for default: '{default_token}'.")


def ask_int_with_min(prompt, min_value=1):
	"""Prompt user for integer >= min_value."""
	while True:
		raw = input(prompt).strip()
		try:
			value = int(raw)
		except ValueError:
			logger.warning(f"Invalid integer input '{raw}'.")
			continue
		if value < min_value:
			logger.warning(f"Input must be >= {min_value}.")
			continue
		return value


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


def ask_choice_with_default(prompt, allowed_choices, default):
	"""Prompt user for a string choice with default value."""
	allowed = {c.lower() for c in allowed_choices}
	default = default.lower()
	while True:
		raw = input(prompt).strip().lower()
		if raw == "":
			return default
		if raw in allowed:
			return raw
		logger.warning(f"Invalid choice '{raw}'. Allowed: {sorted(allowed)}")


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


def summarize_global_bio_cv_pass_across_cutoffs(biocv_by_key, condition_defs, cutoff_values):
	"""Preview global BioRep CV-QC pass counts across candidate cutoffs."""
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
	"""Detect protein identifier column (prefer Accession over Protein Group)."""
	preferred = ["Accession", "Accession ", "Protein Group", "ProteinGroup", "Protein_Group"]
	for c in preferred:
		if c in df.columns:
			return c

	# Tolerant fallback for header variants such as ' accession ' or mixed spacing.
	normalized = {re.sub(r"\s+", "", str(c)).lower(): c for c in df.columns}
	if "accession" in normalized:
		return normalized["accession"]
	if "proteingroup" in normalized:
		return normalized["proteingroup"]

	return df.columns[0]


def canonicalize_protein_id(raw_id):
	"""
	Normalize protein identifier text.
	Supports common formats:
	- Q01105|SET_HUMAN -> Q01105
	- sp|Q01105|SET_HUMAN -> Q01105
	"""
	text = str(raw_id).strip()
	if text == "":
		return ""

	parts = [p.strip() for p in text.split("|") if p is not None]
	if len(parts) >= 3 and parts[0].lower() in {"sp", "tr", "up"}:
		return parts[1].strip()
	if len(parts) >= 2:
		# For formats like Q01105|SET_HUMAN keep the first token as accession.
		return parts[0].strip() if parts[0].strip() else parts[1].strip()

	return text


def detect_sample_blocks(df):
	"""
	Detect 'Sample N' grouped columns. Returns ordered sample indices and mapping.
	If headers contain 'Sample <N>' pattern this collects columns for each sample index.
	"""
	pattern = re.compile(r"Sample\s*(\d+)\b", flags=re.I)
	sample_cols = defaultdict(list)
	sample_order = []
	for col in df.columns:
		m = pattern.search(str(col))
		if m:
			idx = int(m.group(1))
			sample_cols[idx].append(col)
			if idx not in sample_order:
				sample_order.append(idx)
	sample_order.sort()
	return sample_order, sample_cols


def extract_reporter_fragment(colname):
	"""Return a short reporter fragment (e.g. 'TMT-126')."""
	m = re.search(r"(126|127|128|129|130|131)", str(colname))
	if m:
		return f"TMT-{m.group(1)}"
	return str(colname)


def draw_quant_overlap(per_sample_sets, sample_names, out_prefix):
	"""Draw overlap plot of quantified proteins across biological replicates."""
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
		plt.title("Quantified proteins overlap by BioRep")
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


def robust_global_pvals(deltas):
	"""Compute robust p-values globally (no CV binning)."""
	aligned = pd.Series(1.0, index=deltas.index)
	vals = pd.Series(deltas).replace([np.inf, -np.inf], np.nan).dropna()
	if vals.empty:
		return aligned

	p15 = np.percentile(vals.values, 15.87)
	p50 = np.percentile(vals.values, 50.0)
	p84 = np.percentile(vals.values, 84.13)
	left_sigma = p50 - p15
	right_sigma = p84 - p50
	left_sigma = left_sigma if left_sigma > 0 else np.nan
	right_sigma = right_sigma if right_sigma > 0 else np.nan

	for idx, val in vals.items():
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
		aligned.at[idx] = pv

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
	"""Plot volcano for robust branch."""
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


def process_single_file(file_path):
	"""
	Process one CSV containing all conditions and 3 biological replicates encoded
	by TMT pair mapping (BR1 129/126, BR2 130/127, BR3 131/128).
	"""
	df = pd.read_csv(file_path)
	protein_col = detect_protein_col(df)
	unique_col = detect_unique_col(df)

	sample_order, sample_cols_map = detect_sample_blocks(df)
	if not sample_order:
		tmt_cols = [c for c in df.columns if re.search(r"(126|127|128|129|130|131)", str(c))]
		if not tmt_cols:
			raise ValueError("No TMT reporter-like columns (126..131) detected.")

		num_conditions = ask_int_with_min(
			"Enter number of conditions represented in this CSV file: ",
			min_value=1,
		)
		if len(tmt_cols) % num_conditions != 0:
			raise ValueError(
				"Reporter columns cannot be evenly split by conditions: "
				f"reporter_cols={len(tmt_cols)}, conditions={num_conditions}"
			)

		channels_per_condition = len(tmt_cols) // num_conditions
		if channels_per_condition < 6:
			raise ValueError(
				f"Each condition needs at least 6 reporter channels, got {channels_per_condition}."
			)

		sample_order = list(range(1, num_conditions + 1))
		sample_cols_map = {
			s: tmt_cols[i * channels_per_condition : (i + 1) * channels_per_condition]
			for i, s in enumerate(sample_order)
		}

	condition_defs = []
	used_keys = set()
	for i, s in enumerate(sample_order, start=1):
		display_name = input(
			f"Please enter condition name for Sample {s} (default Condition_{i}): "
		).strip()
		if not display_name:
			display_name = f"Condition_{i}"

		base_key = make_safe_key(display_name)
		key = base_key
		suffix = 2
		while key in used_keys:
			key = f"{base_key}_{suffix}"
			suffix += 1
		used_keys.add(key)
		condition_defs.append({"name": display_name, "key": key, "sample_idx": s})

	sample_to_cond = {c["sample_idx"]: c for c in condition_defs}
	mapping_str = "; ".join(
		[f"Sample{c['sample_idx']}->{c['name']}[{c['key']}]" for c in condition_defs]
	)
	logger.info(f"Condition mapping: {mapping_str}")

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
			suffix = 2
			while newname in col_map.values():
				newname = f"{cond_key} {repfrag}_{suffix}"
				suffix += 1
			col_map[old] = newname
			cond_block_cols[cond_key].append(newname)

	df = df.rename(columns=col_map)

	cond_reporter_map = {}
	for cond in condition_defs:
		key = cond["key"]
		reporter_map = {}
		for col in cond_block_cols[key]:
			m = re.search(r"(126|127|128|129|130|131)", col)
			if m:
				reporter_map[int(m.group(1))] = col

		missing = [r for r in REQUIRED_REPORTERS if r not in reporter_map]
		if missing:
			raise ValueError(f"Condition '{cond['name']}' missing reporter channels: {missing}")

		cond_reporter_map[key] = reporter_map

	protein_ids = df[protein_col].map(canonicalize_protein_id)
	valid_pid_mask = protein_ids != ""

	use_unique_filter_for_overlap = False
	if unique_col is not None:
		use_unique_filter_for_overlap = ask_yes_no(
			"Apply unique peptide >=2 filter before overlap plot? [y/N]: ",
			default=False,
		)

	unique_mask = pd.Series(True, index=df.index)
	if use_unique_filter_for_overlap:
		unique_mask = df[unique_col].fillna(0) >= 2

	per_biorep_sets = {}
	for br in BR_LABELS:
		numerator, denominator = BR_REPORTER_PAIRS[br]
		br_cols = []
		for cond in condition_defs:
			key = cond["key"]
			br_cols.extend([cond_reporter_map[key][numerator], cond_reporter_map[key][denominator]])

		has_signal = (df[br_cols].fillna(0) > 0).any(axis=1)
		set_mask = has_signal & valid_pid_mask & unique_mask
		per_biorep_sets[br] = set(protein_ids[set_mask])

	overlap_plot = draw_quant_overlap(
		per_biorep_sets,
		list(BR_LABELS),
		f"{os.path.splitext(file_path)[0]}_BioRepOnly",
	)
	if overlap_plot is not None:
		logger.info(f"Quantified overlap plot saved: {overlap_plot}")

	complete_mask = valid_pid_mask.copy()
	for cond in condition_defs:
		key = cond["key"]
		cols = [cond_reporter_map[key][rid] for rid in REQUIRED_REPORTERS]
		vals = df[cols].apply(pd.to_numeric, errors="coerce")
		cond_complete = np.isfinite(vals).all(axis=1) & (vals > 0).all(axis=1)
		complete_mask &= cond_complete

	df_complete = df.loc[complete_mask].copy()
	df_complete["ProteinID"] = protein_ids.loc[complete_mask]

	metrics = pd.DataFrame({"ProteinID": df_complete["ProteinID"]})
	for cond in condition_defs:
		key = cond["key"]
		reporter_map = cond_reporter_map[key]
		for br in BR_LABELS:
			numerator, denominator = BR_REPORTER_PAIRS[br]
			rm = df_complete[reporter_map[numerator]] / df_complete[reporter_map[denominator]]
			metrics[f"Rm_raw__{key}__{br}"] = rm.replace([np.inf, -np.inf], np.nan)

	if metrics["ProteinID"].duplicated().any():
		numeric_cols = [c for c in metrics.columns if c != "ProteinID"]
		metrics = metrics.groupby("ProteinID", as_index=False)[numeric_cols].median()

	metrics = metrics.set_index("ProteinID").sort_index()

	summary = {
		"File": file_path,
		"Protein_column": protein_col,
		"Unique_column": unique_col if unique_col is not None else "",
		"Overlap_unique_ge2_filter": bool(use_unique_filter_for_overlap),
		"Total_rows": int(df.shape[0]),
		"Non_empty_proteinID_count": int(valid_pid_mask.sum()),
		"Quant_complete_all_channels_count": int(complete_mask.sum()),
		"Condition_mapping": mapping_str,
	}

	logger.info(
		"Quant completeness filter (all conditions x all channels complete) retained proteins: "
		f"{summary['Quant_complete_all_channels_count']}"
	)

	return {
		"metrics_df": metrics,
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


def apply_median_of_medians_normalization(metrics_df, condition_defs, vehicle_key):
	"""
	Batch correction in log domain:
	- Compute Vehicle median per biological replicate (BR1-3)
	- Global reference = median(Vehicle medians) in log space
	- Normalize each condition x BioRep by global_ref / median(condition, that BioRep)
	"""
	vehicle_medians = {}
	for br in BR_LABELS:
		vcol = f"Rm_raw__{vehicle_key}__{br}"
		vmed = float(metrics_df[vcol].replace([np.inf, -np.inf], np.nan).dropna().median())
		vehicle_medians[br] = vmed

	vehicle_vals = np.array([vehicle_medians[br] for br in BR_LABELS], dtype=float)
	if np.any(~np.isfinite(vehicle_vals)) or np.any(vehicle_vals <= 0):
		raise ValueError("Vehicle medians must be finite positive values for log-domain normalization.")

	global_vehicle_median = float(np.exp(np.median(np.log(vehicle_vals))))

	norm_factors = {br: {} for br in BR_LABELS}
	for br in BR_LABELS:
		for cond in condition_defs:
			key = cond["key"]
			raw_col = f"Rm_raw__{key}__{br}"
			sample_median = float(metrics_df[raw_col].replace([np.inf, -np.inf], np.nan).dropna().median())

			if (not np.isfinite(sample_median)) or sample_median <= 0:
				cf = np.nan
			else:
				cf = float(np.exp(np.log(global_vehicle_median) - np.log(sample_median)))

			metrics_df[f"adjusted_Rm__{key}__{br}"] = metrics_df[raw_col] * cf
			metrics_df[f"norm_factor__{key}__{br}"] = cf
			norm_factors[br][key] = cf

	for br in BR_LABELS:
		for cond in condition_defs:
			key = cond["key"]
			if key == vehicle_key:
				continue
			metrics_df[f"DeltaRm__{key}__{br}"] = (
				metrics_df[f"adjusted_Rm__{key}__{br}"] - metrics_df[f"adjusted_Rm__{vehicle_key}__{br}"]
			)

	return global_vehicle_median, vehicle_medians, norm_factors


def compute_biorep_cv_by_condition(metrics_df, condition_defs):
	"""Compute biological-replicate CV from raw Rm across BR1/BR2/BR3 for each condition."""
	biocv_by_key = {}
	for cond in condition_defs:
		key = cond["key"]
		cols = [f"Rm_raw__{key}__{br}" for br in BR_LABELS]
		rm_mat = metrics_df[cols].to_numpy(dtype=float)
		valid_n = np.isfinite(rm_mat).sum(axis=1)

		with np.errstate(divide="ignore", invalid="ignore"):
			cv_vals = np.nanstd(rm_mat, axis=1, ddof=1) / np.nanmean(rm_mat, axis=1)

		cv_vals = np.where(valid_n >= len(BR_LABELS), cv_vals, np.nan)
		cv_vals = np.where(np.isfinite(cv_vals), cv_vals, np.nan)
		s = pd.Series(cv_vals, index=metrics_df.index, name=f"BioRep_CV_rawRm__{key}")
		biocv_by_key[key] = s
		metrics_df[s.name] = s

	return biocv_by_key


def plot_biocv_distribution(metrics_df, condition_defs, out_prefix, cv_cutoff=0.4):
	"""Plot BioRep CV distributions per condition."""
	plot_files = []
	for cond in condition_defs:
		key = cond["key"]
		cv_col = f"BioRep_CV_rawRm__{key}"
		cv_series = (
			metrics_df[cv_col]
			.replace([np.inf, -np.inf], np.nan)
			.dropna()
			.sort_values(ascending=False)
			.reset_index(drop=True)
		)
		if cv_series.empty:
			continue

		rank = cv_series.index + 1
		prop_30 = (cv_series <= 0.3).sum() / len(cv_series) * 100
		prop_40 = (cv_series <= 0.4).sum() / len(cv_series) * 100
		prop_cut = (cv_series <= cv_cutoff).sum() / len(cv_series) * 100
		ymax = max(float(cv_series.max()) * 1.05, 0.5, cv_cutoff + 0.08)

		plt.figure(figsize=(6, 4))
		plt.scatter(rank, cv_series.values, color="#4C6A92", alpha=0.72, s=20)
		if (not np.isclose(cv_cutoff, 0.3)) and (not np.isclose(cv_cutoff, 0.4)):
			plt.axhline(cv_cutoff, linestyle="--", linewidth=1, color="#111827")
		plt.ylim(0, ymax)
		plt.text(0.62 * len(cv_series), min(0.305, ymax * 0.95), f"<=30%: {prop_30:.1f}%", fontsize=10)
		plt.text(0.62 * len(cv_series), min(0.405, ymax * 0.95), f"<=40%: {prop_40:.1f}%", fontsize=10)
		if (not np.isclose(cv_cutoff, 0.3)) and (not np.isclose(cv_cutoff, 0.4)):
			plt.text(
				0.62 * len(cv_series),
				min(cv_cutoff + 0.005, ymax * 0.95),
				f"<={cv_cutoff*100:.0f}%: {prop_cut:.1f}%",
				fontsize=10,
			)
		plt.xlabel("Proteins sorted by BioRep CV (high->low)")
		plt.ylabel("BioRep CV of raw Rm")
		plt.title(f"BioRep CV distribution - {cond['name']}")
		plt.tight_layout()

		out_plot = f"{out_prefix}_BioRep_CV_{key}.tiff"
		plt.savefig(out_plot, dpi=300, format="tiff")
		plt.close()
		plot_files.append(out_plot)

	return plot_files


def build_bio_cv_qc_masks(metrics_df, condition_defs, biocv_by_key, cv_cutoff=0.4):
	"""
	Build biological-replicate CV-QC masks.
	Global rule: protein must pass CV cutoff in all conditions.
	"""
	bio_qc_masks = {}
	rows = []

	for cond in condition_defs:
		key = cond["key"]
		cond_cv = biocv_by_key[key].reindex(metrics_df.index)
		finite_mask = np.isfinite(cond_cv)
		mask = finite_mask & (cond_cv <= cv_cutoff)
		mask = mask.fillna(False)

		bio_qc_masks[key] = mask
		metrics_df[f"BioRep_CV_QC_pass__{key}"] = mask
		rows.append(
			{
				"Condition_name": cond["name"],
				"Condition_key": key,
				"CV_cutoff": float(cv_cutoff),
				"Pass_count": int(mask.sum()),
				"Total_count": int(mask.shape[0]),
				"Pass_percent": float(100.0 * mask.sum() / max(mask.shape[0], 1)),
				"Median_BioRep_CV_rawRm": float(np.nanmedian(cond_cv.values)),
			}
		)

	global_mask = pd.Series(True, index=metrics_df.index)
	for cond in condition_defs:
		key = cond["key"]
		global_mask &= bio_qc_masks[key].reindex(metrics_df.index, fill_value=False)

	metrics_df["Global_BioRep_CV_QC_pass"] = global_mask

	total_n = int(len(metrics_df.index))
	global_summary_df = pd.DataFrame(
		[
			{
				"Mask": "Global_BioRep_CV_QC_all_conditions",
				"Pass_count": int(global_mask.sum()),
				"Total_count": total_n,
				"Pass_percent": float(100.0 * global_mask.sum() / max(total_n, 1)),
			}
		]
	)

	return bio_qc_masks, pd.DataFrame(rows), global_mask, global_summary_df


def evaluate_global_normality(metrics_df, exp_conditions, out_prefix, global_cv_qc_mask):
	"""Global decision: all BioRep x Condition DeltaRm distributions must be normal."""
	all_normal = True
	normality_results = {}

	for br in BR_LABELS:
		for cond in exp_conditions:
			key = cond["key"]
			display_name = cond["name"]
			dcol = f"DeltaRm__{key}__{br}"
			panel_name = f"{br}_{display_name}"

			cv_mask = global_cv_qc_mask.reindex(metrics_df.index, fill_value=False)
			deltas_for_test = metrics_df.loc[cv_mask, dcol]

			plot_files = plot_delta_normality_diagnostics(deltas_for_test, panel_name, out_prefix)
			if plot_files:
				logger.info(f"{panel_name}: normality diagnostics saved: {', '.join(plot_files)}")

			is_normal, msg = assess_delta_normality(deltas_for_test)
			normality_results[(br, key)] = {"is_normal": is_normal, "message": msg}
			logger.info(f"{panel_name}: {msg}; normal={is_normal}; BioRep-CV-QC pass n={int(cv_mask.sum())}")
			if not is_normal:
				all_normal = False

	return all_normal, normality_results


def run_global_normal_branch(
	metrics_df,
	exp_conditions,
	out_prefix,
	global_cv_qc_mask,
	bio_cv_cutoff=0.4,
	delta_cutoff=0.1,
):
	"""
	Global-normal branch:
	1) robust p/FDR within each biological replicate and condition
	2) integrate significance across BR1-3:
	   all FDR<0.1 and >=1 FDR<0.05,
	   plus abs(DeltaRm)>threshold and same-sign DeltaRm in all BRs.
	"""
	mode_text = "global robust sigma (no CV bin)"

	for br in BR_LABELS:
		for cond in exp_conditions:
			key = cond["key"]
			display_name = cond["name"]
			dcol = f"DeltaRm__{key}__{br}"
			pcol = f"pval__{key}__{br}"
			fcol = f"FDR__{key}__{br}"

			valid = global_cv_qc_mask.reindex(metrics_df.index, fill_value=False) & np.isfinite(metrics_df[dcol])
			pvals_series = pd.Series(1.0, index=metrics_df.index)
			fdr_series = pd.Series(1.0, index=metrics_df.index)

			if int(valid.sum()) >= 3:
				valid_idx = metrics_df.index[valid]
				deltas_valid = metrics_df.loc[valid_idx, dcol]
				robust_pvals = robust_global_pvals(deltas_valid)

				pvals_valid = np.where(np.isnan(robust_pvals.values), 1.0, robust_pvals.values)
				_, fdrs_valid, _, _ = multipletests(pvals_valid, method="fdr_bh")
				pvals_series.loc[valid_idx] = pvals_valid
				fdr_series.loc[valid_idx] = fdrs_valid

			metrics_df[pcol] = pvals_series.values
			metrics_df[fcol] = fdr_series.values

			volcano = plot_delta_volcano(
				metrics_df.loc[valid, dcol],
				metrics_df.loc[valid, fcol],
				f"{br}_{display_name}",
				out_prefix,
				delta_cutoff=delta_cutoff,
				fdr_cutoff=0.05,
			)
			if volcano is not None:
				logger.info(f"{br}_{display_name}: volcano saved: {volcano}")

	aggregate = {}
	pass_cv_all = global_cv_qc_mask.reindex(metrics_df.index, fill_value=False).values

	for cond in exp_conditions:
		key = cond["key"]
		delta_mat = np.column_stack([metrics_df[f"DeltaRm__{key}__{br}"].values for br in BR_LABELS])
		pval_mat = np.column_stack([metrics_df[f"pval__{key}__{br}"].values for br in BR_LABELS])
		fdr_mat = np.column_stack([metrics_df[f"FDR__{key}__{br}"].values for br in BR_LABELS])
		pval_mat = np.where(np.isfinite(pval_mat), pval_mat, 1.0)
		fdr_mat = np.where(np.isfinite(fdr_mat), fdr_mat, 1.0)

		pass_delta = np.all(np.isfinite(delta_mat) & (np.abs(delta_mat) > delta_cutoff), axis=1)
		pass_same_sign = np.all(delta_mat > 0, axis=1) | np.all(delta_mat < 0, axis=1)
		pass_fdr_all_01 = np.all(fdr_mat < 0.1, axis=1)
		pass_fdr_one_005 = np.any(fdr_mat < 0.05, axis=1)
		sig = pass_cv_all & pass_delta & pass_same_sign & pass_fdr_all_01 & pass_fdr_one_005

		agg_df = pd.DataFrame(index=metrics_df.index)
		agg_df["DeltaRm_median"] = np.nanmedian(delta_mat, axis=1)
		agg_df["Pass_BioRep_CV_QC_all_conditions"] = pass_cv_all
		agg_df["Pass_absDelta_gt_threshold_all_reps"] = pass_delta
		agg_df["Pass_same_sign_all_reps"] = pass_same_sign
		agg_df["Pass_FDR_lt_0.1_all_reps"] = pass_fdr_all_01
		agg_df["Pass_at_least_one_FDR_lt_0.05"] = pass_fdr_one_005
		agg_df["Sig"] = sig
		agg_df["Sig_basis"] = (
			f"global-normal: BioRep-CV<={bio_cv_cutoff:.3g} all conditions + "
			f"{mode_text} + abs(DeltaRm)>{delta_cutoff:.3g} all reps + "
			"same-sign DeltaRm all reps + all FDR<0.1 + >=1 rep FDR<0.05"
		)

		for i, br in enumerate(BR_LABELS):
			agg_df[f"DeltaRm_{br}"] = delta_mat[:, i]
			agg_df[f"pval_{br}"] = pval_mat[:, i]
			agg_df[f"FDR_{br}"] = fdr_mat[:, i]

		aggregate[key] = agg_df

	return aggregate


def run_global_non_normal_branch(
	metrics_df,
	exp_conditions,
	vehicle_key,
	global_cv_qc_mask,
	paired_method="ttest",
	bio_cv_cutoff=0.4,
	delta_cutoff=0.1,
):
	"""
	Global-non-normal branch:
	paired test across BR1/BR2/BR3 adjusted Rm values for each protein.
	BH correction is applied only to proteins that were actually tested.
	"""
	aggregate = {}
	pass_cv_all = global_cv_qc_mask.reindex(metrics_df.index, fill_value=False).values
	protein_index = metrics_df.index

	for cond in exp_conditions:
		key = cond["key"]
		delta_mat = np.column_stack([metrics_df[f"DeltaRm__{key}__{br}"].values for br in BR_LABELS])

		pvals = np.ones(len(protein_index), dtype=float)
		tested_mask = np.zeros(len(protein_index), dtype=bool)
		for i, pid in enumerate(protein_index):
			if not pass_cv_all[i]:
				pvals[i] = 1.0
				continue

			x = np.array([metrics_df.at[pid, f"adjusted_Rm__{key}__{br}"] for br in BR_LABELS], dtype=float)
			y = np.array([metrics_df.at[pid, f"adjusted_Rm__{vehicle_key}__{br}"] for br in BR_LABELS], dtype=float)

			valid = np.isfinite(x) & np.isfinite(y)
			if valid.sum() < 3:
				pvals[i] = 1.0
				continue

			tested_mask[i] = True

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

		fdrs = np.ones(len(protein_index), dtype=float)
		if tested_mask.any():
			_, fdrs_tested, _, _ = multipletests(pvals[tested_mask], method="fdr_bh")
			fdrs[tested_mask] = fdrs_tested

		pass_delta = np.all(np.isfinite(delta_mat) & (np.abs(delta_mat) > delta_cutoff), axis=1)
		pass_same_sign = np.all(delta_mat > 0, axis=1) | np.all(delta_mat < 0, axis=1)
		sig = pass_cv_all & pass_delta & pass_same_sign & (fdrs < 0.05)

		delta_mean = np.nanmean(delta_mat, axis=1)
		delta_mean = np.where(pass_cv_all, delta_mean, np.nan)

		agg_df = pd.DataFrame(index=protein_index)
		agg_df["DeltaRm_mean"] = delta_mean
		agg_df["DeltaRm_display"] = agg_df["DeltaRm_mean"]
		agg_df["DeltaRm_median"] = np.nanmedian(delta_mat, axis=1)
		agg_df["Pass_BioRep_CV_QC_all_conditions"] = pass_cv_all
		agg_df["Pass_absDelta_gt_threshold_all_reps"] = pass_delta
		agg_df["Pass_same_sign_all_reps"] = pass_same_sign
		agg_df["Paired_test"] = paired_method
		agg_df["pval"] = pvals
		agg_df["FDR"] = fdrs
		agg_df["Sig"] = sig
		agg_df["Sig_basis"] = (
			f"global-non-normal: BioRep-CV<={bio_cv_cutoff:.3g} all conditions + "
			f"paired-{paired_method} FDR<0.05 + abs(DeltaRm)>{delta_cutoff:.3g} all reps + "
			"same-sign DeltaRm all reps; DeltaRm_display=mean"
		)
		for i, br in enumerate(BR_LABELS):
			agg_df[f"DeltaRm_{br}"] = delta_mat[:, i]

		aggregate[key] = agg_df

	return aggregate


def build_full_results(metrics_df, exp_conditions, aggregate_by_cond, normality_results, analysis_mode):
	full_df = pd.DataFrame(index=metrics_df.index)
	full_df["analysis_mode_global"] = analysis_mode

	for (br, cond_key), info in normality_results.items():
		full_df[f"normality__{br}__{cond_key}__is_normal"] = bool(info["is_normal"])
		full_df[f"normality__{br}__{cond_key}__message"] = info["message"]

	for col in metrics_df.columns:
		full_df[col] = metrics_df[col]

	for cond in exp_conditions:
		key = cond["key"]
		agg = aggregate_by_cond[key]
		for col in agg.columns:
			full_df[f"AGG__{key}__{col}"] = agg[col].reindex(full_df.index)

	return full_df


def main():
	logger.info(
		"Branch02 BioRep-only mode: one CSV file with 3 biological replicates encoded by "
		"TMT pair mapping (BR1 129/126, BR2 130/127, BR3 131/128)."
	)

	raw_path = input("Please provide path of ONE CSV file for BioRep-only analysis: ").strip()
	file_path = parse_single_csv_path(raw_path)

	if not os.path.exists(file_path):
		raise FileNotFoundError(file_path)

	logger.info("Step 1/3: parse input, build overlap plot, and compute raw Rm matrix.")
	parsed = process_single_file(file_path)
	metrics_df = parsed["metrics_df"]
	condition_defs = parsed["condition_defs"]
	summary_row = parsed["summary"]

	if metrics_df.empty:
		raise ValueError("No proteins passed quant completeness filter.")

	vehicle_key = pick_vehicle_key(condition_defs)
	exp_conditions = [c for c in condition_defs if c["key"] != vehicle_key]
	if not exp_conditions:
		raise ValueError("No experimental condition found after removing Vehicle.")

	global_vehicle_median, vehicle_medians, norm_factors = apply_median_of_medians_normalization(
		metrics_df, condition_defs, vehicle_key
	)
	logger.info(f"Global vehicle median (median of BR vehicle medians): {global_vehicle_median:.6g}")

	for br in BR_LABELS:
		fac_str = ", ".join([f"{c['key']}={norm_factors[br][c['key']]:.6g}" for c in condition_defs])
		logger.info(f"{br} normalization factors: {fac_str}")

	output_prefix = os.path.splitext(file_path)[0] + "_BioRepOnly"

	logger.info(
		"Step 2/3: build BioRep CV-QC preview first (default cutoff=0.4), "
		"then ask for final threshold."
	)
	preview_biocv_by_key = compute_biorep_cv_by_condition(metrics_df, condition_defs)
	preview_bio_cutoff_df = summarize_global_bio_cv_pass_across_cutoffs(
		preview_biocv_by_key,
		condition_defs,
		BIO_CV_PREVIEW_CANDIDATES,
	)
	preview_plot_files = plot_biocv_distribution(
		metrics_df,
		condition_defs,
		output_prefix,
		cv_cutoff=DEFAULT_BIO_CV_QC_CUTOFF,
	)

	logger.info("CV-QC preview table (global BioRep CV across all conditions):")
	logger.info("\n" + format_preview_table(preview_bio_cutoff_df, float_cols=["Bio_CV_cutoff", "Global_pass_percent"]))
	if preview_plot_files:
		logger.info(f"BioRep CV-QC preview plots saved: {', '.join(preview_plot_files)}")

	input("CV-QC preview generated. Please review tables/plots, then press Enter to input final CV-QC cutoff: ")

	bio_cv_qc_cutoff = ask_float_with_default(
		f"Enter biological-replicate CV-QC cutoff for final analysis (default {DEFAULT_BIO_CV_QC_CUTOFF}): ",
		default=DEFAULT_BIO_CV_QC_CUTOFF,
		min_value=0.0,
	)
	logger.info(f"Final BioRep CV-QC cutoff: {bio_cv_qc_cutoff:.3g}")

	_bio_cv_qc_masks, bio_cv_qc_summary_df, global_bio_cv_qc_mask, global_cv_qc_summary_df = build_bio_cv_qc_masks(
		metrics_df,
		condition_defs,
		biocv_by_key=preview_biocv_by_key,
		cv_cutoff=bio_cv_qc_cutoff,
	)
	bio_cv_plot_files = plot_biocv_distribution(
		metrics_df,
		condition_defs,
		output_prefix,
		cv_cutoff=bio_cv_qc_cutoff,
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

	# Branch03-style: differential analysis is performed only on the CV-QC pass subset.
	analysis_metrics_df = metrics_df.loc[global_bio_cv_qc_mask].copy()
	if analysis_metrics_df.empty:
		raise ValueError("No proteins passed global BioRep CV-QC. Differential analysis cannot proceed.")
	analysis_cv_mask = pd.Series(True, index=analysis_metrics_df.index)
	logger.info(
		"Differential analysis input after global BioRep CV-QC subset: "
		f"{analysis_metrics_df.shape[0]} proteins"
	)

	logger.info("Step 3/3: global normality decision and differential analysis.")
	global_normal, normality_results = evaluate_global_normality(
		analysis_metrics_df,
		exp_conditions,
		output_prefix,
		global_cv_qc_mask=analysis_cv_mask,
	)

	if global_normal:
		logger.info("Global normal -> run global robust sigma (no CV binning)")
		aggregate_by_cond = run_global_normal_branch(
			analysis_metrics_df,
			exp_conditions,
			output_prefix,
			global_cv_qc_mask=analysis_cv_mask,
			bio_cv_cutoff=bio_cv_qc_cutoff,
			delta_cutoff=DEFAULT_DELTA_CUTOFF,
		)
		analysis_mode = "global_normal_robust_global"
	else:
		method_in = ask_choice_with_default(
			"Global non-normal detected. Choose paired test [ttest/wilcoxon], Enter for default ttest: ",
			allowed_choices={"ttest", "wilcoxon"},
			default="ttest",
		)
		logger.info(f"Global non-normal -> run paired {method_in} across biological replicates")
		aggregate_by_cond = run_global_non_normal_branch(
			analysis_metrics_df,
			exp_conditions,
			vehicle_key,
			global_cv_qc_mask=analysis_cv_mask,
			paired_method=method_in,
			bio_cv_cutoff=bio_cv_qc_cutoff,
			delta_cutoff=DEFAULT_DELTA_CUTOFF,
		)
		analysis_mode = f"global_non_normal_paired_{method_in}"

	full_df = build_full_results(analysis_metrics_df, exp_conditions, aggregate_by_cond, normality_results, analysis_mode)

	full_outfile = output_prefix + "_full_Rm_results.csv"
	full_df.reset_index().rename(columns={"index": "ProteinID"}).to_csv(full_outfile, index=False)
	logger.info(f"Full results saved: {full_outfile}")

	sig_cols = [f"AGG__{c['key']}__Sig" for c in exp_conditions]
	sig_mask = full_df[sig_cols].fillna(False).any(axis=1)
	diff_df = full_df.loc[sig_mask].copy()
	diff_outfile = output_prefix + "_Differential_proteins.csv"
	diff_df.reset_index().rename(columns={"index": "ProteinID"}).to_csv(diff_outfile, index=False)
	logger.info(f"Differential proteins saved: {diff_outfile}")

	summary_df = pd.DataFrame([summary_row])
	summary_df["Quant_complete_intersection_count"] = int(metrics_df.shape[0])
	summary_df["BioRep_CV_cutoff_used"] = float(bio_cv_qc_cutoff)
	summary_df["Global_BioRep_CV_QC_pass_count"] = int(global_bio_cv_qc_mask.sum())
	summary_df["Differential_analysis_input_count"] = int(analysis_metrics_df.shape[0])
	summary_df["Analysis_mode"] = analysis_mode
	summary_outfile = output_prefix + "_BioRep_summary.csv"
	summary_df.to_csv(summary_outfile, index=False)
	logger.info(f"BioRep summary saved: {summary_outfile}")

	norm_rows = []
	for br in BR_LABELS:
		for cond in condition_defs:
			key = cond["key"]
			raw_col = f"Rm_raw__{key}__{br}"
			norm_rows.append(
				{
					"BioRep": br,
					"Condition_name": cond["name"],
					"Condition_key": key,
					"Condition_median_raw_Rm": float(metrics_df[raw_col].median()),
					"Vehicle_median_this_BioRep": float(vehicle_medians[br]),
					"Global_vehicle_median_of_medians": float(global_vehicle_median),
					"Normalization_factor": float(norm_factors[br][key]),
				}
			)

	norm_df = pd.DataFrame(norm_rows)
	norm_outfile = output_prefix + "_normalization_factors.csv"
	norm_df.to_csv(norm_outfile, index=False)
	logger.info(f"Normalization factors saved: {norm_outfile}")

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
		sig_summary_rows.append(
			{
				"Condition_name": cond["name"],
				"Condition_key": key,
				"Significant_protein_count": sig_count,
			}
		)

	sig_summary_df = pd.DataFrame(sig_summary_rows)
	sig_summary_outfile = output_prefix + "_significance_summary.csv"
	sig_summary_df.to_csv(sig_summary_outfile, index=False)
	logger.info(f"Significance summary saved: {sig_summary_outfile}")

	logger.info("Refined-TPP plus BioRep-only processing completed successfully.")


if __name__ == "__main__":
	main()
