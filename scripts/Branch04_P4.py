#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Refined-TPP plus - No Biological Replicate / No Technical Replicate Pipeline
-----------------------------------------------------------------------------
Created on: Fri Nov 7 11:59:36 2025
Authors: Haoru Song, Chenxin Li, Bin Fu
Affiliation: Professor Haojie Lu's Lab, Fudan University, Shanghai, China.

Description
-----------
This script is dedicated to the design without technical replicates and without
biological replicates.

Design assumptions in this mode:
1) Each condition uses one TMT reporter pair in one run:
   - Pair 129/126, 130/127, 131/128
2) One run can host up to 3 conditions.
3) If condition count > 3, additional run CSV files can be provided.
4) No geometric averaging is used (because there are no technical replicates).
5) No CV-based QC is used (because there are no replicates).

Core analysis logic:
1) Quantification overview and overlap:
   - Count quantified proteins for each condition (2-channel condition, any
	 non-zero signal in the pair).
   - Count condition overlap (Venn for 2-3 conditions, UpSet for >=4 if
	 available).
   - Report condition-wise complete quantification counts (both channels
	 quantified).
2) Global quantification completeness QC:
   proteins must be quantified in all required channels across all configured
   conditions, i.e. all 2N channels (N = condition count).
3) Rm is computed directly from each condition's reporter pair ratio.
4) Normalization is done directly in log-space by aligning each condition to the
   global Vehicle median in log-space.
5) DeltaRm = adjusted_Rm(condition) - adjusted_Rm(Vehicle).
6) Normality is tested for each condition's DeltaRm distribution.
   - If normal: global robust-sigma p/FDR (no CV binning).
   - If non-normal: no inferential statistics; 
     significance is defined by empirical outliers (top 5% and bottom 5% DeltaRm).
	 Thresholds DeltaRm < -0.1 and DeltaRm > 0.1 are retained for visualization reference only.

Unique peptide filter:
- unique>=2 filtering can be enabled/disabled by user.
- If unique-like column is unavailable, the filter is skipped automatically.

Notes
-----
- If multiple runs are used and Vehicle is not present in some runs, DeltaRm for
  those runs is referenced to the global Vehicle median only. This is a
  descriptive workaround rather than rigorous inter-run batch correction.

Dependencies
------------
- Python 3.8+
- pandas, numpy, matplotlib, scipy, statsmodels
- Optional for overlap plot: matplotlib-venn, upsetplot

Installation
------------
	pip install pandas numpy matplotlib scipy statsmodels matplotlib-venn upsetplot
"""

import os
import re
from collections import defaultdict

import warnings

warnings.filterwarnings("ignore", category=FutureWarning, module="matplotlib")

import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logging.getLogger("numexpr").setLevel(logging.WARNING)
logger = logging.getLogger()

try:
	import matplotlib.pyplot as plt
	import numpy as np
	import pandas as pd
	from scipy import stats
	from scipy.stats import norm
	from statsmodels.stats.multitest import multipletests
except Exception as e:
	raise ImportError(
		f"{e}\nPlease install required packages: "
		"pip install pandas numpy matplotlib scipy statsmodels"
	)

HAVE_VENN = True
HAVE_UPSET = True

try:
	from matplotlib_venn import venn2, venn3
except Exception:
	HAVE_VENN = False

try:
	from upsetplot import from_contents, UpSet
except Exception:
	HAVE_UPSET = False


# ------------------------------
# Helper functions
# ------------------------------
def normalize_path(path):
	return os.path.normpath(os.path.abspath(path))


def parse_csv_paths(raw_text):
	"""Parse input paths split by ';' or ','."""
	parts = [p.strip().strip('"').strip("'") for p in re.split(r"[;,]", raw_text) if p.strip()]
	return [normalize_path(p) for p in parts]


def make_safe_key(text):
	"""Convert display name to an ASCII-safe key for columns and filenames."""
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


def detect_global_reporter_cols(df, run_name):
	"""
	Detect one reporter column for each TMT channel (126..131).
	If multiple candidates exist for one reporter channel, ask user to choose.
	"""
	reporter_hits = defaultdict(list)
	pattern = re.compile(r"(126|127|128|129|130|131)")

	for col in df.columns:
		m = pattern.search(str(col))
		if m:
			reporter_hits[int(m.group(1))].append(col)

	selected = {}
	for reporter in [126, 127, 128, 129, 130, 131]:
		candidates = reporter_hits.get(reporter, [])

		if len(candidates) == 1:
			selected[reporter] = candidates[0]
			continue

		if len(candidates) == 0:
			selected[reporter] = None
			logger.warning(f"[{run_name}] TMT-{reporter} column is not detected.")
			continue

		logger.info(f"[{run_name}] Multiple columns detected for TMT-{reporter}:")
		for i, col in enumerate(candidates, start=1):
			logger.info(f"  {i}. {col}")

		raw_choice = input(f"[{run_name}] Select column index for TMT-{reporter} [default 1]: ").strip()
		choice_idx = 1
		if raw_choice != "":
			try:
				choice_idx = int(raw_choice)
			except ValueError:
				logger.warning(f"[{run_name}] Invalid index '{raw_choice}', fallback to 1.")
				choice_idx = 1

		if choice_idx < 1 or choice_idx > len(candidates):
			logger.warning(f"[{run_name}] Index out of range, fallback to 1.")
			choice_idx = 1

		selected[reporter] = candidates[choice_idx - 1]

	return selected


def draw_quant_overlap(per_condition_sets, condition_labels, out_prefix):
	"""Draw overlap plot for quantified proteins across conditions."""
	n = len(condition_labels)
	if n < 2:
		return None

	out_file = None
	if n < 4 and HAVE_VENN:
		plt.figure(figsize=(6, 6))
		if n == 2:
			venn2([per_condition_sets[label] for label in condition_labels], set_labels=condition_labels)
		elif n == 3:
			venn3([per_condition_sets[label] for label in condition_labels], set_labels=condition_labels)
		plt.title("Quantified proteins overlap (2-channel conditions)")
		out_file = f"{out_prefix}_quantified_overlap_venn.tiff"
		plt.savefig(out_file, dpi=300, bbox_inches="tight", format="tiff")
		plt.close()
	elif HAVE_UPSET:
		upset_data = from_contents({label: per_condition_sets[label] for label in condition_labels})
		plt.figure(figsize=(8, 6))
		UpSet(upset_data, show_counts=True).plot()
		out_file = f"{out_prefix}_quantified_overlap_upset.tiff"
		plt.savefig(out_file, dpi=300, bbox_inches="tight", format="tiff")
		plt.close()
	else:
		logger.warning(
			"Overlap plot skipped: matplotlib-venn/upsetplot are not available. "
			"Install optional packages if overlap figure is needed."
		)

	return out_file


def robust_global_sigma_pvals(deltas):
	"""Compute global robust-sigma p-values from one DeltaRm distribution."""
	series = pd.Series(deltas).replace([np.inf, -np.inf], np.nan)
	valid = series.dropna()

	aligned = pd.Series(1.0, index=series.index)
	if valid.empty:
		return aligned

	vals = valid.values
	p15 = np.percentile(vals, 15.87)
	p50 = np.percentile(vals, 50.0)
	p84 = np.percentile(vals, 84.13)

	left_sigma = p50 - p15
	right_sigma = p84 - p50
	left_sigma = left_sigma if left_sigma > 0 else np.nan
	right_sigma = right_sigma if right_sigma > 0 else np.nan

	pvals = pd.Series(1.0, index=valid.index)
	for idx, val in valid.items():
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
		pvals.loc[idx] = pv

	aligned.loc[pvals.index] = pvals.values
	return aligned


def assess_delta_normality(delta_series):
	"""Assess normality for one DeltaRm series and return (is_normal, message)."""
	arr = pd.Series(delta_series).replace([np.inf, -np.inf], np.nan).dropna().values
	n = len(arr)

	if n < 3:
		return False, f"n={n}: too few proteins for normality test, treat as non-normal empirical-outlier mode"

	if n <= 50:
		stat, p = stats.shapiro(arr)
		is_normal = p >= 0.05
		return is_normal, f"Shapiro-Wilk W={stat:.4f}, p={p:.4g}"

	if n <= 300:
		mean_val = np.mean(arr)
		std_val = np.std(arr, ddof=1)
		if not np.isfinite(std_val) or std_val <= 0:
			return False, "KS skipped because std<=0, treat as non-normal empirical-outlier mode"
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
	"""Plot volcano for normal branch DeltaRm analysis."""
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


def plot_delta_distribution_threshold(delta_series, sample_name, out_prefix):
	"""Plot DeltaRm distribution and annotate DeltaRm<-0.1 / DeltaRm>0.1 counts."""
	arr = pd.Series(delta_series).replace([np.inf, -np.inf], np.nan).dropna().values
	if len(arr) == 0:
		return None

	safe_name = make_safe_key(sample_name)
	n_low = int((arr < -0.1).sum())
	n_high = int((arr > 0.1).sum())
	n_total = int(len(arr))

	plt.figure(figsize=(6.6, 4.2))
	plt.hist(arr, bins=35, color="#B0BEC5", alpha=0.85, edgecolor="white", linewidth=0.4)
	plt.axvline(-0.1, color="#2C7FB8", linestyle="--", linewidth=1)
	plt.axvline(0.1, color="#D7301F", linestyle="--", linewidth=1)
	plt.axvline(float(np.median(arr)), color="black", linestyle="-", linewidth=0.8)

	ymax = plt.ylim()[1]
	plt.text(-0.1, ymax * 0.92, f"<-0.1: {n_low}/{n_total}", color="#2C7FB8", fontsize=9, ha="right")
	plt.text(0.1, ymax * 0.85, f">0.1: {n_high}/{n_total}", color="#D7301F", fontsize=9, ha="left")
	plt.xlabel("DeltaRm")
	plt.ylabel("Count")
	plt.title(f"DeltaRm distribution (non-normal descriptive branch) - {sample_name}")
	plt.tight_layout()

	out_file = f"{out_prefix}_DeltaRm_distribution_threshold_{safe_name}.tiff"
	plt.savefig(out_file, dpi=300, format="tiff")
	plt.close()
	return out_file


def empirical_tail_masks(delta_series, tail_frac=0.05):
	"""Return low/high tail masks and tail size for empirical outlier labeling."""
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


def process_run_no_replicates(file_path, run_name, used_condition_keys, require_unique_ge2=True):
	"""
	Process one run CSV in no-replicate mode.
	One condition is mapped to one reporter pair: 129/126, 130/127, or 131/128.
	"""
	df = pd.read_csv(file_path)
	protein_col = detect_protein_col(df)
	unique_col = detect_unique_col(df)

	protein_ids = df[protein_col].astype(str).str.strip()
	non_empty_protein = protein_ids != ""
	if unique_col is None:
		unique_mask = pd.Series(True, index=df.index)
		unique_filter_effective_ge2 = False
		logger.warning(f"[{run_name}] Unique-like column not found; unique>=2 filter is skipped in this run.")
	else:
		if require_unique_ge2:
			unique_mask = df[unique_col].fillna(0) >= 2
			unique_filter_effective_ge2 = True
		else:
			unique_mask = pd.Series(True, index=df.index)
			unique_filter_effective_ge2 = False

	reporter_cols = detect_global_reporter_cols(df, run_name)
	pair_defs = [(129, 126), (130, 127), (131, 128)]

	run_conditions = []
	for num, den in pair_defs:
		num_col = reporter_cols.get(num)
		den_col = reporter_cols.get(den)
		if num_col is None or den_col is None:
			logger.info(f"[{run_name}] Pair {num}/{den} is unavailable due to missing channel column.")
			continue

		cond_name = input(f"[{run_name}] Enter condition name for pair {num}/{den} (blank to skip): ").strip()
		if cond_name == "":
			continue

		base_key = make_safe_key(cond_name)
		key = base_key
		suffix = 2
		while key in used_condition_keys:
			key = f"{base_key}_{suffix}"
			suffix += 1
		if key != base_key:
			logger.warning(f"[{run_name}] Condition key '{base_key}' duplicated; renamed to '{key}'.")

		used_condition_keys.add(key)
		run_conditions.append(
			{
				"name": cond_name,
				"key": key,
				"num": num,
				"den": den,
				"num_col": num_col,
				"den_col": den_col,
				"pair_label": f"{num}/{den}",
				"run_name": run_name,
				"file_path": file_path,
			}
		)

	if not run_conditions:
		raise ValueError(f"[{run_name}] No condition selected. At least one reporter pair is required.")

	quant_sets_by_condition = {}
	complete_sets_by_condition = {}
	condition_any_counts = {}
	condition_complete_counts = {}

	for cond in run_conditions:
		num_signal = df[cond["num_col"]].fillna(0) != 0
		den_signal = df[cond["den_col"]].fillna(0) != 0

		any_signal_mask = (num_signal | den_signal) & unique_mask & non_empty_protein
		complete_mask = (num_signal & den_signal) & unique_mask & non_empty_protein

		any_set = set(protein_ids[any_signal_mask])
		complete_set = set(protein_ids[complete_mask])

		quant_sets_by_condition[cond["key"]] = any_set
		complete_sets_by_condition[cond["key"]] = complete_set
		condition_any_counts[cond["key"]] = int(len(any_set))
		condition_complete_counts[cond["key"]] = int(len(complete_set))

	quant_union_set = set().union(*quant_sets_by_condition.values()) if quant_sets_by_condition else set()
	quant_intersection_set = set.intersection(*quant_sets_by_condition.values()) if quant_sets_by_condition else set()

	strict_mask = unique_mask & non_empty_protein
	for cond in run_conditions:
		strict_mask &= (df[cond["num_col"]].fillna(0) != 0)
		strict_mask &= (df[cond["den_col"]].fillna(0) != 0)

	df_final = df.loc[strict_mask].copy()
	df_final["ProteinID"] = df_final[protein_col].astype(str).str.strip()
	df_final = df_final[df_final["ProteinID"] != ""].copy()

	metrics = pd.DataFrame({"ProteinID": df_final["ProteinID"]})
	for cond in run_conditions:
		key = cond["key"]
		raw_col = f"raw_Rm__{key}"
		log_col = f"log_Rm__{key}"

		rm = (df_final[cond["num_col"]] / df_final[cond["den_col"]]).replace([np.inf, -np.inf], np.nan)
		rm = rm.where(rm > 0, np.nan)

		metrics[raw_col] = rm
		with np.errstate(divide="ignore", invalid="ignore"):
			metrics[log_col] = np.log(rm)

	if metrics["ProteinID"].duplicated().any():
		numeric_cols = [c for c in metrics.columns if c != "ProteinID"]
		metrics = metrics.groupby("ProteinID", as_index=False)[numeric_cols].median()

	metrics = metrics.set_index("ProteinID").sort_index()

	mapping_str = "; ".join([f"{c['name']}:{c['pair_label']}" for c in run_conditions])
	logger.info(f"[{run_name}] strict-filtered proteins: {metrics.shape[0]}; selected conditions: {mapping_str}")

	summary = {
		"Run": run_name,
		"File": file_path,
		"Unique_filter_requested_ge2": bool(require_unique_ge2),
		"Unique_filter_effective_ge2": bool(unique_filter_effective_ge2),
		"Protein_column": protein_col,
		"Unique_column": unique_col,
		"Run_quantified_union_count_any_signal": int(len(quant_union_set)),
		"Run_quantified_intersection_count_all_conditions_any_signal": int(len(quant_intersection_set)),
		"Run_complete_all_conditions_count": int(metrics.shape[0]),
		"Strict_filtered_count": int(metrics.shape[0]),
		"Condition_count_in_run": int(len(run_conditions)),
		"Channel_count_in_run": int(2 * len(run_conditions)),
	}

	return {
		"metrics_df": metrics,
		"conditions": run_conditions,
		"quant_sets_by_condition": quant_sets_by_condition,
		"complete_sets_by_condition": complete_sets_by_condition,
		"condition_any_counts": condition_any_counts,
		"condition_complete_counts": condition_complete_counts,
		"quant_union_set": quant_union_set,
		"quant_intersection_set": quant_intersection_set,
		"summary": summary,
	}


def pick_vehicle_condition_key(all_conditions):
	"""Pick Vehicle condition from configured conditions."""
	for cond in all_conditions:
		if cond["name"].strip().lower() == "vehicle" or cond["key"].strip().lower() == "vehicle":
			return cond["key"]

	user_text = input("Enter Vehicle condition name (display name or key): ").strip().lower()
	for cond in all_conditions:
		if cond["name"].strip().lower() == user_text or cond["key"].strip().lower() == user_text:
			return cond["key"]

	raise ValueError("Vehicle condition was not found in condition mapping.")


# ------------------------------
# Main workflow
# ------------------------------
def main():
	raw_paths = input(
		"Please provide path(s) of run CSV file(s), separated by ';' or ',' (single file is allowed): "
	).strip()
	file_paths = parse_csv_paths(raw_paths)

	if len(file_paths) == 0:
		raise ValueError("At least one CSV file path is required.")

	for fp in file_paths:
		if not os.path.exists(fp):
			raise FileNotFoundError(fp)

	unique_choice = input(
		"Apply unique peptide filter (>=2) for quantification overview and strict filtering? [YES/NO, default YES]: "
	).strip().upper()
	require_unique_ge2 = unique_choice != "NO"
	logger.info(f"Unique filter setting: {'>=2 enabled' if require_unique_ge2 else 'disabled'}")

	used_condition_keys = set()
	run_entries = []
	all_conditions = []
	summary_rows = []

	for i, fp in enumerate(file_paths, start=1):
		run_name = f"Run{i}"
		run_result = process_run_no_replicates(
			fp,
			run_name,
			used_condition_keys,
			require_unique_ge2=require_unique_ge2,
		)
		run_entries.append(run_result)
		all_conditions.extend(run_result["conditions"])
		summary_rows.append(run_result["summary"])

	if not all_conditions:
		raise ValueError("No conditions were configured across runs.")

	base_prefix = os.path.splitext(file_paths[0])[0] + "_NoRep"

	all_condition_quant_sets = {}
	all_condition_complete_sets = {}
	all_condition_labels = []
	all_condition_sets_for_plot = {}

	for entry in run_entries:
		for cond in entry["conditions"]:
			key = cond["key"]
			label = f"{cond['name']} ({cond['run_name']}, {key})"
			all_condition_quant_sets[key] = entry["quant_sets_by_condition"][key]
			all_condition_complete_sets[key] = entry["complete_sets_by_condition"][key]
			all_condition_labels.append(label)
			all_condition_sets_for_plot[label] = entry["quant_sets_by_condition"][key]

	global_quant_union = set().union(*all_condition_quant_sets.values()) if all_condition_quant_sets else set()
	global_quant_intersection = set.intersection(*all_condition_quant_sets.values()) if all_condition_quant_sets else set()

	logger.info(
		f"Quantification overview (any signal per 2-channel condition): "
		f"union={len(global_quant_union)}, intersection(all conditions)={len(global_quant_intersection)}"
	)

	for cond in all_conditions:
		key = cond["key"]
		quant_n = len(all_condition_quant_sets[key])
		complete_n = len(all_condition_complete_sets[key])
		logger.info(
			f"Condition {cond['name']} ({cond['run_name']}, pair {cond['pair_label']}): "
			f"quantified(any signal)={quant_n}, complete(2 channels)={complete_n}"
		)

	global_overlap_plot = draw_quant_overlap(all_condition_sets_for_plot, all_condition_labels, base_prefix)
	if global_overlap_plot is not None:
		logger.info(f"Global condition overlap plot saved: {global_overlap_plot}")

	vehicle_key = pick_vehicle_condition_key(all_conditions)
	exp_conditions = [c for c in all_conditions if c["key"] != vehicle_key]
	if not exp_conditions:
		raise ValueError("No experimental condition found after removing Vehicle.")

	protein_sets = [set(entry["metrics_df"].index) for entry in run_entries]
	common_proteins = sorted(set.intersection(*protein_sets))
	if len(common_proteins) == 0:
		raise ValueError(
			"No common proteins found across runs after strict 2N-channel quantification completeness filtering."
		)

	for entry in run_entries:
		entry["metrics_df"] = entry["metrics_df"].loc[common_proteins].copy()
		entry["condition_keys"] = [c["key"] for c in entry["conditions"]]

	logger.info(
		f"Global quantification completeness passed proteins: {len(common_proteins)} "
		f"(conditions={len(all_conditions)}, required channels=2N={2 * len(all_conditions)})."
	)
	logger.info(
		"Downstream Rm/DeltaRm calculation uses only proteins passing global 2N complete quantification "
		"across all configured conditions."
	)

	vehicle_log_medians = []
	for entry in run_entries:
		if vehicle_key in entry["condition_keys"]:
			v_col = f"log_Rm__{vehicle_key}"
			v_med = float(entry["metrics_df"][v_col].median())
			if np.isfinite(v_med):
				vehicle_log_medians.append(v_med)

	if not vehicle_log_medians:
		raise ValueError("Vehicle medians are unavailable after filtering; cannot normalize in log-space.")

	global_vehicle_log_median = float(np.median(np.array(vehicle_log_medians, dtype=float)))
	global_vehicle_rm_ref = float(np.exp(global_vehicle_log_median))
	logger.info(
		f"Global Vehicle reference in log-space: {global_vehicle_log_median:.6g}; "
		f"equivalent Rm reference: {global_vehicle_rm_ref:.6g}"
	)

	norm_rows = []
	for entry in run_entries:
		run_df = entry["metrics_df"]

		for cond in entry["conditions"]:
			key = cond["key"]
			log_col = f"log_Rm__{key}"
			adj_col = f"adjusted_Rm__{key}"

			cond_log_median = float(run_df[log_col].median())
			if np.isfinite(cond_log_median):
				norm_factor = float(np.exp(global_vehicle_log_median - cond_log_median))
				adjusted_log = run_df[log_col] - cond_log_median + global_vehicle_log_median
				adjusted_rm = np.exp(adjusted_log)
			else:
				norm_factor = np.nan
				adjusted_rm = pd.Series(np.nan, index=run_df.index)

			run_df[adj_col] = adjusted_rm
			run_df[f"norm_factor__{key}"] = norm_factor

			norm_rows.append(
				{
					"Run": cond["run_name"],
					"Condition_name": cond["name"],
					"Condition_key": key,
					"Condition_log_median": cond_log_median,
					"Global_vehicle_log_median": global_vehicle_log_median,
					"Normalization_factor": norm_factor,
				}
			)

		if vehicle_key in entry["condition_keys"]:
			vehicle_profile = run_df[f"adjusted_Rm__{vehicle_key}"]
			entry["vehicle_reference_mode"] = "in_run_vehicle_profile"
		else:
			vehicle_profile = pd.Series(global_vehicle_rm_ref, index=run_df.index)
			entry["vehicle_reference_mode"] = "global_vehicle_median_constant"
			logger.warning(
				f"[{entry['summary']['Run']}] Vehicle is absent in this run; "
				"DeltaRm is referenced to global Vehicle median only."
			)

		for cond in entry["conditions"]:
			key = cond["key"]
			if key == vehicle_key:
				continue
			run_df[f"DeltaRm__{key}"] = run_df[f"adjusted_Rm__{key}"] - vehicle_profile

	run_lookup = {}
	for entry in run_entries:
		for cond in entry["conditions"]:
			run_lookup[cond["key"]] = entry

	full_df = pd.DataFrame(index=common_proteins)
	full_df["analysis_design"] = "no_biorep_no_techrep"

	for cond in all_conditions:
		key = cond["key"]
		run_df = run_lookup[key]["metrics_df"]

		full_df[f"raw_Rm__{key}"] = run_df[f"raw_Rm__{key}"].reindex(common_proteins)
		full_df[f"adjusted_Rm__{key}"] = run_df[f"adjusted_Rm__{key}"].reindex(common_proteins)
		full_df[f"norm_factor__{key}"] = run_df[f"norm_factor__{key}"].reindex(common_proteins)

		if key == vehicle_key:
			full_df[f"DeltaRm__{key}"] = np.nan
		else:
			full_df[f"DeltaRm__{key}"] = run_df[f"DeltaRm__{key}"].reindex(common_proteins)

	for cond in exp_conditions:
		key = cond["key"]
		sample_name = f"{cond['name']} ({cond['run_name']})"
		delta_col = f"DeltaRm__{key}"
		deltas = full_df[delta_col]

		pcol = f"{delta_col}_pval"
		fcol = f"{delta_col}_FDR"
		mode_col = f"{delta_col}_analysis_mode"
		msg_col = f"{delta_col}_normality_message"
		low_col = f"{delta_col}_BelowMinus0_1"
		high_col = f"{delta_col}_Above0_1"
		abs_col = f"{delta_col}_AbsGt0_1"
		emp_low_col = f"{delta_col}_EmpiricalLow5pct"
		emp_high_col = f"{delta_col}_EmpiricalHigh5pct"
		emp_outlier_col = f"{delta_col}_EmpiricalOutlier"
		emp_tailn_col = f"{delta_col}_EmpiricalTailN"
		sig_col = f"{delta_col}_Sig"
		basis_col = f"{delta_col}_Sig_basis"

		normality_plots = plot_delta_normality_diagnostics(deltas, sample_name, base_prefix)
		if normality_plots:
			logger.info(f"{delta_col}: normality diagnostics saved: {', '.join(normality_plots)}")

		is_normal, normality_msg = assess_delta_normality(deltas)
		branch = "normal_robust_sigma" if is_normal else "non_normal_descriptive"
		full_df[mode_col] = branch
		full_df[msg_col] = normality_msg
		logger.info(f"{delta_col}: {normality_msg}; branch={branch}")

		low_mask = (deltas < -0.1).fillna(False)
		high_mask = (deltas > 0.1).fillna(False)
		abs_mask = (deltas.abs() > 0.1).fillna(False)
		emp_low_mask, emp_high_mask, n_tail = empirical_tail_masks(deltas, tail_frac=0.05)
		emp_outlier_mask = emp_low_mask | emp_high_mask

		full_df[low_col] = low_mask
		full_df[high_col] = high_mask
		full_df[abs_col] = abs_mask
		full_df[emp_low_col] = emp_low_mask.reindex(full_df.index, fill_value=False)
		full_df[emp_high_col] = emp_high_mask.reindex(full_df.index, fill_value=False)
		full_df[emp_outlier_col] = emp_outlier_mask.reindex(full_df.index, fill_value=False)
		full_df[emp_tailn_col] = int(n_tail)
		logger.info(
			f"{delta_col}: empirical outliers top/bottom 5% -> "
			f"low={int(emp_low_mask.sum())}, high={int(emp_high_mask.sum())}, tail_n={int(n_tail)}"
		)

		if is_normal:
			pvals_series = robust_global_sigma_pvals(deltas)
			pvals_arr = np.where(np.isnan(pvals_series.values), 1.0, pvals_series.values)
			_, fdrs, _, _ = multipletests(pvals_arr, method="fdr_bh")

			full_df[pcol] = pvals_arr
			full_df[fcol] = fdrs
			full_df[sig_col] = abs_mask & (full_df[fcol] < 0.05)
			full_df[basis_col] = "statistical(|DeltaRm|>0.1 & FDR<0.05; global robust-sigma)"

			volcano_plot = plot_delta_volcano(deltas, full_df[fcol], sample_name, base_prefix)
			if volcano_plot is not None:
				logger.info(f"{delta_col}: volcano plot saved: {volcano_plot}")
		else:
			full_df[pcol] = np.nan
			full_df[fcol] = np.nan
			full_df[sig_col] = full_df[emp_outlier_col]
			full_df[basis_col] = (
				"descriptive(non-normal; no p/FDR; significance=empirical top/bottom 5% outlier; "
				"|DeltaRm|>0.1 is plotting reference only)"
			)

			dist_plot = plot_delta_distribution_threshold(deltas, sample_name, base_prefix)
			if dist_plot is not None:
				logger.info(f"{delta_col}: threshold distribution plot saved: {dist_plot}")

	full_outfile = base_prefix + "_RefinedTPP_plus_full_Rm_results.csv"
	full_df.reset_index().rename(columns={"index": "ProteinID"}).to_csv(full_outfile, index=False)
	logger.info(f"Full Rm results saved: {full_outfile}")

	sig_cols = [f"DeltaRm__{c['key']}_Sig" for c in exp_conditions]
	if sig_cols:
		sig_mask = full_df[sig_cols].fillna(False).any(axis=1)
	else:
		sig_mask = pd.Series(False, index=full_df.index)

	diff_df = full_df.loc[sig_mask].copy()
	diff_outfile = base_prefix + "_Differential_proteins_RefinedTPP_plus_full_matrix.csv"
	diff_df.reset_index().rename(columns={"index": "ProteinID"}).to_csv(diff_outfile, index=False)
	logger.info(f"Differential/flagged protein output saved: {diff_outfile}")

	mapping_rows = []
	for cond in all_conditions:
		mapping_rows.append(
			{
				"Condition_name": cond["name"],
				"Condition_key": cond["key"],
				"Run": cond["run_name"],
				"File": cond["file_path"],
				"Reporter_pair": cond["pair_label"],
				"Numerator_channel": cond["num"],
				"Denominator_channel": cond["den"],
				"Numerator_column": cond["num_col"],
				"Denominator_column": cond["den_col"],
			}
		)
	mapping_df = pd.DataFrame(mapping_rows)
	mapping_outfile = base_prefix + "_condition_mapping.csv"
	mapping_df.to_csv(mapping_outfile, index=False)
	logger.info(f"Condition mapping saved: {mapping_outfile}")

	norm_df = pd.DataFrame(norm_rows)
	norm_outfile = base_prefix + "_normalization_factors.csv"
	norm_df.to_csv(norm_outfile, index=False)
	logger.info(f"Normalization factors saved: {norm_outfile}")

	cond_quant_rows = []
	for cond in all_conditions:
		key = cond["key"]
		cond_quant_rows.append(
			{
				"Condition_name": cond["name"],
				"Condition_key": key,
				"Run": cond["run_name"],
				"Reporter_pair": cond["pair_label"],
				"Quantified_count_any_signal": int(len(all_condition_quant_sets[key])),
				"Complete_count_2channel": int(len(all_condition_complete_sets[key])),
			}
		)
	cond_quant_df = pd.DataFrame(cond_quant_rows)
	cond_quant_outfile = base_prefix + "_quantification_overview_by_condition.csv"
	cond_quant_df.to_csv(cond_quant_outfile, index=False)
	logger.info(f"Condition quantification overview saved: {cond_quant_outfile}")

	global_quant_overview_df = pd.DataFrame(
		[
			{
				"Unique_filter_requested_ge2": bool(require_unique_ge2),
				"Unique_filter_effective_all_runs_ge2": bool(
					all(entry["summary"]["Unique_filter_effective_ge2"] for entry in run_entries)
				),
				"Global_condition_count": int(len(all_conditions)),
				"Global_required_channel_count_2N": int(2 * len(all_conditions)),
				"Global_quantified_union_count_any_signal": int(len(global_quant_union)),
				"Global_quantified_intersection_count_all_conditions_any_signal": int(len(global_quant_intersection)),
				"Global_complete_quantification_count_2N": int(len(common_proteins)),
				"Global_overlap_plot_file": global_overlap_plot if global_overlap_plot is not None else "",
			}
		]
	)
	global_quant_overview_outfile = base_prefix + "_quantification_overview_global.csv"
	global_quant_overview_df.to_csv(global_quant_overview_outfile, index=False)
	logger.info(f"Global quantification overview saved: {global_quant_overview_outfile}")

	summary_df = pd.DataFrame(summary_rows)
	summary_df["Global_common_protein_count"] = int(len(common_proteins))
	summary_df["Global_condition_count"] = int(len(all_conditions))
	summary_df["Global_required_channel_count_2N"] = int(2 * len(all_conditions))
	summary_df["Global_quantified_union_count_any_signal"] = int(len(global_quant_union))
	summary_df["Global_quantified_intersection_count_all_conditions_any_signal"] = int(len(global_quant_intersection))
	summary_df["Unique_filter_requested_ge2"] = bool(require_unique_ge2)
	summary_df["Unique_filter_effective_all_runs_ge2"] = bool(
		all(entry["summary"]["Unique_filter_effective_ge2"] for entry in run_entries)
	)
	summary_df["Vehicle_condition_key"] = vehicle_key
	summary_df["Global_vehicle_log_median"] = float(global_vehicle_log_median)
	summary_df["Global_vehicle_Rm_reference"] = float(global_vehicle_rm_ref)

	runs_without_vehicle = int(
		sum(1 for e in run_entries if e.get("vehicle_reference_mode") == "global_vehicle_median_constant")
	)
	summary_df["Runs_without_vehicle_reference_mode"] = runs_without_vehicle

	summary_outfile = base_prefix + "_NoRep_summary.csv"
	summary_df.to_csv(summary_outfile, index=False)
	logger.info(f"No-replicate summary saved: {summary_outfile}")

	logger.info("Refined-TPP plus no-replicate processing completed successfully.")


if __name__ == "__main__":
	main()

