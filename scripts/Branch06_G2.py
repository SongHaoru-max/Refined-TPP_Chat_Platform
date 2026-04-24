# -*- coding: utf-8 -*-
"""
Full TMT Proteome & N-glycopeptide Quantification QC & Filtering Pipeline
------------------------------------------------------------------------
Created on: Thu Nov 21 2025
Author: Haoru Song, Chenxin Li, Bin Fu
Affiliation: Professor Haojie Lu's Lab, Fudan University, Shanghai, China

Description
-----------
Integrates:
1. Full 6-plex TMT protein quantification QC & filtering
2. N-glycopeptide site-centric analysis
3. Rm calculation, normalization, and visualization
Outputs CSV and TIFF images
"""

import os
import re
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
            if '+0.98' in mods:
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


def paired_one_tailed_pvalue(bulk_vals, glyco_vals, method="wilcoxon", side="right"):
    """Compute one-tailed paired p-value for one site and return (p, n_pairs)."""
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
                method="auto"
            )
        except TypeError:
            _, p_one = stats.wilcoxon(
                glyco_valid,
                bulk_valid,
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


def consensus_significance_from_replicates(delta_matrix, fdr_matrix, delta_threshold=0.1):
    """
    Consensus significant sites across biological replicates.
    Conditions:
    1) same DeltaRm sign across replicates;
    2) DeltaRm > threshold in all replicates;
    3) all replicate FDR < 0.1;
    4) at least one replicate FDR < 0.05.
    """
    delta_matrix = np.asarray(delta_matrix, dtype=float)
    fdr_matrix = np.asarray(fdr_matrix, dtype=float)

    finite_delta = np.isfinite(delta_matrix).all(axis=1)
    finite_fdr = np.isfinite(fdr_matrix).all(axis=1)
    same_sign = ((delta_matrix > 0).all(axis=1) | (delta_matrix < 0).all(axis=1))
    delta_pass = (delta_matrix > delta_threshold).all(axis=1)
    fdr_all_10 = (fdr_matrix < 0.1).all(axis=1)
    fdr_any_05 = (fdr_matrix < 0.05).any(axis=1)

    sig_mask = finite_delta & finite_fdr & same_sign & delta_pass & fdr_all_10 & fdr_any_05
    return sig_mask, same_sign, delta_pass, fdr_all_10, fdr_any_05

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
file = input("Please provide path of TMT quantification protein CSV file: ").strip()
file = normalize_path(file)
if not os.path.exists(file):
    raise FileNotFoundError(file)

work_dir = os.path.dirname(file)
os.chdir(work_dir)

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
        plt.savefig(os.path.splitext(file)[0]+"_quantified_proteins.tiff", dpi=300)
        plt.close()
    else:
        upset_data = from_contents(per_sample_sets)
        plt.figure(figsize=(8,6))
        UpSet(upset_data, show_counts=True).plot()
        plt.savefig(os.path.splitext(file)[0]+"_quantified_proteins_upset.tiff", dpi=300)
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
df_control_complete.to_csv(os.path.splitext(file)[0]+"_control_fully_quantified.csv", index=False)
df_strict.to_csv(os.path.splitext(file)[0]+"_all_samples_fully_quantified.csv", index=False)

print("TMT protein QC completed.")

# ------------------------------
# 2. Rm calculation & normalization of control group
# ------------------------------
ctrl_cols = sample_block_cols[control_sample]
ch126, ch127, ch128, ch129, ch130, ch131 = ctrl_cols
df_control_complete = df_control_complete.copy()


df_control_complete["Rm1"] = df_control_complete[ch129] / df_control_complete[ch127]
df_control_complete["Rm2"] = df_control_complete[ch130] / df_control_complete[ch128]
df_control_complete["Rm3"] = df_control_complete[ch131] / df_control_complete[ch126]

# global median & correction factor
medians = df_control_complete[["Rm1","Rm2","Rm3"]].median()
log_median_center = np.median(np.log(medians))
CF = np.exp(log_median_center - np.log(medians))
for i, col in enumerate(["Rm1","Rm2","Rm3"]):
    df_control_complete[f"{col}_normalized"] = df_control_complete[col] * CF[col]

cf_df = pd.DataFrame({
    "Repeat": ["Rm1","Rm2","Rm3"],
    "Median_Rm": medians.values,
    "Correction_Factor": CF.values
})
cf_df.to_csv(os.path.splitext(file)[0]+"_Rm_correction_factors.csv", index=False)
df_control_complete.to_csv(os.path.splitext(file)[0]+"_control_fully_quantified_with_Rm.csv", index=False)

print("Rm normalization completed for control group.")

# ------------------------------
# 3. N-glycopeptide analysis
# ------------------------------
peptide_csv_file = input("Enter path of TMT quantification peptide CSV file: ").strip()
peptide_csv_file = normalize_path(peptide_csv_file)

fasta_file = input("Enter path to UniProt fasta file: ").strip()
fasta_file = normalize_path(fasta_file)

work_dir = os.path.dirname(peptide_csv_file)
os.chdir(work_dir)

fasta_csv_file = os.path.join(work_dir, "Uniprot_fasta.csv")
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

output_csv_file = os.path.join(work_dir, "N_glycopeptides_withIntensity.csv")
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
tiff_file = os.path.join(work_dir, "N_glycopeptides_intensity.tiff")
fig.savefig(tiff_file, dpi=300, format='tiff')
plt.close(fig)
print(f"Stacked bar chart saved: {tiff_file}")

# ----------------------------
# Protein N-glycosite distribution (all filtered N-glycopeptides)
# ----------------------------
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

output_file = os.path.join(work_dir, "Protein_Nsite_distribution.tiff")
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
site_csv_file = os.path.join(work_dir, "N_glyco_site_intensity.csv")
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
site_complete_csv = os.path.join(work_dir, "N_glyco_site_intensity_completeQuant.csv")
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
out_plot = os.path.join(work_dir, "Protein_Nsite_distribution_completeQuant.tiff")
plt.savefig(out_plot, dpi=300, format='tiff')
plt.close()
print(f"Protein N-glycosite distribution saved: {out_plot}")

print("All processing completed successfully.")

# ----------------------------
# Fully quantified N-glycosites with fully quantified proteins (in control group) filtering
# ==============================
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
plt.savefig("Protein_level_quantification_completeness_venn.tiff", dpi=300, format='tiff')
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
plt.savefig("Protein_Nsite_distribution_completeQuant_intersect.tiff", dpi=300, format='tiff')
plt.close()

# ==============================
# Step 5: Calculate Rm, CV, and ΔRm for intersected proteins/sites
# ==============================
print(site_df_filtered.columns)
print(df_control_filtered.columns)

# ---------------- Complete sites ----------------
site_df_complete = site_df_complete.copy()
# Calculate raw Rm for each site
site_df_complete['Rm_glyco1'] = site_df_complete[ch129_site] / site_df_complete[ch127_site]
site_df_complete['Rm_glyco2'] = site_df_complete[ch130_site] / site_df_complete[ch128_site]
site_df_complete['Rm_glyco3'] = site_df_complete[ch131_site] / site_df_complete[ch126_site]
# Calculate CV for Rm
site_df_complete['CV_Rm_glyco_complete'] = site_df_complete[['Rm_glyco1','Rm_glyco2','Rm_glyco3']].std(axis=1) / \
                             site_df_complete[['Rm_glyco1','Rm_glyco2','Rm_glyco3']].mean(axis=1)

# Save complete site-level table with Rm and CV
complete_outfile = os.path.join(work_dir, f"N_glyco_sites_Rm_CV_{chosen_group}_complete.csv")
site_df_complete.to_csv(complete_outfile, index=False)
print(f"Complete site table with Rm and CV saved: {complete_outfile}")

# ---------------- Filtered sites ----------------
site_df_filtered = site_df_filtered.copy()
# Calculate raw Rm for each site
site_df_filtered['Rm_glyco1'] = site_df_filtered[ch129_site] / site_df_filtered[ch127_site]
site_df_filtered['Rm_glyco2'] = site_df_filtered[ch130_site] / site_df_filtered[ch128_site]
site_df_filtered['Rm_glyco3'] = site_df_filtered[ch131_site] / site_df_filtered[ch126_site]
# Calculate CV for Rm
site_df_filtered['CV_Rm_glyco_filtered'] = site_df_filtered[['Rm_glyco1','Rm_glyco2','Rm_glyco3']].std(axis=1) / \
                             site_df_filtered[['Rm_glyco1','Rm_glyco2','Rm_glyco3']].mean(axis=1)

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
df_control_filtered['CV_Rm_filtered'] = df_control_filtered[['Rm1','Rm2','Rm3']].std(axis=1, ddof=1) / df_control_filtered[['Rm1','Rm2','Rm3']].mean(axis=1)

# Complete proteins
df_control_complete = df_control_complete.copy()
df_control_complete['CV_Rm'] = df_control_complete[['Rm1','Rm2','Rm3']].std(axis=1, ddof=1) / df_control_complete[['Rm1','Rm2','Rm3']].mean(axis=1)

# Save filtered protein-level CSV
new_protein_cols = ['Rm1','Rm2','Rm3','Rm1_normalized','Rm2_normalized','Rm3_normalized','adjusted_Rm_unmod','CV_Rm_filtered']
original_protein_cols = [c for c in df_control_filtered.columns if c not in new_protein_cols]
protein_cols_to_save_filtered = original_protein_cols + new_protein_cols

protein_outfile_filtered = os.path.join(work_dir, "Control_Proteins_Rm_summary_filtered.csv")
df_control_filtered[protein_cols_to_save_filtered].to_csv(protein_outfile_filtered, index=False)
print(f"Filtered protein-level summary saved: {protein_outfile_filtered}")

# Save complete protein-level CSV
new_protein_cols_complete = ['Rm1','Rm2','Rm3','CV_Rm']
original_protein_cols_complete = [c for c in df_control_complete.columns if c not in new_protein_cols_complete]
protein_cols_to_save_complete = original_protein_cols_complete + new_protein_cols_complete

protein_outfile_complete = os.path.join(work_dir, "Control_Proteins_Rm_summary_complete.csv")
df_control_complete[protein_cols_to_save_complete].to_csv(protein_outfile_complete, index=False)
print(f"Complete protein-level summary saved: {protein_outfile_complete}")


# ---------------- Site-level ΔRm ----------------
# Arithmetic mean of normalized site Rm
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

# Save site-level CSV
new_site_cols = ['Rm_glyco1','Rm_glyco2','Rm_glyco3','CV_Rm_glyco_filtered',
                 'Rm_glyco1_norm','Rm_glyco2_norm','Rm_glyco3_norm',
                 'adjusted_Rm_glyco','adjusted_Rm_unmod',
                 'Rm_unmod1_norm','Rm_unmod2_norm','Rm_unmod3_norm',
                 'ΔRm_rep1','ΔRm_rep2','ΔRm_rep3','ΔRm']
original_site_cols = [c for c in site_df_filtered.columns if c not in new_site_cols]
site_cols_to_save = original_site_cols + new_site_cols

site_outfile = os.path.join(work_dir, "N_glyco_sites_with_DeltaRm_geom.csv")
site_df_filtered[site_cols_to_save].to_csv(site_outfile, index=False)
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

    cv_fig_file = os.path.join(work_dir, fname)
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

    cv_fig_file = os.path.join(work_dir, fname)
    plt.savefig(cv_fig_file, dpi=300)
    plt.close()
    print(f"{label} CV distribution TIFF plot saved: {cv_fig_file}")

# ------------------------------
# Step 8: CV-based QC filtering for sites (user threshold)
# ------------------------------
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
cv_filtered_output = os.path.join(work_dir, f"N_glyco_sites_CVfiltered_{int(cv_threshold*100)}pct.csv")
site_df_CV_pass.to_csv(cv_filtered_output, index=False)
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
    qc_plot_file = os.path.join(work_dir, f"Site_vs_Protein_CV_scatter_{int(cv_threshold*100)}pct.tiff")
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
    plt.savefig(f'Delta_Rm_rep{rep_idx}_histogram.tiff', format='tiff', dpi=300)
    plt.close()

    # --- 2. Q-Q plot ---
    plt.figure(figsize=(6,4))
    sm.qqplot(delta_rm_for_plot, line='s')
    plt.title(f'{rep_col} Q-Q Plot')
    plt.tight_layout()
    plt.savefig(f'Delta_Rm_rep{rep_idx}_qqplot.tiff', format='tiff', dpi=300)
    plt.close()

    # --- 3. P-P plot ---
    plt.figure(figsize=(6,4))
    sm.ProbPlot(delta_rm_for_plot).ppplot(line='45')
    plt.title(f'{rep_col} P-P Plot')
    plt.tight_layout()
    plt.savefig(f'Delta_Rm_rep{rep_idx}_ppplot.tiff', format='tiff', dpi=300)
    plt.close()

    # --- 4. Normality test for each biological replicate ---
    is_normal, normality_result, method_used, n_used = assess_distribution_normality(
        delta_rm_for_plot.values,
        alpha=0.05
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
normality_outfile = os.path.join(work_dir, "DeltaRm_replicate_normality_summary.csv")
normality_df.to_csv(normality_outfile, index=False)
print(f"Replicate-level normality summary saved: {normality_outfile}")

delta_rm_not_normal = not bool(normality_df["is_normal"].all())

# --- 5. Suggest test based on replicate-level normality ---
if delta_rm_not_normal:
    bad_reps = normality_df.loc[~normality_df["is_normal"], "replicate"].tolist()
    print(f"Normality gate failed for replicates: {bad_reps}")
    print("At least one biological replicate is NOT normal, use non-parametric paired testing")
else:
    print("All biological replicates are approximately normal, use CV-binned robust sigma testing per replicate")

# --- 6. Debug check: confirm original table unchanged ---
print("Step 9 check: site_df_CV_pass rows =", site_df_CV_pass.shape[0])

# ------------------------------
# Step 10-11: Significant site identification after CV QC (branch by normality)
# ------------------------------

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
            _, rep_fdr, _, _ = multipletests(rep_pvals, method="fdr_bh")
            rep_fdrs[rep_col] = rep_fdr

        fdr_matrix = np.column_stack([rep_fdrs[c] for c in delta_rep_cols])
        sig_mask, _, _, _, _ = consensus_significance_from_replicates(
            delta_matrix=delta_matrix_sorted,
            fdr_matrix=fdr_matrix,
            delta_threshold=0.1
        )
        sig_count = sig_mask.sum()
        sig_site_counts.append(sig_count)

        print(
            f"Bin count {bin_count}: {sig_count} significant sites "
            f"(same sign + ΔRm_rep>0.1 in all reps + FDR<0.1 in all reps + ≥1 rep FDR<0.05)"
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
        
        plot_file = os.path.join(work_dir, "Significant_sites_vs_bin_count.tiff")
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
        _, rep_fdr, _, _ = multipletests(rep_pvals, method="fdr_bh")
        custom_rep_pvals[rep_col] = rep_pvals
        custom_rep_fdrs[rep_col] = rep_fdr

    custom_pval_matrix = np.column_stack([custom_rep_pvals[c] for c in delta_rep_cols])
    custom_fdr_matrix = np.column_stack([custom_rep_fdrs[c] for c in delta_rep_cols])
    sig_mask_custom, same_sign_custom, delta_pass_custom, fdr_all_10_custom, fdr_any_05_custom = (
        consensus_significance_from_replicates(
            delta_matrix=delta_matrix_sorted,
            fdr_matrix=custom_fdr_matrix,
            delta_threshold=0.1
        )
    )
    
    # 6) Add results to the site table
    site_df_final = site_df_sorted.copy()
    site_df_final["Analysis_method"] = "normal_bin_robust_sigma_right_tailed_per_biorep_consensus"
    site_df_final["p_value_rep1"] = custom_rep_pvals["ΔRm_rep1"]
    site_df_final["p_value_rep2"] = custom_rep_pvals["ΔRm_rep2"]
    site_df_final["p_value_rep3"] = custom_rep_pvals["ΔRm_rep3"]
    site_df_final["FDR_rep1"] = custom_rep_fdrs["ΔRm_rep1"]
    site_df_final["FDR_rep2"] = custom_rep_fdrs["ΔRm_rep2"]
    site_df_final["FDR_rep3"] = custom_rep_fdrs["ΔRm_rep3"]
    # Keep aggregate p/FDR columns for compatibility with downstream table handling.
    site_df_final["p_value"] = custom_pval_matrix.max(axis=1)
    site_df_final["FDR"] = custom_fdr_matrix.max(axis=1)
    site_df_final["Same_sign_across_reps"] = same_sign_custom
    site_df_final["Delta_gt_0p1_all_reps"] = delta_pass_custom
    site_df_final["FDR_lt0p1_all_reps"] = fdr_all_10_custom
    site_df_final["FDR_lt0p05_any_rep"] = fdr_any_05_custom
    site_df_final["Significant"] = sig_mask_custom
    
    # 7) Save table
    custom_output_file = os.path.join(work_dir, f"N_glyco_sites_significant_binSize{custom_bin_size}.csv")
    site_df_final.to_csv(custom_output_file, index=False)
    n_sig = sig_mask_custom.sum()
    print("Consensus rule: same-sign ΔRm across replicates + ΔRm_rep>0.1 in all replicates + all FDR<0.1 + at least one FDR<0.05")
    print(f"Custom bin size {custom_bin_size}: {n_sig} significant sites saved to {custom_output_file}")

    # Volcano plots for normal branch: one plot per biological replicate
    for rep_idx, rep_col in enumerate(delta_rep_cols, start=1):
        fdr_col = f"FDR_rep{rep_idx}"
        if fdr_col not in site_df_final.columns:
            continue

        volcano_df_rep = site_df_final[[rep_col, fdr_col]].copy()
        volcano_df_rep = volcano_df_rep.replace([np.inf, -np.inf], np.nan).dropna()
        if volcano_df_rep.shape[0] == 0:
            print(f"Warning: no valid values for normal-branch volcano plot of {rep_col}.")
            continue

        volcano_df_rep["neglog10FDR"] = -np.log10(np.clip(volcano_df_rep[fdr_col].values, 1e-300, None))
        up_mask = (volcano_df_rep[rep_col] > 0.1) & (volcano_df_rep[fdr_col] < 0.05)
        down_mask = (volcano_df_rep[rep_col] < -0.1) & (volcano_df_rep[fdr_col] < 0.05)
        other_mask = ~(up_mask | down_mask)

        plt.figure(figsize=(6,5))
        plt.scatter(volcano_df_rep.loc[other_mask, rep_col], volcano_df_rep.loc[other_mask, "neglog10FDR"], s=14, color="grey", alpha=0.6)
        plt.scatter(volcano_df_rep.loc[up_mask, rep_col], volcano_df_rep.loc[up_mask, "neglog10FDR"], s=18, color="red", alpha=0.85, label="Up & significant")
        plt.scatter(volcano_df_rep.loc[down_mask, rep_col], volcano_df_rep.loc[down_mask, "neglog10FDR"], s=18, color="blue", alpha=0.85, label="Down & significant")
        plt.axvline(0.1, color="red", linestyle="--", linewidth=1)
        plt.axvline(-0.1, color="blue", linestyle="--", linewidth=1)
        plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
        plt.xlabel(rep_col)
        plt.ylabel("-log10(FDR)")
        plt.title(f"Volcano-like plot (normal branch, {rep_col})")
        plt.legend(frameon=False, fontsize=8)
        plt.tight_layout()
        volcano_plot_rep = os.path.join(work_dir, f"N_glyco_sites_volcano_normal_rep{rep_idx}.tiff")
        plt.savefig(volcano_plot_rep, dpi=300, format='tiff')
        plt.close()
        print(f"Normal-branch volcano plot saved: {volcano_plot_rep}")

else:
    print("Step10 strategy: replicate-level normality gate failed -> use paired one-tailed test (Bulk vs N-Glyco) per site")

    site_df_non_normal = site_df_CV_pass.copy().reset_index(drop=True)
    needed_cols = [
        "Rm_unmod1_norm", "Rm_unmod2_norm", "Rm_unmod3_norm",
        "Rm_glyco1_norm", "Rm_glyco2_norm", "Rm_glyco3_norm", "ΔRm"
    ]
    missing_cols = [c for c in needed_cols if c not in site_df_non_normal.columns]
    if missing_cols:
        raise KeyError(f"Required columns missing for non-normal paired branch: {missing_cols}")

    non_normal_method_input = input(
        "ΔRm not normal. Select paired one-tailed test [wilcoxon/paired_t, default=wilcoxon]: "
    ).strip().lower()
    if non_normal_method_input == "":
        non_normal_method = "wilcoxon"
    elif non_normal_method_input in ["wilcoxon", "paired_t"]:
        non_normal_method = non_normal_method_input
    else:
        raise ValueError("Invalid non-normal paired test. Use 'wilcoxon' or 'paired_t'.")

    print(f"Non-normal branch paired test: {non_normal_method} (one-tailed H1: glyco > bulk)")

    pvals = []
    n_pairs_used = []
    for _, row in site_df_non_normal.iterrows():
        bulk_vals = np.array([row["Rm_unmod1_norm"], row["Rm_unmod2_norm"], row["Rm_unmod3_norm"]], dtype=float)
        glyco_vals = np.array([row["Rm_glyco1_norm"], row["Rm_glyco2_norm"], row["Rm_glyco3_norm"]], dtype=float)

        p, n_pairs = paired_one_tailed_pvalue(
            bulk_vals=bulk_vals,
            glyco_vals=glyco_vals,
            method=non_normal_method,
            side="right"
        )
        pvals.append(float(p))
        n_pairs_used.append(int(n_pairs))

    pvals = np.array(pvals)
    _, fdrs, _, _ = multipletests(pvals, method="fdr_bh")

    site_df_final = site_df_non_normal.copy()
    site_df_final["Analysis_method"] = f"non_normal_paired_{non_normal_method}_right_tailed"
    site_df_final["Paired_test"] = non_normal_method
    site_df_final["n_pairs_used"] = n_pairs_used
    site_df_final["p_value"] = pvals
    site_df_final["FDR"] = fdrs
    site_df_final["Significant"] = (site_df_final["ΔRm"] > 0.1) & (site_df_final["FDR"] < 0.05)

    if non_normal_method == "wilcoxon" and len(n_pairs_used) > 0:
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
        non_normal_output_file = os.path.join(work_dir, "N_glyco_sites_significant_wilcoxon_rightTailed.csv")
    else:
        non_normal_output_file = os.path.join(work_dir, "N_glyco_sites_significant_pairedT_rightTailed.csv")
    site_df_final.to_csv(non_normal_output_file, index=False)
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
            bulk_vals = np.array([row["Rm_unmod1_norm"], row["Rm_unmod2_norm"], row["Rm_unmod3_norm"]], dtype=float)
            glyco_vals = np.array([row["Rm_glyco1_norm"], row["Rm_glyco2_norm"], row["Rm_glyco3_norm"]], dtype=float)

            ax.boxplot(
                [bulk_vals], positions=[1], widths=0.5, showfliers=False, patch_artist=True,
                boxprops=dict(facecolor="#4C78A8", alpha=0.35),
                medianprops=dict(color="#4C78A8"), whiskerprops=dict(color="#4C78A8"), capprops=dict(color="#4C78A8")
            )
            ax.boxplot(
                [glyco_vals], positions=[2], widths=0.5, showfliers=False, patch_artist=True,
                boxprops=dict(facecolor="#F58518", alpha=0.35),
                medianprops=dict(color="#F58518"), whiskerprops=dict(color="#F58518"), capprops=dict(color="#F58518")
            )

            x1 = 1 + rng.normal(0, 0.04, size=len(bulk_vals))
            x2 = 2 + rng.normal(0, 0.04, size=len(glyco_vals))
            ax.scatter(x1, bulk_vals, s=16, color="#4C78A8", zorder=3)
            ax.scatter(x2, glyco_vals, s=16, color="#F58518", zorder=3)

            combined = np.concatenate([bulk_vals, glyco_vals])
            y_min = np.nanmin(combined)
            y_max = np.nanmax(combined)
            y_span = (y_max - y_min) if (y_max > y_min) else 1.0
            y_text = y_max + 0.18 * y_span
            ax.plot([1, 2], [y_text, y_text], color="black", linewidth=0.8)
            star = pval_to_stars(row["FDR"])
            ax.text(1.5, y_text + 0.03 * y_span, f"FDR={row['FDR']:.2e} {star}", ha="center", va="bottom", fontsize=7)

            site_label = f"{row['Protein Accession']}|N{int(row['N-Glycosite'])}"
            ax.set_title(site_label, fontsize=7)
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["Bulk", "N-Glyco"], fontsize=7)
            ax.tick_params(axis='y', labelsize=7)
            ax.set_ylabel("Rm value", fontsize=7)

        for j in range(top_n, len(axes)):
            axes[j].axis("off")

        plt.tight_layout()
        top25_plot = os.path.join(work_dir, f"Top25_sites_Rm_boxplot_{non_normal_method}.tiff")
        plt.savefig(top25_plot, dpi=300, format='tiff')
        plt.close()
        print(f"Top-25 site boxplot saved: {top25_plot}")

    # Volcano-like plot for all sites
    volcano_df = site_df_final.copy()
    print("Non-normal volcano uses mean ΔRm across biological replicates (same ΔRm definition as site-level table).")
    volcano_df["neglog10FDR"] = -np.log10(np.clip(volcano_df["FDR"].values, 1e-300, None))
    up_mask = (volcano_df["ΔRm"] > 0.1) & (volcano_df["FDR"] < 0.05)
    down_mask = (volcano_df["ΔRm"] < -0.1) & (volcano_df["FDR"] < 0.05)
    other_mask = ~(up_mask | down_mask)

    plt.figure(figsize=(6,5))
    plt.scatter(volcano_df.loc[other_mask, "ΔRm"], volcano_df.loc[other_mask, "neglog10FDR"], s=14, color="grey", alpha=0.6)
    plt.scatter(volcano_df.loc[up_mask, "ΔRm"], volcano_df.loc[up_mask, "neglog10FDR"], s=18, color="red", alpha=0.85, label="Up & significant")
    plt.scatter(volcano_df.loc[down_mask, "ΔRm"], volcano_df.loc[down_mask, "neglog10FDR"], s=18, color="blue", alpha=0.85, label="Down & significant")
    plt.axvline(0.1, color="red", linestyle="--", linewidth=1)
    plt.axvline(-0.1, color="blue", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
    plt.xlabel("ΔRm")
    plt.ylabel("-log10(FDR)")
    method_name = "Wilcoxon" if non_normal_method == "wilcoxon" else "paired t-test"
    plt.title(f"Volcano-like plot ({method_name} + BH)")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    volcano_plot = os.path.join(work_dir, f"N_glyco_sites_volcano_{non_normal_method}.tiff")
    plt.savefig(volcano_plot, dpi=300, format='tiff')
    plt.close()
    print(f"Volcano-like plot saved: {volcano_plot}")
