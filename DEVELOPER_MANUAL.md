# 📘 Developer Manual — Chat-RefinedTPP Platform

> This manual is developer-facing and intentionally avoids repeating conceptual content from '[README.md](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/blob/main/README.md)'.
> It focuses on executable setup, local deployment, branch execution, and maintenance procedures.

---

## Table of Contents
- [1. Scope and Audience](#1-scope-and-audience)
- [2. Usage Modes](#2-usage-modes)
- [3. Environment Setup (Python / Conda)](#3-environment-setup-python--conda)
- [4. Repository Setup](#4-repository-setup)
- [5. Unified Execution Flow (All 8 Branches)](#5-unified-execution-flow-all-8-branches)
- [6. 8-Branch Difference Matrix](#6-8-branch-difference-matrix)
- [7. Branch-Specific Run Notes (B1–B8)](#7-branch-specific-run-notes-b1b8)
- [8. Local Agent Deployment (Web-like Experience)](#8-local-agent-deployment-web-like-experience)
- [9. Hosted Backend Operations (Maintainer-only)](#9-hosted-backend-operations-maintainer-only)
- [10. Data Contract & Validation Checklist](#10-data-contract--validation-checklist)
- [11. Testing & QA](#11-testing--qa)
- [12. Troubleshooting Runbook](#12-troubleshooting-runbook)
- [13. Contribution Workflow](#13-contribution-workflow)
- [14. Changelog](#14-changelog)

---

## 1. Scope and Audience

### 1.1 Purpose
This document provides developer procedures for:
- local environment setup
- local script execution (all 8 branches)
- local Agent deployment
- hosted backend maintenance

### 1.2 Intended readers
- Lab developers
- Platform maintainers
- Contributors extending pipelines

---

## 2. Usage Modes

### Mode A — Script-only Local Analysis
Run branch scripts directly with local input files.  
**Best for:** analysis users who do not need Agent/web orchestration.

### Mode B — Local Agent Deployment
Deploy backend + agent (+ optional frontend) locally for web-like interaction.  
**Best for:** labs needing private interactive workflow.

### Mode C — Hosted Web Operations (Maintainer-only)
Maintain production-like backend services for end users.  
**Best for:** core maintainers only.

---

## 3. Environment Setup (Python / Conda)

> Recommended Python version: **3.10–3.11** for stable compatibility.
> Conda is recommended for easier dependency management, but venv is also supported.
> Verified on maintainer machine: Python 3.11.5 (Anaconda), IPython 8.27.0

### 3.1 Option A (Recommended): Conda / Anaconda

#### 3.1.1 Install Conda

Please install **Miniconda** (lightweight, recommended) or **Anaconda** (full distribution):

- Miniconda (official): https://docs.conda.io/en/latest/miniconda.html  
- Anaconda (official): https://www.anaconda.com/download

**OS notes (one-line guidance):**
- **Windows**: Download the `.exe` installer. **Do not check “Add to PATH”** (to avoid conflicts). Finish installation, then open **Anaconda Prompt** (or **Miniconda Prompt**) from the Start menu.
- **macOS**: Download the correct installer for Intel/Apple Silicon, install, then restart terminal and run `conda --version`.
- **Linux**: Download the shell installer (`.sh`), run `bash <installer>.sh`, then run `source ~/.bashrc` (or `~/.zshrc`).

**Verify installation:**
```bash
conda --version
```

#### 3.1.2 Create environment

Create an isolated Conda environment named `refinedtpp` with Python 3.11, then activate it:

```bash
conda create -n refinedtpp python=3.11 -y
conda activate refinedtpp
python --version
```

#### 3.1.3 Install dependencies

Install project dependencies inside the activated environment:

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

Optional quick local setup (temporary, for rapid iteration before updating pinned deps):

```bash
pip install pandas numpy seaborn scipy matplotlib statsmodels matplotlib-venn upsetplot
```

#### 3.1.4 (Optional) Freeze environment

Export current environment for reproducibility:

```bash
conda env export --no-builds > environment.yml
pip freeze > requirements-lock.txt
```

#### 3.1.5 Dependency files in this repository

- [Download `requirements.txt`](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/raw/main/requirements.txt): standard dependency list for most users.
- [Download `requirements-lock.txt`](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/raw/main/requirements-lock.txt): fully pinned environment snapshot for strict reproducibility.
- [Download `environment.yml`](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/raw/main/environment.yml): Conda environment definition file.

Recommended usage:
- Standard install: `pip install -r requirements.txt`
- Strict reproducibility: `pip install -r requirements-lock.txt`
- Conda recreate: `conda env create -f environment.yml`

### 3.2 Option B: venv

If you do not use Conda, use Python built-in venv:

```bash
python -m venv .venv
# Windows
# .venv\Scripts\activate
# macOS/Linux
# source .venv/bin/activate
pip install -r requirements.txt
```

---

## 4. Repository Setup

### 4.1 Clone repositories
> Team baseline repository:

```bash
git clone https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform.git
cd Refined-TPP_Chat_Platform
# The default branch is 'Development', but ensure you are on it:
git checkout Development
```

> Optional personal draft repository (private, principal maintainer only):
> `Chat-RefinedTPP-Python-Draft`  
> Note: This repo is not required for collaborators and is not the source of truth.

### 4.2 Branch policy

- Default branch: `Development` (**all active development, module intergration, and debugging occur here**)
- Stable branch: `main` (Reserved for verified milestones. Code is merged from `Development` to `main` via **Pull Request (PR)** only)
- Feature branches: `feature/<topic>` (for specific tasks, branched from `Development`, merged back via PR)
- Hotfix branches: `hotfix/<topic>` (used for urgent fixes on `main`, then back-merge to both `main` and `Development`)
- Draft scripts repo: private personal use only, no team dependency
- Release/Tag policy:
  - Official versions are tagged only on the `main` branch after successful PR merges.
  - Version format: `vX.Y.Z` (e.g., `v0.1.0`)

### 4.3 Suggested workspace layout

To prevent large Proteomic DB-search and processed result files from being accidentally uploaded to GitHub, follow this layout:

```text
workspace/
├── Refined-TPP_Chat_Platform/    # source of truth
├── data/                         # local runtime data (not versioned)
│   ├── input/
│   ├── output/
│   └── reference/
└── logs/                         # local logs (not versioned)
```

### 4.4 Source-of-truth statement

The `Refined-TPP_Chat_Platform` repository is the authoritative source for the project.

1. **No direct pushes to** `main`: All changes must go through `Development`.

2. **Integration**: Once a module (e.g., a specific Python pipeline branch) is stable, a PR should be opened to merge it into `main` for versioning.

3. The private draft repository is for prototyping only and must not be referenced in production code.

---

## 5. Unified Execution Flow (All 8 Branches)

1. Activate environment 
2. Select branch script (or pipeline entry)
3. Prepare required input files
4. Run scripts and provide interactive parameters
5. Validate outputs (`CSV/TIFF` + `logs`)
6. Archive result package

### 5.1 Generic run command
```bash
python <BranchXX_P/GX>.py
# Replace "X" with the specific number of the scrpit you need
```
### 5.2 Generic input classes
- Protein quantification table
- Peptide quantification table (N-glycosite branches)
- Reference FASTA (N-glycosite branches)
- Branch-specific metadata and thresholds

---

## 6. 8-Branch Difference Matrix

| Branch | Pipeline | Omics | Tech Rep | Bio Rep | Required Inputs | Key Variable Inputs | Main Outputs |
|---|---|---|---:|---:|---|---|---|
| 01 | P1 | Bulk proteome | ✅ | ✅ | [To fill] | [To fill] | [To fill] |
| 02 | P2 | Bulk proteome | ❌ | ✅ | [To fill] | [To fill] | [To fill] |
| 03 | P3 | Bulk proteome | ✅ | ❌ | [To fill] | [To fill] | [To fill] |
| 04 | P4 | Bulk proteome | ❌ | ❌ | [To fill] | [To fill] | [To fill] |
| 05 | G1 | N-glycosite | ✅ | ✅ | [To fill] | [To fill] | [To fill] |
| 06 | G2 | N-glycosite | ❌ | ✅ | Protein CSV + Peptide CSV + FASTA | CV threshold / motif window / significance settings | QC CSV + ΔRm CSV + significance CSV + TIFF |
| 07 | G3 | N-glycosite | ✅ | ❌ | [To fill] | [To fill] | [To fill] |
| 08 | G4 | N-glycosite | ❌ | ❌ | [To fill] | [To fill] | [To fill] |

> Maintenance rule: keep shared steps in Sections 3–5; only branch differences belong here.

---

## 7. Branch-Specific Run Notes (B1–B8)

### 7.1 Branch 01 (P1)
- Script: [To fill]
- Special inputs: [To fill]
- Special outputs: [To fill]
- Notes: [To fill]

### 7.2 Branch 02 (P2)
- Script: [To fill]
- Special inputs: [To fill]
- Special outputs: [To fill]
- Notes: [To fill]

### 7.3 Branch 03 (P3)
- Script: [To fill]
- Special inputs: [To fill]
- Special outputs: [To fill]
- Notes: [To fill]

### 7.4 Branch 04 (P4)
- Script: [To fill]
- Special inputs: [To fill]
- Special outputs: [To fill]
- Notes: [To fill]

### 7.5 Branch 05 (G1)
- Script: [To fill]
- Special inputs: [To fill]
- Special outputs: [To fill]
- Notes: [To fill]

### 7.6 Branch 06 (G2)
- Script: `Branch06_G2_update.py`
- Typical inputs:
  - TMT protein quant CSV
  - TMT peptide quant CSV
  - UniProt FASTA
- Typical interactive parameters:
  - sample naming / control mapping
  - motif window length (odd)
  - CV threshold
  - significance branch settings
- Typical outputs:
  - protein/site QC tables
  - ΔRm-related tables
  - significance tables
  - TIFF visualizations
- Notes:
  - Ensure accession format consistency between peptide CSV and FASTA
  - Ensure selected group names match intensity-column prefixes exactly

### 7.7 Branch 07 (G3)
- Script: [To fill]
- Special inputs: [To fill]
- Special outputs: [To fill]
- Notes: [To fill]

### 7.8 Branch 08 (G4)
- Script: [To fill]
- Special inputs: [To fill]
- Special outputs: [To fill]
- Notes: [To fill]

---

## 8. Local Agent Deployment (Web-like Experience)

### 8.1 Goal
Provide local interactive analysis similar to hosted web usage.

### 8.2 Service components
- frontend (optional for pure API use)
- backend API
- agent process
- pipeline execution layer

### 8.3 Environment variables
```bash
# Example only; adjust to actual implementation
OPENAI_API_KEY=...
BACKEND_PORT=...
FRONTEND_PORT=...
DATA_ROOT=...
LOG_LEVEL=INFO
```

### 8.4 Startup order
1. Start backend API  
2. Start agent service  
3. Start frontend (if used)  
4. Run health checks  
5. Execute one end-to-end sample task

### 8.5 Health checks
- API `/health` reachable
- Agent responds to routing request
- Pipeline job can start and finish successfully

---

## 9. Hosted Backend Operations (Maintainer-only)

### 9.1 Operational boundary
- End users: web upload + interaction only
- Maintainer: backend deployment, monitoring, incident response

### 9.2 Deployment checklist
- environment variables configured
- secrets managed securely
- storage path/permissions validated
- logging/monitoring enabled
- rollback strategy prepared

### 9.3 Release flow
- [To fill] staging validation
- [To fill] production rollout
- [To fill] rollback condition and command

---

## 10. Data Contract & Validation Checklist

### 10.1 Input schema checklist
- required columns present
- accession format consistent
- numeric intensity columns parseable
- missing/zero handling policy clear

### 10.2 Pre-run validation
- [To fill] automated validator script path
- [To fill] manual fallback checklist

### 10.3 Output acceptance
- expected CSV artifacts generated
- expected TIFF artifacts generated
- no fatal warnings in logs

---

## 11. Testing & QA

### 11.1 Unit tests
```bash
# Example
pytest -q
```

### 11.2 Integration tests
- branch run smoke test (all 8)
- agent-to-pipeline invocation test
- output schema regression test

### 11.3 Golden dataset regression
- [To fill] dataset location
- [To fill] acceptance thresholds

---

## 12. Troubleshooting Runbook

### 12.1 Environment errors
- dependency conflicts → recreate clean env
- package import errors → verify active env and reinstall

### 12.2 Data/input errors
- file not found → verify absolute paths
- empty glycosite output → check FASTA mapping and accession formats
- missing columns → compare against schema checklist

### 12.3 Service errors
- backend unavailable → verify port/process
- agent timeout → inspect agent logs and upstream API key status
- frontend can’t fetch results → check backend URL/CORS and job status endpoint

---

## 13. Contribution Workflow

1. Create feature branch  
2. Implement + test  
3. Update docs (README link + manual section)  
4. Open PR with:
   - change summary
   - validation evidence
   - backward compatibility note

---

## 14. Changelog

- **v0.1.0** — Initial single-file developer manual skeleton
- **v0.1.1** — [To fill]