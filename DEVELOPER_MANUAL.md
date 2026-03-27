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
