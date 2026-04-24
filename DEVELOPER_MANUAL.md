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

### 3.3 Node.js / Frontend Setup (Optional)
If you are developing or running the local Agent UI (Mode B), you need Node.js.
1. Download and install Node.js (LTS version recommended: v18 or v20) from [nodejs.org](https://nodejs.org/).
2. Verify installation:
   ```bash
   node -v
   npm -v
   ```
3. Install frontend dependencies:
  ```bash
  cd frontend
  npm install
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

- Default branch: `Development` (**all active development, module integration, and debugging occur here**)
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
# Replace "X" with the specific number of the script you need
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
Provide local interactive analysis similar to hosted web usage, driven by a Large Language Model (LLM) agent orchestrating the 8 analytical branches.

### 8.2 Tech Stack & Service Components
- **API Layer**: **FastAPI**
  > Chosen for native `async` support to efficiently handle long-polling LLM I/O, and native `Pydantic` integration for strict bioinformatic parameter validation.
- **Agent Framework**: **LangChain**
  > Chosen for its robust Tool Calling mechanism (mapping the 8 Python scripts into executable tools), Expression Language (LCEL), and built-in multi-turn conversational memory.
- **Execution Layer**: Subprocess or module imports triggering `BranchXX_P/GX.py`.
- **Frontend** (Optional): [To fill, e.g., React/Vue or Streamlit].

### 8.3 Architectural Principle: Decoupling Control Flow from Data Flow
To prevent LLM context window limits and control API costs (Token explosion), the Agent operates strictly as an orchestrator, **never reading raw data files**.

**The strict separation rule (v0.1.0):**
1. **File Upload**: Raw CSV/FASTA files are stored directly to the local disk/workspace by the FastAPI backend. The LLM is bypassed.
2. **Routing & Tool Calling (LLM)**: The Agent only receives **metadata** (e.g., file paths, user prompts, extracted parameters like CV thresholds). It uses this metadata to decide which Python branch to execute.
3. **Execution (Python)**: The underlying `BranchXX.py` scripts read the large files from disk using `pandas` and perform the heavy lifting locally.
4. **Summary Generation (LLM)**: To generate the final biological report, the Python script must output a condensed `summary.json` (e.g., number of significant proteins, top 5 changes). The Agent reads *only* this JSON, never the resulting full CSVs.

*Anti-pattern to avoid*: Never pass `dataframe.to_string()` or raw CSV contents into the LangChain prompt.

### 8.4 Backend Directory Structure Reference
To maintain separation of concerns, the backend strictly follows this layout:

```text
backend/
├── [main.py](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/tree/Development/backend/main.py)          # FastAPI entry point and API routers
├── [schemas/](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/tree/Development/backend/schemas)         # Pydantic models (Data Contracts for G1-G4, P1-P4)
├── [agent/](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/tree/Development/backend/agent)           # LangChain core logic
│   ├── tools.py       # 8 branches wrapped as LangChain Tools
│   ├── prompts.py     # System prompts and few-shot routing examples
│   └── workflow.py    # Conversation state and execution flow
└── [core/](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/tree/Development/backend/core)            # Global config and environment loading
```

### 8.5 Environment variables
Create a `.env` file in the `backend/` directory:

```bash
# .env.example
OPENAI_API_KEY=sk-...
OPENAI_BASE_URL=https://api.openai.com/v1
AGENT_MODEL_NAME=gpt-4o
BACKEND_PORT=8000
WORKSPACE_DIR=../data/workspace
LOG_LEVEL=INFO
```

### 8.6 Startup Order & Commands
1. Start backend API & Agent service:

```bash
cd backend
# Uvicorn is the ASGI server for FastAPI
uvicorn main:app --reload --port 8000
```
 
2. Start frontend (if used):

```bash
cd frontend
npm run dev
```

3. Run health checks:
- API `/health` reachable (Check Swagger UI)
- Agent responds to routing request
- Pipeline job can start and finish successfully

5. Execute one end-to-end sample task:
- Upload sample CSVs to trigger a simple P2 or G2 job.
- Verify that the Agent correctly extracts parameters, the Python script executes without errors, and the `summary.json` is successfully parsed by the LLM into a final report.

---

## 9. Hosted Backend Operations (Maintainer-only)

### 9.1 Operational boundary
- End users: web upload + interaction only
- Maintainer: backend deployment, LLM API key management, monitoring, incident response.

### 9.2 Deployment checklist
- [ ] **ASGI Server**: Use `gunicorn` with `uvicorn` workers for FastAPI production deployment (e.g., `gunicorn main:app -k uvicorn.workers.UvicornWorker`)
- [ ] **State Management**: Ensure external memory (e.g., Redis or SQLite) is configured for LangChain to persist user sessions across server restarts.
- [ ] **Timeout Configs**: Extend reverse proxy (Nginx/Traefik) timeouts. LLM reasoning and pipelines can run for several minutes.
- [ ] **Secrets Managed:**: API keys managed securely (do NOT hardcode `OPENAI_API_KEY`).
- [ ] **Storage**: Verify `WORKSPACE_DIR` read/write permissions for large CSV/TIFF operations.

### 9.3 Release flow
- [To fill] Staging validation (Automated Agent routing tests)
- [To fill] Production rollout
- [To fill] Rollback condition and command

---

## 10. Data Contract & Validation Checklist

> **Architecture Note**: All Agent-to-Pipeline parameter parsing MUST be strictly typed using **FastAPI / Pydantic** schemas. This prevents LLM "hallucinations" from passing invalid thresholds (e.g., passing a string instead of a float for CV threshold) to the pipeline.

### 10.1 Input Schema Checklist
- Reference implementation: [`backend/schemas/`](https://github.com/SongHaoru-max/Refined-TPP_Chat_Platform/tree/Development/backend/schemas)
- Pydantic models must enforce:
  - required columns present
  - accession format consistent
  - numeric intensity columns parseable (e.g., `cv_threshold: float = Field(ge=0.0, le=1.0)`)
  - missing/zero handling policy clear

### 10.2 Pre-run validation
- [To fill] Automated validator script path (e.g., `backend/core/validator.py`)
- [To fill] Manual fallback checklist

### 10.3 Output Acceptance
- expected CSV artifacts generated
- expected TIFF artifacts generated
- no fatal warnings in logs
- `summary.json` generated for LLM report consumption

---

## 11. Testing & QA

### 11.1 Unit tests
```bash
# Run tests for core bioinformatic functions and schema validation
pytest backend/tests/ -q
```

### 11.2 Integration tests
- **Branch run smoke test (all 8)**: Ensure pure Python scripts run end-to-end without LLM.
- **Agent-to-Pipeline invocation test**: Mock LLM responses to ensure FastAPI triggers the correct Tool.
- **Output schema regression test**: Validate that the generated `summary.json` format remains consistent.

### 11.3 Golden dataset regression
- [To fill] Dataset location
- [To fill] Acceptance thresholds

---

## 12. Troubleshooting Runbook

### 12.1 Environment Errors
- **Dependency conflicts** → recreate clean env from `requirements-lock.txt`.
- **Package import errors** → verify active env (`conda activate refinedtpp`) and reinstall.

### 12.2 Data/Input Errors
- **File not found** → verify absolute paths in `WORKSPACE_DIR`.
- **Empty glycosite output** → check FASTA mapping and accession formats.
- **Missing columns** → compare against schema checklist.

### 12.3 Service & Agent Errors
- **Backend Unavailable** → Verify Uvicorn/Gunicorn port and process.
- **Agent JSON Parsing Error** → If LLM hallucination breaks parameter extraction, enable `LOG_LEVEL=DEBUG` to inspect the raw LLM string. Ensure LangChain's Output Parser is catching and retrying formatting errors.
- **Agent Branch Misrouting** → If the Agent selects the wrong branch (e.g., routes G2 to P1), verify the Tool descriptions in `agent/tools.py` and update the few-shot examples in `prompts.py`.
- **Context Window Exceeded** → Ensure the Agent is NOT reading raw multi-megabyte Proteome CSVs. Refer to the *Decoupling Control Flow from Data Flow* principle in Section 8.
- **Frontend Fetch Failure** → Check backend URL/CORS and job status endpoint.

---

## 13. Contribution Workflow

1. Create feature branch from `Development`.
2. Implement Python logic + update LangChain Tools/Pydantic schemas. 
3. **Format and Lint**: Ensure your code meets PEP8 standards. We recommend using `black` for formatting and `mypy` for static type checking:
   ```bash
   black .
   mypy backend/
   ```
3. Add unit tests for new schemas/functions.
4. Update docs (README link + Manual section). 
5. Open PR with:
   - change summary
   - validation evidence
   - backward compatibility note

---

## 14. Changelog

- **v0.1.0** — Initial single-file developer manual skeleton
- **v0.1.1** — [To fill]