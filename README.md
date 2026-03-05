# Chatt-RefinedTPP Platform
AI Agent–Driven Processing and Analysis Platform for Raw Refined-TPP Quantitative Proteomics


# 🤖 Chat-RefinedTPP: AI-Driven Proteomics Analytics Platform

<div align="center">
  <img src="assets/Logo.png" width="120" height="120" alt="RefinedTPP Logo">
  <br />
  <p align="center">
    <b>Next-Generation Refined-TPP quantitative data Preprocessing & QC & ΔRm Value Analysis Pipeline</b>
    <br />
    <i>Empowering high-precision quantitative proteomics with robust AI-agent collaboration.</i>
  </p>

  [![License](https://img.shields.io/badge/License-Academic%20Only-indigo.svg)](https://opensource.org/licenses/MIT)
  [![Python](https://img.shields.io/badge/Python-3.8+-blue.svg?logo=python&logoColor=white)](https://www.python.org/)
  [![Lab](https://img.shields.io/badge/Lab-Haojie%20Lu%20Lab-emerald.svg)]([http://homepage.fudan.edu.cn/haojielu/](https://chemistry.fudan.edu.cn/89/1f/c45942a690463/page.htm)
  [![Venue](https://img.shields.io/badge/JACS-2025-orange.svg)]([https://pubs.acs.org/journal/jacsat](https://pubs.acs.org/doi/10.1021/jacs.5c08065)
</div>

---

## ⚡ The 8-Branch Decision Matrix

Navigation through complex omics data is simplified via our **Triple-Layer Logic Gate**. The platform automatically routes your data through one of eight specialized branches based on your experimental design.

```mermaid
graph TD
    Start((<b>SUBMIT DATA</b>)) --> Type{<b>Omics Mode</b>}
    
    Type -- "Proteomics (Whole)" --> P_Bio{Bio Reps?}
    Type -- "N-Glyco / Mito" --> G_Bio{Bio Reps?}
    
    %% Proteomics Tree
    P_Bio -- Yes --> P_Tech_Y{Tech Reps?}
    P_Bio -- No --> P_Tech_N{Tech Reps?}
    
    P_Tech_Y -- Yes --> B1[<center><b>Branch 1</b><br/>Full Design</center>]
    P_Tech_Y -- No --> B2[<center><b>Branch 2</b><br/>Bio-Driven</center>]
    P_Tech_N -- Yes --> B3[<center><b>Branch 3</b><br/>Tech-Validated<br/><i>(TMT-6plex Core)</i></center>]
    P_Tech_N -- No --> B4[<center><b>Branch 4</b><br/>Rapid Screen</center>]
    
    %% Glyco Tree
    G_Bio -- Yes --> G_Tech_Y{Tech Reps?}
    G_Bio -- No --> G_Tech_N{Tech Reps?}
    
    G_Tech_Y -- Yes --> B5[<center><b>Branch 5</b><br/>Glyco-Mapping</center>]
    G_Tech_Y -- No --> B6[<center><b>Branch 6</b><br/>Site-Diff</center>]
    G_Tech_N -- Yes --> B7[<center><b>Branch 7</b><br/>Stability Eval</center>]
    G_Tech_N -- No --> B8[<center><b>Branch 8</b><br/>ID Only</center>]

    %% Aesthetics
    style B3 fill:#e0e7ff,stroke:#4f46e5,stroke-width:2px
    style Start fill:#f8fafc,stroke:#334155
    style Type fill:#fdf2f8,stroke:#db2777
