# Chat-RefinedTPP Platform
🤖 AI Agent–Driven Processing and Analysis Platform for Raw Refined-TPP Quantitative Proteomics

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
  [![Lab](https://img.shields.io/badge/Lab-Haojie%20Lu%20Lab-emerald.svg)](https://chemistry.fudan.edu.cn/89/1f/c45942a690463/page.htm)
  [![Venue](https://img.shields.io/badge/JACS-2025-orange.svg)](https://pubs.acs.org/doi/10.1021/jacs.5c08065)
</div>

---

## 🧭 Decision Matrix (8-Branch Architecture)

Our platform automatically routes your data through a **Triple-Layer Logic Gate**. Select your branch based on omics type and replicate design.

```mermaid
graph TD
    %% Node Definitions
    Start((<b>START</b>)) --> Type{<b>Omics Mode</b>}
    
    Type -- "Proteomics" --> P_Bio{Bio Reps?}
    Type -- "N-Glyco/Mito" --> G_Bio{Bio Reps?}
    
    %% Proteomics Sub-tree
    P_Bio -- Yes --> P_Tech_Y{Tech Reps?}
    P_Bio -- No --> P_Tech_N{Tech Reps?}
    
    P_Tech_Y -- Yes --> B1["<b>Branch 1</b><br>Full Design"]
    P_Tech_Y -- No  --> B2["<b>Branch 2</b><br>Bio-Driven"]
    P_Tech_N -- Yes --> B3["<b>Branch 3</b><br>Tech-Validated"]
    P_Tech_N -- No  --> B4["<b>Branch 4</b><br>Rapid Screen"]
    
    %% Glyco Sub-tree
    G_Bio -- Yes --> G_Tech_Y{Tech Reps?}
    G_Bio -- No --> G_Tech_N{Tech Reps?}
    
    G_Tech_Y -- Yes --> B5["<b>Branch 5</b><br>Glyco-Mapping"]
    G_Tech_Y -- No  --> B6["<b>Branch 6</b><br>Site-Diff"]
    G_Tech_N -- Yes --> B7["<b>Branch 7</b><br>Stability Eval"]
    G_Tech_N -- No  --> B8["<b>Branch 8</b><br>ID Only"]

    %% --- Styles ---
    classDef default fill:#fff,stroke:#cbd5e1,stroke-width:1px,color:#334155,font-family:arial;
    classDef decision fill:#f8fafc,stroke:#94a3b8,stroke-width:2px;
    classDef branch fill:#ffffff,stroke:#e2e8f0,stroke-width:1px;
    classDef highlight fill:#eef2ff,stroke:#6366f1,stroke-width:2px,color:#4338ca;

    class Type,P_Bio,G_Bio,P_Tech_Y,P_Tech_N,G_Tech_Y,G_Tech_N decision;
    class B1,B2,B4,B5,B6,B7,B8 branch;
    class B3 highlight;
    style Start fill:#1e293b,color:#fff,stroke:#0f172a
