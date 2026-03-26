## 💡 Description
(Write your description here...)

## 🔍 Type of Change
- [ ] 🐛 Bug fix (non-breaking change which fixes an issue)
- [ ] ✨ New feature (non-breaking change which adds functionality)
- [ ] 💥 Breaking change (fix or feature that would cause existing multi-branch logic to fail)
- [ ] 📚 Documentation or UI update

## 🧪 Pipeline Stability & Validation
**Data Processing Checks:**
- [ ] I have tested the changes using **bulk proteome** datasets (if applicable).
- [ ] I have tested the changes using **N-glycosite** datasets (if applicable).
    - *Crucial Check: If modifying N-glycosylation site retrieval, I have ensured the logic correctly processes both the database search result files and `.fasta` sequence files simultaneously.*
- [ ] I have verified that the **CV (Coefficient of Variation) quality control** step remains in the correct execution sequence (it must occur *after* filtering for unique peptides, but *before* differential analysis).

**General Checks:**
- [ ] The code follows the existing style guidelines of the project.
- [ ] I have updated the documentation (e.g., Markdown files, Mermaid decision trees) if necessary.

## 📝 Additional Notes
(Write your notes here...)
