# Evolutionary Rate Variation at Predicted PTM Sites Highlights Localized Host-Associated Signals in Influenza A Virus

**Teufel AI, Hasan MI, Smyth DS** — *Genome Biology and Evolution*, 2026

This repository contains all scripts, input data, and results for the analysis described in the manuscript. 

---

## Overview

We used [MusiteDeep](https://www.musite.net/) to predict post-translational modification (PTM) sites across Influenza A Virus (IAV) proteins from H1N1, H5N1, and H7N9. PTM states were treated as discrete evolutionary characters and analyzed using Bayesian phylogenetic models in [RevBayes](https://revbayes.github.io/) to estimate site-specific and branch-specific evolutionary rates. Protein structures were predicted with AlphaFold and visualized in ChimeraX.

---

## Repository Structure

```
Flu_PTM_github_reorganized/
├── scripts/
│   ├── pipeline/        # Python pipeline: sequence retrieval → alignment → PTM prediction → RevBayes input
│   ├── revbayes/        # RevBayes (.Rev) scripts for each protein (final converged runs)
│   ├── analysis/        # R scripts for rate analysis, host-shift tests, phylogenetic signal
│   └── figures/         # R scripts used to generate manuscript figures
├── data/
│   ├── raw_sequences/   # Original FASTA sequences per strain (H1N1, H5N1, H7N9), organized by segment
│   ├── alignments/      # MAFFT-aligned FASTA files per protein
│   ├── trees/           # IQ-TREE maximum likelihood trees (.nex, .treefile) per protein
│   └── ptm/
│       ├── predictions/ # MusiteDeep PTM predictions per protein (*_ptm_predictions.txt)
│       ├── multistate/  # Discrete character matrices for RevBayes (*_multistate.txt)
│       └── mappings/    # Position mapping files (alignment ↔ protein coordinates), including paper fast/slow site definitions
├── results/
│   ├── site_analysis/   # Per-protein site rate results (CSVs, PDFs), plus combined summaries
│   ├── host_shift/      # Host-shift branch rate comparisons per protein
│   ├── convergence_plots/ # MCMC trace plots for all proteins (ESS diagnostics)
│   └── summary_figures/ # Combined summary figures and tables
├── structures/
│   ├── pdb/             # AlphaFold-predicted structures (best-ranked model, consensus sequence input)
│   └── chimerax_scripts/ # ChimeraX (.cxc/.cxs) session files and annotation scripts
├── extra_analyses/
│   ├── ecological_analysis/ # PTM variation across geography, host type, and time (not in manuscript)
│   └── lambda_analysis/     # Pagel's lambda phylogenetic signal analysis (not in manuscript
```
---
