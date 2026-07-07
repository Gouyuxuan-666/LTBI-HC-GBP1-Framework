# LTBI-HC-GBP1 Computational Framework

[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![R](https://img.shields.io/badge/R-4.6.0-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.13-yellow)](https://www.python.org/)
[![GROMACS](https://img.shields.io/badge/GROMACS-2020.6-orange)](https://www.gromacs.org/)

**Author**: [Gou Yuxuan](https://orcid.org/0009-0000-4213-0274), Xinjiang Second Medical College

## Key Findings

| Evidence Layer | Module | Key Result |
|---------------|--------|------------|
| ML Feature Ranking | M5 | GBP1 #1 SHAP feature across 113 models (best AUC=0.99) |
| Co-expression Network | M3 | Hub gene in IFN-responsive module |
| Immune Infiltration | M4 | Correlated with M1 macrophage abundance (ρ=0.71) |
| Molecular Docking | M6 | GBP1-Resveratrol -7.6 kcal/mol |
| MD Simulation (50ns) | M7 | Complex stable, RMSD=0.335nm, vdW-driven binding |
| MM-GBSA | M8 | ΔE_vdw=-168 kJ/mol, no H-bonds |
| Virtual Knockout | M9 | 5,001 genes perturbed upon GBP1 ablation |
| Gene Regulatory Network | M10 | STAT1→GBP1 dominant edge (importance=0.953) |
| Bayesian Deconvolution | M11 | GBP1 localized to alveolar macrophages |
| Population-Scale Census | M12 | Confirmed across >61M cells |

## Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────────┐
│                    TIER 1: Transcriptomic Foundation              │
│  M1: Data Integration → M2: DEG+Enrichment → M3: WGCNA → M4: CIBS│
├──────────────────────────────────────────────────────────────────┤
│                  TIER 2: ML & Structural Biology                  │
│   M5: 113 ML + SHAP → M6: Docking → M7: MD 50ns → M8: MMPBSA     │
├──────────────────────────────────────────────────────────────────┤
│                TIER 3: Single-Cell & Network Analysis             │
│        M9: Virtual KO → M10: GRN → M11: BayesPrism               │
├──────────────────────────────────────────────────────────────────┤
│                 TIER 4: Population-Scale Validation               │
│                    M12: CELLxGENE Census                          │
└──────────────────────────────────────────────────────────────────┘
```

## Repository Files

| File | Module | Description |
|------|--------|-------------|
| `01_download_and_merge.R` | M1 | 10 GEO dataset download, probe-to-gene mapping, ComBat batch correction |
| `M6_docking.py` | M6 | AutoDock Vina molecular docking (GBP1 vs 4 anti-TB drugs + resveratrol) |
| `M7_md.py` | M7 | GROMACS full MD pipeline (amber99sb-ildn, 50ns, GPU-accelerated) |
| `M7_plot.py` | M7 | Standalone MD QC plotting (RMSD/RMSF/Rg/SASA/H-bonds) |
| `M8_mmpbsa.py` | M8 | Manual MM-GBSA scheme (mdrun -rerun + energygrps decomposition) |
| `gen_ligand_itp.py` | M7 | GAFF2 ligand topology generation (Kabsch alignment + AM1-BCC charges) |
| `M10_visualize.py` | M10 | GRNBoost2 gene regulatory network visualization |
| `M12_census.py` | M12 | CELLxGENE Census REST API cross-dataset validation (61M+ cells) |

## Requirements

### R (v4.6.0+)

```r
install.packages(c("GEOquery", "limma", "sva", "WGCNA", "clusterProfiler",
    "fgsea", "GSVA", "caret", "fastshap", "shapviz"))
# CIBERSORT: https://cibersort.stanford.edu/
```

### Python (v3.13+)

```bash
pip install numpy matplotlib seaborn networkx arboreto scanpy cellxgene-census
```

### Molecular Simulation

| Software | Version | Notes |
|----------|---------|-------|
| GROMACS | 2020.6 | GPU-accelerated (GTX 1080 tested) |
| AutoDock Vina | 1.2.0 | |

## Quick Start

```bash
# 1. Clone
git clone https://github.com/Gouyuxuan-666/LTBI-HC-GBP1-Framework.git
cd LTBI-HC-GBP1-Framework

# 2. Data acquisition
Rscript 01_download_and_merge.R

# 3. Molecular docking (M6)
python M6_docking.py

# 4. MD simulation (M7) — requires GROMACS + GPU
python M7_md.py --prod    # 50ns production
python M7_plot.py         # QC analysis + figures

# 5. MMPBSA (M8)
python M8_mmpbsa.py

# 6. GRN visualization (M10)
python M10_visualize.py

# 7. Census validation (M12)
python M12_census.py
```

## References

1. Ye GM, Shen SJ, Zhang B, et al. Expression and clinical significance of GBP1 in pulmonary tuberculosis. *Anhui Med Univ*. 2023; 58(2): 214-218. DOI: [10.19405/j.cnki.issn1000-1492.2023.02.007](https://doi.org/10.19405/j.cnki.issn1000-1492.2023.02.007)

2. Ye GM. Effects of Sp110 targeting JAK-STAT signaling pathway on macrophage function [D]. Shihezi University, 2023. (Advisor: Prof. Wu Jiang-dong)

## License

MIT License — see [LICENSE](LICENSE)
