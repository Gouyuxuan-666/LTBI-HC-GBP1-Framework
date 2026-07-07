# LTBI-HC GBP1 Framework

A 12-module reproducible computational framework integrating multi-cohort transcriptomics, ensemble machine learning, molecular dynamics simulation, and single-cell analysis to identify GBP1 as a tuberculosis biomarker.

## Pipeline Overview

| Module | Method | Key Output |
|--------|--------|------------|
| M1 | 10 GEO datasets, ComBat harmonization | Merged expression matrix |
| M2 | limma DEG + GO/KEGG/GSEA enrichment | 347 consensus DEGs |
| M3 | WGCNA co-expression network | IFN-responsive turquoise module |
| M4 | CIBERSORT + ssGSEA immune infiltration | GBP1-M1 macrophage correlation |
| M5 | 113 ML models + SHAP interpretation | Lasso+Stepglm AUC=0.99, GBP1 top SHAP |
| M6 | AutoDock Vina molecular docking | GBP1-Resveratrol -7.6 kcal/mol |
| M7 | GROMACS 50ns MD simulation | RMSD=0.335nm, 0 H-bonds, vdW-driven |
| M8 | MM-GBSA binding free energy | ΔE_vdw=-168 kJ/mol, ΔE_elec=-8 kJ/mol |
| M9 | scTenifoldKnk virtual knockout | 5,001 genes perturbed |
| M10 | GRNBoost2 regulatory network | STAT1→GBP1 r=0.953, p=5.7×10⁻⁹⁴ |
| M11 | BayesPrism deconvolution | Alveolar macrophage localization |
| M12 | CELLxGENE Census validation | 61M cells cross-dataset confirmation |

## Key Finding

**GBP1 is the most robust TB biomarker across every analytical layer** — from ML feature importance to structural biology to network inference — with convergent evidence from 12 independent modules.

The IFN-γ → STAT1 → GBP1 axis is the central regulatory module: STAT1 is the master transcriptional regulator of GBP1 (r=0.953), consistent with known JAK-STAT signaling in TB immunity.

## Repository Structure

```
├── 01_download_and_merge.R     # M1: Data acquisition & ComBat harmonization
├── M6_docking.py               # M6: AutoDock Vina molecular docking
├── M7_md.py                    # M7: GROMACS MD simulation (50ns, amber99sb-ildn)
├── M7_plot.py                  # M7: standalone QC analysis & figure generation
├── M8_mmpbsa.py                # M8: Manual MM-GBSA (mdrun -rerun + energygrps)
├── gen_ligand_itp.py           # GAFF2 ligand topology generation
├── M10_visualize.py            # M10: GRNBoost2 network visualization
├── M12_census.py               # M12: CELLxGENE Census cross-validation
├── M6_M8_报错日志.md           # Detailed error log (10 issues resolved)
└── LTBI_HC_论文写作计划.docx   # Paper writing plan (Chinese)
```

## Requirements

- **R 4.6.0+**: GEOquery, limma, sva, WGCNA, clusterProfiler, CIBERSORT, caret, fastshap
- **Python 3.13+**: numpy, matplotlib, seaborn, networkx, arboreto, scanpy, cellxgene-census
- **GROMACS 2020.6**: GPU-accelerated (NVIDIA GTX 1080 tested)
- **AutoDock Vina 1.2.0**

## Citation

Manuscript in preparation. Target journal: *Briefings in Bioinformatics* (IF 9.5).

## License

MIT
