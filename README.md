# A 12-Module Computational Framework for Tuberculosis Biomarker Discovery

[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![R](https://img.shields.io/badge/R-4.6.0-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.13-yellow)](https://www.python.org/)
[![GROMACS](https://img.shields.io/badge/GROMACS-2020.6-orange)](https://www.gromacs.org/)
[![Status](https://img.shields.io/badge/Status-Manuscript%20in%20Preparation-lightgrey)]()

## 概述 / Overview

本研究构建了一个**完全可复现的12模块计算框架**，系统整合多队列转录组学、集成机器学习（113种模型）、分子动力学模拟（50ns）、单细胞分析和基因调控网络推断，最终从10个GEO数据集中一致地鉴定出 **GBP1（鸟苷酸结合蛋白1）** 作为结核病（TB）最稳健的宿主生物标志物和潜在药物靶点。

We present a **fully reproducible 12-module computational framework** that systematically integrates multi-cohort transcriptomics (10 GEO datasets, >800 samples), systematic ML benchmarking (113 models across 20 algorithm families), all-atom MD simulation (50 ns), MM-GBSA binding free energy calculation, single-cell virtual gene knockout, genome-wide gene regulatory network inference, and population-scale validation. All modules converge on **GBP1 (Guanylate-Binding Protein 1)** as the most robust biomarker for distinguishing latent tuberculosis infection (LTBI) from healthy controls.

---

## 核心发现 / Key Finding

> **GBP1 is the single most robust TB biomarker** — confirmed by 12 independent analytical layers with convergent evidence.

```
            IFN-γ
              ↓
           STAT1 ──(r=0.953, p=5.7×10⁻⁹⁴)→ GBP1
              ↓                                ↓
          IRF9, STAT2                    5,001 genes perturbed
          OAS1, MX1, ISG15               MALAT1, MT-ND2/3
```

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

---

## 管线架构 / Pipeline Architecture

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

---

## 仓库文件 / Repository Structure

| 文件 | 模块 | 说明 |
|------|------|------|
| `01_download_and_merge.R` | M1 | 10 GEO数据集下载、探针→基因映射、ComBat批次校正 |
| `M6_docking.py` | M6 | AutoDock Vina分子对接（GBP1 vs 4种抗TB药+白藜芦醇） |
| `M7_md.py` | M7 | GROMACS全流程MD模拟（amber99sb-ildn力场，50ns，GPU加速） |
| `M7_plot.py` | M7 | 独立MD质控出图脚本（RMSD/RMSF/Rg/SASA/H键） |
| `M8_mmpbsa.py` | M8 | 手动MM-GBSA方案（mdrun -rerun + energygrps能量分解） |
| `gen_ligand_itp.py` | M7 | GAFF2配体拓扑文件生成（Kabsch对齐+AM1-BCC电荷） |
| `M10_visualize.py` | M10 | GRNBoost2基因调控网络可视化（STAT1→GBP1网络图） |
| `M12_census.py` | M12 | CELLxGENE Census REST API 61M细胞跨数据集验证 |
| `M6_M8_报错日志.md` | — | M6-M8详细报错记录（10个issue及解决方案） |
| `LTBI_HC_论文写作计划.docx` | — | 论文写作计划书（含10期刊投稿路线图） |

---

## 环境要求 / Requirements

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
| GROMACS | 2020.6 | GPU-accelerated (GTX 1080 tested, 37.3 ns/day) |
| AutoDock Vina | 1.2.0 | |

---

## 快速开始 / Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/Gouyuxuan-666/-LTBI-HC-GBP1-Framework.git
cd -LTBI-HC-GBP1-Framework

# 2. Data acquisition & ML pipeline
Rscript 01_download_and_merge.R
# ... (M2-M5 R scripts, available upon request)

# 3. Molecular docking (M6)
python M6_docking.py

# 4. MD simulation (M7) — requires GROMACS + GPU
python M7_md.py --prod    # 50ns production
python M7_plot.py --hbond  # QC analysis + figures

# 5. MMPBSA (M8)
python M8_mmpbsa.py

# 6. GRN visualization (M10)
python M10_visualize.py

# 7. Census validation (M12) — requires internet
python M12_census.py
```

---

## 论文 / Manuscript

**Title**: A Computational Framework Integrating Multi-Cohort Transcriptomics, Ensemble Machine Learning, and Molecular Dynamics Identifies GBP1 as a Tuberculosis Biomarker

**Status**: Manuscript in preparation (v2 complete, figure embedding in progress)

**Target Journal**: *Briefings in Bioinformatics* (Oxford, IF=9.5, JCR Q1)

**Writing Plan**: See `LTBI_HC_论文写作计划.docx` (10-journal submission roadmap, 4-tier strategy)

---

## 致谢 / Acknowledgments

- 湿实验数据 (IHC/WB/qRT-PCR): 叶国敏 (2023) — GBP1/SP110 巨噬细胞功能验证
- MD计算资源: NVIDIA GTX 1080 GPU
- 部分分析模块: 阿里云 PAI-DSW 平台

---

## 许可 / License

MIT License — 详见 [LICENSE](LICENSE)

---

## 联系 / Contact

GitHub: [@Gouyuxuan-666](https://github.com/Gouyuxuan-666)

Research Direction: 感染与免疫 (Infection & Immunity)
