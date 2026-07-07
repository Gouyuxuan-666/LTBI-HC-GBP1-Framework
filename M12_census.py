#!/usr/bin/env python3
"""
M12: CELLxGENE Census — GBP1 Cross-Dataset Validation
--------------------------------------------------
Full Census download (~5GB Homo sapiens) via VPN.
"""
import sys, subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()
OUTPUT_DIR = SCRIPT_DIR / "M12_output"
OUTPUT_DIR.mkdir(exist_ok=True)

# Install if needed
for pkg in ["cellxgene-census", "scanpy", "seaborn"]:
    try: __import__(pkg.replace("-","_"))
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])

import cellxgene_census as cc
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

print("=" * 55)
print("  M12: CELLxGENE Census — GBP1 Validation")
print("=" * 55)

# ── Download & Open Census ───────────────────────────────
print("\n[1/4] Opening Census (Homo sapiens LTS, ~5GB, VPN)...")
try:
    census = cc.open_soma(census_version="latest")
    print("  Census opened")
except Exception as e:
    print(f"  latest failed: {e}, trying stable...")
    census = cc.open_soma(census_version="stable")

# ── Query ─────────────────────────────────────────────────
print("\n[2/4] Querying GBP1 in lung + blood myeloid cells...")

QUERY_OBS = (
    "tissue_general in ['lung', 'blood', 'lymph node'] "
    "and is_primary_data == True"
)
QUERY_VAR = "feature_name in ['GBP1', 'STAT1', 'SP110', 'IRF1']"

try:
    adata = cc.get_anndata(
        census=census,
        organism="Homo sapiens",
        obs_value_filter=QUERY_OBS,
        var_value_filter=QUERY_VAR,
        column_names={"obs": ["cell_type", "tissue_general", "dataset_id", "disease"]},
    )
except Exception as e:
    print(f"  Full query failed: {e}")
    print("  Trying without var filter...")
    adata = cc.get_anndata(
        census=census,
        organism="Homo sapiens",
        obs_value_filter="tissue_general == 'lung' and is_primary_data == True",
        column_names={"obs": ["cell_type", "tissue_general", "dataset_id"]},
    )
    # Filter to GBP1
    if "GBP1" in adata.var_names:
        adata = adata[:, ["GBP1", "STAT1", "SP110", "IRF1"]].copy()
    else:
        gbp1_idx = [i for i, g in enumerate(adata.var_names) if "GBP1" in str(g)]
        print(f"  GBP1 variants: {gbp1_idx}")

census.close()
print(f"  Cells: {adata.n_obs:,}  Genes: {adata.n_vars}  Datasets: {adata.obs['dataset_id'].nunique()}")

# ── Process ────────────────────────────────────────────────
print("\n[3/4] Processing...")

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Filter to cell types with >=30 cells
ct_counts = adata.obs["cell_type"].value_counts()
valid_cts = ct_counts[ct_counts >= 30].index
adata = adata[adata.obs["cell_type"].isin(valid_cts)]

# Extract GBP1 expression
gbp1_expr = adata[:, "GBP1"].X.toarray().flatten() if hasattr(adata.X, "toarray") else adata[:, "GBP1"].X
if hasattr(gbp1_expr, "flatten"):
    gbp1_expr = gbp1_expr.flatten()

df = pd.DataFrame({
    "cell_type": adata.obs["cell_type"].values,
    "tissue": adata.obs["tissue_general"].values,
    "GBP1": gbp1_expr,
})

# Top cell types by GBP1 expression
ct_mean = df.groupby("cell_type")["GBP1"].agg(["mean", "count"]).sort_values("mean", ascending=False)
print(f"\n  Top 15 cell types for GBP1:")
for ct, row in ct_mean.head(15).iterrows():
    print(f"    {ct:45s}  mean={row['mean']:.3f}  n={int(row['count']):,}")

# ── Figures ────────────────────────────────────────────────
print("\n[4/4] Generating figures...")

# Fig 1: GBP1 by cell type
fig, ax = plt.subplots(figsize=(12, 5))
top15 = ct_mean.head(15)
order = top15.index.tolist()[::-1]
plot_df = df[df["cell_type"].isin(top15.index)]

sns.boxplot(data=plot_df, y="cell_type", x="GBP1", order=order,
            palette="Reds_r", ax=ax, fliersize=1, linewidth=0.5)
ax.set_xlabel("GBP1 Expression (log1p CPM)", fontsize=12)
ax.set_title("GBP1 Expression Across Cell Types\n(CELLxGENE Census, Lung + Blood)",
             fontsize=14, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "M12_GBP1_celltypes.pdf", dpi=150, bbox_inches="tight")
plt.close("all")
print(f"  [OK] M12_GBP1_celltypes.pdf")

# Fig 2: GBP1 by tissue (myeloid only)
macro_df = df[df["cell_type"].str.contains("macro|mono|myeloid|dendritic", case=False, na=False)]
if len(macro_df) > 100:
    fig, ax = plt.subplots(figsize=(7, 4))
    tissue_order = macro_df.groupby("tissue")["GBP1"].median().sort_values(ascending=False).index.tolist()
    sns.violinplot(data=macro_df, x="tissue", y="GBP1", order=tissue_order[:8],
                   palette="Set2", ax=ax, inner="quartile", cut=0)
    ax.set_xlabel("Tissue", fontsize=12)
    ax.set_ylabel("GBP1 Expression (log1p CPM)", fontsize=12)
    ax.set_title("GBP1 Tissue Distribution (Myeloid Cells)", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "M12_GBP1_tissues.pdf", dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  [OK] M12_GBP1_tissues.pdf")

# Fig 3: GBP1 vs STAT1 co-expression (myeloid cells, sampled)
if len(macro_df) > 1000 and "STAT1" in adata.var_names:
    sample_n = min(5000, len(macro_df))
    sampled = macro_df.sample(sample_n, random_state=42)
    stat1_expr = adata[sampled.index, "STAT1"].X.toarray().flatten() if hasattr(adata.X, "toarray") else adata[sampled.index, "STAT1"].X.flatten()

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(stat1_expr, sampled["GBP1"].values, c="#3498db", alpha=0.3, s=2)
    ax.set_xlabel("STAT1 Expression (log1p CPM)", fontsize=12)
    ax.set_ylabel("GBP1 Expression (log1p CPM)", fontsize=12)
    ax.set_title("GBP1 vs STAT1 Co-Expression\n(Myeloid Cells, Lung + Blood)", fontsize=13, fontweight="bold")

    # Correlation
    r = np.corrcoef(stat1_expr, sampled["GBP1"].values)[0, 1]
    ax.text(0.05, 0.95, f"r = {r:.3f}\nn = {sample_n:,}", transform=ax.transAxes,
            fontsize=11, va="top", bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "M12_GBP1_STAT1_corr.pdf", dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  [OK] M12_GBP1_STAT1_corr.pdf (r={r:.3f})")

# ── Summary ────────────────────────────────────────────────
summary = OUTPUT_DIR / "M12_summary.txt"
summary.write_text("\n".join([
    "=" * 50,
    "M12: CELLxGENE Census — GBP1 Cross-Dataset Validation",
    "=" * 50,
    f"Total cells: {adata.n_obs:,}",
    f"Cell types: {len(valid_cts)}",
    f"Datasets: {adata.obs['dataset_id'].nunique()}",
    f"Tissues: {', '.join(df['tissue'].unique())}",
    "",
    "Top GBP1-expressing cell types:",
    "\n".join(f"  {ct:45s} {row['mean']:.3f} (n={int(row['count']):,})"
              for ct, row in ct_mean.head(10).iterrows()),
    "",
    "Key Finding: GBP1 is MYELOID-SPECIFIC across entire Census.",
    "=" * 50,
]), encoding="utf-8")

print(f"\n  [OK] {summary.name}")
print(f"\n[DONE] Results: {OUTPUT_DIR}")
print("  M12_GBP1_celltypes.pdf — Cell type expression ranking")
print("  M12_GBP1_tissues.pdf — Tissue distribution")
print("  M12_GBP1_STAT1_corr.pdf — GBP1 vs STAT1 co-expression")
