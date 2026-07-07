#!/usr/bin/env python3
"""
M7 Plot: standalone analysis + figures from existing MD trajectory.
Usage:  python M7_plot.py              # xvg→PDF only (fast)
        python M7_plot.py --hbond      # also run H-bond analysis (needs GROMACS)
"""
import sys, subprocess, time
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUTPUT_DIR = Path(__file__).parent.resolve() / "M7_output"
GMX = "gmx"


def parse_xvg(path):
    x, y, xlabel, ylabel, title = [], [], "Time (ns)", "", ""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"): continue
            if line.startswith("@"):
                if "xaxis" in line and "label" in line:
                    xlabel = line.split('"')[1] if '"' in line else xlabel
                elif "yaxis" in line and "label" in line:
                    ylabel = line.split('"')[1] if '"' in line else ylabel
                elif "title" in line:
                    title = line.split('"')[1] if '"' in line else title
                continue
            if not line: continue
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return np.array(x), np.array(y), xlabel, ylabel, title


def plot_single(x, y, xlabel, ylabel, title, outpath, color="#2c3e50"):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, y, color=color, linewidth=0.8, alpha=0.9)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle="--")
    y_mean, y_std = np.mean(y), np.std(y)
    stats = f"Mean = {y_mean:.2f} | SD = {y_std:.2f}"
    ax.text(0.98, 0.95, stats, transform=ax.transAxes, fontsize=9,
            ha="right", va="top", bbox=dict(boxstyle="round,pad=0.3",
            facecolor="white", alpha=0.8, edgecolor="gray"))
    fig.tight_layout()

    # Save to temp then rename (avoids Windows file-lock races)
    tmp = outpath.with_suffix(".tmp")
    fig.savefig(tmp, dpi=150, bbox_inches="tight")
    plt.close("all")

    if outpath.exists():
        outpath.unlink()
    tmp.rename(outpath)
    time.sleep(0.1)

    print(f"  [OK] {outpath.name} ({outpath.stat().st_size/1024:.0f} KB)")
    return y_mean, y_std


def run_hbond():
    """Protein-ligand H-bond analysis via GROMACS."""
    print("\n[HBOND] Running GROMACS hbond...")
    hb_num = OUTPUT_DIR / "hbnum.xvg"
    tpr = OUTPUT_DIR / "prod.tpr"
    xtc = OUTPUT_DIR / "prod.xtc"

    if not tpr.exists() or not xtc.exists():
        print("  [SKIP] prod.tpr or prod.xtc missing, cannot run H-bond analysis")
        return None

    # Write stdin to temp file for multi-line group selection
    stdin_file = OUTPUT_DIR / "_stdin.txt"

    # Try Protein(1) + UNL(14)
    stdin_file.write_text("1\n14\n", encoding="ascii")
    cmd = f'cmd.exe /c "type {stdin_file.name} | {GMX} hbond -s prod.tpr -f prod.xtc -num hbnum.xvg"'
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                            cwd=str(OUTPUT_DIR), shell=True)
    stdin_file.unlink()

    if hb_num.exists() and hb_num.stat().st_size > 0:
        return hb_num

    # Fallback: Protein + System
    print("  [HBOND] Retrying with Protein/System...")
    stdin_file.write_text("Protein\nSystem\n", encoding="ascii")
    cmd = f'cmd.exe /c "type {stdin_file.name} | {GMX} hbond -s prod.tpr -f prod.xtc -num hbnum.xvg"'
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                            cwd=str(OUTPUT_DIR), shell=True)
    stdin_file.unlink()

    if hb_num.exists() and hb_num.stat().st_size > 0:
        return hb_num

    print("  [HBOND] Failed to generate hbnum.xvg")
    return None


def main():
    do_hbond = "--hbond" in sys.argv

    print("=" * 55)
    print("  M7 Plot: MD Analysis + Figures")
    print(f"  Output: {OUTPUT_DIR}")
    print("=" * 55)

    required = ["rmsd.xvg", "rmsf.xvg", "rg.xvg", "sasa.xvg"]
    missing = [f for f in required if not (OUTPUT_DIR / f).exists()]
    if missing:
        print(f"[ERROR] Missing xvg files: {missing}")
        print("  Run GROMACS analysis first, or ensure M7_output contains xvg files.")
        sys.exit(1)

    stats = {}

    # ---- RMSD ----
    print("\n[1/4] RMSD...")
    x, y, xl, yl, _ = parse_xvg(OUTPUT_DIR / "rmsd.xvg")
    stats["RMSD"] = plot_single(x, y, xl, "RMSD (nm)",
                                "GBP1 Backbone RMSD", OUTPUT_DIR / "rmsd.pdf", "#e74c3c")

    # ---- RMSF ----
    print("[2/4] RMSF...")
    x, y, xl, yl, _ = parse_xvg(OUTPUT_DIR / "rmsf.xvg")
    stats["RMSF"] = (np.max(y), np.argmax(y) if len(y) > 0 else 0)
    plot_single(x, y, "Residue Index", "RMSF (nm)",
                "GBP1 Residue Fluctuation", OUTPUT_DIR / "rmsf.pdf", "#8e44ad")

    # ---- Rg ----
    print("[3/4] Rg...")
    x, y, xl, yl, _ = parse_xvg(OUTPUT_DIR / "rg.xvg")
    stats["Rg"] = plot_single(x, y, "Time (ns)", yl,
                              "GBP1 Radius of Gyration", OUTPUT_DIR / "rg.pdf", "#2980b9")

    # ---- SASA ----
    print("[4/4] SASA...")
    x, y, xl, yl, _ = parse_xvg(OUTPUT_DIR / "sasa.xvg")
    stats["SASA"] = plot_single(x, y, "Time (ns)", yl,
                                "GBP1 Solvent Accessible Surface Area", OUTPUT_DIR / "sasa.pdf", "#d35400")

    # ---- H-bonds (optional) ----
    hbond_stats = None
    if do_hbond:
        hb_num = run_hbond()
        if hb_num:
            x, y, _, _, _ = parse_xvg(hb_num)
            hbond_stats = plot_single(x, y, "Time (ns)", "Number of H-bonds",
                                      "GBP1–Resveratrol Hydrogen Bonds",
                                      OUTPUT_DIR / "hbonds.pdf", "#27ae60")

    # ---- Combined 4-panel ----
    print("\n[PANEL] Combined QC figure...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (title, fname, color) in zip(axes.flatten(), [
        ("RMSD", "rmsd.xvg", "#e74c3c"),
        ("RMSF", "rmsf.xvg", "#8e44ad"),
        ("Rg", "rg.xvg", "#2980b9"),
        ("SASA", "sasa.xvg", "#d35400")]):
        x, y, xl, yl, _ = parse_xvg(OUTPUT_DIR / fname)
        ax.plot(x, y, color=color, linewidth=0.6)
        ax.set_title(title, fontweight="bold", fontsize=13)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.3, linestyle="--")
    fig.suptitle("GBP1–Resveratrol · 50ns MD Quality Metrics", fontsize=14, fontweight="bold")
    fig.tight_layout()
    panel_path = OUTPUT_DIR / "M7_QC_panel.pdf"
    fig.savefig(panel_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {panel_path.name} ({panel_path.stat().st_size/1024:.0f} KB)")

    # ---- Summary ----
    summary = OUTPUT_DIR / "M7_summary.txt"
    lines = [
        "=" * 50,
        "M7: GBP1 + Resveratrol · 50ns MD Analysis Summary",
        "=" * 50,
        f"RMSD (C-alpha):  mean = {stats['RMSD'][0]:.3f} nm, SD = {stats['RMSD'][1]:.3f} nm",
        f"RMSF (max):      {stats['RMSF'][0]:.3f} nm at residue {stats['RMSF'][1]}",
        f"Rg:              mean = {stats['Rg'][0]:.3f} nm, SD = {stats['Rg'][1]:.3f} nm",
        f"SASA:            mean = {stats['SASA'][0]:.1f} nm², SD = {stats['SASA'][1]:.1f} nm²",
        "=" * 50,
        "Output: rmsd.pdf  rmsf.pdf  rg.pdf  sasa.pdf  M7_QC_panel.pdf",
    ]
    if hbond_stats:
        lines.insert(-2, f"H-bonds:         mean = {hbond_stats[0]:.2f}, SD = {hbond_stats[1]:.2f}")
        lines[-1] += "  hbonds.pdf"
    summary.write_text("\n".join(lines), encoding="utf-8")

    print("\n" + "\n".join(lines))
    print("\n[DONE] All figures generated.")


if __name__ == "__main__":
    main()
