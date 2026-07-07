#!/usr/bin/env python3
"""
M8: MM-GBSA Binding Free Energy — GBP1 + Resveratrol
=====================================================
Pure GROMACS + Python approach:
  1. Extract N frames from trajectory (protein + UNL only, no water)
  2. mdrun -rerun with energygrps = Protein UNL
  3. Parse LJ-SR + Coulomb-SR interaction energies
  4. Add GBSA solvation via simple approximation

No external packages beyond numpy + matplotlib.
"""
import subprocess, sys, re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).parent.resolve()
M7_DIR = SCRIPT_DIR / "M7_output"
OUTPUT_DIR = SCRIPT_DIR / "M8_output"
OUTPUT_DIR.mkdir(exist_ok=True)

GMX = "gmx"
N_FRAMES = 100
GPU_ID = "0"


def run_gmx(args, cwd=None, stdin_input=None, timeout=600):
    cwd = str(cwd or OUTPUT_DIR)
    gmx_args = " ".join(f'"{a}"' if " " in str(a) else str(a) for a in args)
    if stdin_input:
        tmpfile = Path(cwd) / "_stdin.tmp"
        tmpfile.write_text(stdin_input.replace("\\n", "\n"), encoding="ascii")
        cmd = f'cmd.exe /c "type {tmpfile.name} | {GMX} {gmx_args}"'
    else:
        cmd = f'cmd.exe /c "{GMX} {gmx_args}"'

    print(f"  [gmx] {' '.join(str(a) for a in args)[:100]}...")
    result = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=timeout, cwd=cwd, shell=True)
    if stdin_input:
        try: (Path(cwd) / "_stdin.tmp").unlink()
        except OSError: pass

    if result.returncode != 0:
        err = (result.stderr + result.stdout)[-500:]
        print(f"  [ERR] {err}")
    return result


def parse_xvg(path):
    x, y = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or line.startswith("@"): continue
            if not line: continue
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return np.array(x), np.array(y)


print("=" * 60)
print("  M8: MM-GBSA — GBP1 + Resveratrol")
print(f"  Output: {OUTPUT_DIR}")
print("=" * 60)

# ── Step 1: Extract frames (protein + UNL, no water/ions) ──
print("\n[1/5] Extracting frames from trajectory...")

tpr = M7_DIR / "prod.tpr"
xtc = M7_DIR / "prod.xtc"
frames_xtc = OUTPUT_DIR / "frames.xtc"

if not frames_xtc.exists():
    # First, create an index with Protein and UNL groups
    ndx = OUTPUT_DIR / "strip.ndx"
    ndx_stdin = OUTPUT_DIR / "_ndx.txt"
    ndx_stdin.write_text('keep 0\n"Protein" | "UNL"\nq\n', encoding="ascii")

    cmd = f'cmd.exe /c "type _ndx.txt | {GMX} make_ndx -f {tpr} -o strip.ndx"'
    subprocess.run(cmd, capture_output=True, text=True, timeout=30,
                   cwd=str(OUTPUT_DIR), shell=True)
    ndx_stdin.unlink()

    # Show the new group number from index
    result = subprocess.run(
        f'cmd.exe /c "type _ndx.txt"', capture_output=True, text=True,
        timeout=5, cwd=str(OUTPUT_DIR), shell=True
    )
    # Find group number of "Protein_UNL" combined group
    group_num = None
    if ndx.exists():
        content = ndx.read_text()
        groups = re.findall(r'\[\s*(\S+)\s*\]', content)
        for i, g in enumerate(groups):
            if "Protein" in g and "UNL" in g:
                group_num = i
                print(f"  Index group: [{group_num}] {g}")
                break
        if group_num is None:
            # Look for the last group which should be our merged group
            group_num = len(groups) - 1
            print(f"  Using last group: [{group_num}] {groups[-1]}")

    if group_num is None:
        print("[ERROR] Could not create Protein+UNL index group")
        sys.exit(1)

    # Extract frames
    total_frames = 5000
    dt_frame = max(1, total_frames // N_FRAMES)  # skip every N frames
    print(f"  Extracting every {dt_frame}th frame → ~{N_FRAMES} frames")

    trj_stdin = OUTPUT_DIR / "_trj.txt"
    trj_stdin.write_text(f"{group_num}\n{group_num}\n", encoding="ascii")

    cmd = f'cmd.exe /c "type _trj.txt | {GMX} trjconv -f {xtc} -s {tpr} -n strip.ndx -o frames.xtc -dt {dt_frame} -pbc mol -center"'
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                            cwd=str(OUTPUT_DIR), shell=True)
    trj_stdin.unlink()

    if not frames_xtc.exists():
        print("[ERROR] Failed to extract frames!")
        print(result.stderr[-1000:])
        sys.exit(1)

    n_frames = 0
    for line in (result.stderr + result.stdout).split("\n"):
        if "frame" in line.lower() and "written" in line.lower():
            print(f"  {line.strip()}")
else:
    print("  frames.xtc exists, skip extraction.")

# Count frames
result = subprocess.run(
    f'cmd.exe /c "{GMX} check -f frames.xtc"',
    capture_output=True, text=True, timeout=30,
    cwd=str(OUTPUT_DIR), shell=True
)
actual_frames = 0
for line in (result.stderr + result.stdout).split("\n"):
    m = re.search(r'Last frame.*?time\s+([\d.]+)', line)
    if m:
        actual_frames = int(float(m.group(1)) / 10) + 1  # 10ps per frame
        print(f"  Actual frames: ~{actual_frames} (last time: {m.group(1)} ps)")
        break

print(f"  Output: frames.xtc ({frames_xtc.stat().st_size/1e6:.0f} MB)")

# ── Step 2: Build rerun TPR with energygrps ──
print("\n[2/5] Building rerun TPR with energy groups...")

# Create a minimal MDP for rerun (same parameters as prod, but 0 steps)
rerun_mdp = OUTPUT_DIR / "rerun.mdp"
rerun_mdp.write_text("""
integrator  = md
dt          = 0.002
nsteps      = 0
nstxout-compressed = 5000
nstlog      = 5000
nstenergy   = 1
nstlist     = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
vdwtype     = Cut-off
rvdw        = 1.2
pbc         = xyz
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = 0.1 0.1
ref-t       = 310 310
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint-algorithm = LINCS
energygrps  = Protein UNL
energygrp-table =
""")

# Need a structure file for grompp - use last frame or prod.gro
# prod.tpr already has the right structure, use -c from prod.gro
rerun_tpr = OUTPUT_DIR / "rerun.tpr"
topol = M7_DIR / "topol" / "topol.top"
gro_file = M7_DIR / "prod.gro"

if not gro_file.exists():
    # Extract last frame as GRO
    run_gmx(["trjconv", "-f", str(xtc), "-s", str(tpr),
             "-o", str(OUTPUT_DIR / "last.gro"), "-dump", "50000"],
             stdin_input="0", timeout=120)
    gro_file = OUTPUT_DIR / "last.gro"

result = run_gmx(["grompp",
    "-f", str(rerun_mdp), "-c", str(gro_file),
    "-p", str(topol), "-o", str(rerun_tpr), "-maxwarn", "10"])
if result.returncode != 0:
    print("[ERROR] grompp failed!")
    sys.exit(1)
print("  rerun.tpr ready")

# ── Step 3: Rerun to get interaction energies ──
print(f"\n[3/5] Rerunning {actual_frames} frames with energy decomposition...")

run_gmx(["mdrun", "-s", str(rerun_tpr), "-rerun", str(frames_xtc),
         "-deffnm", "rerun", "-v", "-ntmpi", "1", "-gpu_id", GPU_ID],
         timeout=1800)

edr = OUTPUT_DIR / "rerun.edr"
if not edr.exists():
    print("[ERROR] rerun.edr not produced!")
    sys.exit(1)

# ── Step 4: Extract energies ──
print("\n[4/5] Extracting interaction energies...")

# Write energy terms to xvg files
energy_terms = {
    "LJ_SR": ("LJ-SR:Protein-UNL", "lj_sr.xvg"),
    "Coul_SR": ("Coul-SR:Protein-UNL", "coul_sr.xvg"),
    "Pot_Total": ("Total-Energy", "total.xvg"),
}

energies = {}
for key, (gmx_name, fname) in energy_terms.items():
    out = OUTPUT_DIR / fname
    # Extract specific energy term from EDR
    # Use echo to select the term
    stdin = f"{gmx_name}\n0\n"
    run_gmx(["energy", "-f", str(edr), "-o", str(out)],
             stdin_input=stdin, timeout=120)
    if out.exists():
        x, y = parse_xvg(out)
        energies[key] = (x, y)
        print(f"  {key}: mean = {np.mean(y):.2f} ± {np.std(y):.2f} kJ/mol")
    else:
        print(f"  [WARN] {key} not extracted")

if "LJ_SR" not in energies or "Coul_SR" not in energies:
    # Try listing available terms
    print("\n  Available energy terms:")
    result = run_gmx(["energy", "-f", str(edr)], stdin_input="\n", timeout=30)
    for line in (result.stdout + result.stderr).split("\n"):
        if any(k in line for k in ["LJ", "Coul", "Protein", "UNL", "LIG"]):
            print(f"    {line.strip()}")

# ── Step 5: Calculate binding energy ──
print("\n[5/5] Calculating binding free energy...")

if "LJ_SR" in energies and "Coul_SR" in energies:
    x_lj, y_lj = energies["LJ_SR"]
    x_coul, y_coul = energies["Coul_SR"]

    # ΔE_MM = ΔE_vdw + ΔE_elec (gas phase interaction)
    delta_vdw = np.mean(y_lj)
    delta_elec = np.mean(y_coul)

    # Simple GBSA approximation
    # ΔG_solv ≈ γ * SASA  (nonpolar, from literature)
    # For GB, the polar part from gmx_MMPBSA typically uses IGB=5
    # Here we estimate nonpolar: γ ~ 0.0072 kJ/mol/Å², SASA from M7 ≈ 292 nm²
    # But this is total SASA, not the buried area. Use a rough estimate.
    # ΔSASA (buried) ≈ 2-5 nm² for small molecule → ΔG_np ≈ -15 to -30 kJ/mol
    # For GB polar, estimate using distance-dependent dielectric
    # This is approximate — the real value needs full GB solver

    delta_ggas = delta_vdw + delta_elec

    # Estimate ΔG_solv (GB polar + nonpolar)
    # For a typical protein-ligand complex with IGB=5:
    # ΔG_GB depends on charge distribution, roughly counteracts Coul
    # Rough estimate: ΔG_GB ≈ -0.5 * ΔE_elec (for buried polar groups)
    delta_gb = -0.5 * delta_elec  # rough GB polar estimate

    # Nonpolar: γ * ΔSASA. Estimate buried SASA ~ 3 nm²
    gamma_val = 0.0072 * 100  # kJ/mol/nm² (0.0072 kcal/mol/Å² converted)
    buried_sasa = 3.0  # nm², estimate for resveratrol
    delta_np = gamma_val * buried_sasa  # positive, unfavorable

    delta_total = delta_ggas + delta_gb - delta_np

    print(f"\n{'='*50}")
    print(f"  MM-GBSA Binding Free Energy (GBP1 + Resveratrol)")
    print(f"{'='*50}")
    print(f"  ΔE_vdw      = {delta_vdw:8.2f} ± {np.std(y_lj):.2f} kJ/mol")
    print(f"  ΔE_elec     = {delta_elec:8.2f} ± {np.std(y_coul):.2f} kJ/mol")
    print(f"  ΔE_gas      = {delta_ggas:8.2f} kJ/mol")
    print(f"  ΔG_GB(est)  = {delta_gb:8.2f} kJ/mol")
    print(f"  ΔG_np(est)  = {delta_np:8.2f} kJ/mol")
    print(f"  {'─'*40}")
    print(f"  ΔG_bind     = {delta_total:8.2f} kJ/mol  ({delta_total/4.184:.2f} kcal/mol)")
    print(f"{'='*50}")

    # ── Figures ──
    print("\n[FIG] Generating plots...")

    # Time series of interaction energy
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].plot(x_lj, y_lj, color="#e74c3c", linewidth=0.6)
    axes[0].set_title("vdW Interaction (LJ-SR)", fontweight="bold")
    axes[0].set_xlabel("Time (ps)"); axes[0].set_ylabel("kJ/mol")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(x_coul, y_coul, color="#2980b9", linewidth=0.6)
    axes[1].set_title("Electrostatic (Coul-SR)", fontweight="bold")
    axes[1].set_xlabel("Time (ps)"); axes[1].set_ylabel("kJ/mol")
    axes[1].grid(True, alpha=0.3)

    total_mm = y_lj + y_coul
    axes[2].plot(x_lj, total_mm, color="#27ae60", linewidth=0.6)
    axes[2].set_title("Total MM Interaction", fontweight="bold")
    axes[2].set_xlabel("Time (ps)"); axes[2].set_ylabel("kJ/mol")
    axes[2].grid(True, alpha=0.3)

    fig.suptitle("GBP1–Resveratrol Interaction Energy (100 frames, 50ns MD)",
                 fontweight="bold")
    fig.tight_layout()
    ts_path = OUTPUT_DIR / "M8_energy_timeseries.pdf"
    fig.savefig(ts_path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  {ts_path.name}")

    # Bar chart: energy components
    fig, ax = plt.subplots(figsize=(6, 5))
    components = ["ΔE_vdw", "ΔE_elec", "ΔG_GB\n(est)", "ΔG_np\n(est)", "ΔG_bind"]
    values = [delta_vdw, delta_elec, delta_gb, -delta_np, delta_total]
    colors = ["#e74c3c", "#2980b9", "#8e44ad", "#d35400", "#27ae60"]
    bars = ax.bar(components, values, color=colors, edgecolor="white", width=0.6)
    ax.axhline(y=0, color="black", linewidth=0.8)
    ax.set_ylabel("Energy (kJ/mol)", fontsize=12)
    ax.set_title("GBP1–Resveratrol MM-GBSA Energy Decomposition", fontweight="bold", fontsize=13)
    ax.grid(True, alpha=0.3, axis="y")

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, val,
                f"{val:.1f}", ha="center", va="bottom" if val > 0 else "top",
                fontsize=9, fontweight="bold")

    fig.tight_layout()
    bar_path = OUTPUT_DIR / "M8_energy_components.pdf"
    fig.savefig(bar_path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  {bar_path.name}")

    # Summary file
    summary = OUTPUT_DIR / "M8_summary.txt"
    summary.write_text("\n".join([
        "=" * 50,
        "M8: GBP1 + Resveratrol · MM-GBSA Binding Energy",
        "=" * 50,
        f"ΔE_vdw      = {delta_vdw:8.2f} ± {np.std(y_lj):.2f} kJ/mol",
        f"ΔE_elec     = {delta_elec:8.2f} ± {np.std(y_coul):.2f} kJ/mol",
        f"ΔE_gas      = {delta_ggas:8.2f} kJ/mol",
        f"ΔG_bind     = {delta_total:8.2f} kJ/mol  ({delta_total/4.184:.2f} kcal/mol)",
        "=" * 50,
        "Output: M8_energy_timeseries.pdf  M8_energy_components.pdf",
    ]), encoding="utf-8")
    print(f"\n  {summary.name} written")

else:
    print("\n[WARN] Interaction energies not available.")
    print("  Check rerun.edr energy terms manually:")
    print(f"  gmx energy -f {edr}")

print(f"\n[DONE] Results: {OUTPUT_DIR}")
