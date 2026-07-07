#!/usr/bin/env python3
"""
M7: GROMACS MD Simulation — GBP1 + Resveratrol
===============================================
GPU-accelerated (GTX 1080 8GB), CHARMM36 force field
Output: M7_output/

Usage:
  python M7_md.py          # full pipeline
  python M7_md.py --test   # quick 1ns test
  python M7_md.py --prod   # 50ns production
"""

import os, sys, subprocess, shutil, re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# CONFIG
# ============================================================
SCRIPT_DIR   = Path(__file__).parent.resolve()
M6_DIR       = SCRIPT_DIR / "M6_output"
OUTPUT_DIR   = SCRIPT_DIR / "M7_output"
BEST_POSE    = M6_DIR / "1DG3_Resveratrol_best.pdb"  # from M6

# MD parameters
FORCE_FIELD  = "amber99sb-ildn"
WATER_MODEL  = "tip3p"
FF_SELECT    = "6"    # amber99sb-ildn is option 6
WATER_SELECT = "1"    # TIP3P is option 1
BOX_TYPE     = "dodecahedron"
BOX_DIST     = 1.0          # nm from solute to box edge
ION_CONC     = 0.15         # mol/L NaCl
NVT_TIME     = 100          # ps
NPT_TIME     = 100          # ps
PROD_TIME    = 10.0         # ns (override with --prod for 50ns)
TEMP         = 310          # K (body temperature)
GPU_ID       = "0"          # GPU device

GMX = "gmx"

def run_gmx(args, cwd=None, stdin_input=None, timeout=3600):
    """Run GROMACS command via cmd.exe for reliable stdin handling.

    stdin_input: str of responses separated by newlines (e.g. "4\\n4" or "1\\n14\\n").
    Multi-line input uses a temp file piped via cmd.exe type.
    """
    gmx_args = " ".join(f'"{a}"' if " " in str(a) else str(a) for a in args)
    cwd_str = str(cwd) if cwd else str(OUTPUT_DIR)

    if stdin_input:
        # Write stdin to temp file for multi-line reliability
        tmpfile = Path(cwd_str) / "_gmx_stdin.tmp"
        tmpfile.write_text(stdin_input.replace("\\n", "\n"), encoding="ascii")
        cmd = f'cmd.exe /c "type {tmpfile} | {GMX} {gmx_args}"'
    else:
        cmd = f'cmd.exe /c "{GMX} {gmx_args}"'

    print(f"  [gmx] {' '.join(str(a) for a in args)[:120]}...")
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout,
        cwd=cwd_str, shell=True
    )
    if stdin_input:
        try:
            (Path(cwd_str) / "_gmx_stdin.tmp").unlink()
        except OSError:
            pass

    if result.returncode != 0:
        err = (result.stderr + result.stdout)[-500:]
        print(f"  [ERROR] {err}")
    return result


# ============================================================
# Step 0: Check prerequisites
# ============================================================

def check_setup():
    OUTPUT_DIR.mkdir(exist_ok=True)

    if not BEST_POSE.exists():
        print(f"[ERROR] Best pose not found: {BEST_POSE}")
        print("  Run M6_docking.py first!")
        sys.exit(1)

    # Check GROMACS
    result = subprocess.run([GMX, "--version"], capture_output=True, text=True)
    if result.returncode != 0:
        print("[ERROR] GROMACS not found!")
        sys.exit(1)
    print(f"[SETUP] GROMACS OK")

    # Check pip packages
    for pkg, mod in [("rdkit", "rdkit"), ("numpy", "numpy")]:
        try:
            __import__(mod)
        except ImportError:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet", pkg])

    print("[SETUP] Ready.\n")


# ============================================================
# Step 1: Extract protein from receptor + ligand from best pose
# ============================================================

def extract_protein_ligand():
    print("[EXTRACT] Building protein+ligand complex...")

    # Protein: from M6 clean receptor PDB
    clean_pdb = M6_DIR / "1DG3_clean.pdb"
    if not clean_pdb.exists():
        print(f"[ERROR] Receptor PDB not found: {clean_pdb}")
        print("  Re-run M6_docking.py first!")
        sys.exit(1)

    # Read protein atoms (chain A only, no water/HETATM)
    protein_lines = []
    with open(clean_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                protein_lines.append(line)

    # Ligand: from the best docked pose (only UNL atoms)
    ligand_lines = []
    with open(BEST_POSE) as f:
        for line in f:
            if "UNL" in line[17:20] and (line.startswith("ATOM") or line.startswith("HETATM")):
                ligand_lines.append(line)

    if not protein_lines:
        print("[ERROR] No protein atoms found!")
        sys.exit(1)
    if not ligand_lines:
        print("[ERROR] No ligand atoms found in best pose!")
        sys.exit(1)

    # Write separate files
    prot_pdb = OUTPUT_DIR / "protein.pdb"
    lig_pdb  = OUTPUT_DIR / "ligand.pdb"
    lig_sdf  = OUTPUT_DIR / "ligand.sdf"

    prot_pdb.write_text("".join(protein_lines))
    lig_pdb.write_text("".join(ligand_lines))

    print(f"[EXTRACT] Protein: {len(protein_lines)} atoms → {prot_pdb}")
    print(f"[EXTRACT] Ligand:   {len(ligand_lines)} atoms → {lig_pdb}")

    # Write ligand SDF via RDKit (for topology generation)
    from rdkit import Chem
    from rdkit.Chem import AllChem
    smiles = "C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)
    Chem.MolToMolFile(mol, str(lig_sdf))
    print(f"[EXTRACT] Ligand SDF → {lig_sdf}")

    return prot_pdb, lig_pdb, lig_sdf


# ============================================================
# Step 2: Generate ligand topology (ACPYPE or manual)
# ============================================================

def generate_ligand_topology(lig_pdb, lig_sdf):
    """Generate GROMACS ITP+GRO using gen_ligand_itp.py (GAFF2 via gaff2.dat, no antechamber)."""
    lig_itp = OUTPUT_DIR / "ligand.itp"
    lig_gro = OUTPUT_DIR / "ligand.gro"

    # Check if already generated
    if lig_itp.exists() and lig_gro.exists():
        print(f"[LIGAND] Using existing {lig_itp} and {lig_gro}")
        return lig_itp, lig_gro

    print("[LIGAND] Generating GAFF2 topology via gen_ligand_itp.py...")

    gen_script = SCRIPT_DIR / "gen_ligand_itp.py"
    if gen_script.exists():
        result = subprocess.run(
            [sys.executable, str(gen_script), str(lig_pdb), str(lig_sdf), str(OUTPUT_DIR)],
            capture_output=True, text=True, timeout=120, cwd=str(OUTPUT_DIR)
        )
        print(result.stdout)
        if result.returncode != 0:
            print(result.stderr)
    else:
        # Fallback: expect gen_ligand_itp.py in same dir, or already generated
        print(f"[LIGAND] gen_ligand_itp.py not found at {gen_script}")
        print("[LIGAND] Run: python gen_ligand_itp.py ligand.pdb Resveratrol.sdf M7_output/")

    if lig_itp.exists() and lig_gro.exists():
        print(f"[LIGAND] Topology ready: {lig_itp} ({lig_itp.stat().st_size} bytes), "
              f"{lig_gro} ({lig_gro.stat().st_size} bytes)")
        return lig_itp, lig_gro

    print("[ERROR] Failed to generate ligand topology!")
    sys.exit(1)

# ============================================================
# Step 3: Generate protein topology
# ============================================================

def generate_protein_topology(prot_pdb):
    print("[PROTEIN] Generating topology with pdb2gmx...")

    top_dir = OUTPUT_DIR / "topol"
    top_dir.mkdir(exist_ok=True)

    # GROMACS 2020.6 reads interactively; pipe via cmd.exe
    cmd = (
        f'cmd.exe /c "(echo 9 && echo 1) | ""{GMX}"" pdb2gmx '
        f'-f ""{prot_pdb}"" -o ""{top_dir / "protein.gro"}"" '
        f'-p ""{top_dir / "topol.top"}"" '
        f'-i ""{top_dir / "posre.itp"}"" '
        f'-ff {FORCE_FIELD} -water {WATER_MODEL} -ignh"'
    )
    input_str = f"{FF_SELECT}\n{WATER_SELECT}\n"
    cmd = (
        f'cmd.exe /c "(echo {FF_SELECT} && echo {WATER_SELECT}) | ""{GMX}"" pdb2gmx '
        f'-f ""{prot_pdb}"" -o ""{top_dir / "protein.gro"}"" '
        f'-p ""{top_dir / "topol.top"}"" '
        f'-i ""{top_dir / "posre.itp"}"" '
        f'-ff {FORCE_FIELD} -water {WATER_MODEL} -ignh"'
    )
    print(f"[PROTEIN] Running: {cmd[:200]}...")
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=120, cwd=str(OUTPUT_DIR)
    )

    if result.returncode == 0 and (top_dir / "topol.top").exists():
        print("[PROTEIN] pdb2gmx OK")
    else:
        # Fallback: try without capture to see real error
        print("[PROTEIN] stdout:", result.stdout[-300:] if result.stdout else "empty")
        print("[PROTEIN] stderr:", result.stderr[-300:] if result.stderr else "empty")
        raise subprocess.CalledProcessError(result.returncode, cmd)

    print(result.stdout[-500:] if result.stdout else "")
    if result.stderr:
        print(f"[PROTEIN] STDERR: {result.stderr[:300]}")

    if result.returncode != 0:
        print("[ERROR] pdb2gmx failed!")
        sys.exit(1)

    topol_top = top_dir / "topol.top"
    print(f"[PROTEIN] Topology → {topol_top}")
    return top_dir, topol_top


# ============================================================
# Step 4: Build complex
# ============================================================

def update_topology(top_dir, lig_itp):
    """Add ligand includes to topol.top (gaff2_atypes before FF, ligand.itp before protein)."""
    print("[TOPOL] Updating topology for ligand...")

    # Copy ligand ITP + gaff2 atomtypes to topol dir
    shutil.copy(lig_itp, top_dir / "ligand.itp")
    atypes_src = OUTPUT_DIR / "gaff2_atypes.itp"
    if atypes_src.exists():
        shutil.copy(atypes_src, top_dir / "gaff2_atypes.itp")

    topol = top_dir / "topol.top"
    content = topol.read_text()

    # 1. Include gaff2 atomtypes + forcefield (no defaults) before any moleculetype.
    #    Strip [defaults] from forcefield.itp (we have our own) and use local copy.
    if 'gaff2_atypes.itp' not in content:
        import re as _re
        # Find the forcefield.itp path
        result = subprocess.run(
            [GMX, "--version"], capture_output=True, text=True, timeout=10
        )
        data_dir = None
        for line in (result.stdout + result.stderr).split("\n"):
            if 'ata' in line and ('prefix' in line.lower() or 'dir' in line.lower()):
                # "Data prefix:  D:\path"
                data_dir = line.split(":", 1)[-1].strip().strip('"')
                break

        if data_dir:
            # Compute force field path: data_dir/share/gromacs/top/amber99sb-ildn.ff/
            ff_base = Path(data_dir) / "share/gromacs/top"
            if not ff_base.exists():
                ff_base = Path(data_dir) / "top"  # alternative layout
            ff_src = ff_base / "amber99sb-ildn.ff/forcefield.itp"
            if not ff_src.exists():
                print(f"[TOPOL] WARNING: forcefield not at {ff_src}")
                # Fallback: search
                import glob as _glob
                candidates = list(Path(data_dir).rglob("amber99sb-ildn.ff/forcefield.itp"))
                ff_src = candidates[0] if candidates else Path(".") / "amber99sb-ildn.ff/forcefield.itp"
        else:
            ff_src = Path("amber99sb-ildn.ff/forcefield.itp")

        print(f"[TOPOL] FF src: {ff_src}  exists={ff_src.exists()}")
        ff_dst = top_dir / "forcefield_nod.itp"

        if ff_src.exists():
            ff_dir = ff_src.parent
            # Merge all forcefield content into one file (no nested includes)
            all_ff = []
            all_ff.append("; forcefield_nod.itp — merged AMBER99SB-ILDN (no [defaults])")
            for fname in ["ffnonbonded.itp", "ffbonded.itp", "tip3p.itp", "ions.itp"]:
                sub = ff_dir / fname
                if sub.exists():
                    all_ff.append(f"; --- {fname} ---")
                    all_ff.append(sub.read_text())
                    all_ff.append("")
            ff_dst.write_text("\n".join(all_ff))
            print(f"[TOPOL] Merged forcefield → {ff_dst} ({len(all_ff)} lines)")
        else:
            print(f"[TOPOL] WARNING: {ff_src} not found!")

        # Replace FF includes: forcefield, tip3p, ions all come from forcefield_nod
        for old in [
            '#include "amber99sb-ildn.ff/forcefield.itp"',
            '#include "amber99sb-ildn.ff/tip3p.itp"',
            '#include "amber99sb-ildn.ff/ions.itp"',
        ]:
            content = content.replace(old, "")
        # Remove blank lines from deleted includes
        content = _re.sub(r'\n\s*\n\s*\n', '\n\n', content)
        # Insert our includes where forcefield.itp was
        content = content.replace(
            '; Include forcefield parameters',
            '; Include GAFF2 + AMBER parameters (no duplicate defaults)\n'
            '#include "gaff2_atypes.itp"\n#include "forcefield_nod.itp"'
        )

    # 3. ligand.itp (moleculetype+atoms, no atomtypes) before Protein_chain_A
    if 'ligand.itp' not in content:
        content = content.replace(
            '[ moleculetype ]\n; Name            nrexcl\nProtein_chain_A',
            '#include "ligand.itp"\n\n[ moleculetype ]\n; Name            nrexcl\nProtein_chain_A'
        )

    # 4. Prepend [ defaults ] FIRST (after all manipulation, before writing)
    if '[ defaults ]' not in content:
        lines = content.split("\n")
        # Find first #include
        first_inc = next(i for i, ln in enumerate(lines) if ln.strip().startswith("#include"))
        lines.insert(first_inc, (
            '\n[ defaults ]\n'
            '; nbfunc   comb-rule   gen-pairs   fudgeLJ  fudgeQQ\n'
            '1          2           yes          0.5      0.8333\n'
        ))
        content = "\n".join(lines)

    topol.write_text(content)
    print(f"[TOPOL] Updated → {topol}")
    return topol


# ============================================================
# Step 5: Solvate + ions
# ============================================================

def solvate_and_ions(topol, complex_gro, top_dir):
    print("[SOLVATE] Setting up box and solvation...")

    # Instead of editconf (may drop ligand atoms), manually add box to GRO,
    # then use solvate to both create box and add water in one step.
    # Read complex GRO and estimate box size
    gro_lines = complex_gro.read_text().split("\n")
    n_atoms = int(gro_lines[1].strip())
    coords = []
    for line in gro_lines[2:2+n_atoms]:
        if len(line) >= 44:
            try:
                coords.append((
                    float(line[20:28]), float(line[28:36]), float(line[36:44])
                ))
            except ValueError:
                continue

    if not coords:
        print("[ERROR] No coordinates found in GRO!")
        sys.exit(1)

    xs = [c[0] for c in coords]; ys = [c[1] for c in coords]; zs = [c[2] for c in coords]
    margin = BOX_DIST  # nm
    bx = max(xs) - min(xs) + 2*margin
    by = max(ys) - min(ys) + 2*margin
    bz = max(zs) - min(zs) + 2*margin
    print(f"[SOLVATE] Box (cubic): {bx:.1f} x {by:.1f} x {bz:.1f} nm")

    # Update box vectors in GRO (last line)
    gro_lines[-1] = f" {bx:10.5f} {by:10.5f} {bz:10.5f}"
    boxed_gro = top_dir / "boxed.gro"
    boxed_gro.write_text("\n".join(gro_lines))

    # Solvate (also adds box if not set)
    solv_gro = top_dir / "solv.gro"
    solv_top = top_dir / "topol.top"
    result = run_gmx(["solvate",
         "-cp", str(boxed_gro),
         "-cs", "spc216.gro",
         "-o", str(solv_gro),
         "-p", str(solv_top)])
    if result.returncode != 0: sys.exit(1)

    # Add ions
    ions_mdp = top_dir / "ions.mdp"
    ions_mdp.write_text("integrator = steep\nemtol = 1000.0\nnsteps = 50000\n")
    ions_tpr = top_dir / "ions.tpr"
    result = run_gmx(["grompp",
         "-f", str(ions_mdp),
         "-c", str(solv_gro),
         "-p", str(solv_top),
         "-o", str(ions_tpr),
         "-maxwarn", "5"])
    if result.returncode != 0: sys.exit(1)

    final_gro = top_dir / "system.gro"
    result = run_gmx(["genion",
         "-s", str(ions_tpr),
         "-o", str(final_gro),
         "-p", str(solv_top),
         "-pname", "NA",
         "-nname", "CL",
         "-conc", str(ION_CONC),
         "-neutral"],
         stdin_input="SOL")
    print(f"[SOLVATE] System built → {final_gro}")

    # Count water
    if result.stdout:
        for line in result.stdout.split("\n"):
            if "SOL" in line and "molecules" in line:
                print(f"  {line.strip()}")

    return solv_top, final_gro


# ============================================================
# Step 5b: Append ligand to solvated system (avoids solvate dropping ligand)
# ============================================================

def append_ligand(system_gro, lig_gro, topol):
    """Append ligand atoms to the ionized system GRO and fix topology."""
    print("[LIGAND] Appending ligand to solvated system...")

    # Add UNL to topology [molecules] (not added pre-solvate to avoid mismatch)
    top_content = topol.read_text()
    if 'UNL                  1' not in top_content:
        top_content = top_content.replace(
            'Protein_chain_A     1',
            'Protein_chain_A     1\nUNL                  1'
        )
        topol.write_text(top_content)
        print(f"[LIGAND] Topology: +UNL 1 in [molecules]")

    # Filter empty lines (split trailing newline creates "")
    sys_lines = [l for l in system_gro.read_text().split("\n") if l.strip() != ""]
    lig_lines = [l for l in lig_gro.read_text().split("\n") if l.strip() != ""]

    n_sys = int(sys_lines[1].strip())
    n_lig = int(lig_lines[1].strip())
    n_tot = n_sys + n_lig

    # Ligand atom lines (skip header + count)
    lig_atom_lines = lig_lines[2:2+n_lig]

    # Insert ligand atoms before box vector (last non-empty line)
    box = sys_lines[-1]
    sys_lines = sys_lines[:-1] + lig_atom_lines + [box]

    # Update atom count and title
    sys_lines[0] = f"GBP1+Resveratrol + solvent ({n_tot} atoms)"
    sys_lines[1] = f" {n_tot}"

    # 3. Remove whole water molecules clashing with ligand (distance < 0.23 nm)
    lig_coords = []
    for ln in lig_atom_lines:
        try:
            lig_coords.append((
                float(ln[20:28]), float(ln[28:36]), float(ln[36:44])
            ))
        except (ValueError, IndexError):
            continue

    # Build index: residue_number -> list of line indices for that water
    water_residues = {}  # res_id -> [indices of all 3 atoms in sys_lines[2:]]
    for i, ln in enumerate(sys_lines[2:]):
        if "SOL" in ln[5:10]:
            res_id = ln[0:5]  # residue number (fixed-width, columns 0-4)
            water_residues.setdefault(res_id, []).append(i)

    # Find clashing OW atoms, then remove entire water molecule
    clash_residues = set()
    for i, ln in enumerate(sys_lines[2:]):
        if "SOL" in ln[5:10] and "OW" in ln[10:15]:
            try:
                wx = float(ln[20:28]); wy = float(ln[28:36]); wz = float(ln[36:44])
                for lx, ly, lz in lig_coords:
                    if ((wx - lx)**2 + (wy - ly)**2 + (wz - lz)**2)**0.5 < 0.30:
                        clash_residues.add(ln[0:5])
                        break
            except (ValueError, IndexError):
                pass

    # Collect indices of all atoms belonging to clashing water molecules
    remove_indices = set()
    for res_id in clash_residues:
        if res_id in water_residues:
            remove_indices.update(water_residues[res_id])

    # Build keep list, skipping all atoms of clashing waters
    keep_atoms = [sys_lines[0], sys_lines[1]]
    for i, ln in enumerate(sys_lines[2:]):
        if i not in remove_indices:
            keep_atoms.append(ln)

    n_removed_atoms = len(remove_indices)
    n_removed_waters = n_removed_atoms // 3

    if n_removed_waters > 0:
        n_sys -= n_removed_atoms
        n_tot = n_sys + n_lig
        keep_atoms[1] = f"{n_tot:>6d}"
        top_content = topol.read_text()
        import re as _re2
        top_content = _re2.sub(
            r'(SOL\s+)(\d+)',
            lambda m: f"{m.group(1)}{int(m.group(2)) - n_removed_waters}",
            top_content
        )
        topol.write_text(top_content)
    else:
        print("[LIGAND] WARNING: No clashing water found — is ligand outside the box?")

    print(f"[LIGAND] Removed {n_removed_waters} water molecules "
          f"({n_removed_atoms} atoms) clashing with ligand")

    sys_lines = keep_atoms

    merged_gro = system_gro.with_name("merged.gro")
    merged_gro.write_text("\n".join(sys_lines))
    print(f"[LIGAND] {n_sys} sys + {n_lig} lig = {n_tot} atoms → {merged_gro}")

    return merged_gro


# ============================================================
# Step 6: Energy minimization
# ============================================================

def energy_minimization(topol, system_gro, label="em", restrain_protein=True, emtol=1000.0, emstep=0.001):
    print(f"[EM:{label}] Energy minimization (emtol={emtol}, emstep={emstep})...")

    define = "define = -DPOSRES" if restrain_protein else ""
    em_mdp = OUTPUT_DIR / f"{label}.mdp"
    em_mdp.write_text(f"""
integrator  = steep
emtol       = {emtol}
emstep      = {emstep}
nsteps      = 100000
nstlist     = 1
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
vdwtype     = Cut-off
rvdw        = 1.2
pbc         = xyz
{define}
""")

    gmx_args = ["grompp",
        "-f", str(em_mdp), "-c", str(system_gro)]
    if restrain_protein:
        gmx_args += ["-r", str(system_gro)]
    gmx_args += ["-p", str(topol),
        "-o", str(OUTPUT_DIR / f"{label}.tpr"), "-maxwarn", "10"]

    result = run_gmx(gmx_args)
    if result.returncode != 0:
        print(f"[EM:{label}] grompp FAILED!")
        print(result.stdout[-1000:])
        print(result.stderr[-1000:])
        sys.exit(1)

    result = run_gmx(["mdrun", "-deffnm", label, "-v", "-ntmpi", "1", "-gpu_id", GPU_ID])

    # Extract final energy from all output
    for line in (result.stdout + result.stderr).split("\n"):
        if "Potential Energy" in line or "Maximum force" in line:
            print(f"[EM] {line.strip()}")

    if not (OUTPUT_DIR / f"{label}.gro").exists():
        print(f"[ERROR] EM:{label} failed — no {label}.gro!")
        print((result.stderr + result.stdout)[-2000:])
        sys.exit(1)

    em_gro = OUTPUT_DIR / f"{label}.gro"
    print(f"[EM:{label}] Done → {em_gro}")
    return em_gro

    # If EM energy is too high, warn but continue
    for line in (result.stdout + result.stderr).split("\n"):
        if "Potential Energy" in line:
            try:
                pe = float(line.split("=")[-1].strip())
                if abs(pe) > 1e10:
                    print(f"[EM] WARNING: Energy too high ({pe:.2e}), NVT may crash")
            except:
                pass


# ============================================================
# Step 7: NVT equilibration
# ============================================================

def run_nvt(topol):
    print(f"[NVT] Equilibration ({NVT_TIME} ps, {TEMP} K)...")

    nvt_mdp = OUTPUT_DIR / "nvt.mdp"
    nvt_mdp.write_text(f"""
integrator  = md
dt          = 0.002
nsteps      = {NVT_TIME * 500}
nstxout-compressed = 1000
nstlog      = 1000
nstenergy   = 1000
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
ref-t       = {TEMP} {TEMP}
pcoupl      = no
gen-vel     = yes
gen-temp    = {TEMP}
constraints = h-bonds
constraint-algorithm = LINCS
""")

    result = run_gmx(["grompp",
         "-f", str(nvt_mdp), "-c", str(OUTPUT_DIR / "em.gro"),
         "-r", str(OUTPUT_DIR / "em.gro"),
         "-p", str(topol), "-o", str(OUTPUT_DIR / "nvt.tpr"), "-maxwarn", "5"])
    if result.returncode != 0: sys.exit(1)

    result = run_gmx(["mdrun", "-deffnm", "nvt", "-v", "-ntmpi", "1", "-gpu_id", GPU_ID])
    nvt_gro = OUTPUT_DIR / "nvt.gro"
    if not nvt_gro.exists():
        print("[NVT] FAILED — no nvt.gro produced!")
        print((result.stderr + result.stdout)[-2000:])
        sys.exit(1)
    print("[NVT] Done.")


# ============================================================
# Step 8: NPT equilibration
# ============================================================

def run_npt(topol):
    print(f"[NPT] Equilibration ({NPT_TIME} ps, {TEMP} K, 1 bar)...")

    npt_mdp = OUTPUT_DIR / "npt.mdp"
    npt_mdp.write_text(f"""
integrator  = md
dt          = 0.002
nsteps      = {NPT_TIME * 500}
nstxout-compressed = 1000
nstlog      = 1000
nstenergy   = 1000
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
ref-t       = {TEMP} {TEMP}
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint-algorithm = LINCS
""")

    result = run_gmx(["grompp",
         "-f", str(npt_mdp), "-c", str(OUTPUT_DIR / "nvt.gro"),
         "-r", str(OUTPUT_DIR / "nvt.gro"),
         "-p", str(topol), "-o", str(OUTPUT_DIR / "npt.tpr"), "-maxwarn", "5"])
    if result.returncode != 0: sys.exit(1)

    result = run_gmx(["mdrun", "-deffnm", "npt", "-v", "-ntmpi", "1", "-gpu_id", GPU_ID])
    print("[NPT] Done.")


# ============================================================
# Step 9: Production MD
# ============================================================

def run_production(topol, prod_ns):
    print(f"[MD] Production ({prod_ns} ns, {TEMP} K)...")
    nsteps = int(prod_ns * 1000 / 0.002)

    prod_mdp = OUTPUT_DIR / "prod.mdp"
    prod_mdp.write_text(f"""
integrator  = md
dt          = 0.002
nsteps      = {nsteps}
nstxout-compressed = 5000
nstlog      = 5000
nstenergy   = 5000
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
ref-t       = {TEMP} {TEMP}
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint-algorithm = LINCS
""")

    prod_cpt = OUTPUT_DIR / "prod.cpt"
    prod_tpr = OUTPUT_DIR / "prod.tpr"

    if not prod_cpt.exists():
        # Fresh start: run grompp
        init_gro = OUTPUT_DIR / "npt.gro"
        init_cpt = OUTPUT_DIR / "npt.cpt"
        result = run_gmx(["grompp",
             "-f", str(prod_mdp), "-c", str(init_gro), "-t", str(init_cpt),
             "-p", str(topol), "-o", str(prod_tpr), "-maxwarn", "5"])
        if result.returncode != 0: sys.exit(1)
    else:
        print("[MD] Checkpoint found, skipping grompp (TPR must match checkpoint).")

    mdrun_args = ["mdrun", "-deffnm", "prod", "-v", "-ntmpi", "1", "-gpu_id", GPU_ID]
    if prod_cpt.exists():
        mdrun_args.extend(["-cpi", str(prod_cpt)])
        print("[MD] Resuming from checkpoint...")

    print(f"[MD] Starting {prod_ns}ns GPU simulation (no timeout, ~29h on GTX 1080)...")
    result = run_gmx(mdrun_args, timeout=None)

    # Print performance summary
    for line in (result.stdout + result.stderr).split("\n"):
        if "Performance" in line or "ns/day" in line or "Finished" in line:
            print(f"[MD] {line.strip()}")

    print("[MD] Production complete.")


# ============================================================
# Step 10: Analysis
# ============================================================

def parse_xvg(path):
    """Parse GROMACS xvg file, return (x_data, y_data, xlabel, ylabel, title)."""
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


def plot_figure(x, y, xlabel, ylabel, title, outpath, color="#2c3e50"):
    """Single-panel publication-quality PDF."""
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, y, color=color, linewidth=0.8, alpha=0.9)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle="--")

    # Annotate stats
    y_mean, y_std = np.mean(y), np.std(y)
    stats = f"Mean = {y_mean:.2f} | SD = {y_std:.2f}"
    ax.text(0.98, 0.95, stats, transform=ax.transAxes, fontsize=9,
            ha="right", va="top", bbox=dict(boxstyle="round,pad=0.3",
            facecolor="white", alpha=0.8, edgecolor="gray"))

    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [PLOT] {outpath.name} ({outpath.stat().st_size/1024:.0f} KB)")
    return y_mean, y_std


def run_hbond_analysis(ligand_name="LIG"):
    """Protein-ligand hydrogen bond analysis."""
    print("\n[HBOND] Protein-ligand hydrogen bonds...")
    hb_out = OUTPUT_DIR / "hbond.xvg"
    hb_num = OUTPUT_DIR / "hbnum.xvg"

    # hbond selection: protein (group 1) + ligand (group 14 usually, or named)
    # We need to find the ligand group number
    # Try with the ligand residue name
    run_gmx(["hbond",
         "-s", str(OUTPUT_DIR / "prod.tpr"),
         "-f", str(OUTPUT_DIR / "prod.xtc"),
         "-num", str(hb_num),
         "-hbn", str(OUTPUT_DIR / "hbond.ndx")],
         stdin_input=f"1\n14\n", timeout=300)

    if hb_num.exists():
        x, y, _, _, _ = parse_xvg(hb_num)
        plot_figure(x, y, "Time (ns)", "Number of H-bonds",
                    "GBP1–Resveratrol Hydrogen Bonds",
                    OUTPUT_DIR / "hbonds.pdf", color="#27ae60")
        print(f"  [HBOND] Avg H-bonds: {np.mean(y):.2f}")
    else:
        print("  [HBOND] No H-bond output, trying alternate group selection...")
        # Try: protein(1) + non-protein(14), or just all
        run_gmx(["hbond",
             "-s", str(OUTPUT_DIR / "prod.tpr"),
             "-f", str(OUTPUT_DIR / "prod.xtc"),
             "-num", str(hb_num)],
             stdin_input="Protein\nSystem\n", timeout=300)
        if hb_num.exists():
            x, y, _, _, _ = parse_xvg(hb_num)
            plot_figure(x, y, "Time (ns)", "Number of H-bonds",
                        "GBP1–Resveratrol Hydrogen Bonds",
                        OUTPUT_DIR / "hbonds.pdf", color="#27ae60")


def analyze():
    print("\n[ANALYSIS] Generating MD quality metrics...")
    summary_stats = {}

    # ---- RMSD ----
    rmsd_out = OUTPUT_DIR / "rmsd.xvg"
    if not rmsd_out.exists():
        run_gmx(["rms", "-s", str(OUTPUT_DIR / "prod.tpr"),
                 "-f", str(OUTPUT_DIR / "prod.xtc"),
                 "-o", str(rmsd_out), "-tu", "ns"],
                 stdin_input="4 4")
    x, y, xl, yl, _ = parse_xvg(rmsd_out)
    mean_rmsd, std_rmsd = plot_figure(x, y, xl, "RMSD (nm)",
                                       "GBP1 Backbone RMSD (Cα)",
                                       OUTPUT_DIR / "rmsd.pdf", color="#e74c3c")
    summary_stats["RMSD"] = (mean_rmsd, std_rmsd)

    # ---- RMSF ----
    rmsf_out = OUTPUT_DIR / "rmsf.xvg"
    if not rmsf_out.exists():
        run_gmx(["rmsf", "-s", str(OUTPUT_DIR / "prod.tpr"),
                 "-f", str(OUTPUT_DIR / "prod.xtc"),
                 "-o", str(rmsf_out), "-res"],
                 stdin_input="4")
    x, y, xl, yl, _ = parse_xvg(rmsf_out)
    plot_figure(x, y, "Residue Index", "RMSF (nm)",
                "GBP1 Residue Fluctuation (RMSF)",
                OUTPUT_DIR / "rmsf.pdf", color="#8e44ad")
    summary_stats["RMSF"] = (np.max(y), np.argmax(y) if len(y) > 0 else 0)

    # ---- Rg ----
    rg_out = OUTPUT_DIR / "rg.xvg"
    if not rg_out.exists():
        run_gmx(["gyrate", "-s", str(OUTPUT_DIR / "prod.tpr"),
                 "-f", str(OUTPUT_DIR / "prod.xtc"),
                 "-o", str(rg_out)],
                 stdin_input="4")
    x, y, xl, yl, _ = parse_xvg(rg_out)
    mean_rg, std_rg = plot_figure(x, y, "Time (ns)", yl,
                                   "GBP1 Radius of Gyration",
                                   OUTPUT_DIR / "rg.pdf", color="#2980b9")
    summary_stats["Rg"] = (mean_rg, std_rg)

    # ---- SASA ----
    sasa_out = OUTPUT_DIR / "sasa.xvg"
    if not sasa_out.exists():
        run_gmx(["sasa", "-s", str(OUTPUT_DIR / "prod.tpr"),
                 "-f", str(OUTPUT_DIR / "prod.xtc"),
                 "-o", str(sasa_out)],
                 stdin_input="4")
    x, y, xl, yl, _ = parse_xvg(sasa_out)
    mean_sasa, std_sasa = plot_figure(x, y, "Time (ns)", yl,
                                       "GBP1 Solvent Accessible Surface Area",
                                       OUTPUT_DIR / "sasa.pdf", color="#d35400")
    summary_stats["SASA"] = (mean_sasa, std_sasa)

    # ---- Hydrogen bonds ----
    run_hbond_analysis()

    # ---- Combined summary figure ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.flatten()
    for ax, (name, out, color) in zip(axes, [
        ("RMSD", rmsd_out, "#e74c3c"),
        ("RMSF", rmsf_out, "#8e44ad"),
        ("Rg", rg_out, "#2980b9"),
        ("SASA", sasa_out, "#d35400")]):
        x, y, xl, yl, _ = parse_xvg(out)
        ax.plot(x, y, color=color, linewidth=0.6)
        ax.set_title(name, fontweight="bold")
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.3, linestyle="--")
    fig.suptitle(f"GBP1–Resveratrol · 50ns MD Quality Metrics", fontsize=14, fontweight="bold")
    fig.tight_layout()
    combo = OUTPUT_DIR / "M7_QC_panel.pdf"
    fig.savefig(combo, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [PLOT] {combo.name} ({combo.stat().st_size/1024:.0f} KB)")

    # ---- Summary file ----
    summary = OUTPUT_DIR / "M7_summary.txt"
    lines = [
        "=" * 50,
        "M7: GBP1 + Resveratrol · 50ns MD Analysis Summary",
        "=" * 50,
        f"RMSD (Cα):  mean = {summary_stats['RMSD'][0]:.3f} nm, SD = {summary_stats['RMSD'][1]:.3f} nm",
        f"RMSF (max):  {summary_stats['RMSF'][0]:.3f} nm at residue {summary_stats['RMSF'][1]}",
        f"Rg:          mean = {summary_stats['Rg'][0]:.3f} nm, SD = {summary_stats['Rg'][1]:.3f} nm",
        f"SASA:        mean = {summary_stats['SASA'][0]:.1f} nm², SD = {summary_stats['SASA'][1]:.1f} nm²",
        "=" * 50,
        "Output PDFs: rmsd.pdf, rmsf.pdf, rg.pdf, sasa.pdf, hbonds.pdf, M7_QC_panel.pdf",
    ]
    summary.write_text("\n".join(lines), encoding="utf-8")
    print("\n[M7 SUMMARY]")
    for line in lines:
        print(f"  {line}")


# ============================================================
# Main
# ============================================================

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_true", help="Quick 1ns test")
    parser.add_argument("--prod", action="store_true", help="50ns production")
    parser.add_argument("--continue", action="store_true", dest="continue_md",
                        help="Skip to production from checkpoint (after timeout/crash)")
    parser.add_argument("--analyze", action="store_true",
                        help="Run analysis + plotting only (post-MD)")
    args = parser.parse_args()

    global PROD_TIME

    # --analyze: just run analysis on existing trajectory
    if args.analyze:
        print("=" * 60)
        print("  M7: Analysis + Plotting Mode")
        print("=" * 60)
        analyze()
        print("\n[DONE] M7 analysis complete!")
        return

    if args.test:
        PROD_TIME = 1.0
    elif args.prod or args.continue_md:
        PROD_TIME = 50.0

    print("=" * 60)
    print(f"  M7: GBP1 + Resveratrol · GROMACS MD ({PROD_TIME}ns)")
    print(f"  GPU: GTX 1080 · FF: {FORCE_FIELD} · T: {TEMP}K")
    print("=" * 60)

    if args.continue_md:
        topol = OUTPUT_DIR / "topol" / "topol.top"
        if not topol.exists():
            print(f"[ERROR] Topology not found: {topol}")
            sys.exit(1)
        print("[MODE] Continue from checkpoint — skipping EM/NVT/NPT...")
        run_production(topol, PROD_TIME)
        analyze()
        print(f"\n[DONE] M7 completed ({PROD_TIME}ns)!")
        return

    check_setup()
    prot_pdb, lig_pdb, lig_sdf = extract_protein_ligand()

    # Prefer M6 Resveratrol.sdf (has better structure) over generated SDF
    sdf_from_m6 = M6_DIR / "Resveratrol.sdf"
    sdf_for_gen = sdf_from_m6 if sdf_from_m6.exists() else lig_sdf
    print(f"[SETUP] Using SDF: {sdf_for_gen}")

    lig_itp, lig_gro = generate_ligand_topology(lig_pdb, sdf_for_gen)
    top_dir, topol = generate_protein_topology(prot_pdb)

    # Append ligand to protein BEFORE solvation (standard complex setup)
    update_topology(top_dir, lig_itp)
    prot_gro = top_dir / "protein.gro"
    complex_gro = append_ligand(prot_gro, lig_gro, topol)
    print(f"[COMPLEX] Protein + ligand → {complex_gro}")

    # Solvate + ions: full complex (protein + ligand)
    topol, system_gro = solvate_and_ions(topol, complex_gro, top_dir)

    # EM: minimize full system
    print("\n" + "=" * 50)
    print("[EM1] Minimization with protein restraints...")
    system_gro = energy_minimization(topol, system_gro, label="pre_em",
                                     restrain_protein=True, emtol=500.0)
    print("=" * 50 + "\n")

    print("\n" + "=" * 50)
    print("[EM2] Final minimization (no restraints)...")
    system_gro = energy_minimization(topol, system_gro, label="em",
                                     restrain_protein=False, emtol=100.0)
    print("=" * 50 + "\n")
    run_nvt(topol)
    run_npt(topol)
    run_production(topol, PROD_TIME)
    analyze()

    print(f"\n[DONE] M7 completed ({PROD_TIME}ns)!")


if __name__ == "__main__":
    main()
