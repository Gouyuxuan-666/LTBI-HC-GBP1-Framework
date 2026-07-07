#!/usr/bin/env python3
"""
M6: Molecular Docking — GBP1 (PDB 1DG3) + Resveratrol
=======================================================
Uses meeko (PDBQT prep) + standalone AutoDock Vina .exe

Files needed alongside this script:
  1DG3.pdb   — GBP1 crystal structure (416 KB)
  vina.exe   — AutoDock Vina Windows binary (you download this)

Setup (one time):
  pip install biopython rdkit meeko numpy

Usage:
  python M6_docking.py

Output: M6_output/docking_results.csv + best_pose.pdb + M6_report.txt
"""

import os, sys, subprocess, shutil, csv, re
from pathlib import Path
import numpy as np

# ============================================================
# CONFIG
# ============================================================
SCRIPT_DIR   = Path(__file__).parent.resolve()
PDB_ID       = "1DG3"
PDB_FILE     = SCRIPT_DIR / f"{PDB_ID}.pdb"
VINA_EXE     = SCRIPT_DIR / "vina.exe"
LIGAND_NAME  = "Resveratrol"
LIGAND_SMILES = "C1=CC(=CC=C1/C=C/C2=CC(=CC(=C2)O)O)O"
OUTPUT_DIR   = SCRIPT_DIR / "M6_output"

EXHAUSTIVENESS = 32
NUM_MODES      = 20


# ============================================================
# Step 0: Check prerequisites
# ============================================================

def check_setup():
    OUTPUT_DIR.mkdir(exist_ok=True)

    # Check files
    for f, desc in [(PDB_FILE, "1DG3.pdb"), (VINA_EXE, "vina.exe")]:
        if not f.exists():
            print(f"[ERROR] {desc} not found at {f}")
            print(f"  Put {desc} in the same folder as this script.")
            sys.exit(1)
    print(f"[SETUP] 1DG3.pdb OK ({PDB_FILE.stat().st_size/1024:.0f} KB)")
    print(f"[SETUP] vina.exe OK ({VINA_EXE.stat().st_size/1024:.0f} KB)")

    # Check Python packages
    deps = {"biopython": "Bio", "rdkit": "rdkit", "meeko": "meeko", "numpy": "numpy"}
    missing = []
    for pkg, mod in deps.items():
        try:
            __import__(mod)
        except ImportError:
            missing.append(pkg)
    if missing:
        print(f"[SETUP] Installing: {' '.join(missing)}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "--quiet",
             "-i", "https://pypi.tuna.tsinghua.edu.cn/simple"] + missing
        )
    print("[SETUP] All Python packages OK.\n")


# ============================================================
# Step 1: Prepare receptor (GBP1)
# ============================================================

def prepare_receptor():
    print(f"[RECEPTOR] Reading {PDB_FILE}...")

    from Bio.PDB import PDBParser, PDBIO, Select

    class ChainOnly(Select):
        def accept_chain(self, chain):
            return chain.get_id() == "A"
        def accept_residue(self, residue):
            return residue.get_id()[0] == " "  # skip HETATM/water

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("gbp1", str(PDB_FILE))

    clean_pdb = OUTPUT_DIR / f"{PDB_ID}_clean.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(clean_pdb), ChainOnly())
    print(f"[RECEPTOR] Clean PDB → {clean_pdb}")

    # meeko: PDB → PDBQT
    print("[RECEPTOR] Running meeko (PDB → PDBQT)...")
    pdbqt_path = OUTPUT_DIR / f"{PDB_ID}_receptor.pdbqt"

    result = subprocess.run(
        [sys.executable, "-m", "meeko", "--receptor", str(clean_pdb),
         "-o", str(pdbqt_path)],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        print("[RECEPTOR] meeko failed, using manual PDBQT...")
        write_manual_pdbqt(clean_pdb, pdbqt_path)
    else:
        candidates = list(OUTPUT_DIR.glob("**/*receptor*.pdbqt"))
        if candidates:
            shutil.copy(candidates[0], pdbqt_path)
            print(f"[RECEPTOR] meeko output → {pdbqt_path}")
        else:
            print("[RECEPTOR] WARNING: no PDBQT found, using fallback")
            write_manual_pdbqt(clean_pdb, pdbqt_path)

    return pdbqt_path, clean_pdb


def write_manual_pdbqt(pdb_path, output_path):
    """Build proper PDBQT from PDB atom data."""
    ad_types = {"C": "C", "N": "N", "O": "OA", "S": "SA", "P": "P"}
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                record  = line[0:6]
                serial  = int(line[6:11])
                name    = line[12:16]
                altLoc  = line[16:17]
                resName = line[17:20]
                chain   = line[21:22]
                resNum  = int(line[22:26])
                iCode   = line[26:27]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occupancy = line[54:60].strip() or "1.00"
                tempFactor = line[60:66].strip() or "0.00"

                # Determine element and AutoDock type
                elem = line[76:78].strip()
                if not elem:
                    elem = name[0] if name[0] != " " else name[1:2]
                    if not elem or elem == "H":
                        elem = "C"
                ad = ad_types.get(elem, "C")

                # Build PDBQT line with correct column positions
                # Columns 71-76: partial charge, Columns 77-78: atom type
                pdbqt_line = (
                    f"{record:<6}"                    # 1-6
                    f"{serial:>5}"                    # 7-11
                    f" "                               # 12
                    f"{name:<4}"                       # 13-16
                    f"{altLoc}"                        # 17
                    f"{resName:<3}"                    # 18-20
                    f" "                               # 21
                    f"{chain}"                         # 22
                    f"{resNum:>4}"                     # 23-26
                    f"{iCode}"                         # 27
                    f"   "                             # 28-30
                    f"{x:8.3f}"                        # 31-38
                    f"{y:8.3f}"                        # 39-46
                    f"{z:8.3f}"                        # 47-54
                    f"{float(occupancy):6.2f}"         # 55-60
                    f"{float(tempFactor):6.2f}"        # 61-66
                    f"    "                            # 67-70
                    f" 0.000"                          # 71-76 partial charge
                    f" {ad:<3}"                         # 77-80 atom type
                )
                lines.append(pdbqt_line)
            except Exception as e:
                continue

    output_path.write_text("\n".join(lines))
    print(f"  [manual] {len(lines)} atoms → {output_path}")


# ============================================================
# Step 2: Prepare ligand (Resveratrol)
# ============================================================

def prepare_ligand():
    print(f"[LIGAND] Preparing {LIGAND_NAME} from SMILES...")

    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(LIGAND_SMILES)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)

    sdf_path = OUTPUT_DIR / f"{LIGAND_NAME}.sdf"
    Chem.MolToMolFile(mol, str(sdf_path))

    pdbqt_path = OUTPUT_DIR / f"{LIGAND_NAME}.pdbqt"

    # Try meeko Python API first (gives proper ROOT/BRANCH structure)
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        preparator = MoleculePreparation()
        mol_setup = preparator.prepare(mol)[0]  # returns list
        pdbqt_string, is_ok = PDBQTWriterLegacy.write_string(mol_setup)

        if is_ok:
            pdbqt_path.write_text(pdbqt_string)
            print(f"[LIGAND] meeko API → {pdbqt_path} ({len(pdbqt_string.splitlines())} lines)")
            return pdbqt_path
        else:
            print("[LIGAND] meeko API returned is_ok=False")
    except Exception as e:
        print(f"[LIGAND] meeko API failed: {e}")

    # Fallback: CLI
    print("[LIGAND] Trying meeko CLI...")
    result = subprocess.run(
        [sys.executable, "-m", "meeko", "-i", str(sdf_path), "-o", str(pdbqt_path)],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        candidates = list(OUTPUT_DIR.glob("**/*ligand*.pdbqt"))
        if not candidates:
            candidates = list(OUTPUT_DIR.glob(f"**/*{LIGAND_NAME}*.pdbqt"))
        if candidates:
            shutil.copy(candidates[0], pdbqt_path)
            print(f"[LIGAND] meeko CLI → {pdbqt_path}")
            return pdbqt_path

    # Last resort: manual with torsion tree
    print("[LIGAND] All meeko methods failed, using manual PDBQT with torsion info...")
    write_manual_pdbqt_with_torsions(mol, pdbqt_path)
    return pdbqt_path


def write_manual_pdbqt_with_torsions(mol, output_path):
    """RDKit mol → PDBQT with ROOT/BRANCH/ENDROOT/TORSDOF structure."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    conf = mol.GetConformer()
    ad_types = {"C": "C", "N": "N", "O": "OA", "S": "SA", "P": "P"}
    lines = []

    # Write ROOT
    lines.append("ROOT")

    # Identify rotatable bonds
    rot_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE and not bond.IsInRing():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Skip: terminal groups (OH, CH3), bonds to non-heavy atoms
            if a1.GetAtomicNum() > 1 and a2.GetAtomicNum() > 1:
                n1 = len([x for x in a1.GetNeighbors() if x.GetAtomicNum() > 1])
                n2 = len([x for x in a2.GetNeighbors() if x.GetAtomicNum() > 1])
                if n1 > 1 and n2 > 1:
                    rot_bonds.append((bond.GetIdx(), a1.GetIdx(), a2.GetIdx()))

    # Write all atoms in ROOT (rigid ligand, no branches)
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        elem = atom.GetSymbol()
        ad = ad_types.get(elem, "C")
        # Cols: 1-6 HETATM, 7-11 serial, 12 blank, 13-16 name, 17 altLoc,
        # 18-20 res, 21 blank, 22 chain, 23-26 resNum, 27 iCode, 28-30 blank,
        # 31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 tempF,
        # 67-70 blank, 71-76 charge, 77-80 atom type
        lines.append(
            f"HETATM"                           # 1-6
            f"{i+1:>5}"                         # 7-11
            f" "                                 # 12
            f"{elem:<4}"                         # 13-16
            f" "                                 # 17 altLoc
            f"UNL"                               # 18-20 res
            f" "                                 # 21
            f" "                                 # 22 chain
            f"{1:>4}"                            # 23-26 resNum
            f" "                                 # 27 iCode
            f"   "                               # 28-30
            f"{pos.x:8.3f}"                      # 31-38
            f"{pos.y:8.3f}"                      # 39-46
            f"{pos.z:8.3f}"                      # 47-54
            f"  1.00"                            # 55-60
            f"  0.00"                            # 61-66
            f"    "                              # 67-70
            f" 0.000"                            # 71-76
            f" {ad:<3}"                           # 77-80
        )
    lines.append("ENDROOT")
    lines.append("TORSDOF 0")
    output_path.write_text("\n".join(lines))
    print(f"  [manual+torsions] {mol.GetNumAtoms()} atoms, {len(rot_bonds)} rotatable bonds → {output_path}")


# ============================================================
# Step 3: Define binding pocket
# ============================================================

def define_pocket(clean_pdb_path):
    print("[POCKET] Computing binding box...")
    coords = []
    with open(clean_pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    continue

    arr = np.array(coords)
    center = arr.mean(axis=0)
    size = (arr.max(axis=0) - arr.min(axis=0)) + 8.0

    print(f"[POCKET] Center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    print(f"[POCKET] Size:   ({size[0]:.1f}, {size[1]:.1f}, {size[2]:.1f}) Å")
    return center, size


# ============================================================
# Step 4: Run docking (standalone vina.exe)
# ============================================================

def run_docking(receptor_pdbqt, ligand_pdbqt, center, size):
    print(f"\n[DOCKING] Running {VINA_EXE} ...")

    out_pdbqt = OUTPUT_DIR / f"{PDB_ID}_{LIGAND_NAME}_docked.pdbqt"

    cmd = [
        str(VINA_EXE),
        "--receptor",      str(receptor_pdbqt),
        "--ligand",        str(ligand_pdbqt),
        "--out",           str(out_pdbqt),
        "--center_x",      f"{center[0]:.4f}",
        "--center_y",      f"{center[1]:.4f}",
        "--center_z",      f"{center[2]:.4f}",
        "--size_x",        f"{size[0]:.4f}",
        "--size_y",        f"{size[1]:.4f}",
        "--size_z",        f"{size[2]:.4f}",
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes",      str(NUM_MODES),
    ]

    print(f"[DOCKING] {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    combined = result.stdout + "\n" + result.stderr

    # Save log
    (OUTPUT_DIR / "vina_log.txt").write_text(combined)
    print(f"[DOCKING] Log saved to vina_log.txt")

    # Show condensed output
    for line in result.stdout.strip().split("\n"):
        if "kcal" in line or "mode" in line.lower() or "-----" in line:
            print(f"  {line.strip()}")

    if result.returncode != 0:
        print(f"[DOCKING] ERROR (code {result.returncode}):")
        print(result.stderr[:800])
        return None, None

    if not out_pdbqt.exists():
        print(f"[DOCKING] ERROR: output file not created")
        return None, None

    print(f"[DOCKING] Output → {out_pdbqt} ({out_pdbqt.stat().st_size/1024:.0f} KB)")
    return result.stdout, out_pdbqt


# ============================================================
# Step 5: Parse Vina text output
# ============================================================

def parse_results(stdout, docked_pdbqt):
    print("\n[RESULTS] Parsing scores...")

    results = []
    in_table = False

    for line in stdout.split("\n"):
        line = line.strip()
        if "-----+------------+----------+----------" in line:
            in_table = True
            continue
        if in_table and (line.startswith("Writing") or not line):
            break
        if in_table and line:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    results.append({
                        "mode":    int(parts[0]),
                        "affinity": float(parts[1]),
                        "rmsd_lb":  float(parts[2]),
                        "rmsd_ub":  float(parts[3]),
                    })
                except (ValueError, IndexError):
                    continue

    if not results:
        print("[RESULTS] No docking scores found!")
        return None, None

    # Write CSV
    csv_path = OUTPUT_DIR / "docking_results.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["mode", "affinity", "rmsd_lb", "rmsd_ub"])
        writer.writeheader()
        writer.writerows(results)

    # Best pose → PDB
    best_pdb = extract_best_pose(docked_pdbqt, results[0]["mode"])

    # Print
    best = results[0]
    print(f"  {'='*50}")
    print(f"  BEST: {best['affinity']:.1f} kcal/mol (mode {best['mode']})")
    print(f"  {'='*50}")
    for r in results[:5]:
        print(f"    Mode {r['mode']:2d}: {r['affinity']:7.1f} kcal/mol  RMSD {r['rmsd_lb']:.1f}-{r['rmsd_ub']:.1f}")

    return results, best_pdb


def extract_best_pose(docked_pdbqt, best_mode):
    best_path = OUTPUT_DIR / f"{PDB_ID}_{LIGAND_NAME}_best.pdb"
    models = {}
    current_mode = None
    current_lines = []

    with open(docked_pdbqt) as f:
        for line in f:
            if line.startswith("MODEL"):
                if current_mode is not None:
                    models[current_mode] = "".join(current_lines)
                current_lines = []
                try:
                    current_mode = int(line.split()[1])
                except:
                    current_mode = None
            elif line.startswith("ENDMDL"):
                if current_mode is not None:
                    models[current_mode] = "".join(current_lines)
                current_lines = []
                current_mode = None
            else:
                current_lines.append(line)

    if best_mode in models:
        best_path.write_text(models[best_mode])
        print(f"[RESULTS] Best pose → {best_path}")
        return best_path
    elif models:
        first_mode = min(models.keys())
        best_path.write_text(models[first_mode])
        print(f"[RESULTS] Mode {best_mode} not found, using mode {first_mode} → {best_path}")
        return best_path

    print("[RESULTS] Could not extract best pose")
    return None


# ============================================================
# Step 6: Report
# ============================================================

def write_report(results, best_pdb_path, center, size):
    from datetime import datetime

    report_path = OUTPUT_DIR / "M6_report.txt"
    lines = [
        "=" * 60,
        "M6: Molecular Docking Report",
        "=" * 60,
        f"Date:       {datetime.now():%Y-%m-%d %H:%M}",
        f"Protein:    GBP1 (PDB {PDB_ID})",
        f"Ligand:     {LIGAND_NAME}",
        f"SMILES:     {LIGAND_SMILES}",
        f"Engine:     AutoDock Vina 1.2.5 (standalone)",
        f"Exhaust.:   {EXHAUSTIVENESS}",
        f"",
        f"Docking Box:",
        f"  Center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})",
        f"  Size:   ({size[0]:.1f}, {size[1]:.1f}, {size[2]:.1f}) Å",
        f"",
        f"Docking Results:",
    ]

    if results:
        for r in results[:10]:
            lines.append(f"  Mode {r['mode']:2d}: {r['affinity']:7.1f} kcal/mol  (RMSD {r['rmsd_lb']:.1f}-{r['rmsd_ub']:.1f})")
        lines.append("")
        lines.append(f"Best affinity: {results[0]['affinity']:.1f} kcal/mol")

        a = results[0]["affinity"]
        if a < -9:     interp = "Very strong binding"
        elif a < -7:   interp = "Strong binding — good candidate"
        elif a < -5:   interp = "Moderate binding"
        else:          interp = "Weak binding"
        lines.append(f"Interpretation: {interp}")

    lines += [
        "", "Output:",
        f"  {OUTPUT_DIR / f'{PDB_ID}_receptor.pdbqt'}",
        f"  {OUTPUT_DIR / f'{LIGAND_NAME}.pdbqt'}",
        f"  {OUTPUT_DIR / f'{PDB_ID}_{LIGAND_NAME}_docked.pdbqt'}",
        f"  {OUTPUT_DIR / 'docking_results.csv'}",
        f"  {OUTPUT_DIR / 'vina_log.txt'}",
        f"  {OUTPUT_DIR / 'M6_report.txt'}",
        "=" * 60,
    ]

    text = "\n".join(lines)
    report_path.write_text(text, encoding="utf-8")
    print(f"\n[REPORT] {report_path}")
    print(text)


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 60)
    print(f"  M6: GBP1 + {LIGAND_NAME} · AutoDock Vina")
    print("=" * 60)

    check_setup()
    receptor_pdbqt, clean_pdb = prepare_receptor()
    ligand_pdbqt = prepare_ligand()
    center, size = define_pocket(clean_pdb)
    stdout, out_pdbqt = run_docking(receptor_pdbqt, ligand_pdbqt, center, size)

    if stdout and out_pdbqt:
        results, best_pdb = parse_results(stdout, out_pdbqt)
        if results:
            write_report(results, best_pdb, center, size)
            print(f"\n[DONE] M6 completed! Affinity = {results[0]['affinity']:.1f} kcal/mol")
    else:
        print("\n[FAILED] Check M6_output/vina_log.txt")
        sys.exit(1)


if __name__ == "__main__":
    main()
