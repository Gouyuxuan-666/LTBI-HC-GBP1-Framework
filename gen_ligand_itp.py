#!/usr/bin/env python3
"""
Generate GROMACS ligand topology (ITP + GRO) from docked PDB.


Usage: python gen_ligand_itp.py <ligand.pdb> <ligand.sdf_or_smiles> [output_dir]
"""

import os, sys, subprocess
from pathlib import Path

def find_gaff2():
    import site
    candidates = []
    for sp in site.getsitepackages():
        candidates.append(Path(sp) / "acpype/amber_linux/dat/leap/parm/gaff2.dat")
    candidates.append(
        Path(os.environ.get("LOCALAPPDATA", "")) /
        "Programs/Python/Python313/Lib/site-packages/acpype/amber_linux/dat/leap/parm/gaff2.dat"
    )
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"gaff2.dat not found. Tried: {candidates[:3]}")


def parse_gaff2(path):
    bonds = []
    angles = []
    dihedrals = []
    atomtypes = {}
    masses = {}

    lines = Path(path).read_text(encoding="utf-8", errors="replace").split("\n")

    in_mod4 = False
    for line in lines:
        ls = line.strip()
        if not ls or ls.startswith("!") or ls.startswith("---") or ls.startswith("SOUR") or ls.startswith("Version"):
            continue
        if ls.startswith("MOD4"):
            in_mod4 = True
            continue
        if ls == "END" and in_mod4:
            break

        parts = ls.split()
        if len(parts) < 2:
            continue

        first = parts[0]
        n_dash = first.count("-")

        if in_mod4:
            try:
                atype = parts[0]
                r_half = float(parts[1])
                eps_kcal = float(parts[2])
                sigma_nm = r_half * 2.0 * (2.0 ** (-1.0 / 6.0)) * 0.1
                eps_kj = eps_kcal * 4.184
                if atype not in atomtypes and eps_kj > 0.0001:
                    atomtypes[atype] = (sigma_nm, eps_kj)
            except (ValueError, IndexError):
                pass
            continue

        if n_dash == 0 and len(parts) >= 3:
            if first[0].islower() and len(first) <= 4:
                try:
                    mass = float(parts[1])
                    masses[first] = mass
                except ValueError:
                    pass

        elif n_dash == 1 and len(parts) >= 3:
            types = first.split("-")
            if len(types) == 2:
                try:
                    bonds.append((types[0], types[1], float(parts[1]), float(parts[2])))
                except ValueError:
                    pass

        elif n_dash == 2 and len(parts) >= 4:
            types = first.split("-")
            if len(types) == 3:
                try:
                    angles.append((types[0], types[1], types[2], float(parts[1]), float(parts[2])))
                except ValueError:
                    pass

        elif n_dash == 3 and len(parts) >= 5:
            types = first.split("-")
            if len(types) == 4:
                try:
                    dihedrals.append((
                        types[0], types[1], types[2], types[3],
                        float(parts[1]), float(parts[2]),
                        float(parts[3]), float(parts[4])
                    ))
                except ValueError:
                    pass

    return {"bonds": bonds, "angles": angles, "dihedrals": dihedrals,
            "atomtypes": atomtypes, "masses": masses}


def find_bond(p, a1, a2):
    for b in p["bonds"]:
        if (b[0] == a1 and b[1] == a2) or (b[1] == a1 and b[0] == a2):
            return b[2], b[3]
    print(f"  WARN: bond {a1}-{a2} not found")
    return 200000.0, 0.140


def find_angle(p, a1, a2, a3):
    for ang in p["angles"]:
        if ang[1] == a2:
            if (ang[0] == a1 and ang[2] == a3) or (ang[0] == a3 and ang[2] == a1):
                return ang[3], ang[4]
    print(f"  WARN: angle {a1}-{a2}-{a3} not found")
    return 50.0, 120.0


def find_diheds(p, a1, a2, a3, a4):
    results = []
    for d in p["dihedrals"]:
        if d[1] == a2 and d[2] == a3:
            if (d[0] == a1 and d[3] == a4) or (d[3] == a1 and d[0] == a4):
                results.append((d[4], d[5], d[6], d[7]))
    if not results:
        for d in p["dihedrals"]:
            if d[0] == "X" and d[1] == a2 and d[2] == a3 and d[3] == "X":
                results.append((d[4], d[5], d[6], d[7]))
    return results


def assign_gaff2(mol, atom):
    elem = atom.GetSymbol()
    hyb = atom.GetHybridization()
    is_ring = atom.IsInRing()
    is_aro = atom.GetIsAromatic()
    nbrs = [n for n in atom.GetNeighbors()]
    n_h = sum(1 for n in nbrs if n.GetSymbol() == "H")

    if elem == "C":
        if is_aro:
            return "ca"
        if is_ring:
            return "c3" if hyb.name == "SP3" else "ca"
        if hyb.name == "SP2":
            for nb in nbrs:
                b = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if nb.GetSymbol() == "O" and b.GetBondTypeAsDouble() == 2.0:
                    return "c"
            return "ce"
        if hyb.name == "SP":
            return "c1"
        return "c3"

    if elem == "O":
        if n_h >= 1:
            return "oh"
        if len(nbrs) >= 2:
            return "os"
        for nb in nbrs:
            b = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
            if b.GetBondTypeAsDouble() == 2.0:
                return "o"
        return "os"

    if elem == "H":
        parent = next((n for n in nbrs if n.GetAtomicNum() > 1), None)
        if parent is None:
            return "hc"
        pt = assign_gaff2(mol, parent)
        if pt == "oh":
            return "ho"
        if pt == "ca":
            return "ha"
        if pt in ("ce", "c2", "c"):
            return "h4"
        if pt == "c3":
            return "hc"
        if pt in ("os", "o"):
            return "h1"
        return "ha"

    if elem == "N":
        return "n"
    if elem == "S":
        return "s"
    return "ca"


def main():
    if len(sys.argv) < 3:
        print(f"Usage: python {sys.argv[0]} <ligand.pdb> <ligand.sdf_or_smiles> [output_dir]")
        sys.exit(1)

    lig_pdb = Path(sys.argv[1])
    lig_src = sys.argv[2]
    out_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else lig_pdb.parent

    print(f"PDB:    {lig_pdb}")
    print(f"Source: {lig_src}")
    print(f"Output: {out_dir}")

    # Step 1: Parse gaff2.dat
    gaff2_path = find_gaff2()
    print(f"GAFF2:  {gaff2_path}")
    p = parse_gaff2(gaff2_path)
    print(f"  {len(p['atomtypes'])} atom types, {len(p['bonds'])} bonds, "
          f"{len(p['angles'])} angles, {len(p['dihedrals'])} dihedrals")

    # Step 2: Read molecule
    from rdkit import Chem
    from rdkit.Chem import AllChem

    lig_src_path = Path(lig_src)
    if lig_src_path.suffix in (".sdf", ".mol"):
        suppl = Chem.SDMolSupplier(str(lig_src_path))
        mol = suppl[0]
        if mol is None:
            mols = [m for m in suppl if m is not None]
            mol = mols[0] if mols else None
        src_label = f"SDF ({mol.GetNumAtoms()} atoms)" if mol else "SDF FAILED"
    else:
        mol = Chem.MolFromSmiles(lig_src)
        src_label = f"SMILES ({mol.GetNumAtoms()} atoms)" if mol else "SMILES FAILED"
    print(f"MOL:    {src_label}")

    if mol is None:
        print("[ERROR] Cannot read molecule!")
        sys.exit(1)

    # Add hydrogens (RDKit computes H positions from heavy atom geometry)
    # Do NOT re-embed — preserve SDF/docking 3D coordinates
    mol = Chem.AddHs(mol, addCoords=True)
    n_atoms = mol.GetNumAtoms()
    print(f"Atoms:  {n_atoms} (with H, coords from SDF)")

    # Step 3: Get docked coords from PDB (heavy atoms only)
    pdb_heavy_list = []  # (elem, x, y, z) in PDB file order, nm
    with open(lig_pdb) as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                try:
                    elem = line[76:78].strip() if len(line) > 77 else line[12:16].strip()
                    if elem == "H":
                        continue
                    x = float(line[30:38]) * 0.1  # Å → nm
                    y = float(line[38:46]) * 0.1
                    z = float(line[46:54]) * 0.1
                    pdb_heavy_list.append((elem, x, y, z))
                except Exception:
                    continue
    print(f"PDB:    {len(pdb_heavy_list)} heavy atoms")

    # Step 4: Gasteiger charges via obabel
    mol2_path = out_dir / "ligand.mol2"
    sdf_tmp = out_dir / "_tmp_rdkit.sdf"
    w = Chem.SDWriter(str(sdf_tmp))
    w.write(mol)
    w.close()
    subprocess.run(
        ["obabel", str(sdf_tmp), "-O", str(mol2_path)],
        capture_output=True, text=True
    )
    sdf_tmp.unlink(missing_ok=True)

    mol2_charges = {}
    if mol2_path.exists():
        in_atoms = False
        with open(mol2_path) as f:
            for line in f:
                ls = line.strip()
                if ls.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if ls.startswith("@<TRIPOS>BOND"):
                    break
                if in_atoms and ls:
                    parts = ls.split()
                    if len(parts) >= 9:
                        mol2_charges[int(parts[0])] = float(parts[8])
    print(f"MOL2:   {len(mol2_charges)} charges")

    # Step 5: GAFF2 atom types
    atom_types = []
    for i, atom in enumerate(mol.GetAtoms()):
        atom_types.append(assign_gaff2(mol, atom))

    for i in range(min(5, n_atoms)):
        a = mol.GetAtomWithIdx(i)
        chg = mol2_charges.get(i+1, 0.0)
        print(f"  Atom {i+1:>3}: {a.GetSymbol():<2} type={atom_types[i]:<4} charge={chg:.4f}")
    if n_atoms > 8:
        print(f"  ... ({n_atoms - 8} more)")
        for i in range(n_atoms - 3, n_atoms):
            a = mol.GetAtomWithIdx(i)
            chg = mol2_charges.get(i+1, 0.0)
            print(f"  Atom {i+1:>3}: {a.GetSymbol():<2} type={atom_types[i]:<4} charge={chg:.4f}")

    # Step 6: Build topology
    bonds_list, angles_list, dihedrals_list, impropers_list = [], [], [], []

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        k, r0 = find_bond(p, atom_types[i], atom_types[j])
        bonds_list.append((i, j, k, r0))

    for b1 in mol.GetBonds():
        for b2 in mol.GetBonds():
            if b1.GetIdx() >= b2.GetIdx():
                continue
            shared, a1, a3 = -1, -1, -1
            if b1.GetBeginAtomIdx() == b2.GetBeginAtomIdx():
                shared, a1, a3 = b1.GetBeginAtomIdx(), b1.GetEndAtomIdx(), b2.GetEndAtomIdx()
            elif b1.GetBeginAtomIdx() == b2.GetEndAtomIdx():
                shared, a1, a3 = b1.GetBeginAtomIdx(), b1.GetEndAtomIdx(), b2.GetBeginAtomIdx()
            elif b1.GetEndAtomIdx() == b2.GetBeginAtomIdx():
                shared, a1, a3 = b1.GetEndAtomIdx(), b1.GetBeginAtomIdx(), b2.GetEndAtomIdx()
            elif b1.GetEndAtomIdx() == b2.GetEndAtomIdx():
                shared, a1, a3 = b1.GetEndAtomIdx(), b1.GetBeginAtomIdx(), b2.GetBeginAtomIdx()
            if shared < 0:
                continue
            k, deg = find_angle(p, atom_types[a1], atom_types[shared], atom_types[a3])
            angles_list.append((a1, shared, a3, k, deg))

    for bond in mol.GetBonds():
        b1, b2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        a_atoms = [n.GetIdx() for n in mol.GetAtomWithIdx(b1).GetNeighbors() if n.GetIdx() != b2]
        d_atoms = [n.GetIdx() for n in mol.GetAtomWithIdx(b2).GetNeighbors() if n.GetIdx() != b1]
        for a1 in a_atoms:
            for a4 in d_atoms:
                params = find_diheds(p, atom_types[a1], atom_types[b1], atom_types[b2], atom_types[a4])
                for mult, k, phase, period in params:
                    dihedrals_list.append((a1, b1, b2, a4, max(mult, 1), k, phase, period))

    for i in range(n_atoms):
        a = mol.GetAtomWithIdx(i)
        nb = [n.GetIdx() for n in a.GetNeighbors()]
        if a.GetHybridization().name == "SP2" and len(nb) == 3:
            impropers_list.append((i, nb[0], nb[1], nb[2]))

    # Step 7: Write GRO — Kabsch align SDF geometry to PDB coordinates
    # SDF: correct bond lengths/angles. PDB: correct binding pocket position.
    # Match heavy atoms by element + distance-from-centroid, then Kabsch rotate.
    import numpy as np
    from collections import defaultdict

    conf = mol.GetConformer()

    # PDB heavy atoms (nm)
    pdb_pts = np.array([(x, y, z) for _, x, y, z in pdb_heavy_list])
    pdb_elems = [e for e, _, _, _ in pdb_heavy_list]

    # SDF heavy atoms (Å→nm)
    sdf_pts, sdf_elems, sdf_indices = [], [], []
    for i in range(n_atoms):
        a = mol.GetAtomWithIdx(i)
        if a.GetAtomicNum() > 1:
            pos = conf.GetAtomPosition(i)
            sdf_pts.append((pos.x * 0.1, pos.y * 0.1, pos.z * 0.1))
            sdf_elems.append(a.GetSymbol())
            sdf_indices.append(i)
    sdf_pts = np.array(sdf_pts)

    # Center
    pdb_c = pdb_pts.mean(axis=0)
    sdf_c = sdf_pts.mean(axis=0)
    pdb_ctr = pdb_pts - pdb_c
    sdf_ctr = sdf_pts - sdf_c

    # Match by element + distance-from-centroid (rotation-invariant)
    pdb_by_elem = defaultdict(list)
    for i, e in enumerate(pdb_elems):
        pdb_by_elem[e].append((np.linalg.norm(pdb_ctr[i]), i))
    sdf_by_elem = defaultdict(list)
    for i, e in enumerate(sdf_elems):
        sdf_by_elem[e].append((np.linalg.norm(sdf_ctr[i]), i))

    correspondences = []  # (pdb_idx, sdf_idx)
    for e in set(pdb_by_elem) & set(sdf_by_elem):
        p_sorted = sorted(pdb_by_elem[e])
        s_sorted = sorted(sdf_by_elem[e])
        for (_, pi), (_, si) in zip(p_sorted, s_sorted):
            correspondences.append((pi, si))

    # Kabsch rotation
    if len(correspondences) >= 3:
        P = np.array([pdb_ctr[i] for i, _ in correspondences])
        Q = np.array([sdf_ctr[j] for _, j in correspondences])
        H = Q.T @ P
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        rmsd = np.sqrt(((Q @ R - P) ** 2).mean())
        print(f"GRO:    Kabsch RMSD = {rmsd:.4f} nm ({len(correspondences)} atoms)")
    else:
        R = np.eye(3)
        print(f"GRO:    Kabsch skipped (<3 matches)")

    # Apply to all atoms
    lines = ["Resveratrol (GAFF2)", f" {n_atoms}"]
    for i in range(n_atoms):
        a = mol.GetAtomWithIdx(i)
        elem = a.GetSymbol()
        atom_name = f"{elem}{i+1}"[:4]
        pos = conf.GetAtomPosition(i)
        pt = np.array([pos.x * 0.1, pos.y * 0.1, pos.z * 0.1])
        aligned = R @ (pt - sdf_c) + pdb_c
        lines.append(f"{541:>5}{'UNL':<5}{atom_name:>5}{i+1:>5}"
                     f"{aligned[0]:8.3f}{aligned[1]:8.3f}{aligned[2]:8.3f}")

    lines.append("   1.00000   1.00000   1.00000")
    gro_path = out_dir / "ligand.gro"
    gro_path.write_text("\n".join(lines))
    print(f"\nWrote: {gro_path}")

    # Prefix GAFF2 types to avoid clashes with AMBER force field (both define ca, oh, ha, etc.)
    atom_types = ["G_" + t for t in atom_types]

    # Step 8a: Write gaff2_atypes.itp (just [atomtypes], goes before FF)
    at_lines = []
    at_lines.append("[ atomtypes ]")
    at_lines.append("; name  at.num  mass     charge  ptype  sigma(nm)   epsilon(kJ/mol)")
    atnum_map = {"C": 6, "c": 6, "H": 1, "h": 1, "O": 8, "o": 8, "N": 7, "n": 7, "S": 16, "s": 16}
    written = set()
    for t in atom_types:
        if t in written:
            continue
        written.add(t)
        mass = p["masses"].get(t, 12.01)
        sigma, eps = p["atomtypes"].get(t, (0.35, 0.5))
        atnum = atnum_map.get(t[0], 6)
        at_lines.append(f"  {t:<6}  {atnum:>6}  {mass:>8.3f}  0.000  A  {sigma:>10.5f}  {eps:>10.5f}")
    at_lines.append("")
    atypes_path = out_dir / "gaff2_atypes.itp"
    atypes_path.write_text("\n".join(at_lines))
    print(f"Wrote: {atypes_path}")

    # Step 8b: Write ligand.itp (moleculetype + atoms + bonds, NO atomtypes)
    itp = []
    itp.append("[ moleculetype ]")
    itp.append("; Name  nrexcl")
    itp.append("  UNL    3")
    itp.append("")
    itp.append("[ atoms ]")
    itp.append("; nr  type  resnr  residue  atom  cgnr  charge  mass")
    for i in range(n_atoms):
        a = mol.GetAtomWithIdx(i)
        elem = a.GetSymbol()
        atom_name = f"{elem}{i+1}"[:4]
        chg = mol2_charges.get(i+1, 0.0)
        mass = a.GetMass()
        itp.append(f"  {i+1:>4}  {atom_types[i]:<6}  541  UNL    {atom_name:<5}  {i+1:>4}  {chg:>10.5f}  {mass:>8.3f}")

    itp.append("")
    itp.append("[ bonds ]")
    itp.append("; ai   aj  funct  r(nm)     K(kJ/mol/nm^2)")
    for i, j, k, r0 in bonds_list:
        itp.append(f"  {i+1:>4}  {j+1:>4}  1  {r0*0.1:>8.5f}  {k*418.4:>12.1f}")

    itp.append("")
    itp.append("[ angles ]")
    itp.append("; ai   aj   ak  funct  deg     K(kJ/mol/rad^2)")
    for a1, a2, a3, k, deg in angles_list:
        itp.append(f"  {a1+1:>4}  {a2+1:>4}  {a3+1:>4}  1  {deg:>8.3f}  {k*4.184:>12.3f}")

    if dihedrals_list:
        itp.append("")
        itp.append("[ dihedrals ]")
        itp.append("; ai  aj  ak  al  funct  phase  k(kJ/mol)  period")
        for a1, a2, a3, a4, idivf, k, phase, period in dihedrals_list:
            itp.append(f"  {a1+1:>4}  {a2+1:>4}  {a3+1:>4}  {a4+1:>4}  9  {phase:>8.3f}  {k*4.184/max(idivf,1):>10.3f}  {int(period):>2}")

    if impropers_list:
        itp.append("")
        itp.append("[ dihedrals ]")
        itp.append("; impropers: ai  aj  ak  al  funct  phase  k(kJ/mol)  period")
        for a1, a2, a3, a4 in impropers_list:
            itp.append(f"  {a1+1:>4}  {a2+1:>4}  {a3+1:>4}  {a4+1:>4}  4  {180.0:>8.3f}  {4.184:>10.3f}  {2:>2}")

    itp.append("")

    itp_path = out_dir / "ligand.itp"
    itp_path.write_text("\n".join(itp))
    print(f"Wrote: {itp_path}")

    print(f"\nSummary: {n_atoms} atoms, {len(bonds_list)} bonds, "
          f"{len(angles_list)} angles, {len(dihedrals_list)} dihedrals, "
          f"{len(impropers_list)} impropers")


if __name__ == "__main__":
    main()
