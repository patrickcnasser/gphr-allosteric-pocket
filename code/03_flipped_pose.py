#!/usr/bin/env python3
"""
03_flipped_pose.py

Characterises the decoy pose that every scoring function prefers.

Asks four things:
  (1) Which fragments land in which subpockets when the ligand is reversed?
  (2) Is the pseudo-symmetry a property of the LIGAND (shape self-similarity
      under a 180 deg rotation) or of the POCKET (comparable chemistry at both
      ends of the site)?
  (3) What do the scoring functions plausibly respond to -- buried contact counts,
      hydrophobic contacts, polar pairs -- native vs flipped?
  (4) Is either end of the pocket able to tell the two orientations apart?

INPUTS   ../inputs/7FIH.pdb
         ../results/dock/<config>.sdf   (from run_docking_sweep.sh)
         ../inputs/xtal_lig.sdf
OUTPUT   ../results/flipped_pose.txt
         ../results/flipped_pose_fragment_map.tsv
         ../results/decoy_best.pdb
"""
import os, sys, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gphr_lib import parse_pdb, derive_bonds, write_pdb_ligand
from rmsd_lib import read_sdf as _rs, symmetry_rmsd, mapped_coords, centroid_offset, atom_overlap

HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, '..', 'inputs'); OUT = os.path.join(HERE, '..', 'results')
def op(n): return os.path.join(OUT, n)

FRAG = {'2-methylthio': ['S17', 'C18'],
        'thienopyrimidine core': ['C8', 'C9', 'C10', 'C11', 'C13', 'C15', 'N12', 'N14', 'S16'],
        '5-amino': ['N35'],
        'tert-butyl+carboxamide': ['C1', 'C2', 'C3', 'C4', 'N5', 'C6', 'O7'],
        'pendant phenyl': ['C19', 'C20', 'C21', 'C22', 'C23', 'C24'],
        'anilide linker': ['N25', 'C26', 'O27', 'C28'],
        'morpholine': ['N29', 'C30', 'C31', 'O32', 'C33', 'C34']}

HYDROPHOBIC = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'GLY', 'CYS'}
POLAR = {'SER', 'THR', 'ASN', 'GLN', 'TYR', 'HIS', 'LYS', 'ARG', 'ASP', 'GLU'}


L_res, L_ord, L_het = parse_pdb(os.path.join(IN, '7FIH.pdb'), 'R')
L55 = {n: v for n, v in next(v for k, v in L_het.items() if k[0] == '55Z').items() if v[1] != 'H'}
NAMES = list(L55.keys())
X = np.array([L55[n][0] for n in NAMES])
E = [L55[n][1] for n in NAMES]

# ------------------------------------------------- pick the canonical decoy
REF = _rs(os.path.join(IN, 'xtal_lig.sdf'))[0]
cands = []
for f in sorted(glob.glob(op('dock/*.sdf'))):
    for pi, rec in enumerate(_rs(f), 1):
        r, _na, ok = symmetry_rmsd(rec, REF)
        if not ok:
            continue
        P = mapped_coords(rec, REF)          # pose coords in REFERENCE atom order
        if P is None:
            continue
        cen = centroid_offset(rec, REF)
        ov = atom_overlap(rec, REF)
        cands.append((r, cen, ov, os.path.basename(f), pi, P))
if not cands:
    sys.exit("no docking poses found -- run run_docking_sweep.sh first")
flip = max([c for c in cands if c[1] < 1.5 and c[2] > 0.85], key=lambda c: c[0])
rms, cen, ov, src, pi, D = flip
out = [f"Selected decoy: {src} pose {pi}",
       f"  symmetry-corrected heavy-atom RMSD to crystal pose : {rms:.2f} A",
       f"  centroid offset                                    : {cen:.2f} A",
       f"  fraction of atoms within 2 A of some crystal atom  : {100*ov:.0f}%",
       f"  ({len(cands)} poses examined across {len(set(c[3] for c in cands))} configurations)",
       "  NOTE: coordinates are mapped into the reference atom order by graph isomorphism;",
       "        gnina does not preserve input atom order.", ""]
write_pdb_ligand(op('decoy_best.pdb'), {n: (D[i], E[i]) for i, n in enumerate(NAMES)}, '55Z')

# --------------------------------------------------- (1) fragment relocation
idx = {n: i for i, n in enumerate(NAMES)}
nat_cent = {f: X[[idx[a] for a in ats]].mean(0) for f, ats in FRAG.items()}
dec_cent = {f: D[[idx[a] for a in ats]].mean(0) for f, ats in FRAG.items()}
rows = ["fragment\tdisplacement_A\tlands_nearest_to_native_fragment\tdistance_to_that_centroid"]
out.append("(1) Where each fragment goes in the decoy pose")
out.append(f"{'fragment':<24}{'moves':>8}   lands in the subpocket natively occupied by")
for f in FRAG:
    disp = np.linalg.norm(dec_cent[f] - nat_cent[f])
    near = min(nat_cent, key=lambda g: np.linalg.norm(dec_cent[f] - nat_cent[g]))
    nd = np.linalg.norm(dec_cent[f] - nat_cent[near])
    out.append(f"{f:<24}{disp:>7.1f} A   {near} ({nd:.1f} A)")
    rows.append(f"{f}\t{disp:.2f}\t{near}\t{nd:.2f}")
open(op('flipped_pose_fragment_map.tsv'), 'w').write("\n".join(rows) + "\n")

# ------------------------------------------- (2a) ligand shape self-symmetry
def principal_axes(P):
    C = P - P.mean(0)
    w, v = np.linalg.eigh(C.T @ C)
    return v[:, ::-1], w[::-1]          # columns, largest first


def self_similarity(P, el, axis, centre):
    """RMSD of the point cloud to itself after 180 deg rotation about `axis`,
    matching each rotated atom to the nearest unrotated atom OF THE SAME ELEMENT."""
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    Rr = np.eye(3) + 2 * K @ K            # Rodrigues at theta=pi
    Q = (P - centre) @ Rr.T + centre
    tot, n = 0.0, 0
    for i in range(len(P)):
        same = [j for j in range(len(P)) if el[j] == el[i]]
        d = np.linalg.norm(P[same] - Q[i], axis=1).min()
        tot += d ** 2; n += 1
    return float(np.sqrt(tot / n))


ax, ev = principal_axes(X)
cen0 = X.mean(0)
out.append("\n(2a) Is the LIGAND itself pseudo-symmetric?")
out.append(f"  principal moments (relative): {ev[0]/ev[0]:.2f}, {ev[1]/ev[0]:.2f}, {ev[2]/ev[0]:.2f}"
           f"   -> {'elongated' if ev[1]/ev[0] < 0.5 else 'globular'}")
for i, nm in enumerate(['major', 'middle', 'minor']):
    s = self_similarity(X, E, ax[:, i], cen0)
    out.append(f"  180 deg about {nm:<6} axis: self-RMSD {s:.2f} A")
best_self = min(self_similarity(X, E, ax[:, i], cen0) for i in range(3))
out.append(f"  best self-similarity under any 2-fold: {best_self:.2f} A")
out.append("  interpretation: a genuinely pseudo-symmetric ligand would map onto itself"
           " well under one of these (< ~1.5 A).")

# ----------------------------------------- (2b) pocket chemistry at both ends
long_axis = ax[:, 0]
proj_l = (X - cen0) @ long_axis
end_lo = X[proj_l < np.percentile(proj_l, 33)]
end_hi = X[proj_l > np.percentile(proj_l, 67)]
def env(P, radius=5.0):
    comp = {}
    for k, ats in L_res.items():
        d = min(np.linalg.norm(P - v[0], axis=1).min() for a, v in ats.items() if v[1] != 'H')
        if d < radius:
            comp[k] = d
    nh = sum(1 for k in comp if k[1] in HYDROPHOBIC)
    npz = sum(1 for k in comp if k[1] in POLAR)
    return comp, nh, npz
cl, hl, pl = env(end_lo); ch, hh, ph_ = env(end_hi)
out.append("\n(2b) Is the POCKET pseudo-symmetric? (chemistry at the two ends of the ligand axis)")
out.append(f"  end A ({len(cl)} residues within 5 A): {hl} hydrophobic, {pl} polar "
           f"({100*hl/max(1,hl+pl):.0f}% hydrophobic)")
out.append(f"  end B ({len(ch)} residues within 5 A): {hh} hydrophobic, {ph_} polar "
           f"({100*hh/max(1,hh+ph_):.0f}% hydrophobic)")
out.append(f"  residues shared by both ends: {len(set(cl) & set(ch))}")

# -------------------------------------- (3) what the scoring functions see
def descriptors(P):
    prot = [(k, a, v[0], v[1]) for k, ats in L_res.items() for a, v in ats.items() if v[1] != 'H']
    px = np.array([p[2] for p in prot])
    pel = [p[3] for p in prot]
    pres = [p[0][1] for p in prot]
    nb5 = int(sum((np.linalg.norm(px - x, axis=1) < 5.0).sum() for x in P))
    nb4 = int(sum((np.linalg.norm(px - x, axis=1) < 4.0).sum() for x in P))
    hyd = 0; pol = 0; clash = 0
    for i, x in enumerate(P):
        d = np.linalg.norm(px - x, axis=1)
        for j in np.where(d < 4.5)[0]:
            if E[i] == 'C' and pel[j] == 'C':
                hyd += 1
            if E[i] in 'NO' and pel[j] in 'NO' and d[j] < 3.5:
                pol += 1
        clash += int((d < 2.6).sum())
    return nb5, nb4, hyd, pol, clash


nn = descriptors(X); dd = descriptors(D)
out.append("\n(3) Descriptors the scoring functions are built on -- native vs decoy")
out.append(f"{'descriptor':<34}{'native':>9}{'decoy':>9}{'delta':>9}")
for nm, a, b in zip(["protein neighbours < 5 A", "protein neighbours < 4 A",
                     "apolar C...C contacts < 4.5 A", "polar N/O...N/O pairs < 3.5 A",
                     "steric clashes < 2.6 A"], nn, dd):
    out.append(f"{nm:<34}{a:>9}{b:>9}{b-a:>+9}")
out.append("  The classical terms are dominated by the first three. If the decoy matches or")
out.append("  exceeds the native on those while differing only in polar-pair count, the")
out.append("  functions have no term that distinguishes the orientations.")

# ------------------------------- (4) can any single contact tell them apart?
out.append("\n(4) Contacts unique to one orientation (< 4.0 A in one, absent in the other)")
uniq_n, uniq_d = [], []
for k, ats in L_res.items():
    dn = min(np.linalg.norm(X - v[0], axis=1).min() for a, v in ats.items() if v[1] != 'H')
    dd_ = min(np.linalg.norm(D - v[0], axis=1).min() for a, v in ats.items() if v[1] != 'H')
    if dn < 4.0 and dd_ > 5.0:
        uniq_n.append((k, dn, dd_))
    if dd_ < 4.0 and dn > 5.0:
        uniq_d.append((k, dn, dd_))
out.append(f"  residues contacting ONLY the native pose: {len(uniq_n)}")
for k, a, b in sorted(uniq_n, key=lambda x: x[1]):
    out.append(f"     {k[1]}{k[0]}  native {a:.2f} A, decoy {b:.2f} A")
out.append(f"  residues contacting ONLY the decoy pose: {len(uniq_d)}")
for k, a, b in sorted(uniq_d, key=lambda x: x[2]):
    out.append(f"     {k[1]}{k[0]}  native {a:.2f} A, decoy {b:.2f} A")

open(op('flipped_pose.txt'), 'w').write("\n".join(out) + "\n")
print("\n".join(out))
