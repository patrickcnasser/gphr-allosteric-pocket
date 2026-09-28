#!/usr/bin/env python3
"""
08_manuscript_fragment_table.py

Reproduces the manuscript's per-fragment burial table for Org 43553 (55Z), using the
chemically named partition the text refers to, in BOTH frames:

  - LHCGR      : 55Z in its own receptor (7FIH chain R)
  - FSHR-mapped: 55Z transformed into FSHR (8I2G chain R) by the TM superposition

Script 02 perceives fragments from connectivity alone and therefore reports a coarser
partition: it keeps the thienopyrimidine as one 9-atom ring system and groups the
tert-butyl with the carboxamide it is attached to. The manuscript quotes the finer,
chemically named split, so that split is defined explicitly here.

Both columns are printed side by side and labelled, because an earlier draft of the
manuscript reported the FSHR-mapped column under an LHCGR heading.

INPUTS   ../inputs/7FIH.pdb , ../inputs/8I2G.pdb
OUTPUT   ../results/fragment_table_manuscript.tsv
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gphr_lib import (parse_pdb, sequence, needleman_wunsch, shell, kabsch,
                      neighbour_counts, BLOSUM62)

HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, '..', 'inputs'); OUT = os.path.join(HERE, '..', 'results')

L_res, L_ord, L_het = parse_pdb(os.path.join(IN, '7FIH.pdb'), 'R')
F_res, F_ord, F_het = parse_pdb(os.path.join(IN, '8I2G.pdb'), 'R')
L55 = {n: v for n, v in next(v for k, v in L_het.items() if k[0] == '55Z').items()
       if v[1] != 'H'}

# ---- transmembrane superposition, same construction as script 01 -------------
PAIRS = needleman_wunsch(sequence(L_ord), sequence(F_ord), matrix=BLOSUM62,
                         gap_open=-11, gap_extend=-1)
tm = [(i, j) for i, j in PAIRS if L_ord[i][0] >= 340]
P, Q = [], []
for i, j in tm:
    a = L_res.get(L_ord[i][:2]) or L_res.get(L_ord[i])
    b = F_res.get(F_ord[j][:2]) or F_res.get(F_ord[j])
    ca_a = (a or {}).get('CA'); ca_b = (b or {}).get('CA')
    if ca_a is not None and ca_b is not None:
        P.append(ca_a[0]); Q.append(ca_b[0])
R, t = kabsch(np.array(P), np.array(Q))
M55 = {n: (R @ v[0] + t, v[1]) for n, v in L55.items()}

# ---- the manuscript's chemically named partition -----------------------------
FRAGS = [
    ("2-Methylthio",   ["S17", "C18"]),
    ("tert-Butyl",     ["C1", "C2", "C3", "C4"]),
    ("Pyrimidine",     ["C11", "N12", "C13", "N14"]),
    ("Morpholine",     ["N29", "C30", "C31", "O32", "C33", "C34"]),
    ("Thiophene",      ["C8", "C9", "C10", "C15", "S16"]),
    ("Anilide linker", ["N25", "C26", "O27", "C28"]),
    ("Pendant phenyl", ["C19", "C20", "C21", "C22", "C23", "C24"]),
    ("Carboxamide",    ["C6", "N5", "O7"]),
    ("Amino",          ["N35"]),
]
# Connectivity derived from the deposited coordinates gives a thiophene
# (C8,C9,C10,C15,S16) fused to a pyrimidine (C10,C11,C13,C15,N12,N14) sharing C10 and
# C15. To keep the partition disjoint and complete over all 35 heavy atoms, the two
# fused carbons are counted with the thiophene, so the "Pyrimidine" row covers that
# ring's four exclusive atoms. This is the convention the manuscript table uses.
assigned = [a for _, ats in FRAGS for a in ats]
assert len(assigned) == len(set(assigned)), "partition is not disjoint"
assert set(assigned) == set(L55), "partition does not cover the 35 heavy atoms exactly"
assert len(assigned) == 35, f"expected 35 heavy atoms, got {len(assigned)}"


def mean_neighbours(lig, res, atoms):
    prot = [x for k, ats in res.items() for a, (x, e) in ats.items() if e != 'H']
    prot = np.array(prot)
    out = []
    for a in atoms:
        d = np.linalg.norm(prot - lig[a][0], axis=1)
        out.append(int((d < 5.0).sum()))
    return float(np.mean(out))


rows = []
for name, atoms in FRAGS:
    rows.append((name, len(set(atoms)),
                 mean_neighbours(L55, L_res, atoms),
                 mean_neighbours(M55, F_res, atoms)))

path = os.path.join(OUT, 'fragment_table_manuscript.tsv')
with open(path, 'w') as f:
    f.write("# mean protein heavy atoms within 5 A per ligand atom, by fragment\n")
    f.write("# C10 and C15 are the fused thiophene/pyrimidine atoms, counted with the thiophene only\n")
    f.write("fragment\tn_atoms\tmean_neighbours_5A_LHCGR\tmean_neighbours_5A_FSHR_mapped\n")
    for name, n, a, b in sorted(rows, key=lambda r: -r[2]):
        f.write(f"{name}\t{n}\t{a:.1f}\t{b:.1f}\n")

print(f"{'fragment':16s} {'n':>3s} {'LHCGR':>7s} {'FSHR-mapped':>12s}")
for name, n, a, b in sorted(rows, key=lambda r: -r[2]):
    print(f"{name:16s} {n:3d} {a:7.1f} {b:12.1f}")
print("\nwrote", os.path.relpath(path, HERE))
