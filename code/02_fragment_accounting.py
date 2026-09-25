#!/usr/bin/env python3
"""
02_fragment_accounting.py

Per-fragment burial and contact accounting for BOTH ligands in their own receptors.
Rings are perceived from coordinate-derived connectivity (no hand-drawn partition for
Cpd-21f); acyclic atoms are grouped into the connected components that remain when
ring atoms are removed, then labelled by composition.

INPUTS   ../inputs/7FIH.pdb , ../inputs/8I2G.pdb
OUTPUTS  ../results/fragments_55Z.tsv
         ../results/fragments_O6F.tsv
         ../results/fragment_divergence.txt
"""
import os, sys, numpy as np
from collections import deque
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gphr_lib import (parse_pdb, sequence, needleman_wunsch, shell, derive_bonds,
                      neighbour_counts, BLOSUM62)

HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, '..', 'inputs'); OUT = os.path.join(HERE, '..', 'results')
def op(n): return os.path.join(OUT, n)

L_res, L_ord, L_het = parse_pdb(os.path.join(IN, '7FIH.pdb'), 'R')
F_res, F_ord, F_het = parse_pdb(os.path.join(IN, '8I2G.pdb'), 'R')
L55 = {n: v for n, v in next(v for k, v in L_het.items() if k[0] == '55Z').items() if v[1] != 'H'}
O6F = {n: v for n, v in next(v for k, v in F_het.items() if k[0] == 'O6F').items() if v[1] != 'H'}
PAIRS = needleman_wunsch(sequence(L_ord), sequence(F_ord), matrix=BLOSUM62,
                         gap_open=-11, gap_extend=-1)
AMAP = {L_ord[i][0]: F_ord[j][0] for i, j in PAIRS}
RMAP = {v: k for k, v in AMAP.items()}
LNAME = {k[0]: k[1] for k in L_ord}; FNAME = {k[0]: k[1] for k in F_ord}


def rings(bonds):
    """Smallest set of rings via BFS cycle detection; returns list of atom sets."""
    found = []
    for start in bonds:
        # BFS shortest cycle through `start`
        for nb in bonds[start]:
            q = deque([(nb, [start, nb])])
            seen = {(start, nb)}
            while q:
                cur, path = q.popleft()
                if len(path) > 7:
                    continue
                for nxt in bonds[cur]:
                    if nxt == path[-2]:
                        continue
                    if nxt == start and len(path) >= 5:
                        s = frozenset(path)
                        if s not in found:
                            found.append(s)
                        continue
                    if nxt in path:
                        continue
                    if (cur, nxt) in seen:
                        continue
                    seen.add((cur, nxt))
                    q.append((nxt, path + [nxt]))
    # keep minimal rings (5 and 6 membered), drop supersets that are unions
    found = sorted(set(found), key=len)
    keep = []
    for r in found:
        if len(r) > 6:
            continue
        keep.append(r)
    return keep


def fragmentize(lig):
    bonds, names, E = derive_bonds(lig)
    el = dict(zip(names, E))
    rs = rings(bonds)
    # merge fused rings that share >=2 atoms into ring systems
    systems = []
    for r in rs:
        placed = False
        for s in systems:
            if len(s & r) >= 2:
                s |= set(r); placed = True; break
        if not placed:
            systems.append(set(r))
    merged = True
    while merged:
        merged = False
        for i in range(len(systems)):
            for j in range(i + 1, len(systems)):
                if len(systems[i] & systems[j]) >= 2:
                    systems[i] |= systems[j]; systems.pop(j); merged = True; break
            if merged:
                break
    ringatoms = set().union(*systems) if systems else set()
    frags = {}
    for i, s in enumerate(systems, 1):
        comp = "".join(sorted(set(el[a] for a in s)))
        frags[f"ring{i}_{len(s)}mem_{comp}"] = sorted(s)
    # acyclic components
    rest = [a for a in names if a not in ringatoms]
    seen = set()
    ci = 0
    for a in rest:
        if a in seen:
            continue
        comp = []
        q = deque([a]); seen.add(a)
        while q:
            c = q.popleft(); comp.append(c)
            for nb in bonds[c]:
                if nb in rest and nb not in seen:
                    seen.add(nb); q.append(nb)
        ci += 1
        el_s = "".join(sorted(set(el[x] for x in comp)))
        frags[f"chain{ci}_{len(comp)}at_{el_s}"] = sorted(comp)
    return frags, bonds, el


def report(tag, lig, res, other_name_map, amap, outfile):
    frags, bonds, el = fragmentize(lig)
    total = sum(len(v) for v in frags.values())
    assert total == len(lig) and len(set(sum(frags.values(), []))) == len(lig), (total, len(lig))
    nb = neighbour_counts(lig, res, 5.0)
    sh = shell(res, lig, 4.5)
    lines = [f"# {tag}: {len(lig)} heavy atoms, {len(frags)} fragments, partition verified disjoint+complete"]
    lines.append("fragment\tn_atoms\tatoms\tmean_neighbours_5A\tcontacts_<4.5A\tn_divergent\tdivergent_partners")
    summary = {}
    for fname, ats in sorted(frags.items(), key=lambda x: -np.mean([nb[a] for a in x[1]])):
        mean_nb = float(np.mean([nb[a] for a in ats]))
        hits = []
        for k, d in sh.items():
            dd = min(np.linalg.norm(v[0] - lig[a][0])
                     for an, v in res[k].items() if v[1] != 'H' for a in ats)
            if dd < 4.5:
                partner = amap.get(k[0])
                oname = other_name_map.get(partner)
                div = (oname is not None and oname != k[1])
                hits.append((f"{k[1]}{k[0]}", dd, div, f"{oname}{partner}" if oname else "unaligned"))
        ndiv = sum(1 for h in hits if h[2])
        summary[fname] = (mean_nb, len(hits), ndiv)
        lines.append(f"{fname}\t{len(ats)}\t{','.join(ats)}\t{mean_nb:.1f}\t"
                     f"{';'.join(f'{h[0]}({h[1]:.2f})' for h in sorted(hits, key=lambda x: x[1]))}\t"
                     f"{ndiv}\t{';'.join(h[3] for h in hits if h[2]) or '-'}")
    open(outfile, 'w').write("\n".join(lines) + "\n")
    return summary, frags


print("=== Org 43553 (55Z) in LHCGR ===")
s55, f55 = report("55Z in LHCGR", L55, L_res, FNAME, AMAP, op('fragments_55Z.tsv'))
for k, (n, c, d) in s55.items():
    print(f"  {k:<26} nb/atom {n:5.1f}  contacts {c:2d}  divergent {d}")

print("\n=== Cpd-21f (O6F) in FSHR ===")
sO6, fO6 = report("O6F in FSHR", O6F, F_res, LNAME, RMAP, op('fragments_O6F.tsv'))
for k, (n, c, d) in sO6.items():
    print(f"  {k:<26} nb/atom {n:5.1f}  contacts {c:2d}  divergent {d}")

# the comparison the paper is missing
tot55 = sum(v[1] for v in s55.values()); div55 = sum(v[2] for v in s55.values())
totO6 = sum(v[1] for v in sO6.values()); divO6 = sum(v[2] for v in sO6.values())
txt = ["Fragment-level contact divergence: does the FSHR-selective ligand concentrate",
       "its contacts on divergent positions in a way Org 43553 does not?", "",
       f"Org 43553 in LHCGR : {div55}/{tot55} fragment-contacts are at divergent positions "
       f"({100*div55/tot55:.0f}%)",
       f"Cpd-21f  in FSHR   : {divO6}/{totO6} fragment-contacts are at divergent positions "
       f"({100*divO6/totO6:.0f}%)", ""]
# per-ligand: how many fragments touch NO divergent position
z55 = [k for k, v in s55.items() if v[2] == 0]
zO6 = [k for k, v in sO6.items() if v[2] == 0]
txt.append(f"fragments of Org 43553 touching NO divergent position: {len(z55)}/{len(s55)}")
for k in z55:
    txt.append(f"   {k}")
txt.append(f"fragments of Cpd-21f touching NO divergent position: {len(zO6)}/{len(sO6)}")
for k in zO6:
    txt.append(f"   {k}")
open(op('fragment_divergence.txt'), 'w').write("\n".join(txt) + "\n")
print("\n" + "\n".join(txt))
