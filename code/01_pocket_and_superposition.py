#!/usr/bin/env python3
"""
01_pocket_and_superposition.py

Independent re-derivation of the pocket comparison, the TM superposition and both
reciprocal ligand mappings, written from the deposited coordinates alone.

INPUTS   ../inputs/7FIH.pdb , ../inputs/8I2G.pdb
OUTPUTS  ../results/
           validation.txt              published-distance control
           alignment_LHCGR_FSHR.txt    full pairwise alignment, both scoring schemes
           contact_shells.tsv          per-residue distances, both receptors
           union_positions.tsv         27-position union table
           superposition_matrix.txt    rotation + translation, both directions
           55Z_in_FSHR.pdb             Org 43553 mapped into the FSHR frame
           O6F_in_LHCGR.pdb            Cpd-21f mapped into the LHCGR frame
           mapping_contacts.tsv        recomputed contacts after mapping
           backbone_deviation.tsv      per-residue Calpha deviation profile
           conservation_test.txt       binomial test vs TM background
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gphr_lib import (AA3, parse_pdb, sequence, needleman_wunsch, kabsch, shell,
                      neighbour_counts, write_pdb_ligand, BLOSUM62)

HERE = os.path.dirname(os.path.abspath(__file__))
IN   = os.path.join(HERE, '..', 'inputs')
OUT  = os.path.join(HERE, '..', 'results')
os.makedirs(OUT, exist_ok=True)
def op(n): return os.path.join(OUT, n)

L_res, L_ord, L_het = parse_pdb(os.path.join(IN, '7FIH.pdb'), 'R')
F_res, F_ord, F_het = parse_pdb(os.path.join(IN, '8I2G.pdb'), 'R')
L55 = {n: v for n, v in next(v for k, v in L_het.items() if k[0] == '55Z').items() if v[1] != 'H'}
O6F = {n: v for n, v in next(v for k, v in F_het.items() if k[0] == 'O6F').items() if v[1] != 'H'}

# ---------------------------------------------------------------- validation
PUBLISHED = {350: 3.34, 515: 3.21, 585: 4.48, 589: 3.84, 612: 3.78}
lig_xyz = np.array([v[0] for v in L55.values()])
val = ["Control 1 -- published LHCGR contact distances (Duan et al. 2021)",
       f"{'residue':<10}{'published':>10}{'computed':>10}{'delta':>8}"]
maxdelta = 0.0
for num, pub in PUBLISHED.items():
    k = next(k for k in L_res if k[0] == num)
    d = min(np.linalg.norm(lig_xyz - v[0], axis=1).min()
            for a, v in L_res[k].items() if v[1] != 'H')
    maxdelta = max(maxdelta, abs(d - pub))
    val.append(f"{k[1]}{num:<7}{pub:>10.2f}{d:>10.2f}{d-pub:>8.2f}")
val.append(f"max |delta| = {maxdelta:.4f} A   -> {'PASS' if maxdelta < 0.005 else 'FAIL'}")

# ------------------------------------------------- alignment, both schemes
sa, sb = sequence(L_ord), sequence(F_ord)
schemes = {}
for tag, kw in [('identity(match2/mis-1/gap-8/ext-1)',
                 dict(matrix=None, gap_open=-8, gap_extend=-1, match=2, mismatch=-1)),
                ('BLOSUM62(gap-11/ext-1)',
                 dict(matrix=BLOSUM62, gap_open=-11, gap_extend=-1))]:
    pr = needleman_wunsch(sa, sb, **kw)
    ident = sum(1 for i, j in pr if sa[i] == sb[j])
    schemes[tag] = (pr, ident)

val.append("\nControl 2 -- alignment scheme comparison")
for tag, (pr, ident) in schemes.items():
    amap = {L_ord[i][0]: F_ord[j][0] for i, j in pr}
    y = amap.get(612); f = amap.get(515)
    yname = next((k[1] for k in F_ord if k[0] == y), '?')
    fname = next((k[1] for k in F_ord if k[0] == f), '?')
    val.append(f"  {tag}: {len(pr)} pairs, {100*ident/len(pr):.1f}% identity; "
               f"612->{yname}{y}, 515->{fname}{f}")
PAIRS, _ = schemes['BLOSUM62(gap-11/ext-1)']
pr_id, _ = schemes['identity(match2/mis-1/gap-8/ext-1)']
same = set(PAIRS) == set(pr_id)
val.append(f"  the two schemes give {'IDENTICAL' if same else 'DIFFERENT'} pairings")
AMAP = {L_ord[i][0]: F_ord[j][0] for i, j in PAIRS}
RMAP = {v: k for k, v in AMAP.items()}

with open(op('alignment_LHCGR_FSHR.txt'), 'w') as f:
    f.write("LHCGR(7FIH chain R) vs FSHR(8I2G chain R), Needleman-Wunsch, BLOSUM62, "
            "gap_open -11, gap_extend -1\n")
    f.write(f"{len(PAIRS)} aligned pairs\n\n")
    f.write(f"{'LHCGR':>10} {'FSHR':>10}  identical\n")
    for i, j in PAIRS:
        a, b = L_ord[i], F_ord[j]
        f.write(f"{a[1]}{a[0]:<7} {b[1]}{b[0]:<7}  {'*' if a[1]==b[1] else '.'}\n")

# --------------------------------------------------------- contact shells
L_c45, L_c80 = shell(L_res, L55, 4.5), shell(L_res, L55, 8.0)
F_c45, F_c80 = shell(F_res, O6F, 4.5), shell(F_res, O6F, 8.0)
with open(op('contact_shells.tsv'), 'w') as f:
    f.write("receptor\tresidue\tresnum\tmin_dist_to_own_ligand\tshell\n")
    for tag, sh45, sh80 in (('LHCGR', L_c45, L_c80), ('FSHR', F_c45, F_c80)):
        for k, d in sorted(sh80.items(), key=lambda x: x[1]):
            f.write(f"{tag}\t{k[1]}\t{k[0]}\t{d:.2f}\t{'contact' if k in sh45 else '8A'}\n")

# union over aligned positions
union = []
for k, d in L_c45.items():
    fn = AMAP.get(k[0])
    fk = next((x for x in F_ord if x[0] == fn), None)
    union.append((k, fk, d, F_c45.get(fk)))
for k, d in F_c45.items():
    ln = RMAP.get(k[0])
    if ln is not None and any(u[0][0] == ln for u in union):
        continue
    lk = next((x for x in L_ord if x[0] == ln), None)
    union.append((lk, k, L_c45.get(lk), d))
ident_u = sum(1 for lk, fk, dl, df in union if lk and fk and lk[1] == fk[1])
diverg  = [(lk, fk, dl, df) for lk, fk, dl, df in union if lk and fk and lk[1] != fk[1]]
with open(op('union_positions.tsv'), 'w') as f:
    f.write("LHCGR\tFSHR\td_to_Org43553\td_to_Cpd21f\tstatus\n")
    for lk, fk, dl, df in sorted(union, key=lambda x: (x[0][0] if x[0] else 9999)):
        ln = f"{lk[1]}{lk[0]}" if lk else "-"
        fn = f"{fk[1]}{fk[0]}" if fk else "-"
        st = "identical" if (lk and fk and lk[1] == fk[1]) else ("divergent" if lk and fk else "unaligned")
        f.write(f"{ln}\t{fn}\t{dl if dl is None else f'{dl:.2f}'}\t"
                f"{df if df is None else f'{df:.2f}'}\t{st}\n")

# ------------------------------------------------------------ superposition
ligc = lig_xyz.mean(0)
P, Q, lab = [], [], []
for i, j in PAIRS:
    lk, fk = L_ord[i], F_ord[j]
    if 'CA' not in L_res[lk] or 'CA' not in F_res[fk]:
        continue
    ca = L_res[lk]['CA'][0]
    if lk[0] < 340 or np.linalg.norm(ca - ligc) > 25:
        continue
    P.append(ca); Q.append(F_res[fk]['CA'][0]); lab.append((lk, fk))
P, Q = np.array(P), np.array(Q)
keep = np.ones(len(P), bool)
for _ in range(8):
    R, t = kabsch(P[keep], Q[keep])
    d = np.linalg.norm(P @ R.T + t - Q, axis=1)
    nk = d < 2.0
    if (nk == keep).all():
        break
    keep = nk
R, t = kabsch(P[keep], Q[keep])
d = np.linalg.norm(P @ R.T + t - Q, axis=1)
rmsd = float(np.sqrt((d[keep] ** 2).mean()))
val.append(f"\nControl 3 -- TM superposition: {int(keep.sum())} CA pairs, "
           f"RMSD {rmsd:.2f} A, {int((~keep).sum())} pruned at 2.0 A")

with open(op('superposition_matrix.txt'), 'w') as f:
    f.write("LHCGR(7FIH) -> FSHR(8I2G), Calpha, TM bundle (resnum>=340, CA within 25 A of ligand centroid)\n")
    f.write(f"pairs used: {int(keep.sum())}   RMSD: {rmsd:.4f} A\n\nrotation R (row-major):\n")
    for row in R:
        f.write("  " + "  ".join(f"{v: .8f}" for v in row) + "\n")
    f.write("\ntranslation t:\n  " + "  ".join(f"{v: .8f}" for v in t) + "\n")
    f.write("\napply as: x_FSHR = R @ x_LHCGR + t\ninverse:  x_LHCGR = R.T @ (x_FSHR - t)\n")

# per-residue backbone deviation
with open(op('backbone_deviation.tsv'), 'w') as f:
    f.write("LHCGR\tFSHR\tCA_deviation_A\tused_in_fit\tmin_dist_to_ligand\n")
    for (lk, fk), dev, k_ in sorted(zip(lab, d, keep), key=lambda x: x[0][0][0]):
        dl = L_c80.get(lk)
        f.write(f"{lk[1]}{lk[0]}\t{fk[1]}{fk[0]}\t{dev:.3f}\t{int(k_)}\t"
                f"{'' if dl is None else f'{dl:.2f}'}\n")

# ------------------------------------------------------- reciprocal mapping
M55 = {n: (R @ v[0] + t, v[1]) for n, v in L55.items()}
MO6 = {n: (R.T @ (v[0] - t), v[1]) for n, v in O6F.items()}
rt = max(np.linalg.norm(R.T @ (v[0] - t) - L55[n][0]) for n, v in M55.items())
val.append(f"Control 4 -- round-trip of 55Z through forward+inverse: max error {rt:.4e} A")
tb = ['C1', 'C2', 'C3', 'C4']
bonds_o6 = None
o6_tb = ['C30', 'C31', 'C32', 'C33']
off = np.linalg.norm(np.mean([M55[a][0] for a in tb], 0) - np.mean([O6F[a][0] for a in o6_tb], 0))
val.append(f"Control 5 -- mapped 55Z tert-butyl vs Cpd-21f tert-butyl centroid: {off:.2f} A")
write_pdb_ligand(op('55Z_in_FSHR.pdb'), M55, '55Z')
write_pdb_ligand(op('O6F_in_LHCGR.pdb'), MO6, 'O6F')

with open(op('mapping_contacts.tsv'), 'w') as f:
    f.write("direction\treceptor_residue\tligand_atom\tdistance\tposition_status\n")
    for tag, lig, res, amap_ in (('55Z_into_FSHR', M55, F_res, RMAP),
                                 ('O6F_into_LHCGR', MO6, L_res, AMAP)):
        rows = []
        for k, ats in res.items():
            best = (1e9, None, None)
            for a, (x, e) in ats.items():
                if e == 'H':
                    continue
                for ln, (lx, le) in lig.items():
                    dd = np.linalg.norm(x - lx)
                    if dd < best[0]:
                        best = (dd, a, ln)
            if best[0] < 4.5:
                rows.append((best[0], k, best[1], best[2]))
        for dd, k, a, ln in sorted(rows):
            partner = amap_.get(k[0])
            # the partner number returned by amap_ belongs to the OTHER receptor,
            # so it must be looked up in that receptor's ordered residue list
            other = (next((x[1] for x in (F_ord if tag.endswith('LHCGR') else L_ord)
                           if x[0] == partner), None))
            st = 'unaligned' if other is None else ('divergent' if other != k[1] else 'conserved')
            f.write(f"{tag}\t{k[1]}{k[0]}({a})\t{ln}\t{dd:.2f}\t{st}\n")

# --------------------------------------------- conservation vs TM background
tm_pairs = [(i, j) for i, j in PAIRS if L_ord[i][0] >= 340]
tm_ident = sum(1 for i, j in tm_pairs if L_ord[i][1] == F_ord[j][1])
p0 = tm_ident / len(tm_pairs)
n_u = sum(1 for lk, fk, _, _ in union if lk and fk)
k_u = ident_u
# exact binomial, two-sided, via normal-free enumeration
from math import comb
probs = [comb(n_u, i) * p0 ** i * (1 - p0) ** (n_u - i) for i in range(n_u + 1)]
pk = probs[k_u]
pval = sum(p for p in probs if p <= pk * (1 + 1e-12))
whole_ident = sum(1 for i, j in PAIRS if L_ord[i][1] == F_ord[j][1]) / len(PAIRS)
with open(op('conservation_test.txt'), 'w') as f:
    f.write("Is the allosteric contact shell more conserved than the TM background?\n\n")
    f.write(f"whole aligned chain R      : {100*whole_ident:.1f}% identity over {len(PAIRS)} pairs\n")
    f.write(f"TM background (resnum>=340): {100*p0:.1f}% identity over {len(tm_pairs)} pairs\n")
    f.write(f"contact-shell union        : {100*k_u/n_u:.1f}% identity ({k_u}/{n_u})\n\n")
    f.write(f"exact two-sided binomial, H0: shell identity == TM background ({p0:.4f})\n")
    f.write(f"  n = {n_u}, k = {k_u}, p = {pval:.4f}\n")
    f.write(f"  verdict: {'significantly different' if pval < 0.05 else 'NOT significantly different'}\n")

val.append(f"\nunion positions {n_u}: {k_u} identical, {len(diverg)} divergent "
           f"({100*k_u/n_u:.0f}% conserved)")
val.append("divergent set: " + ", ".join(f"{lk[1]}{lk[0]}/{fk[1]}{fk[0]}"
                                          for lk, fk, _, _ in sorted(diverg, key=lambda x: x[0][0])))
open(op('validation.txt'), 'w').write("\n".join(val) + "\n")
print("\n".join(val))
print(f"\nLHCGR contact shell: {len(L_c45)} residues ({len(L_c80)} within 8 A)")
print(f"FSHR  contact shell: {len(F_c45)} residues ({len(F_c80)} within 8 A)")
print(f"\nwrote outputs to {os.path.abspath(OUT)}")
np.save(op('_state_R.npy'), R); np.save(op('_state_t.npy'), t)
