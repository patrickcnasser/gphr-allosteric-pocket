#!/usr/bin/env python3
"""
07_three_way_tshr.py

Extends the pairwise LHCGR/FSHR comparison to a three-way family analysis by adding
human TSHR.

SCOPE AND ITS LIMIT -- read this before citing the output.
The allosteric pocket is DEFINED here by distance from Org 43553 in 7FIH, i.e. in LHCGR.
Asking what residue each paralog carries at each of those positions is therefore a
SEQUENCE question once the alignment is fixed, and that is what this script answers:
it verifies TSHR Leu570 and Tyr667 and reports the three-way divergence set.
It does NOT compute a TSHR contact shell around a TSHR-bound ligand, which needs a
TSHR structure with an allosteric agonist (e.g. 7XW5 with ML-109). Where that would
change a claim, the script says so.

INPUTS   ../inputs/7FIH.pdb                (human LHCGR, structure = numbering authority)
         ../inputs/8I2G.pdb                (human FSHR)
         TSHR sequence embedded below, UniProt P16473, length verified against the
         database's own stated value from two independent endpoints.
OUTPUT   ../results/three_way_divergence.tsv
         ../results/three_way_summary.txt
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gphr_lib import AA3, parse_pdb, sequence, needleman_wunsch, shell, BLOSUM62

HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, '..', 'inputs'); OUT = os.path.join(HERE, '..', 'results')
def op(n): return os.path.join(OUT, n)

TSHR_P16473 = "".join("""
MRPADLLQLVLLLDLPRDLGGMGCSSPPCECHQEEDFRVTCKDIQRIPSLPPSTQTLKLI
ETHLRTIPSHAFSNLPNISRIYVSIDVTLQQLESHSFYNLSKVTHIEIRNTRNLTYIDPD
ALKELPLLKFLGIFNTGLKMFPDLTKVYSTDIFFILEITDNPYMTSIPVNAFQGLCNETL
TLKLYNNGFTSVQGYAFNGTKLDAVYLNKNKYLTVIDKDAFGGVYSGPSLLDVSQTSVTA
LPSKGLEHLKELIARNTWTLKKLPLSLSFLHLTRADLSYPSHCCAFKNQKKIRGILESLM
CNESSMQSLRQRKSVNALNSPLHQEYEENLGDSIVGYKEKSKFQDTHNNAHYYVFFEEQE
DEIIGFGQELKNPQEETLQAFDSHYDYTICGDSEDMVCTPKSDEFNPCEDIMGYKFLRIV
VWFVSLLALLGNVFVLLILLTSHYKLNVPRFLMCNLAFADFCMGMYLLLIASVDLYTHSE
YYNHAIDWQTGPGCNTAGFFTVFASELSVYTLTVITLERWYAITFAMRLDRKIRLRHACA
IMVGGWVCCFLLALLPLVGISSYAKVSICLPMDTETPLALAYIVFVLTLNIVAFVIVCCC
YVKIYITVRNPQYNPGDKDTKIAKRMAVLIFTDFICMAPISFYALSAILNKPLITVSNSK
ILLVLFYPLNSCANPFLYAIFTKAFQRDVFILLSKFGICKRQAQAYRGQRVPPKNSTDIQ
VQKVTHEMRQGLHNMEDVYELIENSHLTPKKQGQISEEYMQTVL
""".split())
TSHR_EXPECTED_LEN = 764

L_res, L_ord, L_het = parse_pdb(os.path.join(IN, '7FIH.pdb'), 'R')
F_res, F_ord, F_het = parse_pdb(os.path.join(IN, '8I2G.pdb'), 'R')
L55 = {n: v for n, v in next(v for k, v in L_het.items() if k[0] == '55Z').items() if v[1] != 'H'}

out = ["Three-way comparison of the glycoprotein-hormone-receptor allosteric pocket", ""]
out.append(f"TSHR P16473 length: {len(TSHR_P16473)} (UniProt states {TSHR_EXPECTED_LEN}) -> "
           f"{'OK' if len(TSHR_P16473) == TSHR_EXPECTED_LEN else 'MISMATCH'}")
assert len(TSHR_P16473) == TSHR_EXPECTED_LEN

sl = sequence(L_ord)
pairs_F = needleman_wunsch(sl, sequence(F_ord), matrix=BLOSUM62, gap_open=-11, gap_extend=-1)
pairs_T = needleman_wunsch(sl, TSHR_P16473, matrix=BLOSUM62, gap_open=-11, gap_extend=-1)
idF = sum(1 for i, j in pairs_F if sl[i] == sequence(F_ord)[j])
idT = sum(1 for i, j in pairs_T if sl[i] == TSHR_P16473[j])
out.append(f"LHCGR vs FSHR : {len(pairs_F)} pairs, {100*idF/len(pairs_F):.1f}% identity")
out.append(f"LHCGR vs TSHR : {len(pairs_T)} pairs, {100*idT/len(pairs_T):.1f}% identity")

FMAP = {L_ord[i][0]: (F_ord[j][0], F_ord[j][1]) for i, j in pairs_F}
TMAP = {L_ord[i][0]: (j + 1, TSHR_P16473[j]) for i, j in pairs_T}

# ---- the validation that matters: do the literature TSHR numbers fall out?
out.append("\nCONTROL -- published TSHR equivalences, recovered without being supplied:")
ok = True
for lnum, expect_num, expect_aa, label in [(515, 570, 'L', 'ECL2 / the LHCGR-unique position'),
                                           (612, 667, 'Y', 'position 7.42'),
                                           (451, 506, 'E', 'the E3.37 pharmacophore anchor')]:
    tn, ta = TMAP.get(lnum, (None, '?'))
    hit = (tn == expect_num and ta == expect_aa)
    ok &= hit
    lname = next(k[1] for k in L_ord if k[0] == lnum)
    out.append(f"  LHCGR {lname}{lnum} -> TSHR {ta}{tn}   literature says {expect_aa}{expect_num}   "
               f"{'MATCH' if hit else 'MISMATCH'}   ({label})")
out.append(f"  control: {'PASS' if ok else 'FAIL'}")

# ---------------------------------------------------------- three-way table
sh80 = shell(L_res, L55, 8.0)
sh45 = shell(L_res, L55, 4.5)
rows = ["LHCGR\tmin_dist_to_Org43553\tshell\tFSHR\tTSHR\tpattern"]
counts = {'all three identical': 0, 'LHCGR unique': 0, 'FSHR unique': 0,
          'TSHR unique': 0, 'all three differ': 0}
lhcgr_unique_contact = []
for k, d in sorted(sh45.items(), key=lambda x: x[1]):
    pass
for k, d in sorted(sh80.items(), key=lambda x: x[1]):
    la = AA3[k[1]]
    fn, fa3 = FMAP.get(k[0], (None, None))
    fa = AA3.get(fa3) if fa3 else None
    tn, ta = TMAP.get(k[0], (None, None))
    if fa is None or ta is None:
        pattern = 'unaligned'
    elif la == fa == ta:
        pattern = 'all three identical'
    elif la != fa and la != ta and fa == ta:
        pattern = 'LHCGR unique'
    elif fa != la and fa != ta and la == ta:
        pattern = 'FSHR unique'
    elif ta != la and ta != fa and la == fa:
        pattern = 'TSHR unique'
    else:
        pattern = 'all three differ'
    counts[pattern] = counts.get(pattern, 0) + 1
    if pattern == 'LHCGR unique' and k in sh45:
        lhcgr_unique_contact.append((k, d, f"{fa3}{fn}", f"{ta}{tn}"))
    rows.append(f"{k[1]}{k[0]}\t{d:.2f}\t{'contact' if k in sh45 else '8A'}\t"
                f"{fa3}{fn}\t{ta}{tn}\t{pattern}")
open(op('three_way_divergence.tsv'), 'w').write("\n".join(rows) + "\n")

out.append(f"\nPocket defined in LHCGR: {len(sh80)} residues within 8.0 A, {len(sh45)} within 4.5 A")
out.append("\nThree-way pattern over the 8 A shell:")
for kk, vv in sorted(counts.items(), key=lambda x: -x[1]):
    out.append(f"  {kk:<22}{vv:>4}")

out.append("\nPositions unique to LHCGR WITHIN THE 4.5 A CONTACT SHELL "
           "(the claim the Phe515 argument rests on):")
if not lhcgr_unique_contact:
    out.append("  none")
for k, d, f_, t_ in lhcgr_unique_contact:
    out.append(f"  {k[1]}{k[0]}  {d:.2f} A to Org 43553   FSHR {f_}   TSHR {t_}")
out.append(f"  count: {len(lhcgr_unique_contact)}")

# and the contact-shell divergence sets, pairwise
dF = [k for k in sh45 if FMAP.get(k[0]) and AA3[FMAP[k[0]][1]] != AA3[k[1]]]
dT = [k for k in sh45 if TMAP.get(k[0]) and TMAP[k[0]][1] != AA3[k[1]]]
out.append(f"\ncontact-shell divergences vs FSHR: {len(dF)} of {len(sh45)}")
out.append("  " + ", ".join(f"{k[1]}{k[0]}/{FMAP[k[0]][1]}{FMAP[k[0]][0]}"
                            for k in sorted(dF, key=lambda x: x[0])))
out.append(f"contact-shell divergences vs TSHR: {len(dT)} of {len(sh45)}")
out.append("  " + ", ".join(f"{k[1]}{k[0]}/{TMAP[k[0]][1]}{TMAP[k[0]][0]}"
                            for k in sorted(dT, key=lambda x: x[0])))

out.append("\nLIMIT OF THIS ANALYSIS")
out.append("  The pocket is defined by distance to Org 43553 in LHCGR, so these are the")
out.append("  paralog residues at LHCGR's pocket positions. Confirming that TSHR's own")
out.append("  allosteric pocket occupies the same positions, and computing a TSHR contact")
out.append("  shell, requires a TSHR structure with an allosteric agonist bound")
out.append("  (7XW5, ML-109). That file was not reachable from this environment.")
open(op('three_way_summary.txt'), 'w').write("\n".join(out) + "\n")
print("\n".join(out))
