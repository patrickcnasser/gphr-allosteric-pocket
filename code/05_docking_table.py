#!/usr/bin/env python3
"""
05_docking_table.py  -- supplementary docking table, symmetry-corrected RMSD.

NOTE: gnina/openbabel do not preserve input atom order. RMSD is computed through
graph isomorphism plus reference automorphisms (rmsd_lib), NOT row-by-row.

INPUTS   ../inputs/xtal_lig.sdf ; ../results/dock/*.sdf and *.log
OUTPUTS  ../results/docking_supplementary.tsv ; ../results/docking_summary.txt
"""
import os, sys, re, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rmsd_lib import read_sdf, symmetry_rmsd, centroid_offset, atom_overlap
HERE=os.path.dirname(os.path.abspath(__file__))
IN=os.path.join(HERE,'..','inputs'); OUT=os.path.join(HERE,'..','results')
def op(n): return os.path.join(OUT,n)
ref=read_sdf(os.path.join(IN,'xtal_lig.sdf'))[0]

def parse_log(fn):
    txt=open(fn).read() if os.path.exists(fn) else ""
    vals=[float(m.group(2)) for m in (re.match(r'\s*(\d+)\s+(-?\d+\.\d+)',l) for l in txt.split('\n')) if m]
    cnn=[float(v) for v in re.findall(r'CNNscore:\s*([\d.]+)',txt)]
    return vals,cnn

rows=["configuration\tpose_rank\tscore\tcnn_score\trmsd_to_crystal_A\tcentroid_offset_A\tatom_overlap_pct"]
summary=["Supplementary docking table -- all configurations, all poses.",
 "RMSD is symmetry-corrected and order-independent (graph isomorphism + automorphisms).",
 "Bar: the top-ranked pose must lie within 2.0 A of the crystallographic pose.",""]
summary.append(f"{'configuration':<22}{'top RMSD':>10}{'best in set':>13}{'native rank':>13}{'verdict':>9}")
for sdf in sorted(glob.glob(op('dock/*.sdf'))):
    cfg=os.path.basename(sdf)[:-4]; poses=read_sdf(sdf)
    if not poses: continue
    sc,cnn=parse_log(sdf[:-4]+'.log'); rl=[]
    for i,p in enumerate(poses,1):
        r,_,ok=symmetry_rmsd(p,ref)
        if not ok: continue
        cen=centroid_offset(p,ref); ov=atom_overlap(p,ref)
        s = f"{sc[i-1]:.2f}" if i <= len(sc) else ""
        c = f"{cnn[i-1]:.4f}" if i <= len(cnn) else ""
        rows.append(f"{cfg}\t{i}\t{s}\t{c}\t{r:.2f}\t{cen:.2f}\t{100*ov:.0f}")
        rl.append(r)
    if not rl: continue
    best=min(rl); nrank=rl.index(best)+1
    summary.append(f"{cfg:<22}{rl[0]:>9.2f}A{best:>12.2f}A{nrank:>13}{'PASS' if rl[0]<=2.0 else 'FAIL':>9}")
open(op('docking_supplementary.tsv'),'w').write("\n".join(rows)+"\n")
summary+=["","'best in set' is the closest pose the search produced; 'native rank' is its rank",
 "under that configuration's own scoring function."]
open(op('docking_summary.txt'),'w').write("\n".join(summary)+"\n")
print("\n".join(summary))
