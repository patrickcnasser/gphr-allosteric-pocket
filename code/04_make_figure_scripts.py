#!/usr/bin/env python3
"""
04_make_figure_scripts.py

Emits the PyMOL transformation matrix in PyMOL's own 16-element convention and
writes three ready-to-run .pml scripts.

INPUTS   ../results/_state_R.npy , ../results/_state_t.npy  (from script 01)
OUTPUTS  ../results/pymol_matrix.txt
         ../figures/fig1_pockets_superposed.pml
         ../figures/fig2_reciprocal_clash.pml
         ../figures/fig3_flipped_pose.pml
"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results'); FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(FIG, exist_ok=True)
R = np.load(os.path.join(RES, '_state_R.npy')); t = np.load(os.path.join(RES, '_state_t.npy'))

# PyMOL transform_selection: 16 elements, row-major 4x4, translation in col 4.
M = [R[0, 0], R[0, 1], R[0, 2], t[0],
     R[1, 0], R[1, 1], R[1, 2], t[1],
     R[2, 0], R[2, 1], R[2, 2], t[2],
     0.0, 0.0, 0.0, 1.0]
Rt = R.T; ti = -R.T @ t
Minv = [Rt[0, 0], Rt[0, 1], Rt[0, 2], ti[0],
        Rt[1, 0], Rt[1, 1], Rt[1, 2], ti[1],
        Rt[2, 0], Rt[2, 1], Rt[2, 2], ti[2],
        0.0, 0.0, 0.0, 1.0]
fmt = lambda m: "[" + ", ".join(f"{v:.8f}" for v in m) + "]"
with open(os.path.join(RES, 'pymol_matrix.txt'), 'w') as f:
    f.write("LHCGR(7FIH) -> FSHR(8I2G), PyMOL transform_selection convention\n")
    f.write("forward (moves LHCGR onto FSHR):\n  " + fmt(M) + "\n")
    f.write("inverse (moves FSHR onto LHCGR):\n  " + fmt(Minv) + "\n")

HEADER = """# ---------------------------------------------------------------------------
# EDIT THIS LINE ONLY: directory holding 7FIH.pdb, 8I2G.pdb and the results files
set_wd = "/CHANGE/ME/gphr_paper"
# ---------------------------------------------------------------------------
python
import os
from pymol import cmd
wd = set_wd
def P(*a): return os.path.join(wd, *a)
python end
"""

COMMON = """
set ray_shadows, 0
set ray_opaque_background, 0
set antialias, 2
set cartoon_transparency, 0.65
set label_size, 16
set label_color, black
set dash_width, 2.5
set dash_gap, 0.35
set depth_cue, 0
bg_color white
"""

# ---------------------------------------------------------------- figure 1
fig1 = HEADER + f"""
# FIGURE 1 -- the two allosteric pockets superposed, seven divergent positions labelled
# INPUTS: 7FIH.pdb, 8I2G.pdb   (in the directory set above)

load P("7FIH.pdb"), lh
load P("8I2G.pdb"), fs
create lhR, lh and chain R and polymer
create fsR, fs and chain R and polymer
create lig55Z, lh and resn 55Z
create ligO6F, fs and resn O6F
delete lh
delete fs

# apply the transmembrane superposition computed in script 01
# (moves LHCGR and its ligand into the FSHR frame)
python
from pymol import cmd
M = {fmt(M)}
cmd.transform_selection("lhR", M, homogenous=0)
cmd.transform_selection("lig55Z", M, homogenous=0)
python end

hide everything
show cartoon, lhR or fsR
color grey80, lhR
color palecyan, fsR
show sticks, lig55Z or ligO6F
color yellow, lig55Z and elem C
color salmon, ligO6F and elem C
util.cnc("lig55Z or ligO6F")

# the seven divergent contact positions (LHCGR numbering / FSHR numbering)
select div_lh, lhR and resi 515+519+528+531+593+604+612 and not (name C+N+O)
select div_fs, fsR and resi 518+522+531+534+596+607+615 and not (name C+N+O)
show sticks, div_lh or div_fs
color orange, div_lh and elem C
color marine, div_fs and elem C

label lhR and resi 515+519+528+531+593+604+612 and name CA, "%s%s" % (resn, resi)
label fsR and resi 518+522+531+534+596+607+615 and name CA, "%s%s" % (resn, resi)
{COMMON}
orient lig55Z or ligO6F
zoom lig55Z or ligO6F, 5

# ---- SELF-CHECK: is the transform applied with the right matrix convention? ----
# After a correct transform the mapped Org 43553 tert-butyl sits 0.63 A from
# Cpd-21f's own tert-butyl. If this prints something far from 0.63, PyMOL has
# interpreted the matrix transposed: re-run with transpose=1 in transform_selection.
pseudoatom tb_55Z, selection="lig55Z and name C1+C2+C3+C4"
pseudoatom tb_O6F, selection="ligO6F and name C30+C31+C32+C33"
python
from pymol import cmd
d = cmd.get_distance("tb_55Z", "tb_O6F")
print("SELF-CHECK tert-butyl centroid separation: %.2f A (expected 0.63)" % d)
print("TRANSFORM OK" if abs(d - 0.63) < 0.3 else "TRANSFORM WRONG -- try transpose=1")
python end
delete tb_55Z
delete tb_O6F
print("Adjust the view, then: ray 2400, 1800; png fig1.png, dpi=300")
"""

# ---------------------------------------------------------------- figure 2
fig2 = HEADER + f"""
# FIGURE 2 -- reciprocal mapping: Tyr612 clash against the native His615 contact
# INPUTS: 7FIH.pdb, 8I2G.pdb, results/O6F_in_LHCGR.pdb
# Two objects side by side: LEFT  = native FSHR + Cpd-21f (His615, 3.69 A)
#                           RIGHT = LHCGR + mapped Cpd-21f (Tyr612, 2.31 A)

load P("7FIH.pdb"), lh
load P("8I2G.pdb"), fs
load P("results", "O6F_in_LHCGR.pdb"), o6f_mapped

create nativeFSHR, fs and chain R and polymer
create nativeLig,  fs and resn O6F
create mappedLHCGR, lh and chain R and polymer
delete lh
delete fs

# separate the two panels along x so both are visible in one image
translate [38, 0, 0], object=mappedLHCGR, camera=0
translate [38, 0, 0], object=o6f_mapped,  camera=0

hide everything
show cartoon, nativeFSHR or mappedLHCGR
color palecyan, nativeFSHR
color grey80,  mappedLHCGR
show sticks, nativeLig or o6f_mapped
color salmon, nativeLig and elem C
color salmon, o6f_mapped and elem C
util.cnc("nativeLig or o6f_mapped")

select his615, nativeFSHR and resi 615
select tyr612, mappedLHCGR and resi 612
select ser604, mappedLHCGR and resi 604
show sticks, his615 or tyr612 or ser604
color marine, his615 and elem C
color red, tyr612 and elem C
color orange, ser604 and elem C

# the two measurements quoted in the text
distance d_native, nativeFSHR and resi 615 and name CE1, nativeLig and name C01
distance d_clash,  mappedLHCGR and resi 612 and name CE1, o6f_mapped and name C01
distance d_ser,    mappedLHCGR and resi 604 and name OG,  o6f_mapped and name C27
color black, d_native
color red, d_clash
color orange, d_ser
set label_distance_digits, 2

label his615 and name CA, "His615 (FSHR)  3.69 A"
label tyr612 and name CA, "Tyr612 (LHCGR)  2.31 A"

# ---- SELF-CHECK: the two quoted distances must print 3.69 and 2.31 ----
python
from pymol import cmd
a = cmd.get_distance("nativeFSHR and resi 615 and name CE1", "nativeLig and name C01")
b = cmd.get_distance("mappedLHCGR and resi 612 and name CE1", "o6f_mapped and name C01")
print("SELF-CHECK native His615...C01 = %.2f A (expected 3.69)" % a)
print("SELF-CHECK mapped Tyr612...C01 = %.2f A (expected 2.31)" % b)
print("FIGURE OK" if abs(a-3.69) < 0.05 and abs(b-2.31) < 0.05 else "CHECK INPUTS")
python end
{COMMON}
set cartoon_transparency, 0.8
zoom nativeLig or o6f_mapped, 6
print("LEFT: native FSHR/Cpd-21f. RIGHT: Cpd-21f mapped into LHCGR.")
print("Adjust, then: ray 2400, 1400; png fig2.png, dpi=300")
"""

# ---------------------------------------------------------------- figure 3
fig3 = HEADER + """
# FIGURE 3 -- the flipped decoy overlaid on the crystallographic pose
# INPUTS: 7FIH.pdb, results/decoy_best.pdb
# Both occupy the same site; the decoy is the molecule reversed end-for-end.

load P("7FIH.pdb"), lh
load P("results", "decoy_best.pdb"), decoy
create pocket, lh and chain R and polymer
create native, lh and resn 55Z
delete lh

hide everything
show cartoon, pocket
color grey90, pocket
set cartoon_transparency, 0.75

show sticks, native or decoy
color forest, native and elem C
color magenta, decoy and elem C
util.cnc("native or decoy")

# the two ends of the molecule, to make the reversal legible
select tail_nat, native and name N29+C30+C31+O32+C33+C34
select tail_dec, decoy  and name N29+C30+C31+O32+C33+C34
select head_nat, native and name C1+C2+C3+C4+N5+C6+O7
select head_dec, decoy  and name C1+C2+C3+C4+N5+C6+O7
show spheres, tail_nat or tail_dec or head_nat or head_dec
set sphere_scale, 0.30
color green,   tail_nat
color palegreen, head_nat
color deeppurple, tail_dec
color pink,   head_dec

distance swap1, tail_nat, tail_dec, mode=4
distance swap2, head_nat, head_dec, mode=4
color black, swap1
color black, swap2

# the pocket residues that see one orientation but not the other are listed in
# results/flipped_pose.txt -- add them here once you have read that file, e.g.:
# select discriminating, pocket and resi 350+612
# show sticks, discriminating
"""
fig3 += COMMON + """
orient native
zoom native or decoy, 4

# ---- SELF-CHECK: the decoy's tert-butyl should sit where the native morpholine is ----
pseudoatom p_dec_tbu, selection="decoy and name C1+C2+C3+C4+N5+C6+O7"
pseudoatom p_nat_mor, selection="native and name N29+C30+C31+O32+C33+C34"
python
from pymol import cmd
d = cmd.get_distance("p_dec_tbu", "p_nat_mor")
print("SELF-CHECK decoy tert-butyl vs native morpholine centroid: %.2f A (expected 0.5)" % d)
print("HEAD-TO-TAIL EXCHANGE CONFIRMED" if d < 1.5 else "unexpected -- check decoy_best.pdb")
python end
delete p_dec_tbu
delete p_nat_mor
print("Green = crystallographic pose, magenta = decoy. Spheres mark the two ends.")
print("Adjust, then: ray 2400, 1800; png fig3.png, dpi=300")
"""

for name, body in [('fig1_pockets_superposed.pml', fig1),
                   ('fig2_reciprocal_clash.pml', fig2),
                   ('fig3_flipped_pose.pml', fig3)]:
    open(os.path.join(FIG, name), 'w').write(body)
    print("wrote", os.path.join('figures', name))
print("wrote results/pymol_matrix.txt")
