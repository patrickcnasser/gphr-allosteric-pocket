# ---------------------------------------------------------------------------
# EDIT THIS LINE ONLY: directory holding 7FIH.pdb, 8I2G.pdb and the results files
# ---------------------------------------------------------------------------
python
import os
from pymol import cmd
wd = "/Users/patnasrick/Desktop/drugs/paper"
def P(*a): return os.path.join(wd, *a)
cmd.load(P("7FIH.pdb"), "lh")
cmd.load(P("results", "decoy_best.pdb"), "decoy")
cmd.remove("hydro")
python end

# FIGURE 3 -- the flipped decoy overlaid on the crystallographic pose
# INPUTS: 7FIH.pdb, results/decoy_best.pdb
# Both occupy the same site; the decoy is the molecule reversed end-for-end.
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
hide labels, swap1
hide labels, swap2

# the pocket residues that see one orientation but not the other are listed in
# results/flipped_pose.txt -- add them here once you have read that file, e.g.:
# select discriminating, pocket and resi 350+612
# show sticks, discriminating

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

orient native
zoom native or decoy, 4

# ---- SELF-CHECK: the decoy's tert-butyl should sit where the native morpholine is ----
python
from pymol import cmd
cmd.pseudoatom("p_dec_tbu", selection="decoy and name C1+C2+C3+C4+N5+C6+O7")
cmd.pseudoatom("p_nat_mor", selection="native and name N29+C30+C31+O32+C33+C34")
d = cmd.get_distance("p_dec_tbu", "p_nat_mor")
print("SELF-CHECK decoy tert-butyl vs native morpholine centroid: %.2f A (expected 0.5)" % d)
print("HEAD-TO-TAIL EXCHANGE CONFIRMED" if d < 1.5 else "unexpected -- check decoy_best.pdb")
cmd.delete("p_dec_tbu")
cmd.delete("p_nat_mor")
python end
print("Green = crystallographic pose, magenta = decoy. Spheres mark the two ends.")
print("Adjust, then: ray 2400, 1800; png fig3.png, dpi=300")
