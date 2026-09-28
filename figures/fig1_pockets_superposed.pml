# ---------------------------------------------------------------------------
# EDIT THIS LINE ONLY: directory holding 7FIH.pdb, 8I2G.pdb and the results files
set_wd = "/CHANGE/ME/gphr_paper"
# ---------------------------------------------------------------------------
python
import os
from pymol import cmd
wd = set_wd
def P(*a): return os.path.join(wd, *a)
python end

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
M = [0.99611447, -0.08754067, -0.00962212, 12.39545145, 0.08779350, 0.99568417, 0.03008888, -24.81197621, 0.00694659, -0.03081672, 0.99950091, -0.86123302, 0.00000000, 0.00000000, 0.00000000, 1.00000000]
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

set ray_shadows, 0
set ray_opaque_background, 1
set antialias, 2
set cartoon_transparency, 0.65
set label_size, 16
set label_color, black
set dash_width, 2.5
set dash_gap, 0.35
set depth_cue, 0
bg_color white

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
