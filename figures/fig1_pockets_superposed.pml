# ---------------------------------------------------------------------------
# EDIT THIS LINE ONLY: directory holding 7FIH.pdb, 8I2G.pdb and the results files
# ---------------------------------------------------------------------------
python
import os
from pymol import cmd
wd = "/Users/patnasrick/Desktop/drugs/paper"
def P(*a): return os.path.join(wd, *a)
cmd.load(P("7FIH.pdb"), "lh")
cmd.load(P("8I2G.pdb"), "fs")
cmd.remove("hydro")
python end

# FIGURE 1 -- the two allosteric pockets superposed, seven divergent positions labelled
# INPUTS: 7FIH.pdb, 8I2G.pdb   (in the directory set above)
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
cmd.transform_selection("lhR", M, homogenous=1)
cmd.transform_selection("lig55Z", M, homogenous=1)
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
util.cnc("div_lh or div_fs")

python
from pymol import cmd
lh_resi = [519, 528, 531, 593, 604, 612]
fs_resi = [522, 531, 534, 596, 607, 615]
aa = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
      'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
      'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}
for lr, fr in zip(lh_resi, fs_resi):
    ln = aa.get(cmd.get_model("lhR and resi %d and name CA" % lr).atom[0].resn, '?')
    fn = aa.get(cmd.get_model("fsR and resi %d and name CA" % fr).atom[0].resn, '?')
    cmd.label("lhR and resi %d and name CA" % lr, '"%s%d/%s%d"' % (ln, lr, fn, fr))
cmd.set("label_position", [2, 2, 0])
cmd.set("label_position", [4, -3, 0], "lhR and resi 528 and name CA")
python end

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
python
from pymol import cmd
cmd.pseudoatom("tb_55Z", selection="lig55Z and name C1+C2+C3+C4")
cmd.pseudoatom("tb_O6F", selection="ligO6F and name C30+C31+C32+C33")
d = cmd.get_distance("tb_55Z", "tb_O6F")
print("SELF-CHECK tert-butyl centroid separation: %.2f A (expected 0.63)" % d)
print("TRANSFORM OK" if abs(d - 0.63) < 0.3 else "TRANSFORM WRONG -- try transpose=1")
cmd.delete("tb_55Z")
cmd.delete("tb_O6F")
python end
print("Adjust the view, then: ray 2400, 1800; png fig1.png, dpi=300")
