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
cmd.load(P("results", "O6F_in_LHCGR.pdb"), "o6f_mapped")
cmd.remove("hydro")
python end

# FIGURE 2 -- reciprocal mapping: Tyr612 clash against the native His615 contact
# INPUTS: 7FIH.pdb, 8I2G.pdb, results/O6F_in_LHCGR.pdb
# Two objects side by side: LEFT  = native FSHR + Cpd-21f (His615, 3.69 A)
#                           RIGHT = LHCGR + mapped Cpd-21f (Tyr612, 2.31 A)

create nativeFSHR, fs and chain R and polymer
create nativeLig,  fs and resn O6F
create mappedLHCGR, lh and chain R and polymer
delete lh
delete fs

# separate the two panels along x so both are visible in one image
translate [25, 0, 0], object=mappedLHCGR, camera=0
translate [25, 0, 0], object=o6f_mapped,  camera=0

hide everything
select pocket_fs, byres (nativeFSHR within 8 of nativeLig)
select pocket_lh, byres (mappedLHCGR within 8 of o6f_mapped)
show cartoon, pocket_fs or pocket_lh
color palecyan, nativeFSHR
color grey80,  mappedLHCGR
show sticks, nativeLig or o6f_mapped
color salmon, nativeLig and elem C
color salmon, o6f_mapped and elem C
util.cnc("nativeLig or o6f_mapped")
hide everything, hydro

select his615, nativeFSHR and resi 615
select tyr612, mappedLHCGR and resi 612
select ser604, mappedLHCGR and resi 604
show sticks, his615 or tyr612 or ser604
color marine, his615 and elem C
color red, tyr612 and elem C
color orange, ser604 and elem C
util.cnc("his615 or tyr612 or ser604")

# the two measurements quoted in the text
distance d_native, nativeFSHR and resi 615 and name CE1, nativeLig and name C01
distance d_clash,  mappedLHCGR and resi 612 and name CE1, o6f_mapped and name C01
color black, d_native
color red, d_clash
hide labels, d_native
hide labels, d_clash

python
from pymol import cmd
def below(sel, dy=16):
    c = cmd.centerofmass(sel)
    return [c[0], c[1] - dy, c[2]]
p1 = below("nativeLig")
p2 = below("o6f_mapped")
shared_y = min(p1[1], p2[1])
p1[1] = shared_y
p2[1] = shared_y
cmd.pseudoatom("lab_fs", pos=p1, label="His615 (FSHR)  3.69 \xc5")
cmd.pseudoatom("lab_lh", pos=p2, label="Tyr612 (LHCGR)  2.31 \xc5")
python end
hide everything, lab_fs or lab_lh
show label, lab_fs or lab_lh

# ---- SELF-CHECK: the two quoted distances must print 3.69 and 2.31 ----
python
from pymol import cmd
a = cmd.get_distance("nativeFSHR and resi 615 and name CE1", "nativeLig and name C01")
b = cmd.get_distance("mappedLHCGR and resi 612 and name CE1", "o6f_mapped and name C01")
print("SELF-CHECK native His615...C01 = %.2f A (expected 3.69)" % a)
print("SELF-CHECK mapped Tyr612...C01 = %.2f A (expected 2.31)" % b)
print("FIGURE OK" if abs(a-3.69) < 0.05 and abs(b-2.31) < 0.05 else "CHECK INPUTS")
python end

set ray_shadows, 0
set ray_opaque_background, 1
set antialias, 2
set cartoon_transparency, 0.50
set label_size, 16
set label_color, black
set dash_width, 2.5
set dash_gap, 0.35
set depth_cue, 0
bg_color white

zoom nativeLig or o6f_mapped, 6
print("LEFT: native FSHR/Cpd-21f. RIGHT: Cpd-21f mapped into LHCGR.")
print("Adjust, then: ray 2400, 1700; png fig2.png, dpi=300")
