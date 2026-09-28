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

set cartoon_transparency, 0.8
zoom nativeLig or o6f_mapped, 6
print("LEFT: native FSHR/Cpd-21f. RIGHT: Cpd-21f mapped into LHCGR.")
print("Adjust, then: ray 2400, 1400; png fig2.png, dpi=300")
