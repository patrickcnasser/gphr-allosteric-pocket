# gphr-allosteric-pocket

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22969109.svg)](https://doi.org/10.5281/zenodo.22969109)

Structural comparison of the LHCGR and FSHR allosteric binding pockets, using cryo-EM structures of LHCGR with Org 43553 (PDB [7FIH](https://www.rcsb.org/structure/7FIH)) and FSHR with Cpd-21f (PDB [8I2G](https://www.rcsb.org/structure/8I2G)).

This package accompanies:

> Nasser, P. _Asymmetric selectivity at the glycoprotein hormone receptor allosteric site._ Preprint, 2026. \[link forthcoming\]

## What it does

The numbered scripts in `code/` reproduce every structural analysis in the manuscript:

- Sequence alignment (Needleman-Wunsch, BLOSUM62, affine gaps)
- Transmembrane C-alpha superposition (Kabsch algorithm, 187 pairs, 0.64 A RMSD)
- Pocket definition and divergence mapping (7 of 27 contact positions diverge)
- Reciprocal ligand mapping and contact recomputation
- Fragment-level burial and contact accounting
- Three-way comparison with TSHR (sequence-based on the structurally defined pocket)
- Redocking with symmetry-corrected RMSD across eight scoring configurations
- Sampling-independent rescoring control
- Head-to-tail decoy characterisation
- Per-fragment burial table exactly as reported in the manuscript (`08_manuscript_fragment_table.py`)

## Reproducing the results

```bash
./run_all.sh
```

This regenerates every number in the manuscript from the two deposited PDB files (~3 min). Docking requires a [gnina](https://github.com/gnina/gnina) binary and takes ~35 min on 2 cores:

```bash
GNINA=/path/to/gnina ./run_all.sh
```

If gnina is not available, `run_all.sh` re-analyses the shipped docking outputs instead.

Outputs land in `results/` and `figures/`. Nothing else is written.

## Verifying the package

`MANIFEST.md5` lists the md5 of every file. From inside the package directory:

```bash
md5sum -c MANIFEST.md5      # Linux
```

Run it before and after `./run_all.sh`. On Linux x86-64 every file verifies unchanged after a full run. On other platforms (e.g. macOS arm64) three files differ only by floating-point rounding at the 1e-13 level: `results/_state_R.npy`, `results/_state_t.npy` and the round-trip error printed in `results/validation.txt`. No reported value is affected.

See `CHANGELOG.md` for what changed between versions.

## Requirements

| Dependency | Details |
|---|---|
| Python | 3.9 or later |
| Packages | **numpy only** — no scipy, no RDKit, no Biopython |
| Docking (optional) | [gnina](https://github.com/gnina/gnina) v1.1 CPU binary; also provides smina's classical scoring functions, so all eight configurations come from one binary |
| Figures (optional) | [PyMOL](https://pymol.org/) open-source build; the `.pml` scripts use no plugins |

## RMSD computation

RMSD between docked poses and the crystallographic pose is computed by graph isomorphism with enumeration of the ligand's automorphisms (12 for Org 43553), because gnina/Open Babel do not preserve input atom order. Row-order comparison of the same files inflates a 1 A deviation to roughly 9 A. See `code/rmsd_lib.py` for the implementation and its self-tests.

## Layout

```
inputs/       7FIH.pdb, 8I2G.pdb, receptor.pdb, xtal_lig.sdf
code/         analysis scripts (numbered in run order) and shared libraries
results/      all outputs — alignments, contact tables, matrices, docking poses, logs
figures/      PyMOL .pml scripts and rendered PNGs
run_all.sh    regenerates everything from the two PDB files
MANIFEST.md5  md5 of every file in the package
CHANGELOG.md  version history
```

## License

[MIT](LICENSE)
