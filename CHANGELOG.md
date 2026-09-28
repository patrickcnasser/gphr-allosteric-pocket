# Changelog

## 1.1 — 28 September 2026

Three fixes and one addition. No distance, RMSD, alignment or superposition value from
version 1.0 changed; `results/validation.txt` passes identically.

### Fixed: `position_status` in `results/mapping_contacts.tsv`

`code/01_pocket_and_superposition.py` looked the mapped partner residue up in the wrong
receptor's residue list, so conserved positions were reported as divergent. Phe353 in FSHR,
the paralog of LHCGR Phe350, was labelled `divergent`. Corrected counts across both mapping
directions are **34 conserved and 12 divergent**, against 5 and 41 before. Every distance and
atom name in that file is unchanged — only the status column was wrong, and no value reported
in the manuscript depended on it.

### Fixed: non-deterministic ring labels

`code/02_fragment_accounting.py` numbered perceived ring systems in set-iteration order, so
`ring2_6mem_CNO` and `ring3_6mem_CNO` could exchange labels between runs of the same code on
the same input. Ring systems are now sorted by their atom names before numbering. Two
consecutive full runs of `run_all.sh` are now byte-identical across every output file.

### Added: `code/08_manuscript_fragment_table.py`

The manuscript's per-fragment burial table was not regenerable from version 1.0. Script 02
perceives fragments from connectivity alone and reports a coarser partition — it keeps the
thienopyrimidine as one nine-atom ring system and groups the tert-butyl with its carboxamide.
The table in the paper uses the finer, chemically named split, which this script now defines
explicitly and computes in both frames.

It reproduces all eight published rows exactly: 2-methylthio 14.5, tert-butyl 11.0,
pyrimidine 10.2, morpholine 8.8, thiophene 8.2, anilide linker 7.8, pendant phenyl 6.5,
carboxamide 6.3, plus the amino group at 6.0. The LHCGR and FSHR-mapped columns are printed
side by side and labelled, because an earlier draft of the manuscript reported the
FSHR-mapped column under an LHCGR heading.

Connectivity from the deposited coordinates gives a thiophene (C8, C9, C10, C15, S16) fused
to a pyrimidine (C10, C11, C13, C15, N12, N14) sharing C10 and C15. To keep the partition
disjoint and complete over all 35 heavy atoms, the two fused carbons are counted with the
thiophene, so the pyrimidine row covers that ring's four exclusive atoms.

### Figures: opaque backgrounds kept

The three rendered figures (`figures/fig1.png`, `fig2.png`, `fig3.png`) are the opaque-background
renders already published in Zenodo v2 and on GitHub. `04_make_figure_scripts.py` now writes
`set ray_opaque_background, 1`, so regenerated `.pml` scripts match them. Transparent
backgrounds had flattened to black in some journal production pipelines.

### Restored from the GitHub repository

`LICENSE` (MIT), the full `README.md` and the resolved DOI in `PROVENANCE.md` are carried
forward. The citation DOI is now the Zenodo concept DOI, 10.5281/zenodo.22969109, which
always resolves to the newest version.

### Platform note

On Linux x86-64, `run_all.sh` regenerates every output byte-identically. On macOS arm64 three
files differ by floating-point rounding only (at most 4e-13 Å): `results/_state_R.npy`,
`results/_state_t.npy` and one line in `results/validation.txt`.

### Known limitations carried forward

`code/legacy_session_scripts/` is kept for provenance and is **not** part of the reproducible
path. `frag2.py` and `scan.py` both load `/home/claude/super/data.pkl`, which was never part
of this package, and cannot run. Everything either script established is now reproduced by
`08_manuscript_fragment_table.py` and by `rederive_ortholog_control.py` in the project record.
