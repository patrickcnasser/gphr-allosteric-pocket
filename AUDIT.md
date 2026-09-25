# Manuscript audit — what reproduces, what does not, what is wrong

Read section 1 first. It changes the paper.

---

## 1. The methodological claim is wrong, and it was my error

**Manuscript, "The pocket is orientationally ambiguous to every scoring function tested":**
> "Seven scoring methods, six classical and one learned, fail the same control in the same way."
> "Two decoys at approximately 9 Å outscore the crystallographic pose."
> "No affinity estimate, selectivity ratio, or pose-based ranking for this scaffold should be regarded as reliable."

**This is false, and the numbers behind it came from me.** The cause:

> **gnina/openbabel do not preserve input atom order when writing poses.** Comparing pose
> coordinates to the reference row-by-row compares the wrong atoms. For this ligand that
> turns a true ~1 Å RMSD into ~9 Å.

Caught by an internal contradiction: gnina's own log reported that local minimisation of
the crystallographic pose moved it **1.32 Å**, while my row-order comparison of the same
output file said **8.84 Å**. The heavy-atom element sequences confirm the reordering
(`CNCCCCCCCCCCCCNCOCC…` in, `CNCCCCCNSNCCCCCCCNO…` out).

`code/rmsd_lib.py` now computes RMSD through graph isomorphism plus enumeration of the
reference molecule's automorphisms (12 for Org 43553). Controls: self-RMSD 0.0000 Å;
minimised crystal pose 1.32 Å, matching gnina's own figure exactly.

### Corrected results

Redocking, 8 configurations, bar = top-ranked pose within 2.0 Å:

| configuration | top pose | best pose in ensemble | rank of best | verdict |
|---|---|---|---|---|
| ad4_scoring rigid | 10.42 Å | 1.03 Å | 7 | FAIL |
| dkoes_fast rigid | 8.19 Å | 2.17 Å | 6 | FAIL |
| dkoes_scoring rigid | 6.59 Å | 3.47 Å | 3 | FAIL |
| **gnina CNN rescore** | **0.99 Å** | 0.99 Å | 1 | **PASS** |
| vina flex-Lys595 | 10.55 Å | 0.99 Å | 4 | FAIL |
| vina rigid | 10.58 Å | 0.99 Å | 2 | FAIL |
| **vinardo flex-Lys595** | **0.93 Å** | 0.89 Å | 3 | **PASS** |
| vinardo rigid | 3.76 Å | 0.82 Å | 2 | FAIL |

Sampling-independent control (ensemble = minimised crystal pose + 9 docked poses):

| scoring function | top pose RMSD | verdict |
|---|---|---|
| vina | 10.58 Å | FAIL |
| **vinardo** | **0.99 Å** | **PASS** |
| ad4_scoring | 3.77 Å | FAIL |
| dkoes_scoring | 3.77 Å | FAIL |
| dkoes_fast | 5.47 Å | FAIL |
| **gnina CNN** | **0.99 Å** | **PASS** |

**What the section can now honestly say.** Vina, AD4 and the dkoes functions fail badly
and reproducibly. **Vinardo and the gnina CNN recover the crystallographic pose.** The
blanket claim must go, and so must the conclusion that no pose-based ranking at this site
is reliable — the CNN result is the direct counterexample. What survives is narrower and
still worth reporting: the site defeats the most widely used scoring function (Vina) and
several others, while a reweighted classical function and a learned one succeed.

**Also note:** the manuscript's "failure is in ranking, not sampling" **is** supported.
Best-in-ensemble poses are 0.82–3.47 Å across configurations; the search finds the answer.

---

## 2. Reproduced exactly

Re-derived independently from the two PDB files, with code that shares no state with the
session that produced the manuscript.

| manuscript claim | independent result |
|---|---|
| 5 published contact distances, Δ = 0.00 | max deviation **0.0039 Å** |
| contact shell 23 in LHCGR, 24 in FSHR | **23 / 24** |
| 27 union positions, 20 identical, 7 divergent, 74% | **27 / 20 / 7 / 74.1%** |
| divergent set F515/L518, V519/I522, I528/V531, I531/L534, A593/S596, S604/A607, Y612/H615 | **identical set** |
| all 14 distances in that table | **all 14 reproduce to 0.01 Å** |
| 187 Cα pairs, 0.64 Å, 0 pruned | **187 / 0.64 / 0** |
| round-trip error 0.0000 Å | **1.8 × 10⁻¹³ Å** |
| mapped *tert*-butyl 0.63 Å from Cpd-21f's | **0.63 Å** |
| Tyr612 CE1···C01 2.31 Å vs His615 3.69 Å | **2.31 / 3.69 Å** |
| Ser604 OG second-worst at 2.77 Å | **2.77 Å** |
| morpholine contacts 5 residues, all conserved | **Ala349, Phe350, Lys595, Pro597, Val601 — all identical** |
| fragment burial 11.0 / 9.8 / 9.0 / 8.0 | **11.0 / 9.8 / 9.0 / 8.0** |
| backbone 0.48 Å at 7.42 | **0.482 Å** |
| ortholog control: 99.8%, one mismatch at 277 | **99.8%, S277I only** |
| Phe/Leu split across 10 orthologs | **reproduced, all 10** |

The fragment partition was re-derived by automatic ring perception rather than by hand
and returns the same seven fragments.

---

## 3. Errors and inconsistencies short of section 1

**3a. Methods contradicts Results on the alignment.** Methods states BLOSUM62 with gap
open −11, extend −1. That scoring gives **61.2%** identity, not the 60.7% reported — 60.7%
is what an identity matrix with gap −8/−1 gives. The two schemes produce **different
pairings**, though both recover Y612→H615 and F515→L518 and both give the same divergent
set. Pick one, state it, and use its number.

**3b. "74% conserved" is not a finding.** The paper never tests it. Against the correct
background it is unremarkable:

| | identity |
|---|---|
| whole aligned chain R | 61.2% |
| **transmembrane background (resnum ≥ 340)** | **70.0%** |
| contact-shell union | 74.1% (20/27) |

Exact two-sided binomial, n = 27, k = 20, p₀ = 0.700: **p = 0.83.** The allosteric pocket
is **not** more conserved than the transmembrane domain generally. The 60.7% figure the
paper implicitly compares against is whole-receptor identity including the leucine-rich
extracellular domain, which is the wrong comparator.

**3c. The 8 Å shell size for FSHR is not stated.** It is 66 (LHCGR 63).

---

## 4. New analyses

### 4a. The flipped pose, characterised

Selected decoy: 9.55 Å RMSD, centroid offset 0.44 Å, 97% of atoms within 2 Å of some
crystal atom. It is a **reciprocal head-to-tail exchange**, and the fragment mapping is
near-exact:

| fragment | displacement | lands in the subpocket natively occupied by |
|---|---|---|
| *tert*-butyl + carboxamide | 14.5 Å | **morpholine (0.5 Å)** |
| morpholine | 15.2 Å | **tert-butyl + carboxamide (1.7 Å)** |
| pendant phenyl | 4.4 Å | thienopyrimidine core (1.0 Å) |
| thienopyrimidine core | 5.0 Å | anilide linker (1.8 Å) |
| 5-amino | 7.2 Å | 2-methylthio (1.8 Å) |

**19 of 23** high-RMSD, same-centroid decoys across all eight configurations show this
reciprocal exchange (83%).

**The pseudo-symmetry belongs to the pocket, not the ligand.**
- Ligand self-similarity under a 180° rotation about any principal axis: best **1.83 Å** —
  it is *not* pseudo-symmetric (principal moments 1.00 / 0.14 / 0.02, i.e. a rod).
- Pocket: the two ends of the ligand axis are lined by **14 and 12 residues, 64% and 75%
  hydrophobic, sharing only one residue.** Two chemically near-equivalent greasy lobes.

**What the failing functions respond to.** The decoy is *better* on every term they weight:

| descriptor | native | decoy | Δ |
|---|---|---|---|
| protein neighbours < 5 Å | 303 | 322 | **+19** |
| protein neighbours < 4 Å | 52 | 76 | **+24** |
| apolar C···C contacts < 4.5 Å | 60 | 96 | **+36** |
| polar N/O···N/O pairs < 3.5 Å | 3 | 1 | −2 |
| steric clashes < 2.6 Å | 0 | 0 | 0 |

And: **no residue contacts one orientation and not the other** (4.0 Å / 5.0 Å thresholds,
zero in both directions). There is nothing orientation-specific for a contact-counting
function to detect. The decoy is not a scoring artefact — it is a better pose *by those
functions' own criteria*. That is the mechanistic result, and it is consistent with
Vinardo and the CNN succeeding, since both depart from raw contact counting.

### 4b. Cpd-21f fragment accounting — **null result**

Automatic ring perception gives five fragments for Cpd-21f. Does the FSHR-selective
ligand concentrate contacts on divergent positions more than Org 43553 does?

| ligand | fragment-contacts at divergent positions |
|---|---|
| Org 43553 in LHCGR | 10 / 36 (**28%**) |
| Cpd-21f in FSHR | 11 / 37 (**30%**) |

**No difference.** The missing half of the tail argument comes back null at ligand level.

What *is* asymmetric is fragment-specific: **Org 43553's morpholine is the only fragment
in either ligand that touches no divergent position at all** (1 of 7 for Org 43553; 0 of 5
for Cpd-21f). The tail argument should be stated as a property of that one fragment, not
as a general property of the ligand — as written, the paper leaves the stronger reading
available and it is not supported.

### 4c. Backbone deviation profile

187 Cα pairs: mean 0.552 Å, median 0.491 Å, max 1.993 Å. The 63 pocket residues in the fit
average **0.522 Å** — indistinguishable from the global mean, so the superposition hides no
local conformational difference at the site. Tyr612/His615 **0.482 Å**, Phe515/Leu518
**0.182 Å**. The largest deviations (Cys353 1.99 Å, Asp360 1.71 Å, Gly358 1.68 Å,
Tyr359 1.66 Å) are in the TM1–TM2/hinge region, away from the ligand. Full profile in
`results/backbone_deviation.tsv`.

---

## 5. What I could not verify

- **Nothing from the earlier session's docking runs.** The smina outputs behind
  "0.93–1.06 Å", "native ranked third to eighth", "89–100% atom overlap", "centroid offset
  0.2–1.4 Å" and the per-function table do not exist in this session. I regenerated
  equivalents with gnina v1.1 rather than restating them. **The manuscript's specific
  per-function numbers are not supported by any file I can hand you** — use the regenerated
  table instead. Direction of travel agrees (sampling works, ranking often fails); the
  blanket conclusion does not (section 1).
- **All pharmacology.** EC₅₀ values, intrinsic activity 0.8, 82% FSHR efficacy, the
  42/45/110 nM spread, the spare-receptor argument — none is checkable from structures, and
  the primary source (van Koppen 2008) is behind a publisher block from this environment.
- **TSHR Tyr667 and Leu570.** Quoted from literature; no TSHR structure was analysed here.
  Either analyse one or cite explicitly.
- **B-factor claim — verified, and it holds.** All 35 heavy atoms of 55Z carry 72.7 with
  zero variance; chain R has 459 distinct values, σ = 11.5; Cpd-21f has 32 distinct values
  across 33 atoms.

---

## 6. Where the manuscript text would change

1. **Abstract, final sentence** — "including convolutional-network rescoring" must be
   deleted. The CNN passes.
2. **"Orientationally ambiguous" section** — rewrite around the corrected table. Report
   that Vinardo and the CNN succeed. Replace the blanket warning with a narrower one.
   Add the mechanism from 4a: the pocket is two near-equivalent hydrophobic lobes, the
   decoy wins on the functions' own terms, and no residue distinguishes the orientations.
3. **Methods, alignment** — reconcile BLOSUM62 with the 60.7% figure.
4. **"74% conserved" heading and text** — add the binomial test and the TM background.
   The honest framing is that the pocket is *typically* conserved for this domain, which
   makes the seven divergences more interesting, not less.
5. **Tail section** — restrict the claim to the morpholine fragment; add the null
   cross-ligand comparison from 4b.
6. **Add** the backbone-deviation profile (4c) as supplementary support for the
   side-chain-identity argument.
7. **Data availability** — can now name the actual scripts and outputs.

---

## 7. Blast radius of the atom-ordering bug — audited

The bug affects any number derived from a coordinate file **written by gnina or
openbabel**. Nothing else. Audited by input type across every script:

| script | reads gnina/openbabel output? | status |
|---|---|---|
| `01_pocket_and_superposition.py` | **no** — PDB only | unaffected |
| `02_fragment_accounting.py` | **no** — PDB only | unaffected |
| `04_make_figure_scripts.py` | **no** — matrix from `01` | unaffected |
| `gphr_lib.py` | **no** | unaffected |
| `03_flipped_pose.py` | **yes** | **was wrong, rewritten** to map atoms by graph isomorphism |
| `05_docking_table.py` | **yes** | **was wrong, rewritten** to use `rmsd_lib` |

The three specific numbers you asked about are all clean:

- **Superposition (187 pairs, 0.64 Å).** Cα coordinates read from `7FIH.pdb` and
  `8I2G.pdb`. No SDF anywhere in the path.
- **Round-trip control (1.8 × 10⁻¹³ Å).** Applies R then Rᵀ to the 55Z coordinates parsed
  from `7FIH.pdb`. It tests matrix arithmetic, not file I/O, and never leaves memory.
- **Mapped *tert*-butyl 0.63 Å.** Compares 55Z from `7FIH.pdb` (transformed) with O6F from
  `8I2G.pdb`. Both PDB, both HETATM records, deposited atom names.

`inputs/xtal_lig.sdf` is the one SDF in the trusted path, and it is not third-party output:
it was written directly from the 7FIH HETATM coordinates with connectivity derived by
covalent radii. It is the *reference*, so its ordering defines the correct ordering.

Everything in AUDIT §2 ("reproduced exactly") is therefore untouched by the bug, and was in
any case re-derived by code that never opens an SDF.

---

## 8. TSHR — the last unsourced claim is now computed

Human TSHR (UniProt **P16473**, 764 aa, length verified against the database's stated
value) aligned to 7FIH chain R by the same Needleman–Wunsch implementation, BLOSUM62,
gap −11/−1. **550 aligned pairs, 60.7% identity.**

**Control — three published TSHR numbers recovered without being supplied:**

| LHCGR | computed TSHR equivalent | literature | |
|---|---|---|---|
| Phe515 (ECL2) | **Leu570** | Leu570 | MATCH |
| Tyr612 (7.42) | **Tyr667** | Tyr667 | MATCH |
| Glu451 (3.37 anchor) | **Glu506** | Glu506 | MATCH |

**The Phe515 claim is confirmed, and it is tight.** Over the 23-residue contact shell,
**exactly one** position is unique to LHCGR against both paralogs:

> **Phe515**, 3.21 Å from Org 43553 — FSHR Leu518, TSHR Leu570. Count = 1.

Three-way pattern over the 63-residue 8 Å shell:

| pattern | n |
|---|---|
| all three identical | 42 |
| TSHR unique | 7 |
| FSHR unique | 6 |
| all three differ | 5 |
| **LHCGR unique** | **3** |

Contact-shell divergence sets, for the paper:

- **vs FSHR (7 of 23):** F515/L518, V519/I522, I528/V531, I531/L534, A593/S596,
  S604/A607, Y612/H615
- **vs TSHR (7 of 23):** A349/E404, F515/L570, V519/T574, I531/V586, F588/Y643,
  A593/I648, K595/N650

Only **Phe515 and Val519** appear in both lists — which is the family-level statement the
paper wants, and it is now computed rather than cited.

**What this does not cover.** The pocket is defined by distance to Org 43553 *in LHCGR*, so
these are the paralog residues at LHCGR's pocket positions. Confirming that TSHR's own
allosteric site occupies the same positions, and computing a TSHR contact shell around a
bound agonist, needs **7XW5 (TSHR + ML-109)**, which this environment cannot download —
the shell has no route to the PDB. If you attach `7XW5.pdb`, `07_three_way_tshr.py`
extends to a full structural three-way in about ten minutes. Until then the paper should
say the TSHR comparison is sequence-based on a structurally defined LHCGR pocket, which is
honest and still removes the unsourced citation.

**One numbering coincidence worth avoiding in the text:** LHCGR/TSHR identity is 60.7%
under BLOSUM62, and the manuscript currently reports 60.7% for LHCGR/**FSHR** (an
identity-matrix value; BLOSUM62 gives 61.2% there). Two different comparisons would print
the same number for different reasons. Fix §3a and the collision disappears.

---

## 9. Figures — validated, not rendered

**No PyMOL in this environment and no route to install one** (PyPI is blocked), so I could
not render the images. What I did instead:

**Every selection was resolved against the actual files.** All fourteen residue selections
in figure 1 exist and carry the expected identities (7FIH 515 PHE, 519 VAL, 528 ILE,
531 ILE, 593 ALA, 604 SER, 612 TYR; 8I2G 518 LEU, 522 ILE, 531 VAL, 534 LEU, 596 SER,
607 ALA, 615 HIS). Every atom-level selection in figure 2 resolves (His615 CE1, Tyr612 CE1,
Ser604 OG, O6F C01 and C27, present both in 8I2G and in the mapped file). Figure 3's
`decoy_best.pdb` has all 35 heavy atoms with the reference atom names. **No silent empty
selection anywhere** — that is the failure mode that would otherwise produce a
convincing-looking but wrong figure.

**Each script now prints a numeric self-check**, because the one thing I cannot verify
without running PyMOL is whether it interprets the 4×4 matrix the way I laid it out
(translation in the fourth column) or transposed:

- **Fig 1** measures the mapped Org 43553 *tert*-butyl against Cpd-21f's and prints
  `TRANSFORM OK` if it lands at 0.63 Å. If it prints a large number, re-run
  `transform_selection` with `transpose=1`. This is the single most likely thing to go
  wrong, and it now announces itself.
- **Fig 2** prints the two quoted distances; they must read **3.69** and **2.31** Å.
- **Fig 3** measures the decoy's *tert*-butyl against the native morpholine centroid and
  confirms the head-to-tail exchange if it is under 1.5 Å (expected 0.5 Å).

**What each should look like:**

1. **Fig 1** — two transparent cartoons superposed (LHCGR grey, FSHR pale cyan) with two
   sticks ligands almost on top of each other, and fourteen labelled side chains in orange
   and blue clustered around them. The eye should go to Tyr612/His615 and Phe515/Leu518.
2. **Fig 2** — two panels 38 Å apart. Left: FSHR with Cpd-21f and a black dashed line to
   His615 reading 3.69. Right: the same ligand in LHCGR with a **red** dashed line to
   Tyr612 reading 2.31, plus an orange line to Ser604 at 2.77. The point is that the red
   line is visibly shorter and the tyrosine ring is jammed into the ligand.
3. **Fig 3** — one pocket, two ligand copies in the same volume (green native, magenta
   decoy), with spheres marking each molecule's two ends and black dashed lines crossing
   between them. The crossing lines are the figure: they show the molecule is reversed, not
   displaced.

Set `set_wd` at the top of each file to the directory holding `7FIH.pdb`, `8I2G.pdb` and
`results/`, then `pymol figures/fig1_pockets_superposed.pml`.
