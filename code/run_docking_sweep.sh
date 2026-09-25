#!/bin/bash
# Regenerates every docking configuration reported in the paper.
# INPUTS (expected in ../inputs/):  receptor.pdb  xtal_lig.sdf
# REQUIRES: gnina v1.1 binary at $GNINA (contains smina's scoring functions)
# OUTPUT:   ../results/dock/<config>.sdf and <config>.log
set -u
GNINA=${GNINA:-/home/claude/dock/gnina11}
IN=${IN:-/home/claude/paper/inputs}
OUT=${OUT:-/home/claude/paper/results/dock}
mkdir -p "$OUT"
COMMON="-r $IN/receptor.pdb -l $IN/xtal_lig.sdf --autobox_ligand $IN/xtal_lig.sdf --autobox_add 4 \
        --exhaustiveness 16 --num_modes 9 --seed 42 --no_gpu --cpu 2"
for SF in vina vinardo ad4_scoring dkoes_scoring dkoes_fast; do
  echo "=== $SF rigid ==="
  $GNINA $COMMON --scoring $SF --cnn_scoring none -o "$OUT/${SF}_rigid.sdf" > "$OUT/${SF}_rigid.log" 2>&1
done
# Lys595 flexible (LHCGR numbering, chain R)
for SF in vina vinardo; do
  echo "=== $SF flexLys595 ==="
  $GNINA $COMMON --scoring $SF --cnn_scoring none --flexres R:595 \
         -o "$OUT/${SF}_flexK595.sdf" > "$OUT/${SF}_flexK595.log" 2>&1
done
echo "=== gnina CNN rescore ==="
$GNINA $COMMON --cnn_scoring rescore -o "$OUT/gnina_cnn.sdf" > "$OUT/gnina_cnn.log" 2>&1
echo ALLDONE > "$OUT/sweep.status"
