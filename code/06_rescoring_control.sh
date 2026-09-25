#!/bin/bash
# Sampling-independent control: place the crystallographic pose (raw and locally
# minimised) into an ensemble of docked poses and ask each scoring function to rank it.
# INPUTS  ../inputs/receptor.pdb ../inputs/xtal_lig.sdf ../results/dock/vina_rigid.sdf
# OUTPUT  ../results/rescore/<sf>.txt
set -u
GNINA=${GNINA:-/home/claude/dock/gnina11}
IN=${IN:-/home/claude/paper/inputs}; RES=${RES:-/home/claude/paper/results}
mkdir -p "$RES/rescore"
$GNINA -r "$IN/receptor.pdb" -l "$IN/xtal_lig.sdf" --minimize --cnn_scoring none \
       --no_gpu --cpu 2 -o "$RES/rescore/xtal_min.sdf" > "$RES/rescore/minimize.log" 2>&1
head -n -1 "$RES/rescore/xtal_min.sdf" > "$RES/rescore/_a"; echo '$$$$' >> "$RES/rescore/_a"
cat "$RES/rescore/_a" "$RES/dock/vina_rigid.sdf" > "$RES/rescore/ensemble.sdf"; rm -f "$RES/rescore/_a"
for SF in vina vinardo ad4_scoring dkoes_scoring dkoes_fast; do
  $GNINA -r "$IN/receptor.pdb" -l "$RES/rescore/ensemble.sdf" --score_only \
         --scoring $SF --cnn_scoring none --no_gpu --cpu 2 > "$RES/rescore/$SF.txt" 2>&1
done
$GNINA -r "$IN/receptor.pdb" -l "$RES/rescore/ensemble.sdf" --score_only \
       --cnn_scoring rescore --no_gpu --cpu 2 > "$RES/rescore/cnn.txt" 2>&1
echo DONE > "$RES/rescore/status"
