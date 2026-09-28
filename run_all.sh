#!/usr/bin/env bash
# Regenerates every number in the manuscript from the two deposited structures.
#
#   ./run_all.sh                 # structural analyses only (~3 min, no dependencies beyond numpy)
#   GNINA=/path/to/gnina ./run_all.sh   # also redoes all docking (~35 min on 2 cores)
#
# Outputs land in results/ and figures/. Nothing else is written.
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD/code:${PYTHONPATH:-}"

echo "=== 01  pockets, alignment, superposition, reciprocal mappings ==="
python3 code/01_pocket_and_superposition.py

echo "=== 02  per-fragment burial and contacts, both ligands ==="
python3 code/02_fragment_accounting.py

echo "=== 07  three-way family comparison with TSHR ==="
python3 code/07_three_way_tshr.py

echo "=== 08  per-fragment burial table as reported in the manuscript ==="
python3 code/08_manuscript_fragment_table.py

if [ -n "${GNINA:-}" ]; then
  echo "=== docking sweep, 8 configurations ==="
  GNINA="$GNINA" bash code/run_docking_sweep.sh
  echo "=== sampling-independent rescoring control ==="
  GNINA="$GNINA" bash code/06_rescoring_control.sh
  echo "=== 05  supplementary docking table ==="
  python3 code/05_docking_table.py
  echo "=== 03  decoy orientation ==="
  python3 code/03_flipped_pose.py
else
  echo "=== docking skipped (set GNINA=/path/to/gnina to include it) ==="
  if [ -d results/dock ]; then
    echo "    re-analysing the shipped docking outputs instead"
    python3 code/05_docking_table.py
    python3 code/03_flipped_pose.py
  fi
fi

echo "=== 04  PyMOL matrix and figure scripts ==="
python3 code/04_make_figure_scripts.py

echo
echo "done. Controls are in results/validation.txt and results/three_way_summary.txt."
