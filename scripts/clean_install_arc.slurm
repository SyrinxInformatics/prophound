#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Clean-install check on Linux (ARC CPU; internal release checklist, robustness check 5).
# Makes a fresh conda env with only Python, installs the prophound wheel with pip (dependencies
# from PyPI, as a user would), then runs clean_install_check.sh on the example and one
# second-lockbox genome against its 851599 calls. The .out should end "[check] PASS".
# Copy the wheel, clean_install_check.sh, this file and examples/{hin_rd_kw20_region.fna,expected_calls.tsv}
# into /work/$USER/mge_explorer/clean_install, then from there: sbatch clean_install_arc.slurm
#SBATCH -J clean_install
#SBATCH -p compute3
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 0-06:00:00
#SBATCH -o clean_install_%j.out
#SBATCH -e clean_install_%j.err

set -euo pipefail
DIR=/work/$USER/mge_explorer/clean_install
G=GCF_000143845.1_ASM14384v1_genomic
LG=/work/$USER/mge_explorer/raw/lockbox2/$G.fna
LREF=/work/$USER/mge_explorer/data_refactored/analysis/lockbox2/prophound/default/$G/regions.detailed.tsv
cd "$DIR"
for f in prophound-0.1.0-py3-none-any.whl clean_install_check.sh hin_rd_kw20_region.fna expected_calls.tsv "$LG" "$LREF"; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; exit 1; }
done

module load miniconda
rm -rf "$DIR/env"
conda create -y -q -p "$DIR/env" python=3.11 pip
"$DIR/env/bin/pip" install -q "$DIR/prophound-0.1.0-py3-none-any.whl"
"$DIR/env/bin/pip" list 2>/dev/null | grep -i -E "^(prophound|torch|numpy|pandas) "

bash "$DIR/clean_install_check.sh" "$DIR/env/bin/prophound" "$DIR/env/bin/python" \
  "$DIR/hin_rd_kw20_region.fna" "$DIR/expected_calls.tsv" "$LG" "$LREF" "$DIR/run" cpu
