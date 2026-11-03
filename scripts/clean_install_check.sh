#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
#
# Copyright 2025 Jacob Alford <jalford0000@gmail.com>
#                and Benjamin Alford <ssolar5.dev@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Clean-install check (docs/benchmarks.md "Cross-device agreement"): runs an installed prophound,
# from outside the repository, on the example genome and on one second-lockbox genome, and
# compares the calls with the expected ones (the example's expected_calls.tsv; the lockbox
# genome's calls from the ARC V100 run 851599). Devices never agree byte for byte, only in calls.
#
#   bash clean_install_check.sh <prophound executable> <python with pandas> <example.fna> <expected_calls.tsv> \
#        <lockbox genome .fna> <its regions.detailed.tsv from 851599> <output dir> [device]
set -euo pipefail
PH=$1; PY=$2; EX=$3; EXPECTED=$4; LG=$5; LREF=$6; OUT=$7; DEV=${8:-auto}
# resolve relative paths before the cd below (a bare command name is left to PATH)
abs() { case $1 in /*) echo "$1";; */*) echo "$PWD/$1";; *) echo "$1";; esac; }
PH=$(abs "$PH"); PY=$(abs "$PY"); EX=$(abs "$EX"); EXPECTED=$(abs "$EXPECTED"); LG=$(abs "$LG"); LREF=$(abs "$LREF"); OUT=$(abs "$OUT")
for f in "$EX" "$EXPECTED" "$LG" "$LREF"; do [ -s "$f" ] || { echo "[check] FAIL: missing $f" >&2; exit 1; }; done
VERSION=$("$PH" --version) || { echo "[check] FAIL: $PH --version did not run" >&2; exit 1; }
mkdir -p "$OUT"; cd "$OUT"
echo "[check] $VERSION; device $DEV; $(uname -sm)"
"$PH" predict --input "$EX" --outdir "$OUT/example" --device "$DEV" -q --force
"$PY" - "$OUT/example/regions.detailed.tsv" "$EXPECTED" <<'EOF'
import sys, pandas as pd
cols = ["call_id", "contig_id", "start", "end", "boundary", "contig_edge", "seed_support"]
g, e = (pd.read_csv(p, sep="\t")[cols] for p in sys.argv[1:])
pd.testing.assert_frame_equal(g, e); print(f"[check] example: {len(g)} calls, identical to expected_calls.tsv")
EOF
"$PH" predict --input "$LG" --outdir "$OUT/lockbox_one" --device "$DEV" -q --force
"$PY" - "$OUT/lockbox_one/regions.detailed.tsv" "$LREF" <<'EOF'
import sys, pandas as pd
cols = ["contig_id", "start", "end", "boundary", "seed_support"]
m, c = (pd.read_csv(p, sep="\t")[cols].sort_values(cols).reset_index(drop=True) for p in sys.argv[1:])
assert m.equals(c), pd.concat([m.assign(src="here"), c.assign(src="851599")]).drop_duplicates(subset=cols, keep=False).to_string()
print(f"[check] lockbox genome: {len(m)} calls, identical to the ARC V100 run")
EOF
echo "[check] PASS"
