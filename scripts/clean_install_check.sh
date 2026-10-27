#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
# Clean-install check (internal release checklist, robustness check 5): runs an installed prophound,
# from outside the repository, on the example genome and on one second-lockbox genome, and
# compares the calls with the expected ones (the example's expected_calls.tsv; the lockbox
# genome's calls from the ARC V100 run 851599). Devices never agree byte for byte, only in calls.
#
#   bash clean_install_check.sh <prophound executable> <python with pandas> <example.fna> <expected_calls.tsv> \
#        <lockbox genome .fna> <its regions.detailed.tsv from 851599> <output dir> [device]
set -euo pipefail
PH=$1; PY=$2; EX=$3; EXPECTED=$4; LG=$5; LREF=$6; OUT=$7; DEV=${8:-auto}
mkdir -p "$OUT"; cd "$OUT"
echo "[check] $("$PH" --version); device $DEV; $(uname -sm)"
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
