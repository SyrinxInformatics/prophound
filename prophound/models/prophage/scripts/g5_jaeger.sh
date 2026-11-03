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
# G5 fragmented-input comparison (docs/claims-and-evidence.md "G5 preregistration"):
# LOCAL FALLBACK ONLY: superseded 2026-09-30 by scripts/g5_jaeger.slurm (ARC).
# Jaeger 1.26.5 on the complete LAMBDA chromosomes and the four shred conditions,
# locally, both models (57341_1.5M primary, bundled `default` secondary), then
# collects the calls into analysis/g5/jaeger_calls.tsv. Resumable. Inputs are
# checked against analysis/g5/input_md5.tsv first, the same check ARC makes.
# Run from prophound/models/prophage:
#   bash scripts/g5_jaeger.sh
set -uo pipefail
G=analysis/g5
ENV=../../../metadata/g5-comparators
RAW=$G/raw/jaeger
MODELS="57341_1.5M:jaeger_57341_1.5M_fragment:--fsize=3077 default:default:"   # amendment 2026-09-30
CONDS="complete mag_s0 mag_s1 mag_s2 draft_s0"
J="pixi run --manifest-path $ENV/pixi.toml jaeger"

# ---- preflight ------------------------------------------------------------------
fail=0
( cd "$G" && while IFS=$'\t' read -r path want; do
    got=$(md5 -q "$path" 2>/dev/null || md5sum "$path" | cut -d' ' -f1)
    [ "$got" = "$want" ] || { echo "PREFLIGHT FAIL: $path md5 $got != $want" >&2; exit 1; }
  done < input_md5.tsv ) || fail=1
v=$($J --version 2>&1 | grep -o '[0-9]\+\.[0-9]\+\.[0-9]\+' | tail -1)
[ "$v" = "1.26.5" ] || { echo "PREFLIGHT FAIL: jaeger $v, expected 1.26.5" >&2; fail=1; }
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] jaeger $v; all 195 inputs match input_md5.tsv"
# ---------------------------------------------------------------------------------

mkdir -p "$RAW"
for m in $MODELS; do
  IFS=: read -r label model extra <<< "$m"; extra=${extra//=/ }
  for c in $CONDS; do
    src=$([ "$c" = complete ] && echo "$G/fasta" || echo "$G/shreds/$c")
    echo "[$label $c] $(date)"
    while read -r ASM; do
      o="$RAW/$label/$c/$ASM"
      ls "$o"/*/*.log >/dev/null 2>&1 && ls "$o"/*/"$ASM"*.tsv >/dev/null 2>&1 && continue    # resumable
      rm -rf "$o"; s=$(date +%s)
      $J predict -m "$model" -i "$src/$ASM.fna" -o "$o" -p --lc 1500 --cpu $extra --workers 8 -f >"$RAW/last.log" 2>&1 \
        && echo -e "$label\t$c\t$ASM\t$(( $(date +%s) - s ))" >> "$RAW/wall_seconds.tsv" \
        || { echo "FAILED $label $c $ASM" | tee -a "$RAW/failures.txt" >&2; }
    done < "$G/genome_list.txt"
  done
done

pixi run --manifest-path ../../../pyproject.toml python scripts/g5_score.py collect-jaeger --root "$RAW" --out "$G/jaeger_calls.tsv"
cp "$RAW/wall_seconds.tsv" "$G/jaeger_wall_seconds.tsv"
[ -s "$RAW/failures.txt" ] && { echo "WARNING: $(wc -l < "$RAW/failures.txt") failed runs" >&2; exit 2; }
echo "[done] $(date)"
