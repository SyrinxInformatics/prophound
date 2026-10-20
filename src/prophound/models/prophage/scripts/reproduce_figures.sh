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
# Re-run every figure script whose inputs are local, into OUT_DIR, then compare with the
# recorded outputs:
#   scripts/reproduce_figures.sh OUT_DIR && pixi run python scripts/compare_figures.py . OUT_DIR
# Run from the model directory (src/prophound/models/prophage). The commands are the ones that
# produced the figures in docs/release-text-v1.md "Reproduction". Scripts whose inputs live only
# on ARC (offtarget_annotation, small_element_spikein, literature_pilot, eval_containment) are not
# here. docs/tech-debt-audit.md X2.
set -u
O=$1; mkdir -p $O/{lockbox,g5,g5b,lambda,offtarget,integrity}
P="pixi run --quiet python"
run(){ name=$1; shift; if "$@" > $O/$name.log 2>&1; then echo "[ran] $name"; else echo "[FAILED] $name (see $O/$name.log)"; fi; }
run lockbox $P scripts/lockbox_score.py --elements analysis/lockbox/lockbox_elements.csv \
  --genome_list analysis/lockbox/lockbox_genome_list.txt --run default=analysis/lockbox/prophound/default \
  --run exploratory=analysis/lockbox/prophound/exploratory --run conservative=analysis/lockbox/prophound/conservative \
  --out $O/lockbox/results.json
run g5 $P scripts/g5_score.py score --lambda_dir analysis/external_benchmarks/lambda \
  --prophound complete=analysis/lambda/ens --prophound mag_s0=analysis/g5/prophound/mag_s0 \
  --prophound mag_s1=analysis/g5/prophound/mag_s1 --prophound mag_s2=analysis/g5/prophound/mag_s2 \
  --prophound draft_s0=analysis/g5/prophound/draft_s0 \
  --calls analysis/g5/genomad_calls.tsv --calls analysis/g5/jaeger_calls.tsv --shreds analysis/g5/shreds \
  --out $O/g5/results.json
run g5b $P scripts/g5_score.py score --lambda_dir analysis/external_benchmarks/lambda \
  --truth_tsv analysis/g5b/g5b_truth.tsv --genomes_tsv analysis/g5b/g5b_genomes.tsv \
  --prophound complete=analysis/g5b/prophound/complete --prophound mag_s0=analysis/g5b/prophound/mag_s0 \
  --prophound mag_s1=analysis/g5b/prophound/mag_s1 --prophound mag_s2=analysis/g5b/prophound/mag_s2 \
  --prophound draft_s0=analysis/g5b/prophound/draft_s0 \
  --calls analysis/g5b/genomad_calls.tsv --calls analysis/g5b/jaeger_calls.tsv --shreds analysis/g5b/shreds \
  --edge_snap_frac 0.5 --lengths analysis/g5b/contig_lengths.tsv --out $O/g5b/results.json
run lambda $P scripts/lambda_score.py --data analysis/external_benchmarks/lambda --prophound_run analysis/lambda/ens \
  --containment analysis/integrity/containment_lambda.tsv --out $O/lambda/scores.json
run chance $P scripts/lambda_chance.py --data analysis/external_benchmarks/lambda --prophound_run analysis/lambda/ens \
  --out $O/lambda/chance.json
run offtarget $P scripts/offtarget_enrichment_ci.py --annotation analysis/offtarget/annotation_ctrl_elements_only.csv \
  --out $O/offtarget/enrichment_ci.json
run recall42 $P scripts/integrity_checks.py recall --elements scripts/eval_sizestrat_elements.csv \
  --root analysis/multiscale --suffix _3000 --out $O/integrity/recall42.json
run recall18 $P scripts/integrity_checks.py recall --elements scripts/eval_confirm18_elements.csv \
  --root analysis/confirm18 --out $O/integrity/recall18.json
run recallpf $P scripts/integrity_checks.py recall --root analysis/phagefree --out $O/integrity/recall_phagefree.json
mkdir -p $O/seedvote $O/keep_rules
run seedvote18 $P scripts/seed_vote_confirm.py --elements scripts/eval_confirm18_elements.csv --root analysis/confirm18 \
  --out $O/seedvote/confirm18.csv
run seedvotepf $P scripts/seed_vote_confirm.py --root analysis/phagefree --out $O/seedvote/phagefree.csv
run keeprules $P scripts/keep_rule_compare.py $O/keep_rules/results.csv
