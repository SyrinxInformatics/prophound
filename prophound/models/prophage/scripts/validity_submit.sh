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
# Submit the scientific-validity analysis jobs (docs/scientific-validity.md).
#   ssh pke835@arc.utsa.edu 'bash /work/$USER/mge_explorer/data_refactored/scripts/validity_submit.sh'
#
# No conda here on purpose. The `miniconda` module is NOT available on the login
# node ("This module is available on the job nodes" -- docs/environment.md), so
# the env cannot be inspected before submitting; both job scripts load it
# themselves on the compute node. composition_baseline.py imports sklearn at
# module scope so a missing library kills that job in seconds, not after an hour
# of featurization. `sbatch` itself needs no conda.
set -uo pipefail
cd "/work/$USER/mge_explorer/data_refactored" || exit 1
mkdir -p analysis/_slurm_logs analysis/composition_baseline analysis/clades

echo "=== submitting ==="
sbatch scripts/cluster_prophage_elements.slurm
sbatch scripts/composition_baseline.slurm

echo
squeue -u "$USER" -o "%.10i %.22j %.9P %.8T %.10M %.10l %R"
echo
echo "Logs: analysis/_slurm_logs/<jobname>_<jobid>.{out,err}"
echo "Results: analysis/clades/clade_summary.json, analysis/composition_baseline/results.json"
