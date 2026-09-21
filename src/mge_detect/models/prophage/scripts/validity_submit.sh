#!/bin/bash
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
