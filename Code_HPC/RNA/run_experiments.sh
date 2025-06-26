#!/usr/bin/env bash

# Load R module (adjust module name/version as needed)
module load r/4.3.1

# Create a logs directory for job output (optional)
mkdir -p logs

# Loop seeds 1 through 100
for seed in $(seq 1 100); do
  sbatch \
    --mem=100GB \
    --job-name=RNA_Seq_${seed} \
    --output=logs/RNA_Seq_${seed}.out \
    --wrap="Rscript RNA_Seq.R ${seed}"
done
