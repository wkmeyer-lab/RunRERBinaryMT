#!/bin/bash

# First argument is used in job names
j=$1

echo "Current directory: $(pwd)"
echo "Batch script path: $(realpath runScriptSecondPartition.slr)"
file "$(realpath runScriptSecondPartition.slr)"

# Submit main analysis jobs
echo "Submitting: $*"
sbatch -t 3300 --job-name="MainRerMergeEnrichCombineCombine${j}" --dependency=singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunRERAndCorrelation.R "$@"

echo "Submitting: $*"
sbatch -t 3300 --job-name="MainRerMergeEnrichCombineCombine${j}" --dependency=singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunMergedAnalysis.R "$@"

echo "Submitting: $*"
sbatch -t 3300 --job-name="MainRerMergeEnrichCombineCombine${j}" --dependency=singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT EnrichmentAnalysisNew.R "$@"

echo "Submitting: $*"
sbatch -t 3300 --job-name="MainRerMergeEnrichCombineCombine${j}" --dependency=singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT CategoricalCombinePairwiseResults.R "$@"

# Process alternates
for i in $(seq 1 2); do
    echo "Submitting: $*"
    sbatch -t 3300 --job-name="Alternates${j}${i}" --dependency=singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT ProcessAlternates.R "$@" "i=$i"
done

# Run alternate enrichments
for i in $(seq 1 2); do
    echo "Submitting: $*"
    sbatch -t 3300 --job-name="Alternates${j}${i}" --dependency=singleton -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT RunAlternateEnrichments.R "$@" "i=$i"
done

# Get alternate job dependencies
deps=$(squeue -h -o "%A %j" | awk '$2 ~ /Alternates/ {print $1}' | paste -sd:)

# Combine alternates after alternate jobs finish
if [ -n "$deps" ]; then
    echo "Submitting: $*"
    sbatch -t 3300 --job-name="CombineAlts${j}" --dependency="afterok:$deps" -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT CombineAlternates.R "$@"
fi

# Get combine job dependencies
deps2=$(squeue -h -o "%A %j" | awk '$2 ~ /Combine/ {print $1}' | paste -sd:)

# Add alternates to combinations after combine jobs finish
if [ -n "$deps2" ]; then
    echo "Submitting: $*"
    sbatch -t 3300 --job-name="GeneAltsToCombo${j}" --dependency="afterok:$deps2" -n 1 -p rapids runScriptSecondPartition.slr RunRERBinaryMT AddAlternatesToCombinations.R "$@"
fi

