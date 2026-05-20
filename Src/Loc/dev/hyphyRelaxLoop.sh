#!/bin/bash

#full path to your file containing the list of gene names
GENE_FILE="/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/combinedSignificantGenes.txt"

#name of the script Slurm will run for each job.
JOB_SCRIPT="clusterHyphy.slr"

#partition to run on.
PARTITION="rapids"

#foreground value to use.
FOREGROUNDVALUE=2

#maximum number of jobs to run at the same time.
MAX_CONCURRENT_JOBS=100

#base output directory for RERBinaryMT Hyphy runs.
OUTPUT_DIR="/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy"

#count the total number of genes
NUM_GENES=$(wc -l < "$GENE_FILE" | xargs)
if [[ "$NUM_GENES" -eq 0 ]]; then
    echo "Error: Gene file is empty or not found at: $GENE_FILE"
    exit 1
fi

echo "Found $NUM_GENES genes. Submitting jobs to partition '$PARTITION'"
echo "Using '--dependency singleton' to limit to $MAX_CONCURRENT_JOBS concurrent jobs."
echo "This may take a minute to submit all $NUM_GENES jobs..."

#Loop from 1 to NUM_GENES
for (( i=1; i<=$NUM_GENES; i++ )); do

    #Get the gene name for the current line
    CURRENT_GENE=$(sed -n "${i}p" "$GENE_FILE")

    #If gene name is empty, skip
    if [ -z "$CURRENT_GENE" ]; then
        echo "Warning: Empty line at $i, skipping."
        continue
    fi

    #construct the expected output filename.
    #pattern: /.../Hyphy/CategoricalInsVertivoreTree-Hyphy-relax-GENE-Foreground_VALUE.Rds
    EXPECTED_FILE="$OUTPUT_DIR/CategoricalInsVertivoreTree-Hyphy-relax-${CURRENT_GENE}-Foreground_${FOREGROUNDVALUE}.Rds"

    #check if the file already exists.
    if [ -s "$EXPECTED_FILE" ]; then
        echo "Output found for $CURRENT_GENE (Line $i), skipping job submission."
        continue #skip to the next gene in the loop
    fi

    # Calculate the "job slot" (from 0 to 99)
    JOB_SLOT=$(($i % $MAX_CONCURRENT_JOBS))

    #Submit the job
    sbatch \
      -N 1 \
      -p $PARTITION \
      --job-name="hyphy_slot_$JOB_SLOT" \
      --dependency=singleton \
      "$JOB_SCRIPT" \
      "RunRERBinaryMT" \
      "relax" \
      "CategoricalInsVertivoreTree" \
      "$CURRENT_GENE" \
      $FOREGROUNDVALUE \
      "TRUE"

    #Optional: add a small sleep to avoid overwhelming the scheduler
    if (($i % 100 == 0)); then
        echo "Submitted $i / $NUM_GENES jobs..."
        sleep 0.5
    fi

done

echo ""
echo "========================================================"
echo "All $NUM_GENES jobs have been submitted."
echo "Slurm will manage the queue to run max $MAX_CONCURRENT_JOBS at a time."
echo "========================================================"
