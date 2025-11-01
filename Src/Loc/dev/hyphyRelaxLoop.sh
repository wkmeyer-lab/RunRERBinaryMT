#!/bin/bash

#The full path to your file containing the list of gene names.
GENE_FILE="/share/ceph/wym219group/shared/projects/seaverData/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/combinedSignificantGenes.txt"

#The name of the script Slurm will run for each job.
JOB_SCRIPT="clusterHyphy.slr"

#The partition to run on.
PARTITION="rapids"

#The foreground value to use.
FOREGROUNDVALUE

#The maximum number of jobs to run at the same time.
MAX_CONCURRENT_JOBS=100

# Count the total number of genes
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
    
    2. Calculate the "job slot" (from 0 to 99)
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
echo ""
echo "To CANCEL all jobs (if needed):"
echo "  scancel -u $USER -n hyphy_slot_"
echo "========================================================"
echo "Use 'squeue -u $USER' to monitor your jobs."
