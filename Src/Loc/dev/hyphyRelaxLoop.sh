#!/bin/bash

# --- COMMAND LINE ARGUMENTS ---
FILE_PREFIX=""
FOREGROUNDVALUE=""
RATES=""
MAX_CONCURRENT_JOBS=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --tree) FILE_PREFIX="$2"; shift 2 ;;
        --foreground) FOREGROUNDVALUE="$2"; shift 2 ;;
        --rates) RATES="$2"; shift 2 ;;
        --max-jobs) MAX_CONCURRENT_JOBS="$2"; shift 2 ;;
        *) echo "ERROR: Unknown parameter passed: $1"; exit 1 ;;
    esac
done

if [ -z "$FILE_PREFIX" ] || [ -z "$FOREGROUNDVALUE" ] || [ -z "$RATES" ] || [ -z "$MAX_CONCURRENT_JOBS" ]; then
    echo "ERROR: Missing mandatory arguments."
    echo "Usage: bash $0 --tree <phenotype_tree_name> --foreground <foreground_value> --rates <rates_value> --max-jobs <max_jobs>"
    echo "Example: bash $0 --tree ComplexDietCentralAnalysis --foreground 4 --rates 3 --max-jobs 95"
    exit 1
fi

echo "Running master submission script for Phenotype Tree: $FILE_PREFIX | Foreground Category: $FOREGROUNDVALUE | Rates: $RATES | Max Jobs: $MAX_CONCURRENT_JOBS"

# --- CONFIGURATION ---
GENE_FILE="/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Src/Loc/dev/manuscript_genelist.txt"
JOB_SCRIPT="clusterHyphy.sh"
PARTITION="hawkcpu"

# Updated OUTPUT_DIR to nest by Rate
OUTPUT_DIR="/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/${FILE_PREFIX}/Hyphy/Rate_${RATES}"

# Updated LOG_DIR to nest by Rate
LOG_DIR="/share/ceph/wym219group/shared/erroutfiles/seaverErrors/${FILE_PREFIX}/Foreground_${FOREGROUNDVALUE}/Rate_${RATES}"
mkdir -p "$LOG_DIR"
chmod -R 777 "$LOG_DIR"
mkdir -p "$OUTPUT_DIR"
chmod -R 777 "$OUTPUT_DIR"

NUM_GENES=$(wc -l < "$GENE_FILE" | xargs)
if [[ "$NUM_GENES" -eq 0 ]]; then
    echo "Error: Gene file is empty or not found at: $GENE_FILE"
    exit 1
fi

# 1. Create the master filtered list of remaining genes
MASTER_LIST="all_tasks_fg${FOREGROUNDVALUE}_r${RATES}.txt"
> "$MASTER_LIST" 

for (( i=1; i<=$NUM_GENES; i++ )); do
    CURRENT_GENE=$(sed -n "${i}p" "$GENE_FILE")
    [ -z "$CURRENT_GENE" ] && continue

    # Check if the output JSON file for THIS specific rate already exists and is non-empty
    EXPECTED_FILE="$OUTPUT_DIR/${FILE_PREFIX}-Hyphy-relax-${CURRENT_GENE}-Foreground_${FOREGROUNDVALUE}-Rate_${RATES}.json"
    if [ ! -s "$EXPECTED_FILE" ]; then
        echo "$CURRENT_GENE" >> "$MASTER_LIST"
    fi
done

TOTAL_TASKS=$(wc -l < "$MASTER_LIST" | xargs)

if [[ "$TOTAL_TASKS" -eq 0 ]]; then
    echo "All jobs for foreground $FOREGROUNDVALUE and rate $RATES have already been processed."
    rm "$MASTER_LIST"
    exit 0
fi

# 2. Split the master list into chunks using the dynamic max jobs value
CHUNK_PREFIX="${FILE_PREFIX}chunk_fg${FOREGROUNDVALUE}_r${RATES}_"
rm -f ${CHUNK_PREFIX}* 2>/dev/null
split -l "$MAX_CONCURRENT_JOBS" "$MASTER_LIST" "$CHUNK_PREFIX"

FIRST_CHUNK=$(ls ${CHUNK_PREFIX}* 2>/dev/null | head -n 1)
RUN_TASKS=$(wc -l < "$FIRST_CHUNK" | xargs)

echo "Total remaining tasks: $TOTAL_TASKS. Submitting first chunk ($RUN_TASKS tasks)..."

# 3. Create the array wrapper script
WRAPPER_SCRIPT="array_wrapper_fg${FOREGROUNDVALUE}_r${RATES}.slr"
cat << EOF > "$WRAPPER_SCRIPT"
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 48:00:00
#SBATCH -o ${LOG_DIR}/hyphy_%A_%a.out
#SBATCH -e ${LOG_DIR}/hyphy_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-user=aws519@lehigh.edu

TASK_FILE=\$1
FOREGROUNDVALUE=\$2
FILE_PREFIX=\$3
JOB_SCRIPT=\$4
RATES=\$5

# Pull the specific gene mapping directly to this Slurm Array Task ID
CURRENT_GENE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\$TASK_FILE")

if [ ! -z "\$CURRENT_GENE" ]; then
    echo "Task ID \$SLURM_ARRAY_TASK_ID processing gene: \$CURRENT_GENE"
    bash "\$JOB_SCRIPT" \\
        "RunRERBinaryMT" \\
        "relax" \\
        "\${FILE_PREFIX}" \\
        "\${CURRENT_GENE}" \\
        "\${FOREGROUNDVALUE}" \\
        "TRUE" \\
        "Background" \\
        "\${RATES}"
fi
EOF

# 4. Submit exactly the size of the current chunk
sbatch \
  -p "$PARTITION" \
  -J "hyphy_fg${FOREGROUNDVALUE}_r${RATES}" \
  --array=1-${RUN_TASKS} \
  "$WRAPPER_SCRIPT" \
  "$FIRST_CHUNK" \
  "$FOREGROUNDVALUE" \
  "$FILE_PREFIX" \
  "$JOB_SCRIPT" \
  "$RATES"

echo "========================================================"
echo "Submitted array for chunk $FIRST_CHUNK ($RUN_TASKS jobs)."
echo "When this completes, run this script again to process the next chunk."
echo "========================================================"

rm "$MASTER_LIST"
