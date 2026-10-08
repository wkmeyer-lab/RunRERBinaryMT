#!/bin/bash

# Configuration
GENE_LIST="manuscript_genelist.txt"
BASE_DIR="../../../Output"
MAX_JOBS=95

# List of 3 outer directories
MAIN_DIRS=(
    "ComplexDietCentralAnalysisSimplify2"
    "ComplexDietCentralAnalysisSimplifyEqualDrop"
    "ComplexDietCentralAnalysisSimplifyStrictPred"
)

# List of 3 rate subdirectories
RATES=(
    "Rate_2"
    "Rate_3"
    "Rate_4"
)

# Verify base output directory exists before scanning
if [[ ! -d "$BASE_DIR" ]]; then
    echo "Error: Output directory '$BASE_DIR' not found." >&2
    exit 1
fi

# Automatically discover all unique foreground values present in output files
mapfile -t FOREGROUNDS < <(
    find "$BASE_DIR" -type f -name "*-Hyphy-relax-*.json" 2>/dev/null \
    | grep -oE 'Foreground_[a-zA-Z0-9_]+' \
    | sort -u
)

if [[ ${#FOREGROUNDS[@]} -eq 0 ]]; then
    echo "Error: No foreground values detected in '$BASE_DIR'." >&2
    exit 1
fi

# Verify input file exists
if [[ ! -f "$GENE_LIST" ]]; then
    echo "Error: Gene list file '$GENE_LIST' not found in current directory." >&2
    exit 1
fi

declare -A MISSING_COMMANDS
missing_count=0
empty_count=0
total_checked=0

# Loop through each gene and check output files
while IFS= read -r gene || [[ -n "$gene" ]]; do
    gene=$(echo "$gene" | tr -d '\r' | xargs)
    [[ -z "$gene" ]] && continue

    for main_dir in "${MAIN_DIRS[@]}"; do
        for rate in "${RATES[@]}"; do
            for fg in "${FOREGROUNDS[@]}"; do
                filename="${main_dir}-Hyphy-relax-${gene}-${fg}-${rate}.json"
                filepath="${BASE_DIR}/${main_dir}/Hyphy/${rate}/${filename}"

                ((total_checked++))

                if [[ ! -f "$filepath" || ! -s "$filepath" ]]; then
                    if [[ ! -f "$filepath" ]]; then
                        ((missing_count++))
                    else
                        ((empty_count++))
                    fi

                    # Extract raw numerical/string values (e.g., Foreground_4 -> 4, Rate_2 -> 2)
                    fg_val="${fg#Foreground_}"
                    rate_val="${rate#Rate_}"

                    # Store unique run command
                    cmd_key="${main_dir}|${fg_val}|${rate_val}"
                    MISSING_COMMANDS["$cmd_key"]="bash hyphyRelaxLoop.sh --tree ${main_dir} --foreground ${fg_val} --rates ${rate_val} --max-jobs ${MAX_JOBS}"
                fi
            done
        done
    done
done < "$GENE_LIST"

# Output audit summary to stderr so stdout remains clean executable code
echo "Audit complete: $total_checked files checked." >&2
echo "Missing: $missing_count | Empty: $empty_count" >&2
echo "--------------------------------------------------" >&2

# Output the generated bash script header and execution commands to stdout
echo "#!/bin/bash"
echo "# Auto-generated execution script for missing HyPhy RELAX runs"
echo ""

if [[ ${#MISSING_COMMANDS[@]} -eq 0 ]]; then
    echo "# All expected files exist and are non-empty. No runs needed."
else
    for key in $(printf '%s\n' "${!MISSING_COMMANDS[@]}" | sort); do
        echo "${MISSING_COMMANDS[$key]}"
    done
fi
