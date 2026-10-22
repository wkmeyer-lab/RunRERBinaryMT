#!/bin/bash

####################################################################################################################
####    Initial script written by Michael Tene and Alexander Seaver                                             ####
####    Select sections have been revised utilizing Google Gemini 3.6, manually reviewed by Alexander Seaver    ####
####################################################################################################################

# --- ARGUMENTS & SETUP ---
PROJECTNAME=$1
TOOLNAME=$2
FILEPREFIX=$3
GENENAME=$4
FOREGROUNDPATHVALUE=$5
BINARIZEVALUE=$6
BACKGROUNDVALUE=$7
RATES=$8

job_id="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-"NOT_A_SLURM_JOB"}}_${SLURM_ARRAY_TASK_ID:-"NA"}"
CURRENT_STEP="Initialization"

SHARED_HYPHY_DIR="/share/ceph/wym219group/shared/projects/seaverProjects/${PROJECTNAME}/Output/${FILEPREFIX}/Hyphy"
BASE_HYPHY_DIR="${SHARED_HYPHY_DIR}/Rate_${RATES}"

MASTER_ERROR_LOG="${BASE_HYPHY_DIR}/MASTER_JOB_FAILURE_LOG.txt"
MISSING_GENES_FILE="${BASE_HYPHY_DIR}/MISSING_ALIGNMENTS.txt"

error_trap() {
  local exit_code=$?
  if [ $exit_code -ne 0 ]; then
    echo "Job $job_id for $GENENAME failed at step: $CURRENT_STEP (Exit code: $exit_code)"
    ( flock 200; echo "Job $job_id: $GENENAME - ERROR during $CURRENT_STEP" >> "$MASTER_ERROR_LOG" ) 200>"$MASTER_ERROR_LOG.lock"
  fi
}
trap error_trap EXIT

# --- MODULE ENVIRONMENT SETUP ---

# 1. Load HyPhy and R
module load hyphy
module load r

# --- DEBUG & FIND FASTA ---
CURRENT_STEP="Finding Fasta"
ALIGNMENTPATH="/share/ceph/wym219group/shared/projects/seaverProjects/${PROJECTNAME}/Data/masked_genes"
TRANSLATEDPATH="/share/ceph/wym219group/shared/data/fixed_filtered_translated_alignments"

# Sanitize GENENAME
GENENAME=$(echo "$GENENAME" | tr -d '\r\n ')

echo "=================== DEBUG SANITY CHECK ==================="
echo "[1] Raw GENENAME:          '${GENENAME}'"
echo "[2] Check ALIGNMENTPATH:   '${ALIGNMENTPATH}'"
echo "[3] Check TRANSLATEDPATH:  '${TRANSLATEDPATH}'"

# Test lookup in translated directory
TRANSLATED_MATCH=$(ls "${TRANSLATEDPATH}"/*."${GENENAME}".aa.filt.fa 2>/dev/null | head -n 1)
echo "[4] TRANSLATED_MATCH result: '${TRANSLATED_MATCH}'"

if [ -z "$TRANSLATED_MATCH" ]; then
    echo "--> FAIL AT [4]: No file in '${TRANSLATEDPATH}' matched pattern '* .${GENENAME}.aa.filt.fa'"
    exit 1
fi

FASTAFILEBASE=$(basename "$TRANSLATED_MATCH" .aa.filt.fa)
FASTAFILE="${FASTAFILEBASE}.filt.fa"
FASTAFILEPATH="${ALIGNMENTPATH}/${FASTAFILE}"

echo "[5] Extracted FASTAFILEBASE: '${FASTAFILEBASE}'"
echo "[6] Computed FASTAFILE:      '${FASTAFILE}'"
echo "[7] Final FASTAFILEPATH:     '${FASTAFILEPATH}'"
echo "=========================================================="

if [ ! -f "$FASTAFILEPATH" ]; then
    echo "--> FAIL AT [7]: Computed path '${FASTAFILEPATH}' does not exist on disk!"
    exit 1
fi

# --- GENERATE HYPHY INPUT FILE ---
CURRENT_STEP="Generating HyPhy Input (hyphyFileMake.R)"
cd /share/ceph/wym219group/shared/projects/seaverProjects/$PROJECTNAME

# Point the bash script to the main shared directory where R saves the file
HYPHYINPUTFILENAME="${SHARED_HYPHY_DIR}/${FILEPREFIX}${GENENAME}_foreground_${FOREGROUNDPATHVALUE}_HyphyInputFile.fna"
PHENOTYPETREELOCATION="Output/${FILEPREFIX}/${FILEPREFIX}CategoricalTree.rds"

if [ -s "$HYPHYINPUTFILENAME" ] ; then
    echo "Valid Input File already exists for foreground category ${FOREGROUNDPATHVALUE} in shared directory."
else
    echo "Generating Input file for foreground category ${FOREGROUNDPATHVALUE}"
    rm -f "$HYPHYINPUTFILENAME"

    R --no-save --args \
      g=$GENENAME \
      a="$FASTAFILEPATH" \
      r=$FILEPREFIX \
      f=$FOREGROUNDPATHVALUE \
      b=$BINARIZEVALUE \
      p=$PHENOTYPETREELOCATION \
      < Src/Reu/hyphyFileMake.R

    if [ $? -ne 0 ] || [ ! -s "$HYPHYINPUTFILENAME" ]; then
      echo "ERROR: hyphyFileMake.R failed or generated an empty file for ${GENENAME}."
      exit 1
    fi
fi

# --- RUN HYPHY ---
CURRENT_STEP="Running HyPhy ($TOOLNAME)"
cd "$BASE_HYPHY_DIR"
OUTPUTFILENAME="${FILEPREFIX}-Hyphy-${TOOLNAME}-${GENENAME}-Foreground_${FOREGROUNDPATHVALUE}-Rate_${RATES}.json"

if [[ "$TOOLNAME" == "relax" ]] ; then
  BINARIZEVALUE=$(echo "$BINARIZEVALUE" | tr '[:upper:]' '[:lower:]')
  
  # Injecting the $RATES parameter into the hyphy call
  if [[ "$BINARIZEVALUE" == "t" || "$BINARIZEVALUE" == "true" ]] ; then
    hyphy $TOOLNAME --alignment "$HYPHYINPUTFILENAME" --test $FOREGROUNDPATHVALUE --reference Background --output $OUTPUTFILENAME --starting-points 10 --grid-size 250 --rates $RATES
  else
    hyphy $TOOLNAME --alignment "$HYPHYINPUTFILENAME" --test $FOREGROUNDPATHVALUE --reference $BACKGROUNDVALUE --output $OUTPUTFILENAME --rates $RATES
  fi
else
  hyphy $TOOLNAME --alignment "$HYPHYINPUTFILENAME" --branches $FOREGROUNDPATHVALUE --output $OUTPUTFILENAME
fi

if [ $? -ne 0 ]; then
  echo "ERROR: HyPhy command failed. Exiting."
  exit 1
fi
