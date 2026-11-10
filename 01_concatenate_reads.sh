#!/usr/bin/env bash
# =============================================================================
# 01_concatenate_reads.sh
# Concatenate split ONT fastq.gz files into one file per sample,
# correctly named using the sample sheet.
#
# Reads from:  sample_sheet.csv
# Writes to:   $OUTPUT_DIR/{sample}.fastq.gz  (27 files total)
#
# Excludes fastq_fail directories entirely — only fastq_pass is used.
#
# Submit on Setonix:
#   sbatch 01_concatenate_reads.sh
# =============================================================================
#SBATCH --job-name=concat_reads
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/pawsey1018/sloveys/logs/concat_%j.log
#SBATCH --error=/scratch/pawsey1018/sloveys/logs/concat_%j.err
#SBATCH --account=pawsey1018
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# PATHS — edit if your directory structure differs
# -----------------------------------------------------------------------------
BASE_DIR="/scratch/pawsey1018/sloveys/SL_Feacal_Metagenomic_230626"
SAMPLE_SHEET="/scratch/pawsey1018/sloveys/sample_sheet.csv"
OUTPUT_DIR="/scratch/pawsey1018/sloveys/raw_reads"
LOG_DIR="/scratch/pawsey1018/sloveys/logs"

# Flow cell directory names (from system check)
FC1_DIR="${BASE_DIR}/SL_Feacal_Metagenomic_Barcodes_1-13/20260623_1427_P2S-02700-A_PBM52240_be441166/fastq_pass"
FC2_DIR="${BASE_DIR}/SL_Feacal_Metagenomic_Barcodes_14-3/20260623_1427_P2S-02700-B_PBM52099_885d82de/fastq_pass"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

echo "============================================================"
echo " ONT read concatenation — Pawsey Setonix"
echo " Started: $(date)"
echo " Output:  ${OUTPUT_DIR}"
echo "============================================================"

# -----------------------------------------------------------------------------
# Read sample sheet line by line (skip header)
# Columns: sample,individual,timepoint,flowcell,flowcell_id,barcode,abx_status
# -----------------------------------------------------------------------------
TOTAL=0
SUCCESS=0
FAILED=0

while IFS=',' read -r sample individual timepoint flowcell flowcell_id barcode abx_status; do

    # Skip header line
    [[ "$sample" == "sample" ]] && continue

    echo ""
    echo "------------------------------------------------------------"
    echo "Processing: ${sample} | ${flowcell} | ${barcode}"

    # ------------------------------------------------------------------
    # Determine source directory based on flow cell
    # ------------------------------------------------------------------
    if [[ "$flowcell" == "flowcell1" ]]; then
        SOURCE_DIR="${FC1_DIR}/${barcode}"
    elif [[ "$flowcell" == "flowcell2" ]]; then
        SOURCE_DIR="${FC2_DIR}/${barcode}"
    else
        echo "  ERROR: Unknown flowcell value '${flowcell}' for sample ${sample} — skipping"
        (( FAILED++ )) || true
        (( TOTAL++ )) || true
        continue
    fi

    # ------------------------------------------------------------------
    # Check source directory exists
    # ------------------------------------------------------------------
    if [[ ! -d "${SOURCE_DIR}" ]]; then
        echo "  ERROR: Source directory not found: ${SOURCE_DIR}"
        echo "  Skipping ${sample}"
        (( FAILED++ )) || true
        (( TOTAL++ )) || true
        continue
    fi

    # ------------------------------------------------------------------
    # Count input files
    # ------------------------------------------------------------------
    N_FILES=$(find "${SOURCE_DIR}" -name "*.fastq.gz" | wc -l)
    echo "  Found ${N_FILES} split fastq.gz files in ${SOURCE_DIR}"

    if [[ "$N_FILES" -eq 0 ]]; then
        echo "  ERROR: No fastq.gz files found for ${sample} — skipping"
        (( FAILED++ )) || true
        (( TOTAL++ )) || true
        continue
    fi

    # ------------------------------------------------------------------
    # Concatenate all split files into one output file
    # Using cat on .gz files directly concatenates them correctly —
    # the resulting file is valid gzip and can be read by all tools.
    # ------------------------------------------------------------------
    OUTPUT_FILE="${OUTPUT_DIR}/${sample}.fastq.gz"

    echo "  Concatenating to: ${OUTPUT_FILE}"

    # Sort files numerically by split index before concatenating
    # (ensures read order is consistent, though order doesn't affect analysis)
    find "${SOURCE_DIR}" -name "*.fastq.gz" \
        | sort -t_ -k7 -n \
        | xargs cat > "${OUTPUT_FILE}"

    # ------------------------------------------------------------------
    # Verify output file was created and is non-empty
    # ------------------------------------------------------------------
    if [[ -s "${OUTPUT_FILE}" ]]; then
        SIZE=$(du -h "${OUTPUT_FILE}" | cut -f1)
        echo "  SUCCESS: ${sample}.fastq.gz created (${SIZE})"
        (( SUCCESS++ )) || true
    else
        echo "  ERROR: Output file is empty for ${sample}"
        rm -f "${OUTPUT_FILE}"
        (( FAILED++ )) || true
    fi

    (( TOTAL++ )) || true

done < "${SAMPLE_SHEET}"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo " Concatenation complete"
echo " Completed: $(date)"
echo " Total samples processed: ${TOTAL}"
echo " Successful:              ${SUCCESS}"
echo " Failed:                  ${FAILED}"
echo "============================================================"

# -----------------------------------------------------------------------------
# Verify all 27 expected output files exist
# -----------------------------------------------------------------------------
echo ""
echo "Output files created:"
ls -lh "${OUTPUT_DIR}"/*.fastq.gz 2>/dev/null || echo "  No output files found — check logs above"

EXPECTED=27
ACTUAL=$(ls "${OUTPUT_DIR}"/*.fastq.gz 2>/dev/null | wc -l)

if [[ "$ACTUAL" -eq "$EXPECTED" ]]; then
    echo ""
    echo "All ${EXPECTED} sample files accounted for."
else
    echo ""
    echo "WARNING: Expected ${EXPECTED} files but found ${ACTUAL}."
    echo "Check the log above for errors."
fi

# -----------------------------------------------------------------------------
# Print final file list with sizes for records
# -----------------------------------------------------------------------------
echo ""
echo "Final file inventory:"
echo "--------------------------------------------------------------"
printf "%-15s %8s\n" "SAMPLE" "SIZE"
echo "--------------------------------------------------------------"
for f in "${OUTPUT_DIR}"/*.fastq.gz; do
    fname=$(basename "$f" .fastq.gz)
    fsize=$(du -h "$f" | cut -f1)
    printf "%-15s %8s\n" "$fname" "$fsize"
done
echo "--------------------------------------------------------------"
