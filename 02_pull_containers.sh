#!/usr/bin/env bash
# =============================================================================
# 02_pull_containers.sh
# Downloads all Singularity containers needed for the virome pipeline.
# Pulls from depot.galaxyproject.org (Galaxy/BioContainers mirror) —
# faster and more reliable than pulling directly from quay.io on HPC systems.
#
# Submit on Setonix:
#   sbatch 02_pull_containers.sh
#
# Expected runtime: 1–2 hours depending on download speeds.
# All containers are saved to $CONTAINER_DIR as .sif files.
# =============================================================================
#SBATCH --job-name=pull_containers
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/pawsey1018/sloveys/logs/pull_containers_%j.log
#SBATCH --error=/scratch/pawsey1018/sloveys/logs/pull_containers_%j.err
#SBATCH --account=pawsey1018
# =============================================================================

set -uo pipefail
# Note: deliberately NOT using set -e here so a single failed container pull
# does not abort the entire job — failures are reported in the summary instead.

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------
CONTAINER_DIR="/scratch/pawsey1018/sloveys/containers"
CACHE_DIR="/scratch/pawsey1018/sloveys/singularity_cache"
LOG_DIR="/scratch/pawsey1018/sloveys/logs"

mkdir -p "${CONTAINER_DIR}" "${CACHE_DIR}" "${LOG_DIR}"

# Point Singularity cache to scratch (home quota is too small)
export SINGULARITY_CACHEDIR="${CACHE_DIR}"
export SINGULARITY_TMPDIR="${CACHE_DIR}/tmp"
mkdir -p "${SINGULARITY_TMPDIR}"

# Load Singularity module
module load singularity/4.1.0-slurm

echo "============================================================"
echo " Pulling Singularity containers for virome pipeline"
echo " Destination: ${CONTAINER_DIR}"
echo " Started: $(date)"
echo "============================================================"

# -----------------------------------------------------------------------------
# Helper function — pulls a container only if not already present
# Usage: pull_container <output_name> <source_url>
# -----------------------------------------------------------------------------
pull_container() {
    local name="$1"
    local url="$2"
    local outfile="${CONTAINER_DIR}/${name}.sif"

    if [[ -f "${outfile}" ]]; then
        echo "  [SKIP] ${name}.sif already exists"
        return 0
    fi

    echo ""
    echo "  Pulling: ${name}"
    echo "  Source:  ${url}"

    singularity pull \
        --dir "${CONTAINER_DIR}" \
        --name "${name}.sif" \
        "${url}" \
        2>&1

    if [[ -f "${outfile}" ]]; then
        local size
        size=$(du -h "${outfile}" | cut -f1)
        echo "  Done: ${name}.sif (${size})"
    else
        echo "  ERROR: Failed to pull ${name} — check tag and retry"
        # Do NOT return 1 — continue pulling remaining containers
    fi
}

# =============================================================================
# CONTAINERS — grouped by pipeline stage
# All pulled from Galaxy BioContainers depot (Singularity-native format)
# =============================================================================

echo ""
echo "--- Stage 1: QC and preprocessing ---"

pull_container "nanoplot" \
    "docker://quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0"

pull_container "chopper" \
    "docker://quay.io/biocontainers/chopper:0.9.0--hdcf5f25_1"

pull_container "hostile" \
    "docker://quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0"

echo ""
echo "--- Stage 2: Assembly ---"

pull_container "flye" \
    "docker://quay.io/biocontainers/flye:2.9.5--py39hdf45acc_1"

echo ""
echo "--- Stage 3: Bacterial profiling ---"

pull_container "kraken2" \
    "docker://quay.io/biocontainers/kraken2:2.1.3--pl5321h9f5acd7_0"

pull_container "bracken" \
    "docker://quay.io/biocontainers/bracken:2.9--py39hc16433a_0"

echo ""
echo "--- Stage 4: Viral identification ---"

pull_container "vibrant" \
    "docker://quay.io/biocontainers/vibrant:1.2.1--hdfd78af_1"

pull_container "genomad" \
    "docker://quay.io/biocontainers/genomad:1.8.0--pyhdfd78af_0"

echo ""
echo "--- Stage 5: Quality assessment ---"

pull_container "checkv" \
    "docker://quay.io/biocontainers/checkv:1.0.3--pyhdfd78af_0"

echo ""
echo "--- Stage 6: Clustering ---"

pull_container "cdhit" \
    "docker://quay.io/biocontainers/cd-hit:4.8.1--h5b5514e_7"

echo ""
echo "--- Stage 7: Read mapping and abundance ---"

pull_container "minimap2" \
    "docker://quay.io/biocontainers/minimap2:2.28--he4a0461_0"

pull_container "samtools" \
    "docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0"

pull_container "coverm" \
    "docker://quay.io/biocontainers/coverm:0.7.0--h9ee0642_1"

echo ""
echo "--- Stage 8: Functional annotation ---"

pull_container "phold" \
    "docker://quay.io/biocontainers/phold:0.1.17--pyhdfd78af_0"

echo ""
echo "--- Stage 9: Phage-host prediction ---"

pull_container "cherry" \
    "docker://quay.io/biocontainers/cherry:1.0--pyhdfd78af_0"

echo ""
echo "--- Stage 10: Python utilities (custom scripts) ---"

pull_container "python_biopython" \
    "docker://quay.io/biocontainers/python:3.11--1"

# =============================================================================
# VERIFICATION — list all containers and their sizes
# =============================================================================
echo ""
echo "============================================================"
echo " Container pull complete"
echo " Completed: $(date)"
echo "============================================================"
echo ""
echo "Containers available:"
echo "--------------------------------------------------------------"
printf "%-30s %8s\n" "CONTAINER" "SIZE"
echo "--------------------------------------------------------------"

EXPECTED_CONTAINERS=(
    "nanoplot"
    "chopper"
    "hostile"
    "flye"
    "kraken2"
    "bracken"
    "vibrant"
    "genomad"
    "checkv"
    "cdhit"
    "minimap2"
    "samtools"
    "coverm"
    "phold"
    "cherry"
    "python_biopython"
)

MISSING=0
for name in "${EXPECTED_CONTAINERS[@]}"; do
    sif="${CONTAINER_DIR}/${name}.sif"
    if [[ -f "${sif}" ]]; then
        size=$(du -h "${sif}" | cut -f1)
        printf "%-30s %8s  ✓\n" "${name}" "${size}"
    else
        printf "%-30s %8s  ✗ MISSING\n" "${name}" "—"
        (( MISSING++ )) || true
    fi
done

echo "--------------------------------------------------------------"
if [[ "${MISSING}" -eq 0 ]]; then
    echo "All ${#EXPECTED_CONTAINERS[@]} containers present."
else
    echo "WARNING: ${MISSING} container(s) missing — check log for errors."
    echo "Re-run this script to retry failed pulls (existing files are skipped)."
fi
echo "============================================================"
