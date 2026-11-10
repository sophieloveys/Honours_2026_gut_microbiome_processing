#!/usr/bin/env bash
# =============================================================================
# 03_download_databases.sh
# Downloads all reference databases required by the virome pipeline.
#
# Databases downloaded:
#   - Hostile: human T2T-CHM13 + HLA reference (for decontamination)
#   - Kraken2: Standard database (~100GB — check if shared copy exists first)
#   - CheckV: Viral genome quality database
#   - geNomad: Viral identification database
#   - VIBRANT: HMM databases
#   - Phold: Phage gene annotation database
#
# Submit on Setonix:
#   sbatch 03_download_databases.sh
#
# NOTE: This job may take 3–6 hours due to Kraken2 database size (~100GB).
# Check with Pawsey support whether a shared Kraken2 database already exists
# before running — if it does, comment out the Kraken2 section and update
# config.yaml with the shared path instead.
# =============================================================================
#SBATCH --job-name=download_dbs
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/pawsey1018/sloveys/logs/download_dbs_%j.log
#SBATCH --error=/scratch/pawsey1018/sloveys/logs/download_dbs_%j.err
#SBATCH --account=pawsey1018
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------
DB_DIR="/scratch/pawsey1018/sloveys/databases"
CONTAINER_DIR="/scratch/pawsey1018/sloveys/containers"
LOG_DIR="/scratch/pawsey1018/sloveys/logs"

mkdir -p "${DB_DIR}" "${LOG_DIR}"

export SINGULARITY_CACHEDIR="/scratch/pawsey1018/sloveys/singularity_cache"

module load singularity/4.1.0-slurm

echo "============================================================"
echo " Downloading pipeline databases"
echo " Destination: ${DB_DIR}"
echo " Started: $(date)"
echo "============================================================"

# =============================================================================
# 1. Hostile — human reference index for decontamination
# =============================================================================
echo ""
echo "[1/5] Hostile human reference (human-t2t-hla)..."

HOSTILE_DIR="${DB_DIR}/hostile"
mkdir -p "${HOSTILE_DIR}"

if [[ -f "${HOSTILE_DIR}/human-t2t-hla.1.bt2l" ]] || \
   [[ -f "${HOSTILE_DIR}/human-t2t-hla.idx" ]]; then
    echo "  Hostile reference already exists — skipping."
else
    singularity exec \
        --bind "${HOSTILE_DIR}:/data" \
        "${CONTAINER_DIR}/hostile.sif" \
        hostile fetch \
            --out-dir /data \
        2>&1
    echo "  Hostile reference downloaded."
fi

# =============================================================================
# 2. CheckV database
# =============================================================================
echo ""
echo "[2/5] CheckV database..."

CHECKV_DIR="${DB_DIR}/checkv"
mkdir -p "${CHECKV_DIR}"

if [[ -d "${CHECKV_DIR}/checkv-db-v1.5" ]]; then
    echo "  CheckV database already exists — skipping."
else
    singularity exec \
        --bind "${CHECKV_DIR}:/data" \
        "${CONTAINER_DIR}/checkv.sif" \
        checkv download_database /data \
        2>&1
    echo "  CheckV database downloaded."
fi

# =============================================================================
# 3. geNomad database
# =============================================================================
echo ""
echo "[3/5] geNomad database..."

GENOMAD_DIR="${DB_DIR}/genomad"
mkdir -p "${GENOMAD_DIR}"

if [[ -d "${GENOMAD_DIR}/genomad_db" ]]; then
    echo "  geNomad database already exists — skipping."
else
    singularity exec \
        --bind "${GENOMAD_DIR}:/data" \
        "${CONTAINER_DIR}/genomad.sif" \
        genomad download-database /data \
        2>&1
    echo "  geNomad database downloaded."
fi

# =============================================================================
# 4. VIBRANT databases
# =============================================================================
echo ""
echo "[4/5] VIBRANT databases..."

VIBRANT_DIR="${DB_DIR}/vibrant"
mkdir -p "${VIBRANT_DIR}"

if [[ -d "${VIBRANT_DIR}/databases" ]]; then
    echo "  VIBRANT databases already exist — skipping."
else
    singularity exec \
        --bind "${VIBRANT_DIR}:/data" \
        "${CONTAINER_DIR}/vibrant.sif" \
        bash -c "download-db.sh /data" \
        2>&1
    echo "  VIBRANT databases downloaded."
fi

# =============================================================================
# 5. Phold database
# =============================================================================
echo ""
echo "[5/5] Phold database..."

PHOLD_DIR="${DB_DIR}/phold"
mkdir -p "${PHOLD_DIR}"

if [[ -d "${PHOLD_DIR}/db" ]]; then
    echo "  Phold database already exists — skipping."
else
    singularity exec \
        --bind "${PHOLD_DIR}:/data" \
        "${CONTAINER_DIR}/phold.sif" \
        phold install --db_dir /data \
        2>&1
    echo "  Phold database downloaded."
fi

# =============================================================================
# 6. PhaBOX2 database
# =============================================================================
echo ""
echo "[6/7] PhaBOX2 database (host prediction, taxonomy, lifestyle)..."

PHABOX_DIR="${DB_DIR}/phabox"
mkdir -p "${PHABOX_DIR}"

if [[ -d "${PHABOX_DIR}/database" ]]; then
    echo "  PhaBOX2 database already exists — skipping."
else
    singularity exec \
        --bind "${PHABOX_DIR}:/data" \
        "${CONTAINER_DIR}/phabox.sif" \
        phabox2 download \
            --dbdir /data \
        2>&1
    echo "  PhaBOX2 database downloaded."
fi

# =============================================================================
# 7. Kraken2 standard database (~100GB)
# NOTE: Check with Pawsey if a shared copy already exists at:
#       /scratch/references/kraken2/ or similar
#       If so, comment out this entire section and update config.yaml
# =============================================================================
echo ""
echo "[6/6] Kraken2 standard database (~100GB — this will take a while)..."
echo "      If a shared Kraken2 DB exists on Pawsey, Ctrl+C now and update"
echo "      databases.kraken2_db in config.yaml instead."

KRAKEN2_DIR="${DB_DIR}/kraken2/standard"
mkdir -p "${KRAKEN2_DIR}"

if [[ -f "${KRAKEN2_DIR}/taxo.k2d" ]]; then
    echo "  Kraken2 database already exists — skipping."
else
    singularity exec \
        --bind "${KRAKEN2_DIR}:/data" \
        "${CONTAINER_DIR}/kraken2.sif" \
        kraken2-build \
            --standard \
            --db /data \
            --threads 8 \
        2>&1
    echo "  Kraken2 database built."
fi

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo " Database download complete"
echo " Completed: $(date)"
echo "============================================================"
echo ""
echo "Database inventory:"
echo "--------------------------------------------------------------"

check_dir() {
    local label="$1"
    local path="$2"
    if [[ -d "${path}" ]]; then
        local size
        size=$(du -sh "${path}" 2>/dev/null | cut -f1)
        printf "  %-20s %8s  ✓  %s\n" "${label}" "${size}" "${path}"
    else
        printf "  %-20s %8s  ✗  MISSING\n" "${label}" "—"
    fi
}

check_dir "Hostile"   "${DB_DIR}/hostile"
check_dir "CheckV"    "${DB_DIR}/checkv"
check_dir "geNomad"   "${DB_DIR}/genomad"
check_dir "VIBRANT"   "${DB_DIR}/vibrant"
check_dir "Phold"     "${DB_DIR}/phold"
check_dir "PhaBOX2"   "${DB_DIR}/phabox"
check_dir "Kraken2"   "${DB_DIR}/kraken2/standard"

echo "--------------------------------------------------------------"
echo ""
echo "NEXT STEP: Update paths in config.yaml if any database locations differ."
echo "============================================================"
