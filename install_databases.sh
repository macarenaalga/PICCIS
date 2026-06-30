#!/bin/bash
# ================================================================
#  PICCIS v2.0 - Database setup
#  Run once after install.sh:
#    conda activate plasmidos_env
#    bash install_databases.sh
#
#  - Skips any database that already exists
#  - Saves all paths to piccis.conf so the pipeline never asks again
# ================================================================

set -e

DB_DIR="$HOME/databases/piccis"
CONF_FILE="$(cd "$(dirname "$0")" && pwd)/piccis.conf"

echo ""
echo "========================================================"
echo "  PICCIS v2.0 - Database setup"
echo "  All databases will be saved to: $DB_DIR"
echo "========================================================"

mkdir -p "$DB_DIR"

# ── Helper ────────────────────────────────────────────────────
already_exists() {
    # Returns 0 (true) if $1 path exists and is non-empty
    [[ -e "$1" ]] && [[ "$(ls -A "$1" 2>/dev/null)" ]]
}

# ── 1. Abricate ──────────────────────────────────────────────
echo ""
echo "[1/6] Abricate databases (resfinder, card, vfdb)..."
if conda run -n abricate_env abricate --list 2>/dev/null | grep -q "resfinder"; then
    echo "      Already set up, skipping."
else
    conda run -n abricate_env abricate --setupdb
    echo "      Done."
fi

# ── 2. Bakta ─────────────────────────────────────────────────
echo ""
echo "[2/6] Bakta database (~1.3 GB, light version)..."
BAKTA_DB="$DB_DIR/bakta_db"
BAKTA_DB_LIGHT="$BAKTA_DB/db-light"
if already_exists "$BAKTA_DB_LIGHT"; then
    echo "      Already exists at $BAKTA_DB_LIGHT, skipping."
else
    mkdir -p "$BAKTA_DB"
    echo "      Downloading via bakta_db (compatible version auto-detected)..."
    conda run -n bakta_env bakta_db download \
        --output "$BAKTA_DB" \
        --type light
    # bakta_db download deja el tarball sin extraer si falla la extracción interna
    if [[ -f "$BAKTA_DB/db-light.tar.xz" ]] && [[ ! -d "$BAKTA_DB_LIGHT" ]]; then
        echo "      Extrayendo manualmente..."
        tar -xJf "$BAKTA_DB/db-light.tar.xz" -C "$BAKTA_DB/"
        rm "$BAKTA_DB/db-light.tar.xz"
    fi
    echo "      Done. Path: $BAKTA_DB_LIGHT"
fi

# Inicializar AMRFinderPlus (requerido por Bakta, solo una vez)
echo ""
echo "      Initializing AMRFinderPlus database (required by Bakta)..."
if already_exists "$BAKTA_DB_LIGHT/amrfinderplus-db"; then
    echo "      AMRFinderPlus DB already initialized, skipping."
else
    conda run -n bakta_env amrfinder_update \
        --force_update \
        --database "$BAKTA_DB_LIGHT/amrfinderplus-db"
    echo "      AMRFinderPlus DB initialized."
fi

# ── 3. Platon ────────────────────────────────────────────────
echo ""
echo "[3/6] Platon database (~1.8 GB)..."
PLATON_DB="$DB_DIR/platon_db"

if already_exists "$PLATON_DB"; then
    echo "      Already exists at $PLATON_DB, skipping."
else
    mkdir -p "$PLATON_DB"

    # Try the native platon-db command first (cleaner)
    if conda run -n platon_env platon-db --action download \
            --db "$PLATON_DB" 2>/dev/null; then
        echo "      Downloaded with platon-db command."

    # Fallback: direct download from Zenodo
    else
        echo "      Falling back to Zenodo download..."
        wget --show-progress -q \
            "https://zenodo.org/record/4066768/files/db.tar.gz" \
            -O "$PLATON_DB/db.tar.gz"
        tar -xzf "$PLATON_DB/db.tar.gz" -C "$PLATON_DB/"
        rm "$PLATON_DB/db.tar.gz"

        # Zenodo extracts into a subfolder called 'db'
        # Move contents up if needed
        if [[ -d "$PLATON_DB/db" ]]; then
            mv "$PLATON_DB/db"/* "$PLATON_DB/"
            rmdir "$PLATON_DB/db"
        fi
    fi
    echo "      Done. Path: $PLATON_DB"
fi

# ── 4. geNomad ───────────────────────────────────────────────
# Reemplaza a PlasFlow en la v2.0. 'genomad download-database <dir>'
# crea <dir>/genomad_db con todos los archivos del modelo.
echo ""
echo "[4/6] geNomad database (~1.6 GB)..."
GENOMAD_DB_PARENT="$DB_DIR/genomad_db"
GENOMAD_DB="$GENOMAD_DB_PARENT/genomad_db"
if already_exists "$GENOMAD_DB"; then
    echo "      Already exists at $GENOMAD_DB, skipping."
else
    mkdir -p "$GENOMAD_DB_PARENT"
    conda run -n genomad_env genomad download-database "$GENOMAD_DB_PARENT"
    echo "      Done. Path: $GENOMAD_DB"
fi

# ── 5. PlasmidFinder ─────────────────────────────────────────
echo ""
echo "[5/6] PlasmidFinder database..."
if conda run -n plasmidos_env plasmidfinder.py --help 2>&1 | grep -q "database"; then
    # PlasmidFinder downloads its DB automatically on first run
    # We trigger it here so it's ready
    TMP_FA=$(mktemp /tmp/test_XXXXXX.fasta)
    echo ">test" > "$TMP_FA"
    echo "ATCGATCGATCG" >> "$TMP_FA"
    conda run -n plasmidos_env plasmidfinder.py \
        -i "$TMP_FA" -o /tmp/pf_test 2>/dev/null || true
    rm -f "$TMP_FA"
    rm -rf /tmp/pf_test
    echo "      PlasmidFinder database initialized."
fi

# ── 6. EggNOG-mapper (optional) ──────────────────────────────
echo ""
echo "[6/6] EggNOG-mapper database..."
echo "      Opciones:"
echo "        1) Completa  (~14 GB) — eggnog.db + eggnog_proteins.dmnd (recomendado)"
echo "        2) Solo anotación (~6 GB) — solo eggnog.db, sin búsqueda diamond"
echo "        3) Omitir             — EggNOG se saltea (sin anotación COG/GO)"
echo ""
EGGNOG_DB="$DB_DIR/eggnog_db"
if already_exists "$EGGNOG_DB"; then
    echo "      Ya existe en $EGGNOG_DB, saltando."
else
    read -p "      Opción [1/2/3]: " egg_opt
    if [[ "$egg_opt" == "1" ]]; then
        mkdir -p "$EGGNOG_DB"
        conda run -n eggnog_env download_eggnog_data.py -y \
            --data_dir "$EGGNOG_DB"
        echo "      Done. Path: $EGGNOG_DB"
    elif [[ "$egg_opt" == "2" ]]; then
        mkdir -p "$EGGNOG_DB"
        conda run -n eggnog_env download_eggnog_data.py -y \
            --data_dir "$EGGNOG_DB" -D
        echo "      Done. Path: $EGGNOG_DB"
    else
        echo "      Omitido. EggNOG se saltea en el pipeline."
        EGGNOG_DB=""
    fi
fi

# ── Guardar rutas en piccis.conf ─────────────────────────────
echo ""
echo "Saving database paths to: $CONF_FILE"

cat > "$CONF_FILE" << EOF
# PICCIS v2.0 - Database configuration
# Generated automatically by install_databases.sh
# Edit this file if you move your databases.

BAKTA_DB=$BAKTA_DB_LIGHT
PLATON_DB=$PLATON_DB
GENOMAD_DB=$GENOMAD_DB
EGGNOG_DB=$EGGNOG_DB
EOF

echo ""
echo "========================================================"
echo "  Setup complete. Database paths:"
echo ""
echo "  Bakta DB    : $BAKTA_DB_LIGHT"
echo "  Platon DB   : $PLATON_DB"
echo "  geNomad DB  : $GENOMAD_DB"
[[ -n "$EGGNOG_DB" ]] && echo "  EggNOG DB   : $EGGNOG_DB"
echo ""
echo "  Paths saved to piccis.conf"
echo "  The pipeline will read them automatically."
echo "========================================================"
