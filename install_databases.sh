#!/bin/bash
# ================================================================
#  PICCIS v1.0 - Database setup
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

bakta_db_ok() {
    # Returns 0 (true) only if the Bakta DB looks complete.
    # A non-empty folder is not enough: an interrupted download or
    # extraction leaves a partial folder that would be skipped forever.
    local d="$1" f stem tipo s
    [[ -s "$d/version.json" ]] || return 1
    # Índices HMMER obligatorios (los 4 archivos de cada uno)
    for base in antifam pfam; do
        for ext in h3f h3i h3m h3p; do
            [[ -s "$d/$base.$ext" ]] || return 1
        done
    done
    # Cualquier otro índice HMMER (.h3?) o Infernal (.i1?) presente
    # también debe tener sus 4 archivos
    for f in "$d"/*.h3? "$d"/*.i1?; do
        [[ -e "$f" ]] || continue
        stem="${f%.*}"; tipo="${f##*.}"; tipo="${tipo:0:2}"
        for s in f i m p; do
            [[ -s "$stem.$tipo$s" ]] || return 1
        done
    done
    return 0
}

# ── 1. Abricate ──────────────────────────────────────────────
echo ""
echo "[1/5] Abricate databases (resfinder, card, vfdb, plasmidfinder)..."
if conda run -n abricate_env abricate --list 2>/dev/null | grep -q "resfinder"; then
    echo "      Already set up, skipping."
else
    conda run -n abricate_env abricate --setupdb
    echo "      Done."
fi

# ── 2. Bakta ─────────────────────────────────────────────────
echo ""
echo "[2/5] Bakta database (~1.3 GB, light version)..."
BAKTA_DB="$DB_DIR/bakta_db"
BAKTA_DB_LIGHT="$BAKTA_DB/db-light"
if bakta_db_ok "$BAKTA_DB_LIGHT"; then
    echo "      Already exists and is complete at $BAKTA_DB_LIGHT, skipping."
else
    if [[ -e "$BAKTA_DB_LIGHT" ]]; then
        echo "      ⚠  Incomplete Bakta DB found at $BAKTA_DB_LIGHT"
        echo "         (interrupted download/extraction). Removing and downloading again..."
        rm -rf "$BAKTA_DB_LIGHT"
    fi
    rm -f "$BAKTA_DB/db-light.tar.xz"          # tarball parcial de un intento anterior
    mkdir -p "$BAKTA_DB"
    echo "      Downloading via bakta_db (compatible version auto-detected)..."
    # Sin '|| true' set -e cortaría el script antes del rescate manual de abajo
    conda run -n bakta_env bakta_db download \
        --output "$BAKTA_DB" \
        --type light || echo "      ⚠  bakta_db download returned an error, checking files..."
    # bakta_db download deja el tarball sin extraer si falla la extracción interna
    if [[ -f "$BAKTA_DB/db-light.tar.xz" ]] && ! bakta_db_ok "$BAKTA_DB_LIGHT"; then
        echo "      Extrayendo manualmente..."
        rm -rf "$BAKTA_DB_LIGHT"
        tar -xJf "$BAKTA_DB/db-light.tar.xz" -C "$BAKTA_DB/"
        rm "$BAKTA_DB/db-light.tar.xz"
    fi
    if ! bakta_db_ok "$BAKTA_DB_LIGHT"; then
        echo "      ✖  The Bakta DB is still incomplete at $BAKTA_DB_LIGHT."
        echo "         Check disk space and internet connection, then run this script again."
        exit 1
    fi
    echo "      Done. Path: $BAKTA_DB_LIGHT"
fi

# Inicializar AMRFinderPlus (requerido por Bakta, solo una vez).
# 'latest' es el enlace que AMRFinder usa para encontrar la versión
# actual; si falta, la base no sirve aunque la carpeta tenga archivos.
echo ""
echo "      Initializing AMRFinderPlus database (required by Bakta)..."
if [[ -e "$BAKTA_DB_LIGHT/amrfinderplus-db/latest" ]]; then
    echo "      AMRFinderPlus DB already initialized, skipping."
else
    conda run -n bakta_env amrfinder_update \
        --force_update \
        --database "$BAKTA_DB_LIGHT/amrfinderplus-db"
    if [[ ! -e "$BAKTA_DB_LIGHT/amrfinderplus-db/latest" ]]; then
        echo "      ✖  AMRFinderPlus DB could not be initialized."
        exit 1
    fi
    echo "      AMRFinderPlus DB initialized."
fi

# ── 3. Platon ────────────────────────────────────────────────
echo ""
echo "[3/5] Platon database (~1.8 GB)..."
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
echo "[4/5] geNomad database (~1.6 GB)..."
GENOMAD_DB_PARENT="$DB_DIR/genomad_db"
GENOMAD_DB="$GENOMAD_DB_PARENT/genomad_db"
if already_exists "$GENOMAD_DB"; then
    echo "      Already exists at $GENOMAD_DB, skipping."
else
    mkdir -p "$GENOMAD_DB_PARENT"
    conda run -n genomad_env genomad download-database "$GENOMAD_DB_PARENT"
    echo "      Done. Path: $GENOMAD_DB"
fi

# ── 5. EggNOG-mapper (optional) ──────────────────────────────
# EggNOG-mapper no tiene modo en línea: siempre necesita la base local
# completa (eggnog.db para anotar, eggnog.taxa.db para la taxonomía y
# eggnog_proteins.dmnd para la búsqueda con diamond).
# download_eggnog_data.py NO se usa: en versiones anteriores a la 2.1.13
# apunta a eggnogdb.embl.de, dominio que ya no existe. Los archivos se
# bajan directo de eggnog5.embl.de, en la versión que pide el programa.
echo ""
echo "[5/5] EggNOG-mapper database..."
EGGNOG_DB="$DB_DIR/eggnog_db"
# La versión de la base se toma del EggNOG-mapper instalado, así siempre
# coincide con el programa aunque se instale la última versión.
EGG_DB_VER=$(conda run -n eggnog_env python -c \
    "from eggnogmapper.version import __DB_VERSION__; print(__DB_VERSION__)" 2>/dev/null | tr -d '[:space:]')
if [[ -z "$EGG_DB_VER" ]]; then
    EGG_DB_VER="5.0.2"
    echo "      ⚠  Could not read the DB version from eggnog-mapper; using $EGG_DB_VER."
fi
EGG_URL="http://eggnog5.embl.de/download/emapperdb-${EGG_DB_VER}"
EGG_MIN_GB=60          # espacio libre recomendado (descarga + descompresión)

eggnog_db_ok() {
    [[ -s "$1/eggnog.db" ]] && [[ -s "$1/eggnog.taxa.db" ]] && [[ -s "$1/eggnog_proteins.dmnd" ]]
}

# Baja un archivo .gz / .tar.gz y lo descomprime de forma segura:
#  - wget -c retoma descargas cortadas
#  - gzip -t verifica el archivo antes de descomprimir
#  - se descomprime a un temporal y se renombra al final, así un corte
#    a mitad de camino nunca deja un archivo final incompleto
egg_get() {
    local gz="$1" final="$2"
    [[ -s "$final" ]] && return 0
    echo "      → $gz"
    wget -c -q --show-progress "$EGG_URL/$gz" || { echo "      ✖  Download failed: $gz"; return 1; }
    if ! gzip -t "$gz" 2>/dev/null; then
        echo "      ✖  $gz is corrupted; removed. Run the script again to re-download it."
        rm -f "$gz"; return 1
    fi
    if [[ "$gz" == *.tar.gz ]]; then
        rm -rf .extract_tmp && mkdir .extract_tmp
        tar -zxf "$gz" -C .extract_tmp && mv .extract_tmp/* . && rm -rf .extract_tmp "$gz" \
            || { rm -rf .extract_tmp; echo "      ✖  Could not extract $gz"; return 1; }
    else
        gunzip -c "$gz" > "$final.part" && mv "$final.part" "$final" && rm -f "$gz" \
            || { rm -f "$final.part"; echo "      ✖  Could not decompress $gz (disk full?)"; return 1; }
    fi
}

if eggnog_db_ok "$EGGNOG_DB"; then
    echo "      Already exists and is complete at $EGGNOG_DB, skipping."
else
    if [[ -e "$EGGNOG_DB" ]]; then
        echo "      ⚠  Incomplete EggNOG DB at $EGGNOG_DB"
        echo "         (missing eggnog.db, eggnog.taxa.db or eggnog_proteins.dmnd)."
    fi
    echo "      Opciones:"
    echo "        1) Descargar/completar (necesita ~${EGG_MIN_GB} GB libres) — eggnog.db + taxa + diamond"
    echo "        2) Omitir — EggNOG se saltea en el pipeline (sin anotación COG/GO)"
    echo ""
    read -p "      Opción [1/2]: " egg_opt
    if [[ "$egg_opt" == "1" ]]; then
        mkdir -p "$EGGNOG_DB"
        libre_gb=$(df -BG --output=avail "$EGGNOG_DB" | tail -1 | tr -dc '0-9')
        if [[ -n "$libre_gb" && "$libre_gb" -lt "$EGG_MIN_GB" ]]; then
            echo "      ⚠  Only ${libre_gb} GB free in $(dirname "$EGGNOG_DB"); ~${EGG_MIN_GB} GB recommended."
            read -p "         Continue anyway? [s/N]: " seguir
            [[ "$seguir" =~ ^[sSyY]$ ]] || { echo "      Cancelled. Free some space and run the script again."; exit 1; }
        fi
        (
            cd "$EGGNOG_DB"
            find . -maxdepth 1 -name '*.gz' -size 0 -delete     # restos vacíos de intentos fallidos
            ok=0
            egg_get eggnog.taxa.tar.gz     eggnog.taxa.db       || ok=1
            egg_get eggnog_proteins.dmnd.gz eggnog_proteins.dmnd || ok=1
            egg_get eggnog.db.gz           eggnog.db            || ok=1
            exit $ok
        ) || true
        if ! eggnog_db_ok "$EGGNOG_DB"; then
            echo "      ✖  The EggNOG DB is still incomplete at $EGGNOG_DB."
            echo "         Run this script again: partial downloads resume where they stopped."
            exit 1
        fi
        echo "      Done. Path: $EGGNOG_DB  (eggNOG DB $EGG_DB_VER)"
        echo "      Source: $EGG_URL (downloaded $(date +%Y-%m-%d))"
        ram_gb=$(free -g | awk '/^(Mem|Memoria):/ {print $2}')
        if [[ -n "$ram_gb" ]]; then
            echo "      ℹ  Each EggNOG process needs ~4-5 GB of RAM. This machine has ${ram_gb} GB:"
            echo "         use --cores $(( ram_gb / 5 > 0 ? ram_gb / 5 : 1 )) or fewer when running the pipeline."
        fi
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
