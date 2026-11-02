#!/bin/bash
# ================================================================
#  PICCIS v2.0 - Verification script
#  Ejecutar desde el directorio de PICCIS:
#    bash check_piccis.sh
# ================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONF_FILE="$SCRIPT_DIR/piccis.conf"
OK=0
WARN=0
FAIL=0

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

ok()   { echo -e "  ${GREEN}[OK]${NC}      $1"; ((OK++)); }
warn() { echo -e "  ${YELLOW}[WARN]${NC}    $1"; ((WARN++)); }
fail() { echo -e "  ${RED}[MISSING]${NC} $1"; ((FAIL++)); }

check_env_exists() {
    conda env list 2>/dev/null | grep -qE "^$1[[:space:]]"
}

check_tool_in_env() {
    conda run -n "$1" which "$2" &>/dev/null
}

echo ""
echo "========================================================"
echo "  PICCIS v2.0 - System verification"
echo "========================================================"

# ── 1. Entornos conda ────────────────────────────────────────────
echo ""
echo "── Entornos conda ───────────────────────────────────────"

ALL_ENVS=(plasmidos_env spades_env unicycler_env platon_env mob_env \
          panaroo_env bakta_env abricate_env genomad_env tani_env eggnog_env)

for env in "${ALL_ENVS[@]}"; do
    if check_env_exists "$env"; then
        ok "$env"
    elif [[ "$env" == "plasmidos_env" ]]; then
        fail "$env  →  bash install.sh"
    else
        fail "$env  →  mamba env create -f envs/${env}.yml"
    fi
done

# ── 2. Herramientas base en plasmidos_env ────────────────────────
echo ""
echo "── Herramientas base en plasmidos_env ───────────────────"

if check_env_exists "plasmidos_env"; then
    for tool in blastn perl git; do
        if check_tool_in_env "plasmidos_env" "$tool"; then
            ok "$tool  (plasmidos_env)"
        else
            fail "$tool no está en plasmidos_env  →  bash install.sh"
        fi
    done
else
    fail "plasmidos_env no existe — no se pueden comprobar sus herramientas"
fi

# ── 3. Herramientas en entornos separados ───────────────────────
echo ""
echo "── Herramientas en entornos separados ───────────────────"

declare -A TOOL_ENV=(
    ["spades.py"]="spades_env"
    ["unicycler"]="unicycler_env"
    ["platon"]="platon_env"
    ["mob_recon"]="mob_env"
    ["mob_typer"]="mob_env"
    ["panaroo"]="panaroo_env"
    ["bakta"]="bakta_env"
    ["abricate"]="abricate_env"
    ["genomad"]="genomad_env"
    ["emapper.py"]="eggnog_env"
    ["Rscript"]="tani_env"
)

for tool in "${!TOOL_ENV[@]}"; do
    env="${TOOL_ENV[$tool]}"
    if check_env_exists "$env"; then
        if check_tool_in_env "$env" "$tool"; then
            ok "$tool  ($env)"
        else
            fail "$tool no encontrado en $env  →  mamba env create -f envs/${env}.yml"
        fi
    else
        fail "$tool  →  $env no instalado"
    fi
done

# ── 4. Librerías Python en plasmidos_env ─────────────────────────
echo ""
echo "── Librerías Python en plasmidos_env ────────────────────"

if check_env_exists "plasmidos_env"; then
    if conda run -n plasmidos_env python -c \
        "import Bio, pandas, numpy, matplotlib, seaborn, scipy, sklearn" &>/dev/null; then
        ok "biopython, pandas, numpy, matplotlib, seaborn, scipy, scikit-learn"
    else
        fail "Faltan librerías Python obligatorias  →  pip install -r requirements.txt"
    fi
    if conda run -n plasmidos_env python -c \
        "import geopandas, geopy, distinctipy" &>/dev/null; then
        ok "geopandas, geopy, distinctipy  (mapa geográfico habilitado)"
    else
        warn "geopandas/geopy/distinctipy faltan  →  el gráfico planisferio se omite"
    fi
else
    fail "plasmidos_env no existe — no se pueden comprobar librerías Python"
fi

# ── 5. Paquetes R en tani_env ────────────────────────────────────
echo ""
echo "── Paquetes R en tani_env ───────────────────────────────"

if check_env_exists "tani_env"; then
    declare -A R_PKGS=(
        ["ape"]="CRAN"
        ["phangorn"]="CRAN"
        ["ggplot2"]="CRAN"
        ["reshape2"]="CRAN"
        ["ggtree"]="Bioconductor"
        ["treeio"]="Bioconductor"
    )
    for pkg in "${!R_PKGS[@]}"; do
        if conda run -n tani_env Rscript -e "library($pkg)" &>/dev/null; then
            ok "r-$pkg  (${R_PKGS[$pkg]})"
        else
            fail "r-$pkg no instalado  →  conda run -n tani_env Rscript install_r_packages_tani.R"
        fi
    done
else
    fail "tani_env no existe — no se pueden comprobar paquetes R"
fi

# ── 6. Repositorio tANI_tool ─────────────────────────────────────
echo ""
echo "── Repositorio tANI_tool ────────────────────────────────"

TANI_DIR="$SCRIPT_DIR/tANI_tool"
if [[ -f "$TANI_DIR/tANI_tool.pl" ]]; then
    ok "tANI_tool.pl  →  $TANI_DIR"
else
    fail "tANI_tool.pl no encontrado  →  git clone https://github.com/sophiagosselin/tANI_tool.git"
fi
if [[ -f "$TANI_DIR/buildtree_w_support.R" ]]; then
    ok "buildtree_w_support.R  →  $TANI_DIR"
else
    fail "buildtree_w_support.R no encontrado en $TANI_DIR"
fi

# ── 7. Bases de datos ────────────────────────────────────────────
echo ""
echo "── Bases de datos ───────────────────────────────────────"

if [[ -f "$CONF_FILE" ]]; then
    ok "piccis.conf encontrado: $CONF_FILE"
    source "$CONF_FILE"
else
    fail "piccis.conf no encontrado  →  bash install_databases.sh"
fi

# Bakta
if [[ -n "$BAKTA_DB" ]] && [[ -d "$BAKTA_DB" ]] && [[ "$(ls -A "$BAKTA_DB" 2>/dev/null)" ]]; then
    ok "Bakta DB  ($(du -sh "$BAKTA_DB" 2>/dev/null | cut -f1))  →  $BAKTA_DB"
    if [[ -d "$BAKTA_DB/amrfinderplus-db" ]]; then
        ok "AMRFinderPlus DB (requerida por Bakta)"
    else
        warn "AMRFinderPlus DB no inicializada  →  bakta puede fallar; ver install_databases.sh"
    fi
else
    fail "Bakta DB no encontrada  →  ver install_databases.sh"
fi

# Platon (mps.dmnd es el archivo clave)
if [[ -n "$PLATON_DB" ]] && [[ -f "$PLATON_DB/mps.dmnd" ]]; then
    ok "Platon DB  ($(du -sh "$PLATON_DB" 2>/dev/null | cut -f1))  →  $PLATON_DB"
else
    fail "Platon DB no encontrada (falta mps.dmnd)  →  ver install_databases.sh"
fi

# geNomad
if [[ -n "$GENOMAD_DB" ]] && [[ -d "$GENOMAD_DB" ]] && [[ "$(ls -A "$GENOMAD_DB" 2>/dev/null)" ]]; then
    ok "geNomad DB  ($(du -sh "$GENOMAD_DB" 2>/dev/null | cut -f1))  →  $GENOMAD_DB"
else
    fail "geNomad DB no encontrada  →  conda run -n genomad_env genomad download-database <dir>"
fi

# Abricate
if check_env_exists "abricate_env"; then
    if conda run -n abricate_env abricate --list 2>/dev/null | grep -q "resfinder"; then
        ok "Abricate DBs  (resfinder, card, vfdb)"
    else
        fail "Abricate DBs no configuradas  →  conda run -n abricate_env abricate --setupdb"
    fi
fi

# EggNOG (opcional → modo remoto si falta)
if [[ -n "$EGGNOG_DB" ]] && [[ -d "$EGGNOG_DB" ]] && [[ "$(ls -A "$EGGNOG_DB" 2>/dev/null)" ]]; then
    ok "EggNOG DB  ($(du -sh "$EGGNOG_DB" 2>/dev/null | cut -f1))  →  $EGGNOG_DB  (modo local)"
else
    warn "EggNOG DB no encontrada  →  modo remoto (requiere internet)"
fi

# ── 8. Archivos del pipeline ─────────────────────────────────────
echo ""
echo "── Archivos del pipeline ────────────────────────────────"

for f in piccis_pipeline.py requirements.txt install.sh \
         install_databases.sh install_r_packages_tani.R; do
    if [[ -f "$SCRIPT_DIR/$f" ]]; then
        ok "$f"
    else
        fail "$f no encontrado en $SCRIPT_DIR"
    fi
done

# plasmidos_env se crea desde install.sh, no desde un .yml
for yml in "${ALL_ENVS[@]}"; do
    [[ "$yml" == "plasmidos_env" ]] && continue
    if [[ -f "$SCRIPT_DIR/envs/${yml}.yml" ]]; then
        ok "envs/${yml}.yml"
    else
        fail "envs/${yml}.yml no encontrado"
    fi
done

# ── Resumen ──────────────────────────────────────────────────────
echo ""
echo "========================================================"
echo -e "  ${GREEN}OK: $OK${NC}   ${YELLOW}WARN: $WARN${NC}   ${RED}MISSING: $FAIL${NC}"
echo ""

if [[ $FAIL -eq 0 && $WARN -eq 0 ]]; then
    echo -e "  ${GREEN}✔ Todo listo. Corré el pipeline con:${NC}"
    echo "    conda activate plasmidos_env"
    echo "    python piccis_pipeline.py --cores 2"
elif [[ $FAIL -eq 0 ]]; then
    echo -e "  ${YELLOW}⚠ Listo con advertencias. El pipeline funciona, revisá los WARN.${NC}"
    echo "    conda activate plasmidos_env"
    echo "    python piccis_pipeline.py --cores 2"
else
    echo -e "  ${RED}✖ Resolvé los items MISSING antes de correr el pipeline.${NC}"
    echo "    Mirá los comandos de instalación de arriba."
fi
echo "========================================================"
echo ""

[[ $FAIL -eq 0 ]] && exit 0 || exit 1
