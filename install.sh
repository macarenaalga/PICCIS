#!/bin/bash
# ================================================================
#  PICCIS v2.0 - Installation script
#  Run from the repository root:
#    bash install.sh
# ================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo ""
echo "========================================================"
echo "  PICCIS v2.0 - Environment setup"
echo "========================================================"
echo ""

# ── Check mamba ──────────────────────────────────────────────
if ! command -v mamba &> /dev/null; then
    echo "mamba not found. Installing into base..."
    conda install -n base -c conda-forge mamba -y
fi

# ── 1. plasmidos_env ─────────────────────────────────────────
# Entorno principal: desde aquí se ejecuta el pipeline.
echo ""
echo "[1/11] Creating plasmidos_env (Python 3.10)..."
if conda env list | grep -qE "^plasmidos_env[[:space:]]"; then
    echo "       Already exists, skipping."
else
    conda create -n plasmidos_env python=3.10 -y
    conda run -n plasmidos_env pip install -r "$SCRIPT_DIR/requirements.txt"
    mamba install -n plasmidos_env -c bioconda -c conda-forge \
        "blast=2.12" plasmidfinder perl git -y
fi
echo "[1/11] Done."

# ── 2-11. Entornos separados ──────────────────────────────────
# NOTA v2.0: plasflow_env fue reemplazado por genomad_env.
ENVS=(spades_env unicycler_env platon_env mob_env panaroo_env
      bakta_env abricate_env genomad_env tani_env eggnog_env)
N=2
for env in "${ENVS[@]}"; do
    echo ""
    echo "[$N/11] Creating $env..."
    if conda env list | grep -qE "^${env}[[:space:]]"; then
        echo "       Already exists, skipping."
    else
        mamba env create -f "$SCRIPT_DIR/envs/${env}.yml"
    fi
    echo "[$N/11] Done."
    ((N++))
done

# ── Verificación rápida ──────────────────────────────────────
echo ""
echo "========================================================"
echo "  Verifying installations..."
echo ""

ok()   { echo "  [OK]      $1 ($2)"; }
fail() { echo "  [MISSING] $1 in $2 — install manually"; }

check_env() {
    conda run -n "$1" which "$2" &>/dev/null && ok "$2" "$1" || fail "$2" "$1"
}

check_env plasmidos_env  blastn
check_env plasmidos_env  plasmidfinder.py
check_env plasmidos_env  perl
check_env plasmidos_env  git
check_env spades_env     spades.py
check_env unicycler_env  unicycler
check_env platon_env     platon
check_env mob_env        mob_recon
check_env mob_env        mob_typer
check_env panaroo_env    panaroo
check_env bakta_env      bakta
check_env abricate_env   abricate
check_env genomad_env    genomad
check_env tani_env       Rscript
check_env eggnog_env     emapper.py

# ── Instalar paquetes R para tani_env ────────────────────────
echo ""
echo "[ tani_env ] Instalando paquetes R (ape, phangorn, ggtree, OpenMx...)"
conda run -n tani_env Rscript "$SCRIPT_DIR/install_r_packages_tani.R"

# ── Clonar tANI_tool si no existe ────────────────────────────
TANI_DIR="$SCRIPT_DIR/tANI_tool"
if [ ! -f "$TANI_DIR/tANI_tool.pl" ]; then
    echo "[ tANI ] Clonando repositorio tANI_tool..."
    git clone https://github.com/sophiagosselin/tANI_tool.git "$TANI_DIR"
    echo "✔  tANI_tool → $TANI_DIR"
else
    echo "[ tANI ] tANI_tool ya existe → $TANI_DIR"
fi

echo ""
echo "========================================================"
echo "  Installation complete. Next step:"
echo "    conda activate plasmidos_env"
echo "    bash install_databases.sh"
echo "========================================================"
