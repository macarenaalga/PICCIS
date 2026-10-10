#!/bin/bash
# ================================================================
#  PICCIS v1.0 - Installation script
#  Run from the repository root:
#    bash install.sh
# ================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo ""
echo "========================================================"
echo "  PICCIS v1.0 - Environment setup"
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
    # Mismos canales que los envs/*.yml: conda-forge primero, bioconda después,
    # sin 'defaults' (mezclarlo con conda-forge genera conflictos de versiones).
    conda create -n plasmidos_env -c conda-forge --override-channels python=3.10 -y
    conda run -n plasmidos_env pip install -r "$SCRIPT_DIR/requirements.txt"
    mamba install -n plasmidos_env -c conda-forge -c bioconda --override-channels \
        "blast=2.12" perl git -y
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

# PilerCR (usado por Bakta para CRISPR): el binario de bioconda puede estar
# compilado con instrucciones de CPU que procesadores más viejos no tienen.
# 'which' lo encuentra igual, así que hay que ejecutarlo para saber si sirve.
# Un código >= 128 significa que el proceso murió por una señal
# (132 = Illegal instruction, 139 = Segmentation fault).
set +e
conda run -n bakta_env pilercr -options &>/dev/null
rc=$?
set -e
if [ "$rc" -ge 128 ]; then
    echo "  [WARNING] pilercr in bakta_env crashes on this CPU (exit code $rc)."
    echo "            Bakta will fail unless it runs with --skip-crispr."
else
    ok "pilercr" "bakta_env"
fi

# MMseqs2 (usado por geNomad): el binario de bioconda está compilado con AVX2.
# En procesadores sin AVX2 se cae con "Illegal instruction" y geNomad falla en
# todas las muestras. Si pasa, se reemplaza por la compilación oficial para
# SSE4.1 de la MISMA versión (la original queda como mmseqs.avx2_original).
set +e
conda run -n genomad_env mmseqs version &>/dev/null
rc=$?
set -e
if [ "$rc" -ge 128 ]; then
    echo "  [WARNING] mmseqs in genomad_env crashes on this CPU (exit code $rc)."
    if grep -qw sse4_1 /proc/cpuinfo; then
        GENV=$(conda env list | awk '$1=="genomad_env" {print $NF}')
        MMBIN="$GENV/bin/mmseqs"
        MMVER=$(conda list -n genomad_env '^mmseqs2$' | awk '$1=="mmseqs2" {print $2}')
        MMTAG="${MMVER/./-}"                      # 18.8cc5c → 18-8cc5c (etiqueta de GitHub)
        TMPMM=$(mktemp -d)
        echo "            Installing the SSE4.1 build of MMseqs2 ${MMVER:-latest}..."
        set +e
        wget -q -O "$TMPMM/mm.tar.gz" \
            "https://github.com/soedinglab/MMseqs2/releases/download/${MMTAG}/mmseqs-linux-sse41.tar.gz" \
          || wget -q -O "$TMPMM/mm.tar.gz" "https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz"
        tar -xzf "$TMPMM/mm.tar.gz" -C "$TMPMM" \
          && [ -f "$TMPMM/mmseqs/bin/mmseqs" ] \
          && { [ -f "$MMBIN.avx2_original" ] || cp "$MMBIN" "$MMBIN.avx2_original"; } \
          && cp "$TMPMM/mmseqs/bin/mmseqs" "$MMBIN"
        conda run -n genomad_env mmseqs version &>/dev/null
        rc=$?
        set -e
        rm -rf "$TMPMM"
        if [ "$rc" -eq 0 ]; then
            ok "mmseqs (SSE4.1 build)" "genomad_env"
        else
            echo "  [MISSING] mmseqs still fails; geNomad will not work on this computer."
        fi
    else
        echo "            This CPU has neither AVX2 nor SSE4.1: geNomad cannot run here."
    fi
else
    ok "mmseqs" "genomad_env"
fi

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
