#!/usr/bin/env Rscript
# ================================================================
#  PICCIS v2.0 - Paquetes R para tani_env
#  Ejecutado por install.sh:
#    conda run -n tani_env Rscript install_r_packages_tani.R
#
#  Instala los paquetes R que NO vienen por conda en tani_env.yml
#  y que necesita tANI_tool (buildtree_w_support.R) para construir
#  y graficar el árbol filogenético.
#
#  Paquetes ya provistos por conda (no se reinstalan aquí):
#    r-mass, r-ggplot2, r-reshape2, r-gdtools, r-systemfonts, r-aplot
# ================================================================

# Mirror de CRAN (evita que pregunte de forma interactiva)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Número de cores para compilar
ncpus <- max(1L, parallel::detectCores() - 1L)

# ── Helper: instala solo si falta ───────────────────────────
install_if_missing <- function(pkg, installer) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("  → Instalando %s ...", pkg))
        installer(pkg)
    } else {
        message(sprintf("  ✓ %s ya está instalado", pkg))
    }
}

# ── 1. Paquetes de CRAN ─────────────────────────────────────
# ape, phangorn, OpenMx: requeridos por tANI, no vienen por conda.
# stringr, Matrix: requeridos por tANI; normalmente ya presentes
#   (stringr de arrastre por reshape2, Matrix por r-base), se incluyen
#   como red de seguridad por si la resolución de dependencias cambia.
message("\n[1/2] Paquetes de CRAN (ape, phangorn, OpenMx, stringr, Matrix)...")
cran_pkgs <- c("ape", "phangorn", "OpenMx", "stringr", "Matrix")
for (pkg in cran_pkgs) {
    install_if_missing(pkg, function(p)
        install.packages(p, Ncpus = ncpus, quiet = TRUE))
}

# ── 2. Paquetes de Bioconductor ─────────────────────────────
message("\n[2/2] Paquetes de Bioconductor (treeio, tidytree, ggtree)...")
install_if_missing("BiocManager", function(p)
    install.packages(p, Ncpus = ncpus, quiet = TRUE))

bioc_pkgs <- c("treeio", "tidytree", "ggtree")
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("  → Instalando %s (Bioconductor) ...", pkg))
        BiocManager::install(pkg, update = FALSE, ask = FALSE, Ncpus = ncpus)
    } else {
        message(sprintf("  ✓ %s ya está instalado", pkg))
    }
}

# ── Verificación final ──────────────────────────────────────
message("\n── Verificación ────────────────────────────────────────")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
faltan <- character(0)
for (pkg in all_pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("  [OK]      %s", pkg))
    } else {
        message(sprintf("  [MISSING] %s", pkg))
        faltan <- c(faltan, pkg)
    }
}

if (length(faltan) > 0) {
    stop(sprintf("No se pudieron instalar: %s", paste(faltan, collapse = ", ")))
} else {
    message("\n✔ Todos los paquetes R de tani_env están listos.\n")
}
