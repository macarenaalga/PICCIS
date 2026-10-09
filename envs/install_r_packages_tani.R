# Script de instalación de paquetes R para tani_env
# Se llama automáticamente desde install.sh
#
# Los paquetes de CRAN vienen precompilados desde tani_env.yml (conda-forge);
# este script solo instala lo que falte y los paquetes de Bioconductor.
# Al final verifica todo y termina con error si falta algo, así install.sh
# no sigue como si la instalación hubiera funcionado.

repos <- "https://cloud.r-project.org"
options(repos = c(CRAN = repos))
ncpus <- max(1L, parallel::detectCores() - 1L)

instalar_si_falta <- function(pkg, instalar) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        message("  ✓ ", pkg, " ya está instalado")
    } else {
        message("  → Instalando ", pkg, " ...")
        try(instalar(pkg), silent = FALSE)
    }
}

# [1/2] CRAN
message("\n[1/2] Paquetes de CRAN...")
pkgs_cran <- c("ape", "phangorn", "MASS", "ggplot2", "reshape2",
               "OpenMx", "stringr", "Matrix", "codetools")
for (pkg in pkgs_cran) {
    instalar_si_falta(pkg, function(p) install.packages(p, Ncpus = ncpus, quiet = TRUE))
}
# OpenMx desde el repo oficial si falló CRAN
if (!requireNamespace("OpenMx", quietly = TRUE)) {
    message("  → Reintentando OpenMx desde su repositorio oficial ...")
    try(install.packages("OpenMx", repos = "https://openmx.ssri.psu.edu/packages/"))
}

# [2/2] Bioconductor
message("\n[2/2] Paquetes de Bioconductor (treeio, tidytree, ggtree)...")
instalar_si_falta("BiocManager", function(p) install.packages(p, quiet = TRUE))
bioc_pkgs <- c("treeio", "tidytree", "ggtree")
for (pkg in bioc_pkgs) {
    instalar_si_falta(pkg, function(p)
        BiocManager::install(p, ask = FALSE, update = FALSE, Ncpus = ncpus))
}

# Verificación
message("\n── Verificación ────────────────────────────────────────")
todos <- c(pkgs_cran, bioc_pkgs)
faltan <- character(0)
for (pkg in todos) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("  [OK]      %s", pkg))
    } else {
        message(sprintf("  [MISSING] %s", pkg))
        faltan <- c(faltan, pkg)
    }
}

if (length(faltan) > 0) {
    message("\n✖ Faltan paquetes de R: ", paste(faltan, collapse = ", "))
    message("  Revisá los errores de arriba y volvé a correr:")
    message("  conda run -n tani_env Rscript install_r_packages_tani.R")
    quit(status = 1)
}
message("\n✔ Todos los paquetes R de tANI están instalados.")
