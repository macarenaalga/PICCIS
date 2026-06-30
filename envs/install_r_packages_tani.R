# Script de instalación de paquetes R para tani_env
# Se llama automáticamente desde install.sh

repos <- "https://cran.r-project.org"

# Paquetes CRAN
pkgs_cran <- c("ape", "phangorn", "MASS", "ggplot2", "reshape2", "OpenMx")

for (pkg in pkgs_cran) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        message(paste("Instalando", pkg, "..."))
        install.packages(pkg, repos=repos)
    } else {
        message(paste(pkg, "ya instalado."))
    }
}

# OpenMx desde repo oficial si falla CRAN
if (!requireNamespace("OpenMx", quietly=TRUE)) {
    install.packages("OpenMx", repos="https://openmx.ssri.psu.edu/packages/")
}

# BiocManager + paquetes Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos=repos)
}

bioc_pkgs <- c("treeio", "tidytree", "ggtree")
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        message(paste("Instalando", pkg, "desde Bioconductor..."))
        BiocManager::install(pkg, ask=FALSE, update=FALSE)
    } else {
        message(paste(pkg, "ya instalado."))
    }
}

message("✔ Todos los paquetes R de tANI instalados correctamente.")
