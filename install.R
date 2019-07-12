packages <- c("dplyr", "reshape2", "ggplot2", "HiddenMarkov",
              "seqHMM", "doParallel", "RColorBrewer", "BiocManager", "bamboo",
              "reticulate", "DECIPHER")

for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

biomartPackages <- c("biomaRt","org.Hs.eg.db", "org.Dm.eg.db",
                     "org.Mm.eg.db", "org.Ce.eg.db", "GO.db", "clusterProfiler")
for (package in biomartPackages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}