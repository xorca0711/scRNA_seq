# Workspace-local packages only. Run with .tools/R-portable/app/bin/Rscript.exe.
args <- commandArgs(trailingOnly=TRUE)
here <- if (length(args)) normalizePath(args[[1]], winslash="/") else normalizePath("analysis/corrections/statistics", winslash="/")
lib <- file.path(here, ".tools", "R-library")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=300)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", lib=lib, type="binary")
BiocManager::install(c("limma", "edgeR"), version="3.23", lib=lib, ask=FALSE, update=FALSE, type="binary")
writeLines(capture.output(sessionInfo()), file.path(here, "R-session.txt"))
stopifnot(requireNamespace("limma", quietly=TRUE), requireNamespace("edgeR", quietly=TRUE))
expected <- c(limma="3.68.5",edgeR="4.10.5")
for (package in names(expected)) {
  if (as.character(packageVersion(package)) != expected[[package]]) {
    stop("Version differs from the recorded correction: ",package," ",packageVersion(package),
         "; retrieve the recorded binary before claiming exact reproduction")
  }
}
