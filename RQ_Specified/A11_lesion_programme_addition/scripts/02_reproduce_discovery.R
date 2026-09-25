# Must pass before any new Kim module scoring. Reuses the original TMM scope.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
root <- normalizePath(args[1],winslash='/'); data_root <- normalizePath(args[2],winslash='/')
.libPaths(c(file.path(data_root,'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR))
paper <- 'Research Article/gate2_C3_yu_lee_choi_min_2026'
out <- file.path(root,'RQ_Specified/A11_lesion_programme_addition/tables/test_v2')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
target <- file.path(out,'instrument_check.tsv')
if(file.exists(target)) stop('Refusing to overwrite discovery reproduction')
m <- read.csv(file.path(root,paper,'trials/u5_human_niche/unc20_pooled_triad_units.csv'))
x <- as.matrix(read.csv(gzfile(file.path(data_root,paper,'cache/u5_human_niche/unc20_pooled_triad_counts.csv.gz')),
                      row.names=1,check.names=FALSE))
stopifnot(identical(colnames(x),m$unit_id),all(x>=0),all(x==round(x)),all(colSums(x)==m$full_library_sum))
genes <- read.delim(file.path(root,paper,'trials/u6_human_specificity/human_module_genes.tsv'))
module <- 'HPCS_without_ADI_or_operational_markers'
gs <- intersect(genes$gene[genes$module==module],rownames(x))
expected <- read.csv(file.path(root,paper,'trials/u6_human_specificity/paired_program_values.csv'))
rows <- list()
for(label in c('__broad__','AT2')) {
  take <- m$comp=='AT2' & m$label==label & m$cells>=50
  mm <- m[take,,drop=FALSE]
  y <- calcNormFactors(DGEList(counts=x[,take,drop=FALSE]),method='TMM')
  scores <- colMeans(cpm(y,log=TRUE,prior.count=1)[gs,,drop=FALSE])
  e <- expected[expected$module==module & expected$config=='unc20_pooled' &
     expected$label==label & expected$cell_floor==50 & expected$case=='LUAD' & expected$reference=='normal',]
  stopifnot(nrow(e)==23, !anyDuplicated(e$patient))
  for(j in seq_len(nrow(e))) {
    a <- which(mm$patient==e$patient[j] & mm$histology=='LUAD')
    b <- which(mm$patient==e$patient[j] & mm$histology=='normal')
    stopifnot(length(a)==1,length(b)==1)
    got <- unname(scores[a]-scores[b])
    rows[[length(rows)+1]] <- data.frame(label,patient=e$patient[j],expected=e$difference[j],
      reproduced=got,abs_error=abs(got-e$difference[j]),normalization_units=nrow(mm),assayed_genes=length(gs))
  }
}
z <- do.call(rbind,rows)
write.table(z,target,sep='\t',row.names=FALSE,quote=FALSE)
writeLines(trimws(capture.output(sessionInfo()),which='right'),file.path(out,'instrument_R_session.txt'))
stopifnot(nrow(z)==46, all(z$abs_error<=1e-6))
cat('PASS discovery reproduction:',nrow(z),'pairs, max error',max(z$abs_error),'\n')
