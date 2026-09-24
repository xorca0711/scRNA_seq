root <- normalizePath('.',winslash='/');paper <- file.path(root,'Thesis/gate2_C3_yu_lee_choi_min_2026')
.libPaths(c(file.path(root,'analysis/corrections/statistics/.tools/R-library'),.libPaths()));suppressPackageStartupMessages(library(edgeR))
out <- file.path(paper,'trials/u6_ipf_specificity');cache <- file.path(paper,'cache/u6_ipf_specificity')
genesets <- read.delim(file.path(paper,'trials/u6_human_specificity/human_module_genes.tsv'));sets <- split(genesets$gene,genesets$module);source_n <- tapply(genesets$source_genes,genesets$module,unique);output <- list()
for(cohort in c('GSE136831','GSE135893')) {
 m <- read.csv(file.path(out,paste0(cohort,'_units.csv')));x <- as.matrix(read.csv(gzfile(file.path(cache,paste0(cohort,'_counts.csv.gz'))),row.names=1,check.names=FALSE));stopifnot(identical(colnames(x),m$unit_id),all(colSums(x)==m$full_library_sum))
 for(floor in c(50,30,100)) {
  take <- m$cells>=floor;mm <- m[take,,drop=FALSE];xx <- x[,take,drop=FALSE];y <- suppressMessages(calcNormFactors(DGEList(counts=xx),method='TMM'))
  for(prior in if(floor==50)c(1,.5,2) else 1) {
   lc <- cpm(y,log=TRUE,prior.count=prior)
   for(module in names(sets)) {
    gs <- intersect(sets[[module]],rownames(lc));fraction <- length(gs)/source_n[[module]];if(!length(gs)) next
    output[[length(output)+1]] <- data.frame(cohort,cell_floor=floor,prior_count=prior,mm,module,assayed_source_fraction=fraction,eligible=fraction>=.7,score=colMeans(lc[gs,,drop=FALSE]))
   }
  }
 }
}
write.csv(do.call(rbind,output),file.path(out,'donor_state_scores.csv'),row.names=FALSE);capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
