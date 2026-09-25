args <- commandArgs(trailingOnly=TRUE)
cohort <- args[[1]]
root <- normalizePath('.',winslash='/')
paper <- file.path(root,'Thesis/gate2_C3_yu_lee_choi_min_2026')
.libPaths(c(file.path(root,'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR))
cache <- file.path(root,'analysis/corrections/statistics/cache')
x <- as.matrix(read.csv(gzfile(file.path(cache,paste0(cohort,'_subtypes_counts.csv.gz'))),row.names=1,check.names=FALSE))
m <- read.csv(file.path(cache,paste0(cohort,'_subtypes_units.csv')))
stopifnot(identical(colnames(x),m$unit_id),all(x>=0),all(x==round(x)))
genes <- readLines(file.path(paper,'trials/u5_ipf_compatibility/component_genes.txt'))
output <- list(); norm <- list()
for(comp in sort(unique(m$comp))) {
  for(label in c('__broad__',sort(unique(m$label[m$comp==comp])))) {
    take <- m$comp==comp & (label=='__broad__' | m$label==label)
    mm <- m[take,,drop=FALSE]; xx <- x[,take,drop=FALSE]
    if(label=='__broad__') {
      donors <- unique(mm$donor)
      aggregate <- matrix(0,nrow=nrow(x),ncol=length(donors),dimnames=list(rownames(x),donors))
      metadata <- mm[match(donors,mm$donor),,drop=FALSE]
      for(j in seq_along(donors)) {
        same <- mm$donor==donors[[j]]
        stopifnot(length(unique(mm$disease[same]))==1)
        aggregate[,j] <- rowSums(xx[,same,drop=FALSE])
        metadata$cells[[j]] <- sum(mm$cells[same])
      }
      xx <- aggregate; mm <- metadata
    }
    stopifnot(!anyDuplicated(mm$donor))
    for(floor in c(50,30,100)) {
      take <- mm$cells>=floor & mm$disease %in% c('IPF','control')
      if(sum(take)<2) next
      dm <- mm[take,,drop=FALSE]; y <- calcNormFactors(DGEList(xx[,take,drop=FALSE]),method='TMM')
      stopifnot(all(y$samples$lib.size>0),all(is.finite(y$samples$norm.factors)))
      norm[[length(norm)+1]] <- data.frame(cohort,comp,label,cell_floor=floor,donor=dm$donor,disease=dm$disease,cells=dm$cells,library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors)
      for(prior in c(1,0.5,2)) {
        z <- cpm(y,log=TRUE,prior.count=prior)
        for(g in intersect(genes,rownames(z))) {
          output[[length(output)+1]] <- data.frame(cohort,comp,label,cell_floor=floor,prior_count=prior,donor=dm$donor,disease=dm$disease,cells=dm$cells,gene=g,logCPM=as.numeric(z[g,]),raw_count=as.numeric(y$counts[g,]))
        }
      }
    }
  }
}
dir.create(file.path(paper,'cache/u5_ipf_compatibility'),showWarnings=FALSE,recursive=TRUE)
write.csv(do.call(rbind,output),gzfile(file.path(paper,'cache/u5_ipf_compatibility',paste0(cohort,'_components.csv.gz'))),row.names=FALSE)
write.csv(do.call(rbind,norm),file.path(paper,'trials/u5_ipf_compatibility',paste0(cohort,'_normalization.csv')),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(paper,'cache/u5_ipf_compatibility',paste0(cohort,'_R_session.txt')))
cat(cohort,'TMM component export complete\n')
