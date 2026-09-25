# New review-context analysis; historical correction outputs remain unchanged.
args <- commandArgs(trailingOnly=TRUE)
cohort <- args[[1]]
stopifnot(cohort %in% c('GSE136831','GSE135893'))
root <- normalizePath('.', winslash='/')
paper <- file.path(root,'Thesis/gate2_C3_yu_lee_choi_min_2026')
here <- file.path(paper,'trials/u5_ipf_pathways')
dir.create(here, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(file.path(root,'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
cache <- file.path(root,'analysis/corrections/statistics/cache')
counts <- as.matrix(read.csv(gzfile(file.path(cache,paste0(cohort,'_subtypes_counts.csv.gz'))),row.names=1,check.names=FALSE))
meta <- read.csv(file.path(cache,paste0(cohort,'_subtypes_units.csv')))
stopifnot(identical(colnames(counts),meta$unit_id), all(counts>=0), all(counts==round(counts)))
expected <- read.delim(file.path(paper,'trials/u5_ipf_spec/expected_counts.tsv'))
stopifnot(sum(counts)==expected$total_counts[expected$cohort==cohort])
members <- read.delim(file.path(paper,'trials/u5_ipf_spec/pathway_genes.tsv'))
sets <- split(members$gene,members$set)
mapping <- read.delim(file.path(paper,'trials/u5_ipf_spec/compartment_sets.tsv'))
output <- list(); eligibility <- list(); diagnostics <- list()
start <- proc.time()[[3]]
for(comp in sort(unique(meta$comp))) {
  for(label in c('__broad__',sort(unique(meta$label[meta$comp==comp])))) {
    take <- meta$comp==comp & (label=='__broad__' | meta$label==label)
    m0 <- meta[take,,drop=FALSE]; x0 <- counts[,take,drop=FALSE]
    if(label=='__broad__') {
      keys <- unique(paste(m0$donor,m0$disease,sep='|'))
      m <- m0[match(keys,paste(m0$donor,m0$disease,sep='|')),c('donor','disease','cells'),drop=FALSE]
      x <- matrix(0,nrow=nrow(x0),ncol=length(keys),dimnames=list(rownames(x0),keys))
      for(j in seq_along(keys)) {
        sel <- paste(m0$donor,m0$disease,sep='|')==keys[[j]]
        x[,j] <- rowSums(x0[,sel,drop=FALSE]); m$cells[[j]] <- sum(m0$cells[sel])
      }
    } else {m <- m0; x <- x0}
    stopifnot(!anyDuplicated(m$donor))
    chosen_sets <- unique(mapping$set[mapping$compartment==comp])
    chosen_sets <- sets[chosen_sets]
    for(floor in c(50,30,100)) {
      valid <- m$cells>=floor & m$disease %in% c('control','IPF')
      mm <- m[valid,,drop=FALSE]; xx <- x[,valid,drop=FALSE]
      n0 <- sum(mm$disease=='control'); n1 <- sum(mm$disease=='IPF')
      eligible <- min(n0,n1)>=3
      eligibility[[length(eligibility)+1]] <- data.frame(cohort,compartment=comp,label,cell_floor=floor,n_control=n0,n_IPF=n1,eligible)
      if(!eligible) next
      mm$arm <- as.numeric(mm$disease=='IPF')
      design <- model.matrix(~arm,mm)
      stopifnot(qr(design)$rank==ncol(design),nrow(design)-ncol(design)>=2)
      y <- DGEList(counts=xx)
      keep <- filterByExpr(y,design=design,min.count=10)
      y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
      v <- voom(y,design,plot=FALSE)
      assay_coverage <- vapply(chosen_sets,function(g) mean(unique(g) %in% rownames(counts)),0.0)
      ix <- ids2indices(chosen_sets,rownames(y),remove.empty=FALSE)
      pass <- lengths(ix)>=10 & assay_coverage>=0.7
      diagnostics[[length(diagnostics)+1]] <- data.frame(cohort,compartment=comp,label,cell_floor=floor,genes_tested=nrow(y),residual_df=nrow(design)-ncol(design),eligible_sets=sum(pass))
      for(cor_setting in c('estimated','fixed001')) {
        if(!any(pass)) next
        cor <- if(cor_setting=='estimated') NA_real_ else 0.01
        z <- camera(v,ix[pass],design,contrast=2,inter.gene.cor=cor,allow.neg.cor=FALSE,sort=FALSE)
        # CAMERA includes an estimated-correlation column only for some settings.
        # Normalize the output schema before combining primary and sensitivity fits.
        if(!'Correlation' %in% names(z)) z$Correlation <- cor
        z$set <- rownames(z); rownames(z) <- NULL
        z$cohort <- cohort; z$compartment <- comp; z$label <- label
        z$cell_floor <- floor; z$correlation_setting <- cor_setting
        z$n_control <- n0; z$n_IPF <- n1
        z$assayed_fraction <- assay_coverage[z$set]
        z$mean_member_log2FC <- vapply(z$set,function(nm) {
          g <- intersect(chosen_sets[[nm]],rownames(v$E))
          mean(rowMeans(v$E[g,mm$arm==1,drop=FALSE])-rowMeans(v$E[g,mm$arm==0,drop=FALSE]))
        },0.0)
        output[[length(output)+1]] <- z
      }
      if(floor==50) {
        fit <- eBayes(lmFit(v,design))
        de <- topTable(fit,coef=2,number=Inf,sort.by='none',adjust.method='BH')
        de$gene <- rownames(de); rownames(de) <- NULL
        safe_label <- gsub('[^A-Za-z0-9_]','_',paste(comp,label,sep='_'))
        write.csv(de,gzfile(file.path(here,paste0(cohort,'_',safe_label,'_DE.csv.gz'))),row.names=FALSE)
      }
    }
  }
}
result <- if(length(output)) do.call(rbind,output) else data.frame()
if(nrow(result)) {
  result$FDR_global <- NA_real_
  for(floor in unique(result$cell_floor)) for(cor in unique(result$correlation_setting)) {
    idx <- result$cell_floor==floor & result$correlation_setting==cor
    result$FDR_global[idx] <- p.adjust(result$PValue[idx],method='BH')
  }
}
write.csv(result,file.path(here,paste0(cohort,'_camera.csv')),row.names=FALSE)
write.csv(do.call(rbind,eligibility),file.path(here,paste0(cohort,'_eligibility.csv')),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(here,paste0(cohort,'_design.csv')),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(here,paste0(cohort,'_R_session.txt')))
cat(cohort, 'completed', nrow(result), 'set tests including sensitivities in',round(proc.time()[[3]]-start,2),'seconds\n')
