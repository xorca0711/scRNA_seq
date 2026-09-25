# Early mouse niche associations using independently reviewed broad lineages.
root <- normalizePath('.',winslash='/');paper <- file.path(root,'Research Article/gate2_C3_yu_lee_choi_min_2026')
out <- file.path(paper,'trials/u4_mouse_niche');dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache <- file.path(paper,'cache/u4_mouse_niche');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(root,'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
counts <- as.matrix(read.csv(gzfile(file.path(paper,'cache/u3_lineage_annotation/subtype_counts.csv.gz')),row.names=1,check.names=FALSE))
meta <- read.csv(file.path(paper,'trials/u3_lineage_annotation/subtype_units.csv'))
stopifnot(identical(colnames(counts),meta$unit_id),all(counts==round(counts)),all(counts>=0))
members <- read.delim(file.path(paper,'trials/u4_resources/pathway_genes.tsv'));sets <- split(members$gene,members$set)
assignments <- read.delim(file.path(paper,'trials/u4_resources/compartment_sets.tsv'))
resources <- read.csv(file.path(paper,'trials/u4_resources/mouse_resource_edges.csv'))
genes <- sort(unique(unlist(strsplit(c(resources$ligand,resources$receptor),'_',fixed=TRUE))))
genes <- intersect(genes,rownames(counts))
results <- list();eligibility <- list();components <- list();designs <- list();de_eligibility <- list()
for(comp in c('alveolar','fibroblast','myeloid')) {
  for(label in c('__broad__',sort(unique(meta$label[meta$compartment==comp])))) {
    take <- meta$compartment==comp & (label=='__broad__' | meta$label==label);m0 <- meta[take,,drop=FALSE];x0 <- counts[,take,drop=FALSE]
    if(label=='__broad__') {
      ids <- unique(m0$gsm);m <- m0[match(ids,m0$gsm),c('gsm','treatment','cells'),drop=FALSE]
      x <- matrix(0,nrow=nrow(x0),ncol=length(ids),dimnames=list(rownames(x0),ids))
      for(j in seq_along(ids)) {sel <- m0$gsm==ids[[j]];x[,j] <- rowSums(x0[,sel,drop=FALSE]);m$cells[[j]] <- sum(m0$cells[sel])}
    } else {m <- m0;x <- x0}
    stopifnot(!anyDuplicated(m$gsm))
    for(floor in c(50,30,100)) {
      take <- m$cells>=floor;mm <- m[take,,drop=FALSE];xx <- x[,take,drop=FALSE]
      arm <- as.numeric(mm$treatment!='Control IgG');n0 <- sum(arm==0);n1 <- sum(arm==1);eligible <- min(n0,n1)>=3
      eligibility[[length(eligibility)+1]] <- data.frame(compartment=comp,label,cell_floor=floor,n_control=n0,n_antiIL1B=n1,eligible)
      # Component normalization remains within the complete eligible subtype,
      # using all assayed genes; no candidate-gene-only library normalization.
      if(ncol(xx)>=2) {
        yn <- suppressMessages(calcNormFactors(DGEList(counts=xx),method='TMM'))
        for(prior in c(1,.5,2)) {
          lc <- cpm(yn,log=TRUE,prior.count=prior)[genes,,drop=FALSE]
          for(j in seq_len(ncol(lc))) components[[length(components)+1]] <- data.frame(compartment=comp,label,cell_floor=floor,prior_count=prior,gsm=mm$gsm[[j]],treatment=mm$treatment[[j]],gene=genes,logCPM=lc[,j])
        }
      }
      if(!eligible) next
      design <- model.matrix(~arm);y <- DGEList(counts=xx);keep <- filterByExpr(y,design=design,min.count=10)
      y <- suppressMessages(calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM'));v <- voom(y,design,plot=FALSE)
      chosen <- sets[unique(assignments$set[assignments$compartment==comp])];coverage <- vapply(chosen,function(g)mean(unique(g)%in%rownames(counts)),0.)
      ix <- ids2indices(chosen,rownames(y),remove.empty=FALSE);pass <- lengths(ix)>=10 & coverage>=.7
      designs[[length(designs)+1]] <- data.frame(compartment=comp,label,cell_floor=floor,genes_tested=nrow(y),eligible_sets=sum(pass),residual_df=nrow(design)-ncol(design))
      for(setting in c('estimated','fixed001')) {
        if(!any(pass)) next
        cor <- if(setting=='estimated') NA_real_ else .01
        z <- camera(v,ix[pass],design,contrast=2,inter.gene.cor=cor,allow.neg.cor=FALSE,sort=FALSE)
        if(!'Correlation'%in%names(z)) z$Correlation <- cor
        z$set <- rownames(z);rownames(z) <- NULL;z$compartment <- comp;z$label <- label;z$cell_floor <- floor;z$correlation_setting <- setting;z$n_control <- n0;z$n_antiIL1B <- n1;z$assayed_fraction <- coverage[z$set]
        results[[length(results)+1]] <- z
      }
      if(floor==50) {
        de <- topTable(eBayes(lmFit(v,design)),coef=2,number=Inf,sort.by='none',adjust.method='BH');de$gene <- rownames(de);rownames(de) <- NULL
        safe <- gsub('[^A-Za-z0-9_]','_',paste(comp,label,sep='_'));write.csv(de,gzfile(file.path(cache,paste0(safe,'_DE.csv.gz'))),row.names=FALSE)
        de_eligibility[[length(de_eligibility)+1]] <- data.frame(compartment=comp,label,genes_tested=nrow(de),q05_up=sum(de$adj.P.Val<.05 & de$logFC>0),q05_down=sum(de$adj.P.Val<.05 & de$logFC<0))
      }
    }
    cat(comp,label,'done\n');flush.console()
  }
}
z <- do.call(rbind,results);z$FDR_global_early <- NA_real_;z$FDR_reserving_late_family <- NA_real_
for(floor in unique(z$cell_floor)) for(setting in unique(z$correlation_setting)) {
  take <- z$cell_floor==floor & z$correlation_setting==setting
  z$FDR_global_early[take] <- p.adjust(z$PValue[take],method='BH')
  z$FDR_reserving_late_family[take] <- p.adjust(z$PValue[take],method='BH',n=2*sum(take))
}
write.csv(z,file.path(out,'camera.csv'),row.names=FALSE)
write.csv(do.call(rbind,eligibility),file.path(out,'eligibility.csv'),row.names=FALSE)
write.csv(do.call(rbind,designs),file.path(out,'designs.csv'),row.names=FALSE)
write.csv(do.call(rbind,de_eligibility),file.path(out,'DE_target_eligibility.csv'),row.names=FALSE)
write.csv(do.call(rbind,components),gzfile(file.path(cache,'normalized_components.csv.gz')),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
cat('Completed early mouse pseudobulk analysis\n')
