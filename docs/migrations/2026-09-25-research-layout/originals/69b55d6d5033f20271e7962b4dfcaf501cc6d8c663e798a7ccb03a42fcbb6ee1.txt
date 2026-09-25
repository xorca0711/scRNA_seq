# Patient-blocked human lesion comparisons; no cell-level replication.
root <- normalizePath('.',winslash='/');paper <- file.path(root,'Thesis/gate2_C3_yu_lee_choi_min_2026')
out <- file.path(paper,'trials/u5_human_niche');cache <- file.path(paper,'cache/u5_human_niche')
.libPaths(c(file.path(root,'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
members <- read.delim(file.path(paper,'trials/u5_ipf_spec/pathway_genes.tsv'));sets <- split(members$gene,members$set)
assignments <- read.delim(file.path(paper,'trials/u5_ipf_spec/compartment_sets.tsv'))
edges <- read.csv(file.path(paper,'trials/u5_ipf_compatibility/resource_edges.csv'))
components <- sort(unique(c(unlist(strsplit(c(edges$ligand,edges$receptor),'_',fixed=TRUE)),'IL1RN','IL1R2','SIGIRR','IL1RAP')))
contrasts <- list(c('AAH','normal'),c('AIS','normal'),c('MIA','normal'),c('LUAD','normal'),c('LUAD','AAH'),c('LUAD','AIS'),c('LUAD','MIA'))
program_members <- read.delim(file.path(paper,'trials/u6_human_specificity/human_module_genes.tsv'))
programs <- split(program_members$gene,program_members$module)
program_source_n <- tapply(program_members$source_genes,program_members$module,unique)
program_rows <- list()
results <- list();eligible <- list();target_gate <- list();coverage_rows <- list();normalizations <- list();designs <- list()
configs <- c('unc20_pooled','unc30_pooled','unc20_largest_library')
for(config in configs) {
 m <- read.csv(file.path(out,paste0(config,'_triad_units.csv')));x <- as.matrix(read.csv(gzfile(file.path(cache,paste0(config,'_triad_counts.csv.gz'))),row.names=1,check.names=FALSE))
 stopifnot(identical(colnames(x),m$unit_id),all(x>=0),all(x==round(x)),all(colSums(x)==m$full_library_sum))
 cg <- intersect(components,rownames(x));floors <- if(config=='unc20_pooled') c(50,30,100) else 50
 for(comp in unique(m$comp)) for(label in unique(m$label[m$comp==comp])) for(floor in floors) {
  take <- m$comp==comp & m$label==label & m$cells>=floor;mm <- m[take,,drop=FALSE];xx <- x[,take,drop=FALSE]
  stopifnot(!anyDuplicated(paste(mm$patient,mm$histology)))
  # All retained patient/histology units normalize together within each view;
  # compatibility contrasts later retain only complete source/receiver pairs.
  if(ncol(xx)>=2) {
   yn <- suppressMessages(calcNormFactors(DGEList(counts=xx),method='TMM'))
   priors <- if(config=='unc20_pooled' && floor==50) c(1,.5,2) else 1
   for(prior in priors) {
    lc_all <- cpm(yn,log=TRUE,prior.count=prior);lc <- lc_all[cg,,drop=FALSE]
    for(j in seq_len(ncol(lc))) normalizations[[length(normalizations)+1]] <- data.frame(config,comp,label,cell_floor=floor,prior_count=prior,patient=mm$patient[[j]],histology=mm$histology[[j]],cells=mm$cells[[j]],gene=cg,logCPM=lc[,j])
    if(prior==1 && comp=='AT2') for(module in names(programs)) {
     gs <- intersect(programs[[module]],rownames(lc_all));fraction <- length(gs)/program_source_n[[module]]
     if(!length(gs)) next
     values <- colMeans(lc_all[gs,,drop=FALSE])
     program_rows[[length(program_rows)+1]] <- data.frame(config,comp,label,cell_floor=floor,patient=mm$patient,histology=mm$histology,cells=mm$cells,module,assayed_source_fraction=fraction,eligible=fraction>=.7,mean_logCPM=values)
    }
   }
  }
  for(pair in contrasts) {
   case <- pair[[1]];reference <- pair[[2]];patients <- intersect(mm$patient[mm$histology==case],mm$patient[mm$histology==reference]);ok <- length(patients)>=3
   eligible[[length(eligible)+1]] <- data.frame(config,comp,label,cell_floor=floor,case,reference,complete_patients=length(patients),eligible=ok)
   if(!ok) next
   sel <- mm$patient%in%patients & mm$histology%in%pair;meta <- mm[sel,,drop=FALSE];counts <- xx[,sel,drop=FALSE]
   patient <- factor(meta$patient);arm <- as.numeric(meta$histology==case);design <- model.matrix(~patient+arm)
   stopifnot(qr(design)$rank==ncol(design),nrow(design)==2*length(patients),nrow(design)-ncol(design)==length(patients)-1)
   y <- DGEList(counts=counts);keep <- filterByExpr(y,design=design,min.count=10)
   if(sum(keep)<20) next
   y <- suppressMessages(calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM'));v <- voom(y,design,plot=FALSE)
   chosen <- sets[unique(assignments$set[assignments$compartment==comp])];coverage <- vapply(chosen,function(g)mean(unique(g)%in%rownames(x)),0.)
   ix <- ids2indices(chosen,rownames(y),remove.empty=FALSE);pass <- lengths(ix)>=10 & coverage>=.7
   coverage_rows[[length(coverage_rows)+1]] <- data.frame(config,comp,label,cell_floor=floor,case,reference,set=names(chosen),assayed_fraction=coverage,tested_genes=lengths(ix),eligible=pass)
   designs[[length(designs)+1]] <- data.frame(config,comp,label,cell_floor=floor,case,reference,n_patients=length(patients),genes_tested=nrow(y),residual_df=nrow(design)-ncol(design),patient_ids=paste(sort(patients),collapse=';'))
   for(setting in c('estimated','fixed001')) {
    if(!any(pass)) next
    cor <- if(setting=='estimated') NA_real_ else .01
    z <- camera(v,ix[pass],design,contrast=ncol(design),inter.gene.cor=cor,allow.neg.cor=FALSE,sort=FALSE)
    if(!'Correlation'%in%names(z)) z$Correlation <- cor
    z$set <- rownames(z);rownames(z) <- NULL;z$config <- config;z$comp <- comp;z$label <- label;z$cell_floor <- floor;z$case <- case;z$reference <- reference;z$correlation_setting <- setting;z$n_patients <- length(patients);z$assayed_fraction <- coverage[z$set];results[[length(results)+1]] <- z
   }
   if(config=='unc20_pooled' && floor==50) {
    de <- topTable(eBayes(lmFit(v,design)),coef=ncol(design),number=Inf,sort.by='none',adjust.method='BH');de$gene <- rownames(de);rownames(de) <- NULL
    safe <- gsub('[^A-Za-z0-9_]','_',paste(comp,label,case,reference,sep='_'));write.csv(de,gzfile(file.path(cache,paste0(safe,'_DE.csv.gz'))),row.names=FALSE)
    write.csv(meta,file.path(cache,paste0(safe,'_DE_units.csv')),row.names=FALSE)
    target_gate[[length(target_gate)+1]] <- data.frame(comp,label,case,reference,n_patients=length(patients),genes_tested=nrow(de),q05_up=sum(de$adj.P.Val<.05 & de$logFC>0),q05_down=sum(de$adj.P.Val<.05 & de$logFC<0),file_stem=safe)
   }
  }
 }
 cat(config,'paired analyses complete\n');flush.console()
}
z <- do.call(rbind,results);z$FDR_global <- NA_real_
for(config in unique(z$config)) for(floor in unique(z$cell_floor)) for(setting in unique(z$correlation_setting)) {
 take <- z$config==config & z$cell_floor==floor & z$correlation_setting==setting
 z$FDR_global[take] <- p.adjust(z$PValue[take],method='BH')
}
write.csv(z,file.path(out,'camera.csv'),row.names=FALSE)
write.csv(do.call(rbind,eligible),file.path(out,'paired_eligibility.csv'),row.names=FALSE)
write.csv(do.call(rbind,coverage_rows),file.path(out,'pathway_assay_coverage.csv'),row.names=FALSE)
write.csv(do.call(rbind,designs),file.path(out,'paired_designs.csv'),row.names=FALSE)
write.csv(do.call(rbind,target_gate),file.path(out,'DE_target_eligibility.csv'),row.names=FALSE)
write.csv(do.call(rbind,normalizations),gzfile(file.path(cache,'normalized_components.csv.gz')),row.names=FALSE)
write.csv(do.call(rbind,program_rows),file.path(paper,'trials/u6_human_specificity/patient_histology_scores.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
