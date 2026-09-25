# NicheNet v2 prior Pearson statistic and complete donor leave-one-out refits.
# Uses base-R Pearson arithmetic; does not claim to run all nichenetr metrics.
args <- commandArgs(trailingOnly=TRUE);cohort <- args[[1]]
stopifnot(cohort %in% c('GSE136831','GSE135893'))
root <- normalizePath('.',winslash='/');paper <- file.path(root,'Research Article/gate2_C3_yu_lee_choi_min_2026')
out <- file.path(paper,'trials/u5_ligand_targets',cohort);dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache <- file.path(paper,'cache/u5_ligand_targets')
.libPaths(c(file.path(root,'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
prior <- readRDS(file.path(cache,'ligand_target_matrix_nsga2r_final.rds'))
lr <- readRDS(file.path(cache,'lr_network_human_21122021.rds'))
write.csv(lr,file.path(out,'prior_ligand_receptor.csv'),row.names=FALSE)
stopifnot(is.matrix(prior),!anyDuplicated(rownames(prior)),!anyDuplicated(colnames(prior)),all(is.finite(prior)))
pb <- file.path(root,'analysis/corrections/statistics/cache')
counts <- as.matrix(read.csv(gzfile(file.path(pb,paste0(cohort,'_subtypes_counts.csv.gz'))),row.names=1,check.names=FALSE))
meta <- read.csv(file.path(pb,paste0(cohort,'_subtypes_units.csv')))
stopifnot(identical(colnames(counts),meta$unit_id),all(counts>=0),all(counts==round(counts)))
score_rows <- list();eligibility <- list();target_rows <- list();units <- list();qa <- list()
start <- proc.time()[[3]]
for(comp in sort(unique(meta$comp))) {
  labels <- sort(unique(meta$label[meta$comp==comp]))
  # AT2 has one deposited subtype: avoid counting an identical broad fit twice.
  views <- if(length(labels)==1) labels else c('__broad__',labels)
  for(label in views) {
    take <- meta$comp==comp & (label=='__broad__' | meta$label==label)
    m0 <- meta[take,,drop=FALSE];x0 <- counts[,take,drop=FALSE]
    if(label=='__broad__') {
      keys <- unique(paste(m0$donor,m0$disease,sep='|'))
      m <- m0[match(keys,paste(m0$donor,m0$disease,sep='|')),c('donor','disease','cells'),drop=FALSE]
      x <- matrix(0,nrow=nrow(x0),ncol=length(keys),dimnames=list(rownames(x0),keys))
      for(j in seq_along(keys)) {sel <- paste(m0$donor,m0$disease,sep='|')==keys[[j]];x[,j] <- rowSums(x0[,sel,drop=FALSE]);m$cells[[j]] <- sum(m0$cells[sel])}
    } else {m <- m0;x <- x0}
    ok <- m$cells>=50 & m$disease %in% c('control','IPF');m <- m[ok,,drop=FALSE];x <- x[,ok,drop=FALSE]
    stopifnot(!anyDuplicated(m$donor))
    if(min(sum(m$disease=='control'),sum(m$disease=='IPF'))<3) next
    units[[length(units)+1]] <- data.frame(compartment=comp,label,m[,c('donor','disease','cells')])
    full_eligible <- FALSE
    for(omit in c('__none__',as.character(m$donor))) {
      if(omit!='__none__' && !full_eligible) next
      valid <- m$donor!=omit;mm <- m[valid,,drop=FALSE];xx <- x[,valid,drop=FALSE]
      n0 <- sum(mm$disease=='control');n1 <- sum(mm$disease=='IPF')
      if(min(n0,n1)<3) {
        for(direction in c('up','down')) eligibility[[length(eligibility)+1]] <- data.frame(compartment=comp,label,omitted_donor=omit,direction,n_control=n0,n_IPF=n1,tested_genes=NA,mapped_background=NA,targets=NA,mapped_targets=NA,eligible=FALSE,reason='fewer_than_three_units_per_arm_after_omission')
        next
      }
      mm$arm <- as.numeric(mm$disease=='IPF');design <- model.matrix(~arm,mm)
      y <- DGEList(counts=xx);keep <- filterByExpr(y,design=design,min.count=10)
      y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM');v <- voom(y,design,plot=FALSE)
      de <- topTable(eBayes(lmFit(v,design)),coef=2,number=Inf,sort.by='none',adjust.method='BH');de$gene <- rownames(de)
      if(omit=='__none__') {
        nm <- gsub('[^A-Za-z0-9_]','_',paste(comp,label,sep='_'))
        old <- read.csv(gzfile(file.path(paper,'trials/u5_ipf_pathways',paste0(cohort,'_',nm,'_DE.csv.gz'))))
        stopifnot(identical(old$gene,de$gene),max(abs(old$logFC-de$logFC))<1e-8,max(abs(old$adj.P.Val-de$adj.P.Val))<1e-8)
        qa[[length(qa)+1]] <- data.frame(compartment=comp,label,genes=nrow(de),max_q_difference=max(abs(old$adj.P.Val-de$adj.P.Val)))
      }
      bg <- intersect(de$gene,rownames(prior));pp <- prior[bg,,drop=FALSE]
      for(direction in c('up','down')) {
        target <- de$gene[de$adj.P.Val<.05 & if(direction=='up') de$logFC>0 else de$logFC<0]
        mapped <- intersect(target,bg);eligible <- length(mapped)>=10 && length(mapped)<length(bg)
        eligibility[[length(eligibility)+1]] <- data.frame(compartment=comp,label,omitted_donor=omit,direction,n_control=n0,n_IPF=n1,tested_genes=nrow(de),mapped_background=length(bg),targets=length(target),mapped_targets=length(mapped),eligible,reason=if(eligible) 'eligible' else 'fewer_than_ten_mapped_targets_or_no_background_negatives')
        if(omit=='__none__' && length(mapped)) target_rows[[length(target_rows)+1]] <- data.frame(compartment=comp,label,direction,gene=mapped)
        if(!eligible) next
        if(omit=='__none__') full_eligible <- TRUE
        response <- as.numeric(bg %in% mapped)
        r <- as.numeric(cor(response,pp,method='pearson'))
        # Independent direct centered-vector calculation for a deterministic ligand.
        j <- which(is.finite(r))[[1]];a <- response-mean(response);b <- pp[,j]-mean(pp[,j]);manual <- sum(a*b)/sqrt(sum(a*a)*sum(b*b))
        stopifnot(abs(manual-r[[j]])<1e-12)
        score_rows[[length(score_rows)+1]] <- data.frame(compartment=comp,label,omitted_donor=omit,direction,ligand=colnames(pp),pearson=r)
      }
    }
    cat(cohort,comp,label,'refits done\n');flush.console()
  }
}
write.csv(do.call(rbind,score_rows),gzfile(file.path(out,'all_prior_pearson_scores.csv.gz')),row.names=FALSE)
write.csv(do.call(rbind,eligibility),file.path(out,'target_eligibility.csv'),row.names=FALSE)
write.csv(do.call(rbind,target_rows),file.path(out,'full_fit_mapped_targets.csv'),row.names=FALSE)
write.csv(do.call(rbind,units),file.path(out,'receiver_units.csv'),row.names=FALSE)
write.csv(do.call(rbind,qa),file.path(out,'full_fit_DE_parity.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
cat(cohort,'completed',round(proc.time()[[3]]-start,1),'seconds\n')
