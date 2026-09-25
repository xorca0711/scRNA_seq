# Conditional paired patient target refits and exact NicheNet Pearson statistic.
root <- normalizePath('.',winslash='/');paper <- file.path(root,'Research Article/gate2_C3_yu_lee_choi_min_2026')
out <- file.path(paper,'trials/u5_human_ligand_targets');dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache <- file.path(paper,'cache/u5_human_niche');prior_cache <- file.path(paper,'cache/u5_ligand_targets')
.libPaths(c(file.path(root,'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
prior <- readRDS(file.path(prior_cache,'ligand_target_matrix_nsga2r_final.rds'))
lr <- readRDS(file.path(prior_cache,'lr_network_human_21122021.rds'));write.csv(lr,file.path(out,'prior_ligand_receptor.csv'),row.names=FALSE)
m <- read.csv(file.path(paper,'trials/u5_human_niche/unc20_pooled_triad_units.csv'))
x <- as.matrix(read.csv(gzfile(file.path(cache,'unc20_pooled_triad_counts.csv.gz')),row.names=1,check.names=FALSE));stopifnot(identical(colnames(x),m$unit_id))
gate <- read.csv(file.path(paper,'trials/u5_human_niche/DE_target_eligibility.csv'));results <- list();eligibility <- list();targets <- list();parity <- list();units <- list()
for(i in seq_len(nrow(gate))) {
 g <- gate[i,];take <- m$comp==g$comp & m$label==g$label & m$cells>=50 & m$histology%in%c(g$case,g$reference)
 mm <- m[take,,drop=FALSE];xx <- x[,take,drop=FALSE];patients <- intersect(mm$patient[mm$histology==g$case],mm$patient[mm$histology==g$reference]);take <- mm$patient%in%patients;mm <- mm[take,,drop=FALSE];xx <- xx[,take,drop=FALSE]
 stopifnot(length(patients)==g$n_patients,nrow(mm)==2*length(patients));full_eligible <- FALSE
 units[[length(units)+1]] <- data.frame(comp=g$comp,label=g$label,case=g$case,reference=g$reference,mm[,c('patient','histology','cells')])
 for(omit in c('__none__',as.character(patients))) {
  if(omit!='__none__' && !full_eligible) next
  valid <- as.character(mm$patient)!=omit;meta <- mm[valid,,drop=FALSE];counts <- xx[,valid,drop=FALSE];n <- length(unique(meta$patient))
  if(n<3) {
   for(direction in c('up','down')) eligibility[[length(eligibility)+1]] <- data.frame(comp=g$comp,label=g$label,case=g$case,reference=g$reference,omitted_patient=omit,direction,n_patients=n,mapped_background=NA,mapped_targets=NA,eligible=FALSE,reason='fewer_than_three_complete_patients_after_omission')
   next
  }
  patient <- factor(meta$patient);arm <- as.numeric(meta$histology==g$case);design <- model.matrix(~patient+arm)
  y <- DGEList(counts=counts);keep <- filterByExpr(y,design=design,min.count=10);y <- suppressMessages(calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM'));v <- voom(y,design,plot=FALSE)
  de <- topTable(eBayes(lmFit(v,design)),coef=ncol(design),number=Inf,sort.by='none',adjust.method='BH');de$gene <- rownames(de)
  if(omit=='__none__') {
   old <- read.csv(gzfile(file.path(cache,paste0(g$file_stem,'_DE.csv.gz'))));stopifnot(identical(old$gene,de$gene),max(abs(old$logFC-de$logFC))<1e-8,max(abs(old$adj.P.Val-de$adj.P.Val))<1e-8)
   parity[[length(parity)+1]] <- data.frame(comp=g$comp,label=g$label,case=g$case,reference=g$reference,genes=nrow(de),maximum_q_difference=max(abs(old$adj.P.Val-de$adj.P.Val)))
  }
  bg <- intersect(de$gene,rownames(prior));pp <- prior[bg,,drop=FALSE]
  for(direction in c('up','down')) {
   target <- de$gene[de$adj.P.Val<.05 & if(direction=='up') de$logFC>0 else de$logFC<0];mapped <- intersect(target,bg);ok <- length(mapped)>=10 && length(mapped)<length(bg)
   eligibility[[length(eligibility)+1]] <- data.frame(comp=g$comp,label=g$label,case=g$case,reference=g$reference,omitted_patient=omit,direction,n_patients=n,mapped_background=length(bg),mapped_targets=length(mapped),eligible=ok,reason=if(ok) 'eligible' else 'fewer_than_ten_mapped_targets_or_no_background_negatives')
   if(omit=='__none__' && length(mapped)) targets[[length(targets)+1]] <- data.frame(comp=g$comp,label=g$label,case=g$case,reference=g$reference,direction,gene=mapped)
   if(!ok) next
   if(omit=='__none__') full_eligible <- TRUE
   response <- as.numeric(bg%in%mapped);r <- as.numeric(cor(response,pp,method='pearson'));j <- which(is.finite(r))[[1]];a <- response-mean(response);b <- pp[,j]-mean(pp[,j]);stopifnot(abs(sum(a*b)/sqrt(sum(a*a)*sum(b*b))-r[[j]])<1e-12)
   results[[length(results)+1]] <- data.frame(comp=g$comp,label=g$label,case=g$case,reference=g$reference,omitted_patient=omit,direction,ligand=colnames(pp),pearson=r)
  }
 }
 cat(g$file_stem,'paired target refits complete\n');flush.console()
}
write.csv(do.call(rbind,results),gzfile(file.path(out,'all_prior_pearson_scores.csv.gz')),row.names=FALSE)
write.csv(do.call(rbind,eligibility),file.path(out,'target_eligibility.csv'),row.names=FALSE)
write.csv(do.call(rbind,targets),file.path(out,'full_fit_mapped_targets.csv'),row.names=FALSE)
write.csv(do.call(rbind,parity),file.path(out,'full_fit_DE_parity.csv'),row.names=FALSE)
write.csv(do.call(rbind,units),file.path(out,'receiver_units.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
