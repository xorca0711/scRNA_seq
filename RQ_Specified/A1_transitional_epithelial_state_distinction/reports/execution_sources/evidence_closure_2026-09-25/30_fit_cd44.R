# Paired contrasts and direct interaction; config/closure_analysis_contract.json.
args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==1)
out <- normalizePath(args[1],winslash='/')
.libPaths(c(file.path(getwd(),'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR))
h <- read.delim(file.path(out,'input_hashes.tsv'))
stopifnot(unname(tools::md5sum(file.path(out,'counts_input.tsv.gz')))==h$value[h$key=='counts_md5'],
          unname(tools::md5sum(file.path(out,'sample_manifest.tsv')))==h$value[h$key=='samples_md5'])
m <- read.delim(file.path(out,'sample_manifest.tsv'),check.names=FALSE)
x <- as.matrix(read.delim(gzfile(file.path(out,'counts_input.tsv.gz')),row.names=1,check.names=FALSE))
markers <- read.delim(file.path(out,'marker_mapping.tsv'),colClasses='character')
stopifnot(identical(colnames(x),m$sample_id),all(x>=0),all(x==round(x)),!anyDuplicated(rownames(x)),
          length(unique(m$mouse))==8,all(table(m$mouse,m$sort)==1),all(table(m$genotype,m$sort)==4))
write_tsv <- function(d,n) write.table(d,file.path(out,n),sep='\t',quote=FALSE,row.names=FALSE)
make_design <- function(meta) {
  meta$mouse <- factor(meta$mouse)
  positive <- as.integer(meta$sort=='positive')
  cbind(model.matrix(~0+mouse,meta),positiveWT=positive*(meta$genotype=='WT'),
        positiveMutant=positive*(meta$genotype=='Mutant'))
}
design <- make_design(m)
stopifnot(qr(design)$rank==10,ncol(design)==10,nrow(design)-ncol(design)==6)
y <- DGEList(x)
keep <- filterByExpr(y,group=interaction(m$genotype,m$sort),min.count=10,min.total.count=15)
write_tsv(data.frame(feature_id=rownames(x),retained=keep),'feature_filter.tsv')
y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
fit <- glmQLFit(estimateDisp(y,design,robust=TRUE),design,robust=TRUE)
contrasts <- matrix(0,ncol(design),3,dimnames=list(colnames(design),c('WT_positive_minus_negative','Mutant_positive_minus_negative','Mutant_minus_WT_CD44_interaction')))
contrasts['positiveWT',c(1,3)] <- c(1,-1)
contrasts['positiveMutant',c(2,3)] <- c(1,1)
tables <- lapply(seq_len(3),function(i) {
  z <- topTags(glmQLFTest(fit,contrast=contrasts[,i]),n=Inf,sort.by='none')$table
  z$feature_id <- rownames(z);z$contrast <- colnames(contrasts)[i];z
})
all <- do.call(rbind,tables);all$FDR_all_three <- p.adjust(all$PValue,'BH')
write_tsv(all,'all_gene_effects.tsv')
focus <- merge(markers,all,by='feature_id',all.x=TRUE,sort=FALSE)
write_tsv(focus,'focus_effects.tsv')
write.table(design,file.path(out,'design.tsv'),sep='\t',quote=FALSE,col.names=NA)
write.table(contrasts,file.path(out,'contrasts.tsv'),sep='\t',quote=FALSE,col.names=NA)
write_tsv(cbind(m,raw_library_count=colSums(x),retained_count=colSums(y$counts),
                y$samples[,c('lib.size','norm.factors')]),'normalization.tsv')
logcpm <- cpm(y,log=TRUE,prior.count=2)
ids <- intersect(markers$feature_id,rownames(logcpm))
write_tsv(data.frame(feature_id=ids,logcpm[ids,,drop=FALSE],check.names=FALSE),'marker_logCPM.tsv')
paired <- do.call(rbind,lapply(unique(m$mouse),function(mouse) {
  sel <- which(m$mouse==mouse); pos <- sel[m$sort[sel]=='positive'];neg <- sel[m$sort[sel]=='negative']
  data.frame(feature_id=ids,mouse=mouse,genotype=m$genotype[pos],logCPM_positive=logcpm[ids,pos],
             logCPM_negative=logcpm[ids,neg],difference=logcpm[ids,pos]-logcpm[ids,neg])
}))
write_tsv(merge(markers,paired,by='feature_id',sort=FALSE),'marker_paired_differences.tsv')
vary <- apply(logcpm,1,var);top <- names(sort(vary,decreasing=TRUE))[seq_len(min(2000,length(vary)))]
pc <- prcomp(t(logcpm[top,]),center=TRUE,scale.=FALSE)
write_tsv(cbind(m,pc$x),'sample_PCA.tsv')
write_tsv(data.frame(component=seq_along(pc$sdev),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),'PCA_variance.tsv')
write_tsv(data.frame(feature_id=top),'PCA_features.tsv')
summary <- do.call(rbind,lapply(colnames(contrasts),function(name) {
  z <- all[all$contrast==name,]
  data.frame(contrast=name,mice=8,rank=10,residual_df=6,genes=nrow(z),FDR05=sum(z$FDR<.05),
             joint_FDR05=sum(z$FDR_all_three<.05),joint_up=sum(z$FDR_all_three<.05 & z$logFC>0),
             joint_down=sum(z$FDR_all_three<.05 & z$logFC<0))
}))
write_tsv(summary,'fit_summary.tsv')
# Result-dependent launch rule frozen before the primary fit. Gene universe unchanged.
if (any(focus$FDR_all_three<.05,na.rm=TRUE)) {
  sensitivity <- list()
  for (mouse in unique(m$mouse)) {
    selected <- m$mouse!=mouse;sm <- m[selected,];sd <- make_design(sm)
    stopifnot(qr(sd)$rank==ncol(sd),nrow(sd)-ncol(sd)==5)
    sy <- calcNormFactors(DGEList(y$counts[,selected,drop=FALSE]),method='TMM')
    sf <- glmQLFit(estimateDisp(sy,sd,robust=TRUE),sd,robust=TRUE)
    coeff <- sf$coefficients[ids,,drop=FALSE]/log(2)
    sensitivity[[mouse]] <- data.frame(feature_id=rep(ids,3),omitted_mouse=mouse,
      contrast=rep(colnames(contrasts),each=length(ids)),
      logFC=c(coeff[,'positiveWT'],coeff[,'positiveMutant'],coeff[,'positiveMutant']-coeff[,'positiveWT']),rank=ncol(sd),residual_df=5)
    cat('Omit ',mouse,' completed\n',sep='');flush.console()
  }
  write_tsv(merge(markers,do.call(rbind,sensitivity),by='feature_id',sort=FALSE),'focus_leave_one_pair_out.tsv')
}
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
