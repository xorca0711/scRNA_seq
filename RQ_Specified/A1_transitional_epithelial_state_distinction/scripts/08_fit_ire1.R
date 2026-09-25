# Frozen GSE190821 unpaired mouse-level contrast; do not reuse the paired pilot.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- normalizePath(args[[1]], winslash='/')
.libPaths(c(file.path(getwd(),'analysis/corrections/statistics/.tools/R-library'), .libPaths()))
suppressPackageStartupMessages(library(edgeR))
parameters <- read.delim(file.path(out,'input_contract.tsv'),stringsAsFactors=FALSE)
cfg <- as.list(setNames(parameters$value, parameters$key))
stopifnot(unname(tools::md5sum(cfg$counts_path))==cfg$counts_md5,
          unname(tools::md5sum(cfg$samples_path))==cfg$samples_md5)
meta <- read.delim(cfg$samples_path,check.names=FALSE)
x <- as.matrix(read.delim(cfg$counts_path,row.names=1,check.names=FALSE))
stopifnot(identical(colnames(x),meta$sample_id),!anyDuplicated(rownames(x)),
          !anyDuplicated(meta$mouse),all(x>=0),all(is.finite(x)),all(x==round(x)),
          all(colSums(x)>=as.numeric(cfg$minimum_library_count)),nrow(meta)==10,
          all(table(meta$group)==5),all(meta$group %in% c('Vehicle','KIRA8')))
meta$batch <- factor(meta$batch)
meta$sex <- factor(meta$sex)
meta$group <- factor(meta$group,levels=c('Vehicle','KIRA8'))
design <- model.matrix(as.formula(cfg$primary_model),meta)
sensitivity <- model.matrix(as.formula(cfg$sensitivity_model),meta)
stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design),
          qr(sensitivity)$rank==ncol(sensitivity),'groupKIRA8' %in% colnames(design))
write_tsv <- function(frame,name) write.table(frame,file.path(out,name),sep='\t',quote=FALSE,row.names=FALSE)
write.table(design,file.path(out,'design.tsv'),sep='\t',quote=FALSE,col.names=NA)
write.table(sensitivity,file.path(out,'sensitivity_design.tsv'),sep='\t',quote=FALSE,col.names=NA)
y <- DGEList(x)
keep <- filterByExpr(y,group=meta$group,min.count=as.numeric(cfg$min_count),min.total.count=as.numeric(cfg$min_total_count))
write_tsv(data.frame(feature_id=rownames(x),retained=keep),'feature_filter.tsv')
stopifnot(sum(keep)>1)
y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
fit_model <- function(d) {
  dy <- estimateDisp(y,d,robust=TRUE)
  fit <- glmQLFit(dy,d,robust=TRUE)
  result <- topTags(glmQLFTest(fit,coef=which(colnames(d)=='groupKIRA8')),n=Inf,sort.by='none')$table
  result$feature_id <- rownames(result)
  result
}
tab <- fit_model(design)
secondary <- fit_model(sensitivity)
write_tsv(tab,'gene_effects.tsv')
write_tsv(secondary,'gene_effects_without_sex.tsv')
write_tsv(cbind(meta,y$samples[,c('lib.size','norm.factors')]),'normalization.tsv')
logcpm <- cpm(y,log=TRUE,prior.count=2)
variance <- apply(logcpm,1,var)
features <- names(sort(variance[variance>0],decreasing=TRUE))[seq_len(min(as.numeric(cfg$pca_features),sum(variance>0)))]
pc <- prcomp(t(logcpm[features,,drop=FALSE]),center=TRUE,scale.=FALSE)
write_tsv(cbind(meta,pc$x),'sample_PCA.tsv')
write_tsv(data.frame(component=seq_along(pc$sdev),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),'PCA_variance.tsv')
write_tsv(data.frame(feature_id=features),'PCA_features.tsv')
markers <- read.delim(file.path(out,'marker_mapping.tsv'))
ids <- intersect(markers$feature_id,rownames(logcpm))
write_tsv(data.frame(feature_id=ids,logcpm[ids,,drop=FALSE],check.names=FALSE),'marker_logCPM.tsv')
sets <- read.delim(file.path(out,'gene_set_mapping.tsv'))
sets$retained <- sets$feature_id %in% rownames(y)
write_tsv(sets,'gene_set_mapping_tested.tsv')
set_names <- unique(sets$gene_set)
coverage <- do.call(rbind,lapply(set_names,function(s) {
  g <- sets[sets$gene_set==s,]
  data.frame(gene_set=s,source_symbols=nrow(g),mapped_symbols=sum(g$feature_id!='unmapped'),
             tested_symbols=sum(g$retained),unique_tested_genes=length(unique(g$feature_id[g$retained])),
             fraction_tested=mean(g$retained))
}))
coverage$eligible <- coverage$fraction_tested>=as.numeric(cfg$gene_set_min_fraction) & coverage$unique_tested_genes>=as.numeric(cfg$gene_set_min_genes)
write_tsv(coverage,'gene_set_coverage.tsv')
index <- lapply(coverage$gene_set[coverage$eligible],function(s) which(rownames(y) %in% sets$feature_id[sets$gene_set==s & sets$retained]))
names(index) <- coverage$gene_set[coverage$eligible]
if (length(index)) {
  v <- voom(y,design,plot=FALSE)
  cam <- camera(v,index,design,contrast=which(colnames(design)=='groupKIRA8'),inter.gene.cor=NA,allow.neg.cor=FALSE)
  cam$FDR <- p.adjust(cam$PValue,method='BH')
  cam$gene_set <- rownames(cam)
  write_tsv(cam,'gene_set_camera.tsv')
}
write_tsv(data.frame(status='completed',input_features=nrow(x),tested_features=nrow(tab),
                     mice=10,vehicle=5,KIRA8=5,design_rank=qr(design)$rank,residual_df=nrow(design)-ncol(design),
                     genes_FDR05=sum(tab$FDR<.05),higher_KIRA8_FDR05=sum(tab$FDR<.05 & tab$logFC>0),
                     lower_KIRA8_FDR05=sum(tab$FDR<.05 & tab$logFC<0),
                     sensitivity_logFC_Pearson=cor(tab$logFC,secondary$logFC),
                     primary_direction='KIRA8 minus Vehicle',gene_set_tests=length(index)), 'fit_summary.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
