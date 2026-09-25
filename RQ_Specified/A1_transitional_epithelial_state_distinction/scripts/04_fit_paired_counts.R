# Called by 03_run_paired_counts.py after the verified contract gate passes.
# Source: edgeR User's Guide (Bioconductor); no integrated embeddings as input.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- normalizePath(args[[1]], winslash='/')
local_lib <- file.path(getwd(), 'analysis/corrections/statistics/.tools/R-library')
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages(library(edgeR))
parameters <- read.delim(file.path(out,'input_contract.tsv'),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(parameters$key))
c <- as.list(setNames(parameters$value,parameters$key))
for (key in c('minimum_independent_pairs','min_count','min_total_count')) c[[key]] <- as.numeric(c[[key]])
stopifnot(c$status=='frozen', c$input_scale=='raw_integer_counts')
stopifnot(unname(tools::md5sum(c$counts_path))==c$counts_md5)
stopifnot(unname(tools::md5sum(c$samples_path))==c$samples_md5)
meta <- read.delim(c$samples_path, check.names=FALSE)
x <- as.matrix(read.delim(c$counts_path,row.names=1,check.names=FALSE))
stopifnot(identical(colnames(x),meta$sample_id), !anyDuplicated(rownames(x)),
          all(is.finite(x)),all(x>=0),all(x==round(x)),all(colSums(x)>0))
meta$unit <- factor(meta$biological_unit_id)
meta$group <- factor(meta$group,levels=c(c$reference,c$case))
stopifnot(!anyNA(meta$group),nlevels(meta$unit)>=c$minimum_independent_pairs,
          all(table(meta$unit,meta$group)==1))
design <- model.matrix(~unit+group,meta)
stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
y <- DGEList(x)
keep <- filterByExpr(y,group=meta$group,min.count=c$min_count,min.total.count=c$min_total_count)
write.table(data.frame(feature_id=rownames(x),retained=keep),file.path(out,'feature_filter.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
stopifnot(sum(keep)>1)
y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
y <- estimateDisp(y,design,robust=TRUE)
fit <- glmQLFit(y,design,robust=TRUE)
test <- glmQLFTest(fit,coef=ncol(design))
tab <- topTags(test,n=Inf,sort.by='none')$table
tab$feature_id <- rownames(tab)
write.table(tab,file.path(out,'paired_effects.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(cbind(meta,y$samples[,c('lib.size','norm.factors')]),file.path(out,'normalization.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(design,file.path(out,'design.tsv'),sep='\t',quote=FALSE,col.names=NA)
logcpm <- cpm(y,log=TRUE,prior.count=2)
vary <- apply(logcpm,1,var)
top <- names(sort(vary[vary>0],decreasing=TRUE))[seq_len(min(2000,sum(vary>0)))]
if (length(top)>1) {
  pc <- prcomp(t(logcpm[top,,drop=FALSE]),center=TRUE,scale.=FALSE)
  write.table(cbind(meta,pc$x),file.path(out,'sample_PCA.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(component=seq_along(pc$sdev),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),file.path(out,'PCA_variance.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(feature_id=top),file.path(out,'PCA_features.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
write.table(data.frame(status='completed',tested_features=nrow(tab),independent_pairs=nlevels(meta$unit),
                       direction='case minus reference',multiple_testing='BH within this contrast',
                       coefficient_tested=colnames(design)[ncol(design)],
                       interpretation='Within-study paired association, not causal or taxonomy evidence'),
            file.path(out,'fit_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
