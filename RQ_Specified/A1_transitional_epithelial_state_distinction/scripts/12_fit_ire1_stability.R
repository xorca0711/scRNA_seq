# Sensitivity fits under config/second_batch.md; original primary outputs immutable.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
base <- normalizePath(args[1],winslash='/')
out <- normalizePath(args[2],winslash='/')
.libPaths(c(file.path(getwd(),'analysis/corrections/statistics/.tools/R-library'),.libPaths()))
suppressPackageStartupMessages(library(edgeR))
meta <- read.delim(file.path(base,'tables/ire1/sample_manifest.tsv'),check.names=FALSE)
x <- as.matrix(read.delim(file.path(base,'processed/ire1/counts.tsv'),row.names=1,check.names=FALSE))
primary <- read.delim(file.path(base,'tables/ire1/gene_effects.tsv'))
sets <- read.delim(file.path(base,'tables/ire1/gene_set_mapping_tested.tsv'))
coverage <- read.delim(file.path(base,'tables/ire1/gene_set_coverage.tsv'))
stopifnot(identical(colnames(x),meta$sample_id),!anyDuplicated(meta$mouse),nrow(meta)==10,
          all(table(meta$group)==5),nrow(primary)==14811,all(primary$feature_id %in% rownames(x)))
x <- x[primary$feature_id,,drop=FALSE]
markers <- read.delim(file.path(base,'tables/ire1/marker_mapping.tsv'))
hit <- primary$feature_id[primary$FDR<.05]
focus <- unique(c(markers$feature_id,hit))
write_tsv <- function(d,n) write.table(d,file.path(out,n),sep='\t',quote=FALSE,row.names=FALSE)
results <- list(); diagnostics <- list(); pathways <- list()
do_fit <- function(name,selected,formula,kind,omitted='') {
  m <- droplevels(meta[selected,,drop=FALSE])
  m$group <- factor(m$group,levels=c('Vehicle','KIRA8'))
  m$sex <- factor(m$sex); m$batch <- factor(m$batch)
  design <- model.matrix(formula,m)
  stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design),'groupKIRA8' %in% colnames(design),
            all(table(m$group)>=3))
  counts <- x[,selected,drop=FALSE]
  stopifnot(all(rowSums(counts)>0))
  y <- calcNormFactors(DGEList(counts),method='TMM')
  d <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(d,design,robust=TRUE)
  result <- topTags(glmQLFTest(fit,coef=which(colnames(design)=='groupKIRA8')),n=Inf,sort.by='none')$table
  result$feature_id <- rownames(result); result$run <- name; result$kind <- kind; result$omitted_mouse <- omitted
  stopifnot(identical(result$feature_id,primary$feature_id))
  results[[name]] <<- result
  diagnostics[[name]] <<- data.frame(run=name,kind=kind,omitted_mouse=omitted,n=nrow(m),
      Vehicle=sum(m$group=='Vehicle'),KIRA8=sum(m$group=='KIRA8'),rank=qr(design)$rank,
      residual_df=nrow(design)-ncol(design),genes=nrow(result),genes_FDR05=sum(result$FDR<.05),
      Pearson_to_primary=cor(result$logFC,primary$logFC),
      Spearman_to_primary=cor(result$logFC,primary$logFC,method='spearman'))
  ix <- lapply(coverage$gene_set[coverage$eligible],function(s) which(rownames(y) %in% sets$feature_id[sets$gene_set==s & sets$retained]))
  names(ix) <- coverage$gene_set[coverage$eligible]
  v <- voom(y,design,plot=FALSE)
  cam <- camera(v,ix,design,contrast=which(colnames(design)=='groupKIRA8'),inter.gene.cor=NA,allow.neg.cor=FALSE)
  cam$FDR <- p.adjust(cam$PValue,'BH'); cam$gene_set <- rownames(cam); cam$run <- name; cam$kind <- kind
  pathways[[name]] <<- cam
  write_tsv(cbind(m,y$samples[,c('lib.size','norm.factors')]),paste0(name,'_normalization.tsv'))
  cat(name,' completed: ',nrow(m),' mice, ',sum(result$FDR<.05),' genes FDR<0.05\n',sep=''); flush.console()
}
for (mouse in meta$mouse) do_fit(paste0('omit_',mouse),meta$mouse!=mouse,~batch+sex+group,'leave_one_mouse_out',as.character(mouse))
do_fit('S061',meta$batch=='S061',~sex+group,'within_batch')
# The smaller batch remains descriptive under the existing sample-floor contract.
small <- meta$batch=='S135'; y <- calcNormFactors(DGEList(x[,small,drop=FALSE]),method='TMM')
values <- cpm(y)
group <- meta$group[small]
desc <- data.frame(feature_id=rownames(y),Vehicle_mean_CPM=rowMeans(values[,group=='Vehicle',drop=FALSE]),
                    KIRA8_mean_CPM=rowMeans(values[,group=='KIRA8',drop=FALSE]))
desc$log2_ratio_offset_0.5 <- log2((desc$KIRA8_mean_CPM+.5)/(desc$Vehicle_mean_CPM+.5))
write_tsv(desc,'S135_descriptive.tsv')
all <- do.call(rbind,results)
write_tsv(all,'all_gene_sensitivity.tsv')
write_tsv(do.call(rbind,diagnostics),'fit_diagnostics.tsv')
write_tsv(do.call(rbind,pathways),'pathway_sensitivity.tsv')
write_tsv(all[all$feature_id %in% focus,],'focus_effects.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
