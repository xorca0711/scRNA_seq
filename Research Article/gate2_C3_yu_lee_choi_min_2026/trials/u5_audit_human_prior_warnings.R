# Reconstruct primary prior scores and identify undefined constant-column fits.
root <- normalizePath('.',winslash='/')
paper <- file.path(root,'Research Article/gate2_C3_yu_lee_choi_min_2026')
out <- file.path(paper,'trials/u5_human_ligand_targets')
cache <- file.path(paper,'cache/u5_human_niche')
prior <- readRDS(file.path(paper,'cache/u5_ligand_targets/ligand_target_matrix_nsga2r_final.rds'))
gate <- read.csv(file.path(out,'target_eligibility.csv'))
gate <- gate[gate$omitted_patient=='__none__' & gate$eligible,,drop=FALSE]
targets <- read.csv(file.path(out,'full_fit_mapped_targets.csv'))
scores <- read.csv(gzfile(file.path(out,'all_prior_pearson_scores.csv.gz')))
scores <- scores[scores$omitted_patient=='__none__',,drop=FALSE]
rows <- list()
for(i in seq_len(nrow(gate))) {
 g <- gate[i,]
 same <- function(x) x$comp==g$comp & x$label==g$label & x$case==g$case & x$reference==g$reference & x$direction==g$direction
 target <- targets$gene[same(targets)]
 stem <- gsub('[^A-Za-z0-9_]','_',paste(g$comp,g$label,g$case,g$reference,sep='_'))
 de <- read.csv(gzfile(file.path(cache,paste0(stem,'_DE.csv.gz'))))
 bg <- intersect(de$gene,rownames(prior)); pp <- prior[bg,,drop=FALSE]
 response <- as.numeric(bg%in%target)
 stopifnot(sum(response)==g$mapped_targets,length(bg)==g$mapped_background,sd(response)>0)
 messages <- character()
 r <- as.numeric(withCallingHandlers(cor(response,pp,method='pearson'),warning=function(w) {
  messages <<- c(messages,conditionMessage(w)); invokeRestart('muffleWarning')
 }))
 constant <- vapply(seq_len(ncol(pp)),function(j) all(pp[,j]==pp[1,j]),logical(1))
 stopifnot(identical(!is.finite(r),constant))
 old <- scores[same(scores),]; stopifnot(!anyDuplicated(old$ligand),nrow(old)==ncol(pp))
 old <- old[match(colnames(pp),old$ligand),]
 stopifnot(identical(is.na(old$pearson),!is.finite(r)),max(abs(old$pearson[!constant]-r[!constant]))<1e-12)
 stopifnot(all(messages=='the standard deviation is zero'))
 rows[[length(rows)+1]] <- data.frame(comp=g$comp,label=g$label,case=g$case,reference=g$reference,direction=g$direction,
   mapped_background=length(bg),mapped_targets=sum(response),prior_ligands=length(r),constant_prior_columns=sum(constant),
   nonfinite_exactly_constant=TRUE,maximum_finite_score_difference=max(abs(old$pearson[!constant]-r[!constant])),
   reconstructed_warning_messages=paste(unique(messages),collapse=';'))
}
write.csv(do.call(rbind,rows),file.path(out,'primary_prior_warning_audit.csv'),row.names=FALSE)
cat('Primary prior-score reconstructions passed:',length(rows),'\n')
