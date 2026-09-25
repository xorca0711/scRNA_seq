args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==1)
root <- normalizePath(args[1],winslash='/')
source(file.path(root,'RQ_Specified/A5_A11_shared_component_contract/scripts/paired_inference.R'))
out <- file.path(root,'RQ_Specified/A5_developmental_programme_reuse/tables/test_v1')
if(file.exists(file.path(out,'inference.tsv'))) stop('Refusing to overwrite A5 inference')
z <- read.delim(file.path(out,'paired_differences.tsv'))
coverage <- read.delim(file.path(out,'coverage.tsv'))
tests <- data.frame(test=c('primary','identity_excluded','identity_controls_excluded','resting_reference'),
  module=c('Guo_AT1_AT2_external','Guo_minus_identity','Guo_minus_identity_and_controls','Guo_AT1_AT2_external'),
  reference=c(rep('AT2 activated',3),'AT2'))
rows <- list()
for(j in seq_len(nrow(tests))) {
  d <- z$difference_pp[z$module==tests$module[j] & z$reference==tests$reference[j]]
  pass <- length(d)>=3 && coverage$eligible[coverage$module==tests$module[j]]
  r <- if(pass) mean_test(d) else data.frame(n=length(d),mean=NA,sd=NA,mean_low=NA,mean_high=NA,p_mean=NA)
  rows[[j]] <- cbind(tests[j,],eligible=pass,r)
}
r <- do.call(rbind,rows); r$q_secondary <- NA_real_
r$q_secondary[-1] <- p.adjust(ifelse(is.finite(r$p_mean[-1]),r$p_mean[-1],1),'holm')
r$direction <- ifelse(!r$eligible,'ineligible',ifelse(r$mean_low>0,'positive',ifelse(r$mean_high<0,'negative','unresolved')))
write.table(r,file.path(out,'inference.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
d <- z[z$module=='Guo_AT1_AT2_external' & z$reference=='AT2 activated',]
loo <- data.frame(omitted_mouse=d$sample_id,remaining_mean=sapply(seq_len(nrow(d)),function(j)mean(d$difference_pp[-j])))
write.table(loo,file.path(out,'omission_diagnostics.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
days <- aggregate(difference_pp~day,data=d,FUN=mean)
write.table(days,file.path(out,'day_means.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session.txt'))
print(r[,c('test','n','mean','mean_low','mean_high','p_mean','q_secondary','direction')],row.names=FALSE)
cat('Equal-day mean:',mean(days$difference_pp),'pp\n')
