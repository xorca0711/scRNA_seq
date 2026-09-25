# Retrospective sensitivity/correction. Never overwrites the historical trials.
# Actual limma::camera, not a Python approximation. Estimated residual correlation
# is explicit; this requests CAMERA's residual-df treatment of correlation error.
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) args[[1]] else "mouse"
here <- normalizePath("analysis/corrections/statistics", winslash="/")
root <- normalizePath(".", winslash="/")
.libPaths(c(file.path(here, ".tools/R-library"), .libPaths()))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
options(stringsAsFactors=FALSE)
cache <- file.path(here, "cache")
out <- file.path(here, "tables")
dir.create(out, showWarnings=FALSE)
trials <- file.path(root, "Thesis/gate1_01_niethamer_2025/trials")

read_gmt <- function(path) {
  x <- strsplit(readLines(path), "\t", fixed=TRUE)
  setNames(lapply(x, function(a) unique(a[-c(1,2)])), vapply(x, `[[`, "", 1))
}
read_bundle <- function(name) {
  x <- as.matrix(read.csv(gzfile(file.path(cache,paste0(name,"_counts.csv.gz"))), row.names=1, check.names=FALSE))
  m <- read.csv(file.path(cache,paste0(name,"_units.csv")))
  stopifnot(identical(colnames(x),m$unit_id))
  list(x=x,m=m)
}
collections <- function(species) {
  version <- if (species=="Mm") "2024.1.Mm" else "2024.1.Hs"
  first <- if (species=="Mm") "mh.all.v" else "h.all.v"
  second <- if (species=="Mm") "m5.go.bp.v" else "c5.go.bp.v"
  c(read_gmt(file.path(root,"raw_data/msigdb",paste0(first,version,".symbols.gmt"))),
    read_gmt(file.path(root,"raw_data/msigdb",paste0(second,version,".symbols.gmt"))))
}
all_human <- collections("Hs")
all_mouse <- collections("Mm")
frozen <- read.csv(file.path(trials,"g2_gsea_ipf/g2_replication.csv"))
focus_mouse <- c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT",
  "GOBP_MITOTIC_DNA_REPLICATION","GOBP_CELL_CYCLE_DNA_REPLICATION",
  "GOBP_ARGININE_METABOLIC_PROCESS","GOBP_ORNITHINE_METABOLIC_PROCESS")
design_audit <- list()

fit_camera <- function(x,m,sets,formula=~arm,method="voom_TMM",minsize=15, fixed_genes=NULL,correlation=NA) {
  m <- droplevels(m)
  design <- model.matrix(formula,m)
  if (qr(design)$rank < ncol(design)) stop("rank-deficient design")
  if (nrow(design)-ncol(design)<2) stop("fewer than two residual df")
  k <- which(colnames(design)=="arm")
  stopifnot(length(k)==1)
  keep <- if(is.null(fixed_genes)) rowSums(x>=10)>=3 else rownames(x) %in% fixed_genes
  x <- x[keep,,drop=FALSE]
  if (method=="historical_logCPM") {
    y <- log2(t(t(x)/colSums(x))*1e6+1)
  } else {
    d <- calcNormFactors(DGEList(counts=x),method="TMM")
    y <- voom(d,design=design,plot=FALSE)
  }
  ix <- ids2indices(sets,rownames(x),remove.empty=FALSE)
  sizes <- lengths(ix)
  ix <- ix[sizes>=minsize & sizes<=500]
  if (!length(ix)) stop("no eligible sets")
  z <- camera(y,ix,design=design,contrast=k,inter.gene.cor=correlation,
              allow.neg.cor=FALSE,use.ranks=FALSE,trend.var=FALSE,sort=FALSE)
  if (!"Correlation" %in% names(z)) z$Correlation <- correlation
  z$correlation_setting <- if(is.na(correlation)) "estimated" else paste0("fixed_",correlation)
  z$set <- rownames(z)
  z$n_first <- sum(m$arm==1)
  z$n_second <- sum(m$arm==0)
  z$residual_df <- nrow(design)-ncol(design)
  z$genes_tested <- nrow(x)
  z$method <- method
  z$FDR_family <- p.adjust(z$PValue,method="BH")
  rownames(z) <- NULL
  list(table=z,genes=rownames(x),design=design)
}

audit <- function(label,m,formula) {
  d <- model.matrix(formula,droplevels(m))
  design_audit[[length(design_audit)+1]] <<- data.frame(label=label,
    formula=paste(deparse(formula),collapse=""),n=nrow(d),columns=ncol(d),rank=qr(d)$rank,
    residual_df=nrow(d)-qr(d)$rank,full_rank=qr(d)$rank==ncol(d),
    n_first=sum(m$arm==1),n_second=sum(m$arm==0))
}

mouse <- function(correlation=NA,tag="",only_g1=FALSE) {
  am <- read_bundle("mouse_amac")
  am$m$batch <- factor(am$m$round)
  am$m$sex <- factor(am$m$sex)
  am$m$phase <- ifelse(am$m$day %in% c(6,11,19,25),"active",ifelse(am$m$day %in% c(42,90),"resolution",ifelse(am$m$day==366,"longterm","baseline")))
  am$m$homo <- grepl("Cre/Cre",am$m$genotype,fixed=TRUE)
  wsets <- lapply(all_mouse[focus_mouse],setdiff,y="Mki67")
  wr <- list()
  for (contrast in if(only_g1) character(0) else c("A","B")) {
    phases <- if(contrast=="A") c("active","resolution") else c("resolution","longterm")
    sel <- am$m$cells>=50 & !am$m$homo & am$m$phase %in% phases
    m <- am$m[sel,]
    m$arm <- as.numeric(m$phase==phases[1])
    audit(paste0("W1_",contrast),m,~batch+sex+arm)
    for (method in c("historical_logCPM","voom_TMM")) {
      z <- fit_camera(am$x[,sel,drop=FALSE],m,wsets,~batch+sex+arm,method,minsize=10)$table
      z$contrast <- contrast
      wr[[length(wr)+1]] <- z
      message("W1 ",contrast," ",method,": ",sum(z$FDR_family<0.05)," / ",nrow(z)," hits")
    }
  }
  if(length(wr)) write.csv(do.call(rbind,wr),file.path(out,"w1_reference_camera.csv"),row.names=FALSE)
  my <- read_bundle("mouse_myeloid")
  my$m$batch <- factor(my$m$round)
  my$m$sex <- factor(my$m$sex)
  my$m$homo <- as.numeric(grepl("Cre/Cre",my$m$genotype,fixed=TRUE))
  my$m$arm <- as.numeric(my$m$day %in% c(42,90))
  base <- my$m$cells>=50 & my$m$day %in% c(42,90,366)
  specs <- list(all_unadjusted=list(sel=base,f=~arm),
    all_adjusted=list(sel=base,f=~batch+sex+homo+arm),
    heterozygote_adjusted=list(sel=base & my$m$homo==0,f=~batch+sex+arm),
    december_heterozygote=list(sel=base & my$m$homo==0 & my$m$round=="2022-12-06",f=~sex+arm))
  gr <- list()
  for (label in names(specs)) {
    sp <- specs[[label]]
    m <- my$m[sp$sel,]
    audit(paste0("G1_",label),m,sp$f)
    z <- fit_camera(my$x[,sp$sel,drop=FALSE],m,all_mouse,sp$f,correlation=correlation)$table
    z$specification <- label
    write.csv(z,gzfile(file.path(cache,paste0("g1_",label,tag,"_allsets.csv.gz"))),row.names=FALSE)
    gr[[label]] <- z[z$set %in% focus_mouse,]
    message("G1 ",label,": ",sum(z$FDR_family<0.05)," / ",nrow(z)," hits")
  }
  # Age/time is not separable from phase plus infection round in these contrasts.
  audit("G1_age_time_not_identifiable",my$m[base,],~batch+sex+homo+arm+day)
  male <- base & my$m$homo==0 & my$m$round=="2022-12-06" & my$m$sex=="M"
  audit("G1_same_round_male_only_below_arm_floor",my$m[male,],~arm)
  write.csv(do.call(rbind,gr),file.path(out,paste0("g1_confounded_sensitivity",tag,".csv")),row.names=FALSE)
  if(!only_g1) write.csv(do.call(rbind,design_audit),file.path(out,"mouse_design_audit.csv"),row.names=FALSE)
}

human <- function(cohort,correlation=NA,tag="") {
  d <- read_bundle(cohort)
  full <- list()
  for (comp in unique(d$m$comp)) {
    sel <- d$m$comp==comp
    m <- d$m[sel,]
    m$arm <- as.numeric(m$disease=="IPF")
    candidates <- frozen$set[frozen$compartment==comp]
    sets <- if(cohort=="GSE136831") all_human else all_human[candidates]
    z <- fit_camera(d$x[,sel,drop=FALSE],m,sets,correlation=correlation)$table
    z$cohort <- cohort
    z$compartment <- comp
    z$frozen_candidate <- z$set %in% candidates
    full[[comp]] <- z
    message(cohort," ",comp,": ",sum(z$FDR_family<0.05)," / ",nrow(z)," hits")
  }
  z <- do.call(rbind,full)
  # All tested discovery hypotheses enter discovery BH, not just selected hits.
  # Validation BH covers the entire frozen 254-candidate family across compartments.
  z$FDR_global <- p.adjust(z$PValue,method="BH",n=if(cohort=="GSE136831") nrow(z) else nrow(frozen))
  write.csv(z,gzfile(file.path(cache,paste0(cohort,tag,"_camera_allsets.csv.gz"))),row.names=FALSE)
  candidates <- merge(frozen[,c("compartment","set")],z[z$frozen_candidate,],
                      by=c("compartment","set"),all.x=TRUE,sort=FALSE)
  candidates$eligible <- !is.na(candidates$PValue)
  candidates$cohort <- cohort
  candidates$not_tested_reason <- ifelse(candidates$eligible,"","set intersection outside 15 to 500 genes")
  write.csv(candidates,file.path(out,paste0(cohort,tag,"_camera_candidates.csv")),row.names=FALSE)
}

lodo <- function(cohort) {
  d <- read_bundle(cohort)
  reference <- read.csv(file.path(out,paste0(cohort,"_camera_candidates.csv")))
  all <- list()
  for (comp in unique(d$m$comp)) {
    sel <- d$m$comp==comp
    x <- d$x[,sel,drop=FALSE]
    m <- d$m[sel,]
    m$arm <- as.numeric(m$disease=="IPF")
    candidates <- frozen$set[frozen$compartment==comp]
    sets <- all_human[candidates]
    genes <- rownames(x)[rowSums(x>=10)>=3]
    for (i in seq_len(nrow(m))) {
      z <- fit_camera(x[,-i,drop=FALSE],m[-i,],sets,fixed_genes=genes)$table
      z$cohort <- cohort
      z$compartment <- comp
      z$omitted_donor <- m$donor[i]
      z$below_original_arm_floor <- min(z$n_first[1],z$n_second[1])<3
      ref <- reference[reference$compartment==comp,]
      z$direction_matches_full <- z$Direction==ref$Direction[match(z$set,ref$set)]
      all[[length(all)+1]] <- z
    }
    message(cohort," LODO ",comp," complete")
  }
  z <- do.call(rbind,all)
  write.csv(z,gzfile(file.path(cache,paste0(cohort,"_lodo_all.csv.gz"))),row.names=FALSE)
  groups <- split(seq_len(nrow(z)),paste(z$compartment,z$set,sep="|"))
  summary <- lapply(groups,function(i) data.frame(cohort=cohort,compartment=z$compartment[i[1]],set=z$set[i[1]],
    omissions=length(i),direction_reversals=sum(!z$direction_matches_full[i]),
    min_p=min(z$PValue[i]),max_p=max(z$PValue[i]),max_FDR_candidate_family=max(z$FDR_family[i]),
    below_floor_omissions=sum(z$below_original_arm_floor[i]),
    worst_donor=z$omitted_donor[i[which.max(z$PValue[i])]]))
  write.csv(do.call(rbind,summary),file.path(out,paste0(cohort,"_lodo_summary.csv")),row.names=FALSE)
}

subtypes <- function(cohort) {
  d <- read_bundle(paste0(cohort,"_subtypes"))
  all <- list()
  auditrows <- list()
  for (comp in unique(d$m$comp)) for (label in unique(d$m$label[d$m$comp==comp])) {
    sel <- d$m$comp==comp & d$m$label==label & d$m$cells>=50
    m <- d$m[sel,]
    m$arm <- as.numeric(m$disease=="IPF")
    enough <- min(sum(m$arm==1),sum(m$arm==0))>=3
    auditrows[[length(auditrows)+1]] <- data.frame(cohort=cohort,compartment=comp,label=label,
      n_IPF=sum(m$arm==1),n_control=sum(m$arm==0),tested=enough)
    if (!enough) next
    sets <- all_human[frozen$set[frozen$compartment==comp]]
    z <- fit_camera(d$x[,sel,drop=FALSE],m,sets)$table
    z$cohort <- cohort
    z$compartment <- comp
    z$label <- label
    all[[length(all)+1]] <- z
  }
  write.csv(do.call(rbind,auditrows),file.path(out,paste0(cohort,"_subtype_eligibility.csv")),row.names=FALSE)
  z <- do.call(rbind,all)
  z$FDR_exploratory_subtype_family <- p.adjust(z$PValue,"BH")
  write.csv(z,file.path(out,paste0(cohort,"_subtype_camera.csv")),row.names=FALSE)
}

if(mode=="mouse") mouse()
if(mode %in% c("GSE136831","GSE135893")) human(mode)
if(mode=="lodo") for(cohort in c("GSE136831","GSE135893")) lodo(cohort)
if(mode=="subtypes") for(cohort in c("GSE136831","GSE135893")) subtypes(cohort)
if(mode=="fixed001") {
  for(cohort in c("GSE136831","GSE135893")) human(cohort,correlation=0.01,tag="_fixed001")
  mouse(correlation=0.01,tag="_fixed001",only_g1=TRUE)
}
writeLines(capture.output(sessionInfo()),file.path(here,"R-session.txt"))
