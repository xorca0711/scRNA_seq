# Small paired-unit inference helpers. No cells are treated as replicates.
mean_test <- function(x) {
  stopifnot(length(x) >= 3, all(is.finite(x)))
  z <- t.test(x, mu=0)
  data.frame(n=length(x), mean=mean(x), sd=sd(x), mean_low=z$conf.int[1],
             mean_high=z$conf.int[2], p_mean=z$p.value)
}

location_test <- function(x, margin=.10) {
  m <- mean_test(x)
  warnings <- character()
  z <- tryCatch(withCallingHandlers(
    wilcox.test(x, mu=0, alternative='two.sided', exact=TRUE, conf.int=TRUE,
                conf.level=.95, tol.root=1e-10),
    warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning') }),
    error=function(e) e)
  available <- !inherits(z,'error') && grepl('exact',z$method,ignore.case=TRUE) &&
    is.finite(z$p.value) && length(z$conf.int)==2 && all(is.finite(z$conf.int)) &&
    is.finite(z$estimate) && !any(grepl('cannot compute|not achievable', warnings))
  if (!available) return(cbind(m,data.frame(HL=NA_real_,low=NA_real_,high=NA_real_,
    p_exact=NA_real_,exact_available=FALSE,direction='unavailable',magnitude='unavailable',
    inference_note=paste(c(if(inherits(z,'error')) conditionMessage(z) else z$method,warnings),collapse='; '))))
  direction <- if(z$p.value < .05 && z$estimate > 0) 'positive' else
    if(z$p.value < .05 && z$estimate < 0) 'negative' else 'unresolved'
  magnitude <- if(z$conf.int[1] > margin) 'meaningful_positive_supported' else
    if(z$conf.int[2] < margin) 'positive_margin_ruled_out' else 'unresolved'
  cbind(m,data.frame(HL=unname(z$estimate),low=z$conf.int[1],high=z$conf.int[2],
    p_exact=z$p.value,exact_available=TRUE,direction,magnitude,
    inference_note=paste(c(z$method,warnings),collapse='; ')))
}
