# Native R implementation of the HypoForge core. No Python runtime required.
#' Create a feature specification
#' @param name Feature name.
#' @param op Transformation identifier.
#' @param inputs Source column names.
#' @param params Named parameter list.
#' @return A HypoForge feature specification.
#' @export
hf_spec <- function(name, op, inputs, params=list()) {
  structure(list(name=name, op=op, inputs=inputs, params=params), class='hf_spec')
}

#' Evaluate a feature specification
#' @param spec A specification created with hf_spec.
#' @param data Data frame.
#' @return A vector with one element per row.
#' @export
hf_evaluate <- function(spec, data) {
  x <- data[[spec$inputs[[1]]]]; p <- spec$params; op <- spec$op
  if (is.null(x)) stop('Missing source column')
  if (op=='identity') return(x)
  if (op=='log1p_abs') return(sign(x)*log1p(abs(x)))
  if (op=='sqrt_abs') return(sign(x)*sqrt(abs(x)))
  if (op=='square') return(x*x)
  if (op=='cube') return(x*x*x)
  if (op=='abs') return(abs(x))
  if (op=='reciprocal') return(ifelse(x==0, NA_real_, 1/x))
  if (op=='hill') return(ifelse(x>=0, x^p$n/(p$k^p$n+x^p$n), NA_real_))
  if (op=='is_missing') return(as.integer(is.na(x)))
  if (op=='threshold_gt') return(as.integer(x>p$value))
  if (op=='threshold_lt') return(as.integer(x<p$value))
  if (op=='range') return(as.integer(x>=p$low & x<=p$high))
  if (op=='category_eq') return(as.integer(x==p$value))
  if (op=='weekday') return(as.POSIXlt(as.Date(x))$wday)
  if (op=='month') return(as.POSIXlt(as.Date(x))$mon+1L)
  if (op=='dayofyear') return(as.POSIXlt(as.Date(x))$yday+1L)
  if (op=='is_weekend') return(as.integer(as.POSIXlt(as.Date(x))$wday %in% c(0,6)))
  lag <- function(v, k) {if (k<1 || k>=length(v)) return(rep(NA_real_,length(v))); c(rep(NA_real_,k),head(v,-k))}
  if (op=='lag') return(lag(x,p$periods))
  if (op=='diff') return(x-lag(x,p$periods))
  if (op=='pct_change') {v<-lag(x,p$periods); return(ifelse(v==0,NA_real_,(x-v)/v))}
  if (op=='ewm_mean') {
    z<-lag(x,1); a<-2/(p$span+1); out<-rep(NA_real_,length(z)); prev<-NA_real_
    for (i in seq_along(z)) {if (!is.na(z[i])) prev<-if(is.na(prev)) z[i] else a*z[i]+(1-a)*prev;out[i]<-prev}
    return(out)
  }
  if (startsWith(op,'rolling_')) {
    w<-p$window; z<-lag(x,1); out<-rep(NA_real_,length(z))
    for (i in seq_along(z)) {
      v<-z[max(1,i-w+1):i];v<-v[is.finite(v)]
      if(length(v)<max(2,floor(w/3))) next
      out[i]<-switch(op,rolling_mean=mean(v),rolling_std=sd(v),rolling_min=min(v),
        rolling_max=max(v),rolling_median=median(v),rolling_range=max(v)-min(v))
    }
    return(out)
  }
  y<-data[[spec$inputs[[2]]]]
  if (is.null(y)) stop('Missing second source column')
  switch(op,add=x+y,subtract=x-y,abs_diff=abs(x-y),multiply=x*y,
    ratio=ifelse(y==0,NA_real_,x/y),relative_diff=ifelse(y==0,NA_real_,(x-y)/abs(y)),
    and=as.integer(as.logical(x)&as.logical(y)),stop('Unknown operation'))
}

#' Generate candidate transformations
#' @param data Data frame.
#' @param exclude Column names excluded from generation.
#' @param lags Integer lag sizes.
#' @param windows Integer rolling window sizes.
#' @return Named list of feature specifications.
#' @export
hf_generate <- function(data, exclude=character(), lags=c(1L,7L,30L), windows=c(7L,30L)) {
  out<-list(); add<-function(s) {out[[s$name]] <<- s}
  for (nm in setdiff(names(data),exclude)) {
    x<-data[[nm]]; add(hf_spec(nm,'identity',nm));add(hf_spec(paste0(nm,'__missing'),'is_missing',nm))
    if (is.numeric(x)) {
      for(op in c('log1p_abs','sqrt_abs','square','cube','reciprocal','abs')) add(hf_spec(paste(nm,op,sep='__'),op,nm))
      for(k in lags) for(op in c('lag','diff','pct_change')) add(hf_spec(paste(nm,op,k,sep='__'),op,nm,list(periods=k)))
      for(w in windows) {
        for(op in c('rolling_mean','rolling_std','rolling_min','rolling_max','rolling_median','rolling_range'))
          add(hf_spec(paste(nm,op,w,sep='__'),op,nm,list(window=w)))
        add(hf_spec(paste(nm,'ewm_mean',w,sep='__'),'ewm_mean',nm,list(span=w)))
      }
      for(q in c(.1,.25,.5,.75,.9)) {
        v<-as.numeric(quantile(x,q,na.rm=TRUE))
        for(op in c('threshold_gt','threshold_lt')) add(hf_spec(paste(nm,op,q,sep='__'),op,nm,list(value=v)))
      }
      pos<-x[is.finite(x)&x>0]
      if(length(pos)>0) for(k in unique(as.numeric(quantile(pos,c(.25,.5,.75))))) for(n in c(.5,1,2,3))
        add(hf_spec(paste(nm,'hill',k,n,sep='__'),'hill',nm,list(k=k,n=n)))
    } else if (inherits(x,'Date') || inherits(x,'POSIXt')) {
      for(op in c('weekday','month','dayofyear','is_weekend')) add(hf_spec(paste(nm,op,sep='__'),op,nm))
    } else {
      for(v in head(names(sort(table(x),decreasing=TRUE)),20))
        add(hf_spec(paste(nm,'eq',v,sep='__'),'category_eq',nm,list(value=v)))
    }
  }
  out
}

#' Numeric safety diagnostics
#' @param spec Feature specification.
#' @param data Data frame.
#' @return List with numerical safety and finite fraction.
#' @export
hf_safety <- function(spec,data) {
  x<-suppressWarnings(as.numeric(hf_evaluate(spec,data)))
  fraction<-mean(is.finite(x));v<-x[is.finite(x)]
  list(numerically_safe=length(v)>1 && length(unique(v))>1 && fraction>=.8 && max(abs(v))<1e150,
       finite_fraction=fraction)
}

#' Rank exploratory feature hypotheses by absolute correlation
#' @param data Data frame.
#' @param target Numeric target name.
#' @param top_k Maximum returned features.
#' @return Data frame of ranked candidates. Correlation is NOT out-of-sample evidence.
#' @export
hf_discover <- function(data,target,top_k=50L) {
  if(!is.numeric(data[[target]])) stop('Numeric target required')
  specs<-hf_generate(data,exclude=target)
  rows<-lapply(specs,function(s) {
    safe<-tryCatch(hf_safety(s,data),error=function(e) NULL)
    if(is.null(safe)||!safe$numerically_safe) return(NULL)
    x<-suppressWarnings(as.numeric(hf_evaluate(s,data)))
    z<-is.finite(x)&is.finite(data[[target]])
    if(sum(z)<10) return(NULL)
    data.frame(feature=s$name,score=abs(cor(x[z],data[[target]][z])),stringsAsFactors=FALSE)
  })
  ans<-do.call(rbind,rows);if(is.null(ans)) return(data.frame(feature=character(),score=numeric()))
  head(ans[order(-ans$score),,drop=FALSE],top_k)
}

#' Conservative redundancy screening
#' @param data Data frame.
#' @param ranked Feature names or hf_discover result.
#' @param baseline_features Baseline numeric columns.
#' @param max_features Maximum candidates.
#' @param threshold Maximum absolute correlation.
#' @return List containing selected names and audit table.
#' @export
hf_select <- function(data,ranked,baseline_features=character(),max_features=5L,threshold=.98) {
  specs<-hf_generate(data); names_in<-if(is.data.frame(ranked)) ranked$feature else ranked
  selected<-character(); audit<-list(); base<-lapply(baseline_features,function(nm) data[[nm]])
  for(nm in names_in) {
    s<-specs[[nm]];if(is.null(s)) next
    x<-tryCatch(as.numeric(hf_evaluate(s,data)),error=function(e) NULL)
    if(is.null(x)) next
    correlations<-vapply(c(base,lapply(selected,function(nm2) as.numeric(hf_evaluate(specs[[nm2]],data)))),
      function(y) {ok<-is.finite(x)&is.finite(y);if(sum(ok)<3 || sd(x[ok])==0 || sd(y[ok])==0) return(0);abs(cor(x[ok],y[ok]))},numeric(1))
    mc<-if(length(correlations)) max(correlations) else 0
    reason<-if(mc>=threshold) 'redundant_correlation' else 'selected'
    audit[[length(audit)+1L]]<-data.frame(feature=nm,decision=reason,max_abs_correlation=mc)
    if(reason=='selected') selected<-c(selected,nm)
    if(length(selected)>=max_features) break
  }
  list(selected=selected,audit=if(length(audit)) do.call(rbind,audit) else data.frame())
}
