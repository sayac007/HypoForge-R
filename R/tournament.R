#' Rolling-origin cross-validation splits
#' @param n Number of rows.
#' @param min_train Minimum initial training length.
#' @param horizon Validation length.
#' @param step Step between origins.
#' @param gap Embargo between train and validation.
#' @return List of train/test index lists.
#' @export
hf_rolling_splits <- function(n,min_train,horizon,step=horizon,gap=0L) {
  if(min_train<2 || horizon<1 || step<1 || gap<0) stop('Invalid split parameters')
  ends<-seq(min_train,n-gap-horizon,by=step)
  lapply(ends,function(e) list(train=seq_len(e),test=seq.int(e+gap+1,e+gap+horizon)))
}

#' Compare candidate features to a fixed ridge regression baseline
#' @param data Time-ordered data frame.
#' @param target Numeric target column.
#' @param time Unique time column.
#' @param baseline_features Known-future or static numeric predictors.
#' @param candidates Named list of hf_spec objects.
#' @param availability Named vector mapping each predictor to known_future or static.
#' @param min_train,horizon,step,gap Rolling split configuration.
#' @param lambda Ridge penalty.
#' @return List of summary and per-fold results. NOT an unbiased final test.
#' @export
hf_tournament <- function(data,target,time,baseline_features,candidates,
                          availability, min_train=104L,horizon=13L,step=horizon,gap=0L,lambda=10) {
  if(anyDuplicated(data[[time]]) || anyNA(data[[time]])) stop('Time must be unique and complete')
  if(!is.numeric(data[[target]]) || anyNA(data[[target]])) stop('Numeric complete target required')
  if(!length(baseline_features)) stop('Baseline required')
  data<-data[order(data[[time]]),,drop=FALSE]
  allowed<-function(cols) all(cols %in% names(data)) &&
    all(availability[cols] %in% c('known_future','static')) && !anyNA(availability[cols]) &&
    !any(cols %in% c(target,time))
  if(!allowed(baseline_features)) stop('Baseline input not known at prediction time')
  banned<-c('hill','threshold_gt','threshold_lt','range','category_eq','lag','diff','pct_change',
    'ewm_mean','rolling_mean','rolling_std','rolling_min','rolling_max','rolling_median','rolling_range')
  candidates<-Filter(function(s) !(s$op %in% banned) && allowed(s$inputs),candidates)
  splits<-hf_rolling_splits(nrow(data),min_train,horizon,step,gap)
  if(length(splits)<3) stop('At least three folds required')
  fit_pred<-function(X,y,tr,te) {
    a<-X[tr,,drop=FALSE];b<-X[te,,drop=FALSE]
    for(j in seq_len(ncol(a))) {
      m<-median(a[,j],na.rm=TRUE);if(!is.finite(m)) stop('No finite training data')
      a[!is.finite(a[,j]),j]<-m;b[!is.finite(b[,j]),j]<-m
      center<-mean(a[,j]);sc<-sd(a[,j]);if(!is.finite(sc)||sc<1e-12) sc<-1
      a[,j]<-(a[,j]-center)/sc;b[,j]<-(b[,j]-center)/sc
    }
    A<-cbind(1,a);B<-cbind(1,b);pen<-diag(c(0,rep(lambda,ncol(a))))
    coef<-solve(crossprod(A)+pen,crossprod(A,y[tr]));as.vector(B%*%coef)
  }
  X<-as.matrix(data[,baseline_features,drop=FALSE]);storage.mode(X)<-'double'
  y<-data[[target]];rows<-list()
  for(i in seq_along(splits)) {
    tr<-splits[[i]]$train;te<-splits[[i]]$test
    base<-sqrt(mean((y[te]-fit_pred(X,y,tr,te))^2))
    for(s in candidates) {
      v<-tryCatch(as.numeric(hf_evaluate(s,data)),error=function(e) NULL)
      if(is.null(v) || sum(is.finite(v[tr]))<3 || sd(v[tr],na.rm=TRUE)<1e-12) next
      score<-tryCatch(sqrt(mean((y[te]-fit_pred(cbind(X,v),y,tr,te))^2)),error=function(e) NA_real_)
      if(is.finite(score)) rows[[length(rows)+1L]]<-data.frame(feature=s$name,fold=i,baseline_rmse=base,
        candidate_rmse=score,improvement=base-score)
    }
  }
  folds<-if(length(rows)) do.call(rbind,rows) else data.frame()
  if(!nrow(folds)) return(list(summary=data.frame(),folds=folds))
  summary<-do.call(rbind,lapply(split(folds,folds$feature),function(z)
    data.frame(feature=z$feature[1],mean_improvement=mean(z$improvement),
      positive_fold_fraction=mean(z$improvement>0),folds=nrow(z))))
  summary<-summary[order(-summary$mean_improvement),,drop=FALSE]
  list(summary=summary,folds=folds)
}
