#' Finite causal geometric adstock
#' @param x Numeric ordered exposure series; NA treated as zero.
#' @param alpha Decay in [0,1).
#' @param l_max Maximum lag.
#' @param normalize Normalize by available-history weights.
#' @return Numeric vector.
#' @export
hf_geometric_adstock <- function(x, alpha=.5, l_max=8L, normalize=TRUE) {
  if (!is.numeric(x) || !is.finite(alpha) || alpha < 0 || alpha >= 1 ||
      l_max < 0 || l_max != as.integer(l_max)) stop('Invalid input/parameters')
  x[is.na(x)] <- 0
  w <- alpha^(0:l_max)
  .hf_convolve(x, w, normalize)
}

.hf_convolve <- function(x, weights, normalize=TRUE) {
  n <- length(x); out <- numeric(n)
  if (!n) return(out)
  for (t in seq_len(n)) {
    count <- min(t, length(weights)); w <- weights[seq_len(count)]
    v <- x[seq.int(t, t-count+1L)]
    out[t] <- sum(w*v)/(if (normalize) sum(w) else 1)
  }
  out
}

#' Delayed peaked adstock
#' @param peak Integer peak lag.
#' @inheritParams hf_geometric_adstock
#' @export
hf_delayed_adstock <- function(x, alpha=.5, peak=2L, l_max=12L, normalize=TRUE) {
  if (!is.numeric(x) || alpha <= 0 || alpha >= 1 || peak < 0 ||
      peak != as.integer(peak) || l_max < peak || l_max != as.integer(l_max)) stop('Invalid parameters')
  x[is.na(x)] <- 0
  .hf_convolve(x, alpha^abs(0:l_max-peak), normalize)
}

#' Discrete Weibull survival-difference adstock
#' @param shape Positive Weibull shape.
#' @param scale Positive Weibull scale.
#' @inheritParams hf_geometric_adstock
#' @export
hf_weibull_adstock <- function(x, shape=2, scale=3, l_max=12L, normalize=TRUE) {
  if (!is.numeric(x) || shape <= 0 || scale <= 0 || l_max < 0 || l_max != as.integer(l_max)) stop('Invalid parameters')
  x[is.na(x)] <- 0
  lag <- 0:l_max
  w <- exp(-(lag/scale)^shape)-exp(-((lag+1)/scale)^shape)
  .hf_convolve(x,w,normalize)
}

.hf_nonnegative <- function(x) {
  if (!is.numeric(x) || any(x < 0, na.rm=TRUE)) stop('Nonnegative numeric input required')
}

#' Zero-based logistic saturation
#' @param x Nonnegative numeric vector.
#' @param lam Positive steepness.
#' @export
hf_logistic_saturation <- function(x, lam=1) {
  .hf_nonnegative(x); if (lam <= 0) stop('lam must be positive')
  tanh(lam*x/2)
}

#' Michaelis-Menten saturation
#' @param half_saturation Positive half-saturation.
#' @param maximum Positive maximum.
#' @inheritParams hf_logistic_saturation
#' @export
hf_michaelis_menten <- function(x, half_saturation=1, maximum=1) {
  .hf_nonnegative(x)
  if (half_saturation <= 0 || maximum <= 0) stop('Positive parameters required')
  maximum*x/(half_saturation+x)
}

#' Concave power saturation
#' @param exponent Exponent in (0,1].
#' @inheritParams hf_logistic_saturation
#' @export
hf_power_saturation <- function(x, exponent=.5) {
  .hf_nonnegative(x)
  if (exponent <= 0 || exponent > 1) stop('Invalid exponent')
  x^exponent
}

#' Hyperbolic tangent saturation
#' @param scale Positive scale.
#' @inheritParams hf_logistic_saturation
#' @export
hf_tanh_saturation <- function(x, scale=1) {
  .hf_nonnegative(x); if (scale <= 0) stop('scale must be positive')
  tanh(x/scale)
}

#' Reach times Hill-saturated frequency
#' @param reach Nonnegative reach.
#' @param frequency Nonnegative average frequency.
#' @param half_saturation Positive frequency half-saturation.
#' @param slope Positive Hill slope.
#' @export
hf_reach_frequency <- function(reach, frequency, half_saturation=2, slope=1) {
  .hf_nonnegative(reach); .hf_nonnegative(frequency)
  if (length(reach)!=length(frequency) || half_saturation <= 0 || slope <= 0) stop('Invalid inputs')
  reach*frequency^slope/(half_saturation^slope+frequency^slope)
}
