#' Gregorian Easter Sunday
#' @param year Integer Gregorian year >= 1583.
#' @return Date.
#' @export
hf_easter <- function(year) {
  if (length(year)!=1 || is.na(year) || year < 1583 || year != as.integer(year)) stop('Invalid year')
  a <- year %% 19; b <- year %/% 100; c <- year %% 100
  d <- b %/% 4; e <- b %% 4; f <- (b+8)%/%25; g <- (b-f+1)%/%3
  h <- (19*a+b-d-g+15)%%30; i <- c%/%4; k <- c%%4
  l <- (32+2*e+2*i-h-k)%%7; m <- (a+11*h+22*l)%/%451
  z <- h+l-7*m+114; month <- z%/%31; day <- z%%31+1
  as.Date(sprintf('%04d-%02d-%02d',year,month,day))
}

#' Known-calendar features (daily dates)
#' @param dates Date vector.
#' @param offsets_weeks Signed week offsets from event date.
#' @param extra_events Named list of Date vectors.
#' @return data.frame, one row per input date.
#' @export
hf_calendar_features <- function(dates, offsets_weeks=c(-10,-4,-2,-1,0,1,2,4,10), extra_events=list()) {
  dates <- as.Date(dates)
  if (anyNA(dates)) stop('Missing dates')
  yr <- as.integer(format(dates,'%Y'))
  start <- as.Date(sprintf('%04d-01-01',yr))
  next_start <- as.Date(sprintf('%04d-01-01',yr+1))
  phase <- as.numeric(dates-start)/as.numeric(next_start-start)
  out <- data.frame(annual_sin_1=sin(2*pi*phase),annual_cos_1=cos(2*pi*phase),
                    annual_sin_2=sin(4*pi*phase),annual_cos_2=cos(4*pi*phase),
                    weekday=as.POSIXlt(dates)$wday, month=as.integer(format(dates,'%m')),
                    day_of_year=as.integer(format(dates,'%j')))
  out$weekday <- (out$weekday+6)%%7
  years <- seq(min(yr)-1L,max(yr)+1L)
  easters <- as.Date(vapply(years, function(y) as.character(hf_easter(y)), character(1)))
  events <- c(list(easter=easters,good_friday=easters-2,easter_monday=easters+1,
                   ascension=easters+39,pentecost=easters+49,
                   christmas=as.Date(sprintf('%04d-12-25',years)),
                   new_year=as.Date(sprintf('%04d-01-01',years))),extra_events)
  if (anyDuplicated(names(events))) stop('Duplicate event names')
  for (nm in names(events)) for (offset in offsets_weeks) {
    label <- sprintf('event_%s_w%+d',nm,as.integer(offset))
    out[[label]] <- as.integer(dates %in% (as.Date(events[[nm]])+7L*as.integer(offset)))
  }
  out
}

#' Exploratory repeated-phase residual scan (TRAIN ONLY)
#' @param dates Date vector.
#' @param residuals Numeric training residuals.
#' @param period_days Integer phase period (364 or 365 typical).
#' @param offsets_days Offsets to inspect.
#' @param min_cycles Minimum distinct years.
#' @param min_hits Minimum observations for a phase.
#' @return Descriptive data.frame, not statistical significance.
#' @export
hf_scan_recurring <- function(dates, residuals, period_days=365L,
                              offsets_days=seq(-70,70,7),min_cycles=3L,min_hits=3L) {
  dates <- as.Date(dates)
  if (length(dates)!=length(residuals) || anyNA(dates) || period_days<=0) stop('Invalid input')
  cols <- c('phase_day','offset_days','n','mean_residual','standardized_effect')
  empty <- as.data.frame(setNames(replicate(5,numeric(),simplify=FALSE),cols))
  if (length(unique(format(dates,'%Y')))<min_cycles) return(empty)
  ok <- is.finite(residuals); dates <- dates[ok]; y <- residuals[ok]
  if (length(y)<min_hits || is.na(sd(y)) || sd(y)<=0) return(empty)
  phase <- as.integer(dates-as.Date('2000-01-01'))%%period_days
  rows <- list()
  for (base in 0:(period_days-1L)) for (offset in offsets_days) {
    hit <- phase == (base+offset)%%period_days; n <- sum(hit)
    if (n >= min_hits) {
      effect <- mean(y[hit])-mean(y)
      rows[[length(rows)+1L]] <- c(base,offset,n,effect,effect/sd(y))
    }
  }
  if (!length(rows)) return(empty)
  out <- as.data.frame(do.call(rbind,rows)); names(out) <- cols
  out[order(-abs(out$standardized_effect)),,drop=FALSE]
}
