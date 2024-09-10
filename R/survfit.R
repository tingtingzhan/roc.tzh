
#' @title conf_level.survfit
#' 
#' @param object ..
#' 
#' @returns ..
#' 
#' @export
conf_level.survfit <- function(object) object[['conf.int']]








#' @title Convert \link[survival]{survfit.object} to \link[base]{data.frame}
#' 
#' @description ..
#' 
#' @param x \link[survival]{survfit.object}.  
#' `'survfitcox'` object returned from 
#' `survival:::survfit.coxph` is also supported, as of 2022-09-19
#' 
#' @param ... additional parameters of \link[base]{as.data.frame.list}
#' 
#' @param times (optional) \link[base]{numeric} vector
#' 
#' @param rm01 \link[base]{logical} scalar, whether to remove survival rates (and confidence bands) 
#' being 0 or 1, as some y-axis transformation may yield `NA` or `Inf` values.  Default `FALSE`.
#' 
#' @details 
#' A column named `'strata'`
#' of all-equal values `'all_subjects'` 
#' for \link[survival]{survfit.object} without a strata.
#' 
#' @returns 
#' Function [as.data.frame.survfit()] returns a \link[base]{data.frame}
#' 
#' @seealso 
#' `survival:::summary.survfit`
#' 
#' @export as.data.frame.survfit
#' @export
as.data.frame.survfit <- function(
    x, ..., 
    times,
    rm01 = FALSE
) {
  
  fn <- function(xs, ...) { # `xs` is return from ?survival:::summary.survfit
    if (!length(xs$strata)) xs$strata <- 'all_subjects'
    nm <- c('surv', 'time', 'upper', 'lower', 
            'strata', # 'strata' can be NULL
            'n.censor', 'n.event', 'n.risk') # these should be 'integer'
    as.data.frame.list(xs[nm], ...)
  }
  
  # ?survival:::summary.survfit
  # `strata` will be \link[base]{factor}
  # `censored`, ignored if !missing(times); should the censoring times be included in the output?
  # `extend`, ignored if missing(times); even if there are no subjects left at the end of the specified times
  if (!missing(times)) {
    
    ret <- fn(summary(x, times = times, extend = TRUE), ...) 
    ret$txt <- sprintf(
      fmt = '%.0f%% (95%% CI %.0f%%~%.0f%%) at t=%d', 
      1e2*ret$surv, 1e2*ret$lower, 1e2*ret$upper, ret$time)
    
  } else {
    
    d1 <- fn(summary(x, times = 0, extend = TRUE), ...)
    d2 <- fn(summary(x, censored = TRUE), ...)
    ret <- if (length(x$strata)) {
      ds1 <- split.data.frame(d1, f = d1[['strata']])
      ds2 <- split.data.frame(d2, f = d2[['strata']])
      do.call(rbind, args = .mapply(FUN = rbind, dots = list(ds1, ds2), MoreArgs = NULL))
    } else rbind(d1, d2)
    
  }
  
  if (rm01) {
    id <- (ret$time > 0) & (ret$surv < 1) & (ret$surv > 0) & (ret$upper < 1) & (ret$upper > 0) & (ret$lower < 1) & (ret$lower > 0)
    ret <- ret[id, ]
  } # else do_nothing
  
  if ((fom <- x$call$formula)[[1L]] == '~') { # has true formula
    if (is.symbol(fom[[3L]])) {
      if (is.factor(ret[['strata']])) {
        ptn0 <- glob2rx(paste0(deparse1(fom[[3L]]), '='))
        ptn <- substr(ptn0, start = 2L, stop = nchar(ptn0) - 1L)
        attr(ret[['strata']], which = 'levels') <- gsub(pattern = ptn, replacement = '', x = attr(ret[['strata']], which = 'levels', exact = TRUE))
      } else stop('should not happen?')
    }
  }
  
  return(ret)
}




#' @importFrom stats nobs
#' @export
nobs.survfit <- function(object, ...) sum(object[['n']])




#' @title Layers of Kaplan-Meier Curve of \link[survival]{survfit.object} using \CRANpkg{ggplot2}
#' 
#' @description ..
#' 
#' @param object \link[survival]{survfit.object}
#' 
#' @param times (optional) \link[base]{numeric} scalar or \link[base]{vector},
#' time points where survival rates as well as the confidence intervals are plotted
#' 
#' @param ribbon \link[base]{logical} scalar, default `TRUE`
#' 
#' @param labels (optional) \link[base]{character} \link[base]{vector}
#' 
#' @param ... ..
#' 
#' @note
#' It's very difficulty to grab the time \link[base]{units} of the response of a 
#' \link[survival]{survfit} object.
#' The only way is `eval(object$call$formula[[2L]], envir = eval(object$call$data))`,
#' which is very vulnerable as `object$call$data` might be defined in some other internal functions.
#' Therefore we leave the specification of `xlab` (ideally the time unit) to the end user.
#' 
#' 
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 autolayer aes geom_ribbon geom_step scale_fill_discrete geom_point scale_colour_discrete
#' @importFrom scales percent
#' @importFrom survival survdiff
#' @importFrom stats aggregate.data.frame quantile
#' @importFrom utils glob2rx
#' @export
autolayer.survfit <- function(
    object, 
    times,
    ribbon = TRUE,
    labels,
    ...
) {
  
  fom <- object$call$formula
  if (fom[[1L]] != '~') stop('must preserve the true `formula`')
  
  d <- as.data.frame.survfit(object, ...)
  
  old_labs <- if (length(object$strata)) attr(d[['strata']], which = 'levels', exact = TRUE) #else 'all_subjects'

  if (missing(labels) || !length(labels)) {
    ns <- aggregate.data.frame(d[c('n.event', 'n.censor')], by = d['strata'], FUN = sum)
    sm <- quantile(object, probs = .5, conf.int = TRUE) # survival:::quantile.survfit
    #obj_sum <- summary(object)
    #maxFU <- vapply(split.default(obj_sum$time, f = obj_sum$strata), FUN = max, FUN.VALUE = 0)
    # I have difficulty getting the censor time from 'survfit' object
    maxFU <- max(object$time)
    sm_txt <- ifelse(
      test = !is.na(c(sm$quantile)), 
      yes = sprintf(fmt = 't50 = %.1f (%.1f, %.1f)', c(sm$quantile), c(sm$lower), c(sm$upper)),
      no = sprintf(fmt = 't50 > %.1f', maxFU))
    new_labs <- sprintf(fmt = '%d event; %d censor\n%s\n', 
                        ns[['n.event']], ns[['n.censor']], 
                        sm_txt)
    labels <- if (length(old_labs)) paste0(old_labs, '\n', new_labs) else new_labs
    
  } else if (is.character(labels)) {
    if (!length(labels) || anyNA(labels) || !all(nzchar(labels))) stop('illegal labels')
    if (length(labels) != length(old_labs)) stop('user specified `labels` not match the strata of `survfit` object')
    if (!is.null(names(labels))) stop('user specified `labels` must be unnamed!  Watch the order of default saput first, then write user-specified labels')
    cat(sprintf(fmt = '%s was %s\n', sQuote(labels), sQuote(old_labs)), sep = '')
    
  } else if (isFALSE(labels)) {
    labels <- old_labs
    
  } else stop('illegal `labels`')
  
  d_c <- d[d$n.censor > 0L, , drop = FALSE]
  
  strata_nm <- if (is.symbol(fom[[3L]])) deparse1(fom[[3L]]) # else NULL
  
  ret <- list( 
    geom_step(mapping = aes(x = d$time, y = d$surv, group = d$strata, colour = d$strata)), 
    if (ribbon) geom_ribbon(mapping = aes(x = d$time, ymax = d$upper, ymin = d$lower, group = d$strata, fill = d$strata), alpha = .1),
    geom_point(mapping = aes(x = d_c$time, y = d_c$surv, colour = d_c$strata), shape = 3L),
    if (!missing(times)) {
      yr <- as.data.frame.survfit(x = object, times = times)
      geom_label_repel(mapping = aes(x = yr$time, y = yr$surv, fill = yr$strata, label = yr$txt), colour = 'white', size = 3)
    },
    if (length(labels)) scale_colour_discrete(labels = labels) else scale_colour_discrete(guide = 'none'),
    if (length(labels)) scale_fill_discrete(labels = labels) else scale_fill_discrete(guide = 'none'),
    labs(y = deparse1(fom[[2L]]), colour = strata_nm, fill = strata_nm)
  )

  attr(ret, which = 'data') <- d
  return(ret)
  
}





#' @title Kaplan-Meier Curve of \link[survival]{survfit.object} using \CRANpkg{ggplot2}
#' 
#' @description ..
#' 
#' @param object \link[survival]{survfit.object}
#' 
#' @param ... additional parameters of [autolayer.survfit()]
#' 
#' @seealso 
#' `survival:::plot.survfit` `survival:::quantile.survfit`
#' 
#' @importFrom scales percent
#' @importFrom ggplot2 ggplot scale_y_continuous
#' @export
autoplot.survfit <- function(object, ...) {
  ggplot() + 
    autolayer.survfit(object, ...) + 
    scale_y_continuous(labels = percent)
}










#' @title Sprintf.survfit
#' 
#' @description ..
#' 
#' @param model \link[survival]{survfit.object}
#' 
#' @param ... ..
#' 
#' @examples
#' library(survival)
#' Sprintf.survfit(survfit(Surv(time, status) ~ 1, data = aml))
#' Sprintf.survfit(survfit(Surv(time, status) ~ x, data = aml))
#' 
#' @export Sprintf.survfit
#' @export
Sprintf.survfit <- function(model, ...) {
  # read ?survival::survfit.formula carefully
  # how to tell 'single event' or not ?
  fom <- model$call$formula
  edp <- deparse1(fom[[2L]])
  fmt <- 'Kaplan-Meier estimates and curves of endpoint `%s` are obtained using [R] package [survival]'
  if (identical(fom[[3L]], 1)) {
    # no predictor
    sprintf(fmt = paste0(fmt, '.'), edp)
  } else {
    sprintf(fmt = paste0(fmt, ', by predictor(s) %s.'),
            edp,
            paste0('`', all.vars(fom[[3L]]), '`', collapse = ','))
  }
}




#' @title endpoint.survfit
#' 
#' @param x ..
#' 
#' @param ... ..
#' 
#' @returns ..
#' 
#' @export
endpoint.survfit <- function(x, ...) {
  if (is.symbol(fom <- x$call$formula)) return(as.character.default(fom))
  if (!is.call(fom) || (fom[[1L]] != '~')) stop('illegal formula in call of survfit object: ', deparse1(fom))
  return(deparse1(fom[[2L]]))
}

#' @title model_desc.survfit
#' 
#' @param x ..
#' 
#' @param ... ..
#' 
#' @returns ..
#' 
#' @export
model_desc.survfit <- function(x, ...) {
  if (length(x$strata)) {
    model_desc_survdiff_rho(...)
  } else 'Kaplan-Meier'
}







