

# ?survival:::duplicated.Surv is good (using ?survival:::as.matrix.Surv)
# ?survival:::unique.Surv good



#' @title Time Unit of a \link[survival]{Surv} object
#' 
#' @description ..
#' 
#' @param x \link[survival]{Surv} object
#' 
#' @seealso 
#' \link[base]{units}, `Hmisc:::units.Surv`
#' 
#' @export
units.Surv <- function(x) {
  at <- attributes(x)
  if (length(at$units)) stop('I have retired this usage')
  if (!length(ia <- at$inputAttributes)) return(invisible())
  un <- ia$time$units
  un2 <- ia$time2$units
  if (length(un)) {
    if (length(un2) && (un != un2)) stop('`time` and `time2` in Surv not having same unit')
    return(un)
  }
  if (length(un2)) return(un2)
  return(invisible())
}




#' @title Get Time from \link[survival]{Surv} Object
#' 
#' @description ..
#' 
#' @param x \link[survival]{Surv} objects
#' 
#' @param type \link[base]{character} scalar
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @note 
#' Do not clash with S3 \link[stats]{time}!
#' 
#' @export
getTime.Surv <- function(x, type = c('event', 'censor', 'any'), ...) {
  evt <- as.logical(x[, 'status'])
  id <- switch(match.arg(type), any = TRUE, event = evt, censor = !evt)
  nc <- dim(x)[2L]
  if (nc == 3L) return((x[, 'stop'] - x[, 'start'])[id]) #Not sure if correct!!
  if (nc == 2L) return(x[, 'time'][id])
  stop('Surv object cannot have ncol>3L')
}







#' @title Kaplan-Meier Curve of \link[survival]{Surv} Objects using \CRANpkg{ggplot2}
#' 
#' @description ..
#' 
#' @param object \link[survival]{Surv} object
#' 
#' @param ... ..
#' 
#' @examples 
#' library(survival)
#' autoplot(with(veteran, Surv(time, status)))
#' 
#' @importFrom ggplot2 ggplot autoplot scale_y_continuous
#' @importFrom scales percent
#' @importFrom survival survfit.formula
#' @export
autoplot.Surv <- function(object, ...) {
  nc <- dim(object)[2L]
  cl0 <- match.call()
  if (nc == 2L) {
    sfit <- eval(call('survfit.formula', formula = call('~', cl0$object, 1)))
    p <- ggplot() + 
      autolayer.survfit(sfit, ...) + 
      scale_y_continuous(labels = percent) +
      labs(x = units.Surv(object), y = deparse1(substitute(object)))
    return(p)
  }
  if (nc == 3L) stop('autoplot.Surv not good for col-3 Surv object yet')
  stop('nc > 3L')
}


