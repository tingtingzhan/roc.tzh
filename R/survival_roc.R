


#' @title Survival ROC and Risk Set ROC
#' 
#' @description
#' Convert the returned objects of \link[survivalROC]{survivalROC} and \link[risksetROC]{risksetROC}
#' into `'roc'` objects (as returned by \link[pROC]{roc}), 
#' to take advantages of S3 methods such as \link[pROC]{plot.roc} and [autoplot.roc()].
#' 
#' @param formula ..
#' 
#' @param data ..
#' 
#' @param predict.time,span,... additional parameters of 
#' \link[survivalROC]{survivalROC} and \link[risksetROC]{risksetROC}
#' 
#' @references
#' \url{https://rpubs.com/kaz_yos/survival-auc}
#' 
#' @examples
#' ?survivalROC::survivalROC
#' data(mayo, package = 'survivalROC')
#' library(survival)
#' library(ggplot2)
#' m1 = mayo |> 
#'  within.data.frame(expr = {
#'   time = as.difftime(time, units = 'days')
#'   edp = Surv(time, censor)
#'  }) |>
#'  survival_roc(formula = edp ~ mayoscore4, predict.time = 365)
#' autoplot.roc(m1)
#' \dontrun{autoplot(m1)} # unicode error on devtools::check
#' 
#' library(rmd.tzh); list(
#'   'survival_roc' = m1
#' ) |> render_(file = 'survival_roc')
#' @name survival_roc
#' @importFrom pROC coords
#' @importFrom survivalROC survivalROC
#' @importFrom survival.tzh units.Surv
#' @export
survival_roc <- function(
    formula, data, predict.time, 
    span = .25 * NROW(data)^(-.2), # from examples in ?survivalROC::survivalROC
    ...
) {
  
  yval <- eval(formula[[2L]], envir = data)
  if (!length(unt <- units.Surv(yval))) stop('time-to-event endpoint must carry `unit`')
  
  xval <- eval(formula[[3L]], envir = data)
  
  tmp <- survivalROC(Stime = yval[,1L], status = yval[,2L], marker = xval, predict.time = predict.time, span = span, ...)
  # `tmp` is a 'list'
  # convert `tmp` to class `?pROC::roc`, so that I can use `?pROC:::plot.roc.roc`
  
  ret <- with(tmp, list(
    percent = FALSE,
    sensitivities = TP,
    specificities = 1 - FP,
    thresholds = cut.values,
    auc = structure(AUC, partial.auc = FALSE, class = 'auc')
  ))
  class(ret) <- c('survival_roc', 'roc')
  #attr(ret, which = 'main') <- paste0(deparse1(formula[[2L]]), ' at Time ', predict.time)
  attr(ret, which = 'main') <- paste0(deparse1(formula), ', t = ', predict.time)
  attr(ret, which = 'formula') <- formula
  attr(ret, which = 'data') <- data
  attr(ret, which = 'units') <- unt
  attr(ret, which = 'survivalROC') <- tmp
  attr(ret, which = 'threshold') <- coords(ret, x = 'best', best.method = c('youden'))$threshold # ?pROC:::coords.roc
  return(ret)

}





#' @title md_.survival_roc
#' 
#' @description
#' ..
#' 
#' @param x `'survival_roc'` object
#' 
#' @importFrom rmd.tzh md_
#' @export md_.survival_roc
#' @export
md_.survival_roc <- function(x, xnm, ...) {
  
  fom <- attr(x, which = 'formula', exact = TRUE)
  sroc <- attr(x, which = 'survivalROC', exact = TRUE)
  txt <- sprintf(
    fmt = 'Receiver operating characteristic (ROC) curve of %d-%s `%s` by the predictor `%s` is created using <u>**`R`**</u> package <u>**`survivalROC`**</u>. The Youden\'s index (e.g., the `%s` threshold that maximizes the sensitivity and specificity) is marked. Kaplan-Meier curves stratified by Youden\'s index are provided using <u>**`R`**</u> package <u>**`survival`**</u>.', 
    sroc$predict.time, 
    gsub('s$', replacement = '', attr(x, which = 'units', exact = TRUE)),
    deparse1(fom[[2L]]),
    deparse1(fom[[3L]]), deparse1(fom[[3L]]))
  
  return(c(
    txt,
    '\n',
    '```{r}',
    (attr(x, which = 'fig.height', exact = TRUE) %||% 4) |> sprintf(fmt = '#| fig-height: %.1f'),
    (attr(x, which = 'fig.width', exact = TRUE) %||% 7) |> sprintf(fmt = '#| fig-width: %.1f'),
    sprintf(fmt = '(%s) |> autoplot.survival_roc()', xnm),
    '```'
  ))
  
}








#' @title autoplot.survival_roc
#' 
#' @description
#' ..
#' 
#' @param object `'survival_roc'` object
#' 
#' @param ... ..
#' 
#' @importFrom ggplot2 autoplot ggplot labs
#' @importFrom survival.tzh autolayer.survfit .pval.survdiff
#' @importFrom scales.tzh label_pvalue_sym
#' @export autoplot.survival_roc
#' @export
autoplot.survival_roc <- function(object, ...) {
  attrs <- attributes(object)
  fom <- attrs$formula
  data <- attrs$data
  thres <- attrs$threshold
  x0 <- eval(call('>', fom[[3L]], thres), envir = data)
  x <- structure(x0 + 1L, levels = paste(deparse1(fom[[3L]]), c('>', '<='), sprintf(fmt = '%.2f', thres)), class = 'factor')
  sfit <- eval(call('survfit.formula', formula = eval(call('~', fom[[2L]], quote(x))), data = quote(data)))
  sdiff <- eval(call('survdiff', formula = eval(call('~', fom[[2L]], quote(x))), data = quote(data)))
  ggplot() + autolayer.survfit(sfit) + 
    labs(x = attrs$units, 
         title = sprintf('Threshold determined by %d-%s survival-ROC', attrs$survivalROC$predict.time, gsub('s$', replacement = '', attrs$units)),
         caption = sdiff |>
           .pval.survdiff() |>
           label_pvalue_sym(add_p = TRUE)() |> 
           sprintf(fmt = '%s, Log-rank (unweighted)'),
         colour = NULL, fill = NULL)
}






#' @rdname survival_roc
#' 
#' @param plot must be `FALSE`, to suppress the plot from \link[risksetROC]{risksetROC}
#' 
#' @examples
#' library(ggplot2)
#' data(VA, package = 'MASS')
#' va1 = within(VA, {
#'  cell = factor(cell)
#'  tx = (VA$treat == 1)
#'  status[stime>500] <- 0
#'  stime[stime>500 ] <- 500
#' })
#' library(survival)
#' fit0 = coxph(Surv(stime,status) ~ Karn + cell + tx + age, data = va1, na.action = na.omit)
#' va1$eta = fit0$linear.predictor
#' new = riskset_roc(Surv(stime,status) ~ eta, data = va1, predict.time = 30)
#' autoplot(new)
#' 
#' @importFrom risksetROC risksetROC
#' @export
riskset_roc <- function(
    formula, data, predict.time,
    plot = FALSE,
    ...
) {
  
  yval <- eval(formula[[2L]], envir = data)
  xval <- eval(formula[[3L]], envir = data)

  tmp <- risksetROC(Stime = yval[,1L], status = yval[,2L], marker = xval, predict.time = predict.time, 
                    plot = FALSE, ...)
  # `tmp` is a 'list'
  # convert `tmp` to class `?pROC::roc`, so that I can use `?pROC:::plot.roc.roc`
  
  ret <- with(tmp, list(
    percent = FALSE,
    sensitivities = TP,
    specificities = 1 - FP,
    thresholds = marker,
    auc = structure(AUC, partial.auc = FALSE, class = 'auc')
  ))
  attr(ret, which = 'main') <- paste0('Risk Set of ', deparse1(formula[[2L]]), ' at Time ', predict.time)
  class(ret) <- 'roc'
  return(ret)

}

