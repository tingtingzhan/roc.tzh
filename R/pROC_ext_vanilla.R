

#' @title Helper Files of \link[pROC]{roc} Plot
#' 
#' @description ..
#' 
#' @param formula \link[stats]{formula} of the univariable or multivariable logistic regression model
#' 
#' @param data \link[base]{data.frame}
#' 
#' @param pctR \link[base]{numeric} scalar between 0 and 1, 
#' percentage of training data from input `data` 
#' (stratified sampling from response and non-response).
#' Or use `FALSE` to indicate all subjects are used in both training data and test data. 
#' Default `.6`.
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @references 
#' Inspired by \url{https://daviddalpiaz.github.io/r4sl/logistic-regression.html}
#' 
#' @examples 
#' # ?datasets::infert
#' # using the wrong model just to illustrate the creation of ROC plot
#' 
#' # will change with each run, because of sampling
#' set.seed(12); list(
#'   'with-seed' = get_roc(case ~ spontaneous + induced, data = infert, pctR = .6)
#' ) |> rmd.tzh::render_(file = 'roc1')
#' 
#' # will not change with each run
#' list(
#'   'without-seed' = get_roc(case ~ spontaneous + induced, data = infert, pctR = FALSE)
#' ) |> rmd.tzh::render_(file = 'roc2')
#' 
#' @keywords internal
#' @importFrom pROC roc
#' @export
get_roc <- function(
    formula, 
    data, 
    pctR = .6, 
    ...
) {
  
  edp <- formula[[2L]]
  edp1 <- as.logical(eval(edp, envir = data))
  if (anyNA(edp1)) stop('do not allow missingness in the binary endpoint')
  edp0 <- !edp1
  
  if (isFALSE(pctR)) {
    # all are train and all are test
    # `id0`: training index
    # `id1`: test index
    id0 <- id1 <- TRUE
  } else {
    if (!is.numeric(pctR) || anyNA(pctR) || length(pctR) != 1L) stop('illegal training percentage `pctR`')
    if (pctR <= 0 || pctR >= 1) stop('illegal training percentage `pctR`')
    id0 <- sort.default(c(
      sample(which(edp1), size = floor(sum(edp1) * pctR)),
      sample(which(edp0), size = floor(sum(edp0) * pctR))
    )) # sort is not needed, just for developer's eye
    id1 <- -id0
  }
  
  m0 <- glm(formula, data = data[id0,], family = binomial(link = 'logit')) # training model
  
  trms <- m0[['terms']] |>
    attr(which = 'dataClasses', exact = TRUE)
  if (length(trms) == 2L && trms[2L] == 'numeric') { # single numeric predictor
    X <- eval(formula[[3L]], envir = data)
    rocS <- roc(edp1[id1] ~ X, quiet = TRUE, percent = TRUE) # test ROC
  } else { # must use linear predictor
    probS <- predict.glm(m0, newdata = data[id1,], type = 'link') # predicted test linear-predictor, using training model
    rocS <- roc(edp1[id1] ~ probS, quiet = TRUE, percent = TRUE) # test ROC
  }
  
  attr(rocS, which = 'main') <- if (isFALSE(pctR)) {
    fom <- formula
    edp <- fom[[2L]]
    edp_v <- all.vars(edp)
    if (length(edp_v) == 1L) {
      fom[[2L]] <- as.symbol(edp_v)
    }
    sprintf(fmt = '%s (%d subjects)', deparse1(fom), length(edp1))
  } else {
    sprintf(fmt = 'ROC curve based on %d:%d training:test subjects', length(id0), length(edp1) - length(id0))
  }
  
  return(invisible(rocS))
}



#' @title autoplot.roc
#' 
#' @description
#' ..
#' 
#' 
#' @param object \link[pROC]{roc} object
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' Function `pROC:::plot.roc.roc` does \emph{not} return a \link[ggplot2]{ggplot} figure.
#' 
#' @importFrom pROC plot.roc
#' @importFrom ggplot2 autoplot
#' @export autoplot.roc
#' @export
autoplot.roc <- function(
    object, 
    ...
) {
  plot.roc(
    object, 
    print.thres = 'best', 
    print.auc = TRUE, print.auc.x = .2, print.auc.y = 0,
    main = attr(object, which = 'main', exact = TRUE), 
    ...)
}






#' @title md_.roc
#' 
#' @description
#' ..
#' 
#' @param x \link[pROC]{roc} object
#' 
#' @param xnm ..
#' 
#' @param ... ..
#' 
#' @keywords internal
#' @importFrom methods new
#' @importFrom rmd.tzh md_
#' @importClassesFrom rmd.tzh md_lines
#' @export md_.roc
#' @export
md_.roc <- function(x, xnm, ...) {
  
  fom <- attr(x, which = 'formula', exact = TRUE)
  sroc <- attr(x, which = 'survivalROC', exact = TRUE)
  
  #z1 <- sprintf(
  #  fmt = 'Receiver operating characteristic (ROC) curve of %d-%s `%s` by the predictor `%s` is created using <u>**`R`**</u> package <u>**`survivalROC`**</u>. The Youden\'s index (e.g., the `%s` threshold that maximizes the sensitivity and specificity) is marked. Kaplan-Meier curves stratified by Youden\'s index are provided using <u>**`R`**</u> package <u>**`survival`**</u>.', 
  #  sroc$predict.time, 
  #  gsub('s$', replacement = '', attr(x, which = 'units', exact = TRUE)),
  #  deparse1(fom[[2L]]),
  #  deparse1(fom[[3L]]), deparse1(fom[[3L]])
  #) |> # fix in future!!
  #  new(Class = 'md_lines', package = 'pROC')
  
  z2 <- c(
    '```{r}',
    '#| echo: false',
    (attr(x, which = 'fig.height', exact = TRUE) %||% 4) |> sprintf(fmt = '#| fig-height: %.1f'),
    (attr(x, which = 'fig.width', exact = TRUE) %||% 7) |> sprintf(fmt = '#| fig-width: %.1f'),
    sprintf(fmt = '(%s) |> autoplot.roc()', xnm),
    '```'
  ) |> 
    new(Class = 'md_lines', package = 'pROC')
  
  #return(c(z1, z2)) # ?rmd.tzh::c.md_lines
  return(z2)
  
}




