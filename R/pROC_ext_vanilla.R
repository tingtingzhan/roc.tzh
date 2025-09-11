

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
#' ?infert
#' # using the wrong model just to illustrate the creation of ROC plot
#' 
#' # will change with each run, because of sampling
#' # use ?base::set.seed if needed
#' m0 = get_roc(case ~ spontaneous + induced, data = infert, pctR = .6)
#' \dontrun{autoplot(m0)}
#' 
#' # will not change with each run
#' m = get_roc(case ~ spontaneous + induced, data = infert, pctR = FALSE)
#' \dontrun{autoplot(m)}
#' 
#' @importFrom stats glm family binomial predict.glm
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
    # `idxR`: training index
    # `idxS`: test index
    idxR <- idxS <- TRUE
  } else {
    if (!is.numeric(pctR) || anyNA(pctR) || length(pctR) != 1L) stop('illegal training percentage `pctR`')
    if (pctR <= 0 || pctR >= 1) stop('illegal training percentage `pctR`')
    idxR <- sort.default(c(
      sample(which(edp1), size = floor(sum(edp1) * pctR)),
      sample(which(edp0), size = floor(sum(edp0) * pctR))
    )) # sort is not needed, just for developer's eye
    idxS <- -idxR
  }
  
  dataR <- data[idxR, ] # training data
  dataS <- data[idxS, ] # test data
  
  modR <- glm(formula, data = dataR, family = binomial(link = 'logit')) # training model
  
  trms <- modR[['terms']] |>
    attr(which = 'dataClasses', exact = TRUE)
  if (length(trms) == 2L && trms[2L] == 'numeric') { # single numeric predictor
    X <- eval(formula[[3L]], envir = data)
    rocS <- roc(edp1[idxS] ~ X, quiet = TRUE, percent = TRUE) # test ROC
  } else { # must use linear predictor
    probS <- predict.glm(modR, newdata = dataS, type = 'link') # predicted test linear-predictor, using training model
    rocS <- roc(edp1[idxS] ~ probS, quiet = TRUE, percent = TRUE) # test ROC
  }
  
  attr(rocS, which = 'main') <- if (isFALSE(pctR)) {
    sprintf('%s (%d subjects)', deparse1(formula), length(edp1))
  } else {
    sprintf('ROC curve based on %d:%d training:test subjects', length(idxR), length(edp1) - length(idxR))
  }
  
  return(invisible(rocS))
}



#' @title autoplot.roc
#' 
#' @description
#' ..
#' 
#' 
#' @param object roc
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
