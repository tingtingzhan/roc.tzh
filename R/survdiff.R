
#' @importFrom stats nobs
#' @export
nobs.survdiff <- function(object, ...) sum(object$n)

#' @importFrom stats pchisq
#' @export .pval.survdiff
#' @export
.pval.survdiff <- function(x) {
  pchisq(x$chisq, df = length(x$n) - 1L, lower.tail = FALSE)
}


model_desc_survdiff_rho <- function(rho, ...) {
  # see ?survival::survdiff
  if (missing(rho) || is.null(rho) || (rho == 0)) {
    # see ?survival::survdiff and ?coin::logrank_test
    # 'Logrank [unweighted]' # 'Logrank [Mantel-Cox]' # 'hypergeometric variance'
    'Logrank'
  } else if (rho == 1) {
    'Peto-Gehan-Wilcoxon'
  } else {
    'G-rho family of tests'
  }
}


#' @title model_desc.survdiff
#' 
#' @param x ..
#' 
#' @param ... ..
#' 
#' @returns ..
#' 
#' @export model_desc.survdiff
#' @export
model_desc.survdiff <- function(x, ...) model_desc_survdiff_rho(x$call$rho)


