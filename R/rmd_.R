

#' @title rmd_.survival_roc
#' 
#' @description
#' ..
#' 
#' @param x [survival_roc]
#' 
#' @param xnm ..
#' 
#' @param ... ..
#' 
#' @export rmd_.survival_roc
#' @export
rmd_.survival_roc <- function(x, xnm, ...) {
  h <- attr(x, which = 'fig.height', exact = TRUE) %||% 4
  w <- attr(x, which = 'fig.width', exact = TRUE) %||% 7
  r_figsz <- sprintf(fmt = '```{r results = \'asis\', fig.height = %.1f, fig.width = %.1f}', h, w)
  return(c(
    Sprintf.survival_roc(x),
    r_figsz, 
    sprintf(fmt = 'autoplot.survival_roc(%s)', xnm), 
    '```'#,
    #r_figsz, 
    #sprintf(fmt = 'ggKM(%s)', .obj), # do we (mathematically) have this method?
    # '```'
  ))
}



