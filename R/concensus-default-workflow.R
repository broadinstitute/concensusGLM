#' @title Run concensusGLM analysis
#' @description Do all caluclations on count data to calculate batch effects, dispersions, and final effect sizes and p-values.
#' @param x concensusWorkflow or concensusDataSet.
#' @param grouping Character vector. Combination of data columns to identify unique observations.
#' @param conditions Character vector. Combination of data columns to idenitfy unique conditions.
#' @param scattering Character. Column identifying analytically independent chunks for parallelization.
#' @param ... Other arguments passed to \link{execute}.
#' @return concensusWorkflow or concensusDataSet with a new \code{mean_variance_relationship} and \code{dispersion} attribute.
#' @export
analyze <- function(x, ...) UseMethod('analyze')

#' @rdname analyze
#' @export
analyze.default <- function(x, ...) stop('Can\'t analyze ', class(x))

#' @rdname analyze
#' @export
analyze.concensusDataSet <- function(x,
                                     grouping=c('compound', 'concentration', 'strain'),
                                     conditions=c('compound', 'concentration'),
                                     scattering='strain',
                                     ...) {

  # get rough dispersions
  x <- getRoughDispersions(x, grouping=grouping)

  # get batch effects
  x <- getBatchEffects(x, grouping=grouping)

  # get final dispersions
  x <- getFinalDispersions(x, grouping=grouping)

  # resample untreated
  #x <- resampleNegative(x)

  # final models
  x <- getFinalModel(x, conditions=conditions, grouping=scattering)

  return (x)

}

#' @rdname analyze
#' @export
analyze.concensusWorkflow <- function(x,
                                      grouping=c('compound', 'concentration', 'strain'),
                                      conditions=c('compound', 'concentration'),
                                      scattering='strain',
                                      locality='local', parallel=FALSE,
                                      ...) {

  # split on strain
  x <- scatter(x, scattering)

  x <- workflows::delay(x, analyze, grouping=grouping, conditions=conditions, scattering=scattering)

  x <- execute(x, locality=locality, parallel=parallel, ...)

  # gather back together
  x <- gather(x)

  return (x)

}
