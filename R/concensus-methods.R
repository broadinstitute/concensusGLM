#' @title Clean up ConcensusGLM data
#' @description Clean up ConcensusGLM data. Removes strains with total count across data below \code{threshold} and plates
#' with total counts below 1000 x \code{threshold}
#' @param x data.frame to clean.
#' @param threshold Numeric. Strains below this total count threshold will be discarded.
#' Plates below 1000 x \code{threshold} will be discarded.
#' @param ... Other arguments.
#' @return Cleaned data.frame
#' @export
clean <- function(x, ...) UseMethod('clean')

#' @rdname clean
#' @export
clean.default <- function(x, ...) stop('Can\'t clean', class(x), '\n')

#' @rdname clean
#' @importFrom magrittr %>%
#' @export
clean.data.frame <- function(x, threshold=1000, ...) {

  # remove strains with low usage
  unused_strains <- x %>%
    dplyr::group_by(strain, plate_name) %>%
    dplyr::summarise(sum_count=sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(strain) %>%
    dplyr::summarize(n_plates_absent=sum(sum_count <= threshold),
              sum_count=mean(sum_count[sum_count <= threshold])) %>%
    dplyr::arrange(n_plates_absent) %>%
    dplyr::filter(n_plates_absent > 10) %>%
    get_unique_values('strain')

  println('Removing strains with total counts <', threshold + 1, ':', paste(unused_strains, collapse=', '))

  x <- x %>% filter(! strain %in% unused_strains)

  # remove plates with low usage

  unused_plates <- x %>%
    group_by(plate_name) %>%
    summarise(sum_count=sum(count)) %>%
    ungroup() %>%
    filter(sum_count <= 1000 * threshold) %>%
    get_unique_values('plate_name')

  println('Removing plates with total counts <', 1000 * threshold + 1, ':', paste(unused_plates, collapse=', '))

  x <- x %>% filter(! plate_name %in% unused_plates)

}

#' @title Print ConcensusWorkflow
#' @description Print data of a ConcensusWorkflow to the console
#' @param x ConcensusWorkflow.
#' @param ... Other arguments.
#' @export
print.concensusWorkflow <- function(x, ...) {

  lapply(x$pipelines, function(y) print(y$data))

}

#' @title Print ConcensusDataSet
#' @description Print summary stitistics of a ConcensusDataSet to the console
#' @param x ConcensusDataSet.
#' @param ... Other arguments.
#' @return \code{x}, invisibly.
#' @export
print.concensusDataSet <- function(x, ...) {

  println('ConCensus Dataset: created on', x$date_created, 'from filename', x$data_filename, '\n')

  cat('File prefix:', x$output_prefix, '\n')

  cat('Number of strains:', x$n_strains, '\n')
  cat('Number of plates:', x$n_plates, '\n')
  cat('Number of compounds:', x$n_compounds, '\n')
  cat('Number of lanes:', x$n_lanes, '\n')
  cat('Number of reads:', x$n_reads, '\n')
  cat('\n')

  return ( invisible(x) )

}

#' @title Scatter concensusWorkflow
#' @description Scatter concensusWorkflow based on unique values of column \code{by} in the input data.
#' Allows parallelization of downstream analysis.
#' @param x concensusWorkflow.
#' @param by Character. Column of input data containing values idenitfying analytically independent chunks, e.g. strains.
#' @param ... Otehr arguments.
#' @return Scattered concensusWorkflow
#' @seealso \link{scatter}
#' @export
scatter.concensusWorkflow <- function(x, by, ...) {

  println('Scattering by', by, '...')
  println('Scattering by', by, 'will give', length(get_unique_values(x$pipelines[[1]]$data$data, by)), 'chunks')

  class(x) <- 'workflow'

  new_concensus_workflow <- structure(workflows::scatter(x, elements='data', by=by), class=c('concensusWorkflow', 'workflow'))

  return ( new_concensus_workflow )

}

#' @title Gather concensusWorkflow
#' @description Gather concensusWorkflow back into a single pipeline.
#' @param x concensusWorkflow.
#' @param ... Other arguments passed to \code{workflows::execute}.
#' @return Gathered concensusWorkflow
#' @seealso \link{gather}
#' @export
gather.concensusWorkflow <- function(x, ...) {

  println('Gathering', x$pipelines[[1]]$scattered_by, '...')

  class(x) <- 'workflow'

  new_concensus_workflow <- structure(workflows::gather(x, elements=c('data', 'mean_variance_relationship',
                                                                    'dispersion', 'batch_effect_model',
                                                                    'model_parameters'), ...),
                                      class=c('concensusWorkflow', 'workflow'))

  return ( new_concensus_workflow )

}

#' @title Execute concensusWorkflow
#' @description Execute concensusWorkflow analysis pipeline.
#' @param x concensusWorkflow.
#' @param ... Other arguments passed to \code{workflows::execute}.
#' @return Executed concensusWorkflow
#' @seealso \link{execute}
#' @export
execute.concensusWorkflow <- function(x, locality='local', parallel=FALSE, clobber=FALSE, ...) {

  class(x) <- 'workflow'

  new_concensus_workflow <- structure(workflows::execute(x, locality=locality, parallel=parallel, clobber=clobber, ...),
                                      class=c('concensusWorkflow', 'workflow'))

  return ( new_concensus_workflow )

}

#' @title Write ConcensusDataSet to disk
#' @description Serialize ConcensusDataSet to disk under \code{filename}. A wrapper to \code{write_workflow} from the
#' \code{workflows} package.
#' @param x ConcensusDataSet.
#' @param filename Character. Where to save.
#' @param ... Other arguments passed to \code{write_workflow}.
#' @seealso \link{saveRDS}, \link{write_workflow}
#' @export
write_concensusDataSet <- function(x, filename, clean=TRUE, ...)
  workflows::write_workflow(structure(x, class='workflow'), file=filename, clean=TRUE, ...)

