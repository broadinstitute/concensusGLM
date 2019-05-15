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
clean.data.frame <- function(x, threshold=100, plate_fraction=0.6, ...) {

  n_plates <- length(get_unique_values(x, 'plate_name'))

  # remove strains with low usage
  unused_strains0 <- x %>%
    dplyr::group_by(strain, plate_name) %>%
    dplyr::summarise(sum_count=sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(strain) %>%
    dplyr::summarize(n_plates_absent=sum(sum_count <= threshold),
                     mean_sum_count=mean(sum_count[sum_count <= threshold])) %>%
    dplyr::arrange(n_plates_absent) %>%
    dplyr::filter(n_plates_absent > n_plates * plate_fraction)

  unused_strains <- unused_strains0 %>% get_unique_values('strain')

  if ( length(unused_strains) > 0 ) {

    println('Removing', length(unused_strains), 'strains with total counts <',
            threshold + 1, 'in more than', plate_fraction, 'of the plates (',
            n_plates, ') :')
    print(knitr::kable(unused_strains0))

    x <- x %>% filter(! strain %in% unused_strains)

  }

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

  elements <- intersect(c('data', 'model_parameters', 'resampled', 'mean_variance_relationship'),
                        names(x$pipelines[[1]]$data))

  println("Chunking elements", pyjoin(elements, ', '), '...')

  new_concensus_workflow <- structure(workflows::scatter(x, elements=elements, by=by), class=c('concensusWorkflow', 'workflow'))

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
                                                                    'model_parameters',
                                                                    'resampled', 'selected_cutoff', 'cutoffs',
                                                                    'resampled_roc_data',
                                                                    'resampled_roc_summary'), ...),
                                      class=c('concensusWorkflow', 'workflow'))

  return ( new_concensus_workflow )

}

#' @title Execute concensusWorkflow
#' @description Execute concensusWorkflow analysis pipeline.
#' @param x concensusWorkflow.
#' @param locality Character. \code{"local"} for local machine (the default) or \code{"sge"} for GridEngine.
#' @param parallel Logical. Run in multicore parallel.
#' @param clobber Logical. If \code{TRUE}, start from beginning of workflow; if \code{FALSE}, start from last built target.
#' @param submit_script Character. If \code{locality="sge"}, path to template SGE script
#' @param ... Other arguments passed to \code{workflows::execute}.
#' @return Executed concensusWorkflow
#' @seealso \link{execute}
#' @export
execute.concensusWorkflow <- function(x, locality='local', parallel=FALSE, clobber=FALSE,
                                      submit_script=file.path('exec/sge-template.sh'), ...) {

  class(x) <- 'workflow'

  new_concensus_workflow <- structure(workflows::execute(x, locality=locality, parallel=parallel, clobber=clobber,
                                                         submit_script=submit_script, ...),
                                      class=c('concensusWorkflow', 'workflow'))

  return ( new_concensus_workflow )

}

#' @title Write ConcensusDataSet to disk
#' @description Serialize ConcensusDataSet to disk under \code{filename}. A wrapper to \code{write_workflow} from the
#' \code{workflows} package.
#' @param x ConcensusDataSet.
#' @param filename Character. Where to save.
#' @param input_csv Character. Write.
#' @param output_csv Character. Where to save.
#' @param output_matrix Logical. Save a matrix of strains and conditions. Make sure concentrations (if used for
#' final fit) are common across all compounds or this will make a mess.
#' @param ... Other arguments passed to \code{write_workflow}.
#' @seealso \link{saveRDS}, \link{write_workflow}
#' @export
write_concensusDataSet <- function(x, filename,
                                   clean=TRUE,
                                   input_csv=TRUE, output_csv=TRUE, output_matrix=FALSE,
                                   ...){

  workflows::write_workflow(structure(x, class='workflow'), file=filename, clean=TRUE, ...)

  filename_stub <- strsplit(filename, '.', fixed=TRUE)[[1]]
  filename_stub <- paste(filename_stub[1:(length(filename_stub) - 1)], collapse='.')
  write_plaintext(x, filename=filename_stub, input_csv=input_csv, output_csv=output_csv, output_matrix=output_matrix)


}

#' @title Write ConcensusGLM data as plain text
#' @description Write ConcensusGLM data as plain text
#' @param x concensusWorkflow
#' @param input_csv Character. Write.
#' @param output_csv Character. Where to save.
#' @param output_matrix Logical. Save a matrix of strains and conditions. Make sure concentrations (if used for
#' final fit) are common across all compounds or this will make a mess.
#' @param ... Other arguments.
#' @export
write_plaintext <- function(x, ...) UseMethod('write_plaintext')

#' @rdname write_plaintext
#' @export
write_plaintext.default <- function(x, ...) stop('Can\'t write_plaintext', class(x), '\n')

#' @rdname write_plaintext
#' @importFrom magrittr %>%
#' @export
write_plaintext.concensusWorkflow <- function(x, filename,
                                             input_csv=TRUE, output_csv=TRUE,
                                             output_matrix=TRUE,
                                             ...) {

  if ( input_csv ) {

    input_data <- x$pipelines[[1]]$data$data

    this_filename <- paste(filename, 'cleaned-input', 'csv', sep='.')

    println('Writing cleaned input data as', this_filename)
    readr::write_csv(dplyr::filter(input_data, !grepl('^pseudo_', strain)),
                     this_filename)

    if ( length(grep('^pseudo', get_unique_values(input_data, 'strain'))) ) {

      this_filename <- paste(filename, 'cleaned-input-pseudostrains', 'csv', sep='.')

      println('Writing pseudostrain input data as', this_filename)
      readr::write_csv(dplyr::filter(input_data, grepl('^pseudo_', strain)),
                       this_filename)

    }

  }


  if ( output_csv ) {

    if ( ! 'model_parameters' %in% names(x$pipelines[[1]]$data) )
      stop('Analysis not complete; can\'t write inference CSV' )

    output_data <- x$pipelines[[1]]$data$model_parameters

    this_filename <- paste(filename, 'output', 'csv', sep='.')

    println('Writing inference output data as', this_filename)
    readr::write_csv(dplyr::filter(output_data, !grepl('^pseudo_', strain)),
                     this_filename)

    if ( length(grep('^pseudostrain_', get_unique_values(output_data, 'strain'))) > 0 ) {

      this_filename <- paste(filename, 'output-pseudostrains', 'csv', sep='.')

      println('Writing pseudostrain output data as', this_filename)
      readr::write_csv(dplyr::filter(output_data, grepl('^pseudostrain_', strain)),
                       this_filename)

    }

  }

  if ( output_matrix ) {

    if ( ! 'model_parameters' %in% names(x$pipelines[[1]]$data) )
      stop('Analysis not complete; can\'t write inference matrix')

    output_data <- reshape2::dcast(x$pipelines[[1]]$data$model_parameters,
                                   strain ~ condition_group,  # smooths out conc inconsistencies
                                   value.var='l2fc',
                                   fun.aggregate=mean)  # not ideal but prevents erroring

    this_filename <- paste(filename, 'output-matrix', 'csv', sep='.')

    println('Writing inference output matrix as', this_filename)
    readr::write_csv(dplyr::filter(output_data, !grepl('^pseudostrain_', strain)),
                     this_filename)

    if ( length(grep('^pseudostrain_', get_unique_values(output_data, 'strain'))) > 0 ) {

      this_filename <- paste(filename, 'output-matrix-pseudostrains', 'csv', sep='.')

      println('Writing pseudostrain output data as', this_filename)
      readr::write_csv(dplyr::filter(output_data, grepl('^pseudostrain_', strain)),
                       this_filename)

    }

  }

}

#' @title Invoke a Data Viewer for concensus objects
#' @description Invoke a spreadsheet-style data viewer on a Concensus R object.
#' @param x concensusWorkflow
#' @param scatter_n Numeric. If \code{x} is a \code{concensusWorkflow} which has been scattered, use this chunk.
#' @param element Character. Which element of \code{concensusDataset} to look at. Defaults to the input data, \code{"data"}.
#' @param limit_size Numeric. Limit number of rows passed to \code{View}. For very large data, this can stop
#' your RStudio session from crashing! Default \code{1000}.
#' @param ... Other arguments.
#' @seealso \link{View}
#' @export
View_c <- function(x, ...) UseMethod('View_c')

#' @rdname write_plaintext
#' @export
View_c.default <- function(x, ...) View(x, ...)

#' @rdname write_plaintext
#' @importFrom magrittr %>%
#' @export
View_c.concensusWorkflow <- function(x, scatter_n=1, element='data', limit_size=1000, ...) {

  View_c(x$pipelines[[scatter_n]]$data, element=element, limit_size=limit_size, ...)

}

#' @rdname write_plaintext
#' @importFrom magrittr %>%
#' @export
View_c.concensusDataSet <- function(x, element='data', limit_size=1000, ...) {

  View(head(getElement(x, element), n=limit_size), title=element, ...)

}



