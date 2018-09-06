#' @title Set up a directory tree for ConcensusGLM
#' @description Set up a directory tree for ConcensusGLM
#' @param path Character. Path to root of directory tree
#' @return list of class "directoryTree"
#' @export
directoryTree <- function ( path='.' ) {

  top.dir <- file.path(path, 'concensusGLM-results')
  plot.dir <- 'plots'

  plot.path <- file.path(top.dir, plot.dir)
  plot.subdirs <- c(
    'dispersions', 'heatmaps', 'histograms', 'resampling',
    'ma-plots', 'rocs', 'dose-response-curves', 'nuisance-factors', 'pca',
    'plate-heatmaps', 'chemistry', 'hits'
  )

  csv.dir <- 'tables'
  csv.path <- file.path(top.dir, csv.dir)

  dir.tree <- lapply(plot.subdirs, function (x) {
    if ( ! file.exists(file.path(plot.path, x)) ) {
      println('Creating', file.path(plot.path, x))
      dir.create(file.path(plot.path, x), showWarnings=FALSE, recursive=TRUE)
    }
  })
  dir.create(csv.path, showWarnings=FALSE, recursive=TRUE)

  return.list <- structure(list(
    top=top.dir,
    tables=csv.path,
    #checkpoints=file.path(top.dir, 'checkpoints'),
    plots=plot.path,
    #hist=file.path(plot.path, 'histograms'),
    dispersions=file.path(plot.path, 'dispersions'),
    #heatmaps=file.path(plot.path, 'heatmaps'),
    plate.heatmaps=file.path(plot.path, 'plate-heatmaps'),
    ma.plots=file.path(plot.path, 'ma-plots'),
    #rocs=file.path(plot.path, 'rocs'),
    pca=file.path(plot.path, 'pca'),
    dr.curves=file.path(plot.path, 'dose-response-curves'),
    #nuisance=file.path(plot.path, 'nuisance-factors'),
    #resampling=file.path(plot.path, 'resampling'),
    #chemistry=file.path(plot.path, 'chemistry'),
    hits=file.path(plot.path, 'hits')
  ), class='directoryTree')

  return ( return.list )

}

#' @title Set up analysis container for ConcensusGLM
#' @description Set up analysis container for ConcensusGLM. This is the workhorse function; most users should use
#' \link{concensusDataSetFromFile}.
#' @param checkpoint Logical. Save intermediate results as checkpoints?
#' @param working_directory Character. Path to where you want the analysis output.
#' @param load_checkpoint If \code{"latest"}, reload from latest checkpoint files in \code{working_directory}.
#' If a string pattern matching existing checkpoint, then loads those files. If \code{NULL} (the default) then doesn nothing.
#' @return list of class "concensusDataSet"
#' @seealso \link{workflows::newWorkflow}, \link{concensusDataSetFromFile}
#' @export
newConcensusDataSet <- function(checkpoint=FALSE, working_directory='.', load_checkpoint=NULL, ...) {

  workflow <- workflows::newWorkflow(structure(list(...), class='concensusDataSet'),
                                     checkpoint=checkpoint, working_directory=working_directory,
                                     load_checkpoint=load_checkpoint)
  class(workflow) <- c('concensusWorkflow', class(workflow))

  #print(str(workflow))

  return ( workflow )

}

#' @title Set up analysis container for ConcensusGLM
#' @description Set up analysis container for ConcensusGLM using path to data as input. Performs validity checks on the input.
#' Optionally annotates input with experimental metadata from \code{annotation_filename}.
#' @param data_filename Character. Path to a table of counts.
#' @param annotation_filename Character. Optional. Path to experimental annotations. These will be joined to the input data
#' and used for batch correction.
#' @param output_path Character. Path to directory where you want the analysis output. Default is current working directory.
#' @param controls Named list with elements \code{positive} and \code{negative}.
#' Adds logical columns to data called \code{positive_control} and \code{negative_control}.
#' @param test Logical. Run in test mode?
#' @param checkpoint Logical. Save intermediate results as checkpoints?
#' @param threshold Numeric. Strains below this total count threshold will be discarded.
#' Plates below 1000 x \code{threshold} will be discarded.
#' @return list of class "concensusDataSet"
#' @seealso \link{workflows::newWorkflow}, \link{concensusDataSetFromFile}
#' @import workflows
#' @importFrom magrittr %>%
#' @export
concensusDataSetFromFile <- function(data_filename, annotation_filename=NULL, output_path='.',
                                     controls=NULL, rename=NULL, test=FALSE, checkpoint=FALSE,
                                     threshold=1000,
                                     ...) {

  library('readr')
  library('dplyr')

  println('Loading dataset from', data_filename, '...')
  if ( test ) data_ <- readr::read_csv(data_filename, progress=TRUE, n_max=5e6)
  else        data_ <- readr::read_csv(data_filename, progress=TRUE)

  check_headers(data_)

  original_columns <- names(data_)

  stopifnot(all(c('positive_control', 'negative_control') %in% original_columns) | !is.null(controls))

  data_ <- clean(data_, threshold=threshold)

  if ( ! is.null(annotation_filename) ) {

    println('Loading annotations from', annotation_filename, '...')
    annotations <- read_csv(annotation_filename)

  }

  annotation_columns <- names(annotations)

  if ( ! is.null(annotation_filename) ) data_ <- data_ %>% inner_join(annotations)

  if ( ! is.null(controls) ) data_ <- data_ %>%
    dplyr::mutate(positive_control=grepl(controls$positive, compound),
                  negative_control=grepl(controls$negative, compound))

  output_prefix <- file.path(output_path, strsplit(basename(data_filename), '.csv', fixed=TRUE)[[1]][1])

  concensus_dataset <- newConcensusDataSet(data_filename=data_filename,
                                           annotation_filename=annotation_filename,
                                           date_created=date(),
                                           data=data_,
                                           data_columns=original_columns,
                                           annotation_columns=annotation_columns,
                                           controls=controls,
                                           n_datapoints=nrow(data_),
                                           n_compounds=length(get_unique_values(data_, 'compound')),
                                           n_strains=length(get_unique_values(data_, 'strain')),
                                           n_plates=length(get_unique_values(data_, 'plate_name')),
                                           n_lanes=length(get_unique_values(data_, 'id')),
                                           n_reads=sum(as.numeric(column2vector(data_, 'count'))),
                                           output_prefix=output_prefix,
                                           checkpoint=checkpoint,
                                           working_directory=output_path)

  return ( concensus_dataset )

}
