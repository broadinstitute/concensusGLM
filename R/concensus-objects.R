#' @title Set up a directory tree for ConcensusGLM
#' @description Set up a directory tree for ConcensusGLM
#' @param path Character. Path to root of directory tree
#' @return list of class "directoryTree"
#' @export
directoryTree <- function ( path='.' ) {

  top.dir <- file.path(path, 'concensusGLM-results')

  csv.dir <- 'tables'
  csv.path <- file.path(top.dir, csv.dir)

  dir.create(csv.path, showWarnings=FALSE, recursive=TRUE)

  return.list <- structure(list(
    top=top.dir,
    tables=csv.path
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
#' @details This creates the concensusDataSet object from \code{data_filename} on which downstream analysis is carried out.
#'
#' The input CSV should be found at \code{data_filename}. This CSV should at least have the headers
#' \code{'id', 'compound', 'concentration', 'strain', 'plate_name', 'count', 'well'}, and must have one row per
#' strain-compound-concentration-plate_name-well combination. Together, "id", "plate_name" and "well" define unique experimental
#' samples; "id" refers to sequencing (technical) replicates (if any) and "plate_name" and "well" are biological replicates
#' (recommended). Any given condition should have at least 2 replicates of some kind. If "row" and "column"
#' are present but not "well", then well is constructed by concatenating "row" and "column".
#'
#' Firstly, this function loads a CSV from \code{data_filename}. This may take some time if it is a large file. It then checks
#' that the minimum headers are present.
#'
#' It then checks for either a \code{negative_control} column or a list of control compounds supplied to \code{controls} argument.
#'
#' Under-represented (assumed to be spurious) strains and plates (as defined by \code{threshold}) are removed.
#'
#' Then, pseudo-strains are built if requested (the default). So far, the only pseudo-strain defined is "total", which
#' is the total counts per well of non-spike-in strains.
#'
#' If defined, annotations (like experimental meta-data) are loaded from \code{annotation_filename}. This CSV file needs a
#' column name in common with \code{data_filename} (case sensitive) since the next step is a join on the shared columns.
#' Also, every observation in the input data that you want to keep must be annotated.
#'
#' Finally, it adds a \code{negative_control} column if necessary and a \code{positive_control} column if defined, and checks
#' that there are at least 2 negative control observations.
#' @param data_filename Character. Path to a table of counts.
#' @param annotation_filename Character. Optional. Path to experimental annotations. These will be joined to the input data
#' and used for batch correction.
#' @param output_path Character. Path to directory where you want the analysis output. This is also where
#' checkpoints and logs from cluster execution are stored. Default is current working directory.
#' @param controls Named list of Characters with elements \code{positive} and \code{negative}.
#' Adds logical columns to data called \code{positive_control} and \code{negative_control} based on regular expression
#' matches.
#' @param test Logical. Run in test mode? If so, only reads the first 5 million lines of \code{data_filename}.
#' @param checkpoint Logical. Save intermediate results as checkpoints?
#' @param threshold Numeric. Strains below this total count threshold will be discarded. Default \code{1000}.
#' Plates below 1000 x \code{threshold} will be discarded. Set to \code{0} to skip this.
#' @param spike_in Character. A regular expression to match spike-in controls.
#' @param pseudostrains Logical. Make pseudostrains such as "total", which is the sum of all non-spike-ins.
#' @return list of class "concensusDataSet"
#' @seealso \link{newWorkflow}, \link{pipeline}, \link{concensusDataSetFromFile}
#' @import workflows
#' @importFrom magrittr %>%
#' @export
concensusDataSetFromFile <- function(data_filename, annotation_filename=NULL, output_path='.',
                                     controls=NULL, rename=NULL, test=FALSE, checkpoint=FALSE,
                                     threshold=100, spike_in='^intcon',
                                     pseudostrains=TRUE,
                                     ...) {

  println('Loading dataset from', data_filename, '...')
  if ( test ) data_ <- readr::read_csv(data_filename, progress=TRUE, n_max=5e6, col_types=readr::cols(concentration=readr::col_double()))
  else        data_ <- readr::read_csv(data_filename, progress=TRUE, col_types=readr::cols(concentration=readr::col_double()))
  check_headers(data_, essential_headers=c('compound', 'concentration', 'strain', 'plate_name', 'count'))

  if ( ! 'well' %in% names(data_) & all(c('row', 'column') %in% names(data_)) ) {

    data_ <- data_ %>% dplyr::mutate(well=paste0(row, column))

  } else if ( ! 'well' %in% names(data_) ) {

    stop('Column "well" or both "row" and "column" must be present in data file.\n')

  }

  original_columns <- names(data_)

  if( ! ('negative_control' %in% original_columns) & is.null(controls) ) {

    stop('Negative control column absent or controls not supplied.\n')

  }

  data_ <- clean(data_, threshold=threshold)

  if ( pseudostrains ) {

    println('Building pseudostrains...')

    pseudostrain_total <- data_ %>%
      dplyr::filter(!grepl(spike_in, strain)) %>%
      dplyr::group_by(id, plate_name, well) %>%
      dplyr::summarize(count=sum(count)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(strain='pseudostrain_total') %>%
      dplyr::left_join(data_ %>% dplyr::ungroup() %>% dplyr::select(-count, -strain) %>% dplyr::distinct()) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()

    stopifnot(length(setdiff(names(pseudostrain_total), names(data_))) == 0)
    stopifnot(nrow(pseudostrain_total) == nrow(data_ %>% dplyr::filter(strain == get_unique_values(data_, 'strain')[1])))

    data_ <- dplyr::bind_rows(data_, pseudostrain_total)

  }

  if ( ! is.null(annotation_filename) ) {

    println('Loading annotations from', annotation_filename, '...')
    annotations <- readr::read_csv(annotation_filename)

  }

  annotation_columns <- names(annotations)

  if ( ! is.null(annotation_filename) ) {

    intersection <- intersect(names(data_), annotation_columns)

    if ( length(intersection) == 0 ) stop('No columns in common between data and annotations.\n')

    data_ <- data_ %>% dplyr::inner_join(annotations, by=intersection)

  }

  if ( ! is.null(controls) ) {

    if ( ! 'negative' %in% names(controls) ) stop('List of controls must include an element named "negative".\n')

    println('Assigning negative controls based on compound matching "', controls$negative, '"')
    data_ <- data_ %>%
      dplyr::mutate(negative_control=grepl(controls$negative, compound))

    if ( 'positive' %in% names(controls) ) {

      if ( ! is.null(controls$positive) ) {

        println('Assigning positive controls based on compound matching "', controls$positive, '"')
        data_ <- data_ %>%
          dplyr::mutate(positive_control=grepl(controls$positive, compound))

      }

    }

    # check at least some negative controls are present
    n_negative_controls <- sum(getElement(data_, 'negative_control'))
    n_positive_controls <- sum(getElement(data_, 'positive_control'))
    println('There are', n_negative_controls, 'negative controls')
    println('There are', n_positive_controls, 'positive controls')
    if ( n_negative_controls  < 2 ) {

      stop('Not enough negative controls. Only ', n_negative_controls, 'present; you need at least 2.\n')
    }

  }

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
