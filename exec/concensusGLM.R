#!/usr/bin/env Rscript

println <- function(...) cat(date(), '>', ..., '\n')

main <- function(cl_options) {
  
  library(concensusGLM)

  # data input
  file_inputs <- cl_options[c('data', 'meta')]
  controls <- setNames(cl_options[c('pos', 'neg')], c('positive', 'negative'))

  working_dir <- file.path(dirname(file_inputs$data), 'processed')

  if (!file.exists(working_dir)) {

    println('Creating', working_dir)
    dir.create(working_dir)

  }

  dir_tree <- directoryTree(path=working_dir)

  concensus_data <- concensusDataSetFromFile(data_filename=file_inputs$data,
                                             annotation_filename=file_inputs$meta,
                                             output_path=dir_tree$top,
                                             controls=controls,
                                             
                                             pseudostrains = FALSE,

                                             # don't discard plates and strains in in interactive mode
                                             threshold=ifelse(cl_options[['no-count-threshold']], 0, 1000),

                                             checkpoint=cl_options$checkpoint)

  println('Analysis output will be in', concensus_data$pipelines[[1]]$working_directory)

  concensus_data <- analyze(concensus_data,
                            scattering='strain',
                            conditions=c('compound', 'concentration'),

                            # use on local computer or SGE cluster?
                            locality=ifelse(is.null(cl_options$sge), 'local', 'sge'),
                            # parallelize on local computer?
                            parallel=ifelse(is.null(cl_options$sge), cl_options$parallel, FALSE),

                            # parameters for using SGE
                            submit_script=cl_options$sge,
                            # sets memory reservation
                            mem_multiplier=1000,
                            # conservative run time
                            run_time='8:00:00')

  write_concensusDataSet(concensus_data,
                         paste0(concensus_data$pipelines[[1]]$data$output_prefix, '-concensus-data.rds'))

}

VERSION <- '0.5.0'

'Usage:
concensusGLM.R --help
concensusGLM.R --data <data-file> --neg <negative-control> [--meta <experiment-metadata> --pos <positive-control> --no-count-threshold --checkpoint --parallel --sge <template-script>]

Options:
-h --help                      Show this message and exit
--data <data-file>             Input count table CSV, one line per lane per strain per well
--meta <experiment-metadata>   CSV of handling records
--neg <negative-control>       Name of negative control compound, e.g. DMSO or none or untreated
--pos <positive-control>       Name of positive control compound, e.g. rifampin or BRD-K01507359-001-19-5 or BRD-K01507359
--checkpoint                   Save checkpoints to allow restart in case of failure
--parallel                     Analyze strains in parallel (needs multiple cores)
--sge <template-script>        Use GridEngine cluster with template to analyze strains in parallel
--no-count-threshold           Don\'t discard plates or strains based on low counts
' -> doc

if ( interactive() ) {

  # for testing; otherwise edit to the values you want
  COMMANDARGS <- c('--data tests/input/count-data.csv',
                    '--meta tests/input/annotations.csv',
                    '--neg untreated',
                    '--pos BRD-K01507359')

} else {

  COMMANDARGS <- commandArgs(trailingOnly=TRUE)

}

cl_options <- docopt::docopt(doc,  COMMANDARGS, version=VERSION)

println('ConcensusGLM version', VERSION)

main(cl_options)

println('Done!')
