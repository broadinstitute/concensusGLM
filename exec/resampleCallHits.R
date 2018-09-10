#!/usr/bin/env Rscript

println <- function(...) cat(date(), '>', ..., '\n')

main <- function(cl_options) {

  library(concensusGLM)

  concensus_data <- readRDS(cl_options$data)

  class(concensus_data) <- c('concensusWorkflow', class(concensus_data))

  println('Analysis output will be in', concensus_data$pipelines[[1]]$working_directory)

  concensus_data_res <- scatter(concensus_data, 'strain')

  concensus_data_res <- resample(concensus_data_res,
                                 sample_size=as.numeric(cl_options$`sample-size`),
                                 n_replicates=as.numeric(cl_options$replicates))

  concensus_data_res <- getFinalModel(concensus_data_res,
                                grouping='strain',
                                conditions=c('compound', 'concentration'))

  concensus_data_res <- execute(concensus_data_res,
                                locality=ifelse(is.null(cl_options$sge), 'local', 'sge'),
                                # parallelize on local computer?
                                parallel=ifelse(is.null(cl_options$sge), cl_options$parallel, FALSE),

                                # parameters for using SGE
                                submit_script=cl_options$sge,
                                # sets memory reservation
                                mem_multiplier=1000,
                                # conservative run time
                                run_time='8:00:00')

  concensus_data_res <- gather(concensus_data_res)

  write_concensusDataSet(concensus_data_res,
                         paste0(concensus_data$pipelines[[1]]$data$output_prefix, '-concensus-data-resampled.rds'))

  concensus_data <- callHits(concensus_data, concensus_data_res,
                             false_discovery_rate=as.numeric(cl_options$fdr))


}

VERSION <- '0.5.0'

'Usage:
resampleCallHits.R --help
resampleCallHits.R --data <data-file> --sample-size <sample-size> --replicates <n-replicates> --fdr <false-discovery-rate> [--checkpoint --parallel --sge <template-script> ]

Options:
-h --help                      Show this message and exit
--data <data-file>             Input count table CSV, one line per lane per strain per well
--sample-size <sample-size>    CSV of handling records
--replicates <n-replicates>    Name of negative control compound, e.g. DMSO or none or untreated
--fdr <false-discovery-rate>   Name of positive control compound, e.g. rifampin or BRD-K01507359-001-19-5 or BRD-K01507359
--checkpoint                   Save checkpoints to allow restart in case of failure
--parallel                     Analyze strains in parallel (needs multiple cores)
--sge <template-script>        Use GridEngine cluster with template to analyze strains in parallel
' -> doc

if ( interactive() ) {

  # for testing; otherwise edit to the values you want
  COMMANDARGS <- c('--data tests/input/count-data-concensus-data.rds',
                   '--sample-size 500',
                   '--replicates 2',
                   '--fdr 0.05')

} else {

  COMMANDARGS <- commandArgs(trailingOnly=TRUE)

}

cl_options <- docopt::docopt(doc,  COMMANDARGS, version=VERSION)

println('ConcensusGLM version', VERSION)

main(cl_options)

println('Done!')
