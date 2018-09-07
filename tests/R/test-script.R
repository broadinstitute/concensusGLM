library(concensusGLM)

if ( interactive() ) {

  file_inputs <- list(data=file.path('tests/input/count-data.csv'),
                      annotation=file.path('tests/input/annotations.csv'))
  controls    <- list(positive='BRD-K01507359',
                      negative='untreated')

} else {

  file_inputs <- get_file_inputs(commandArgs())
  controls <- get_controls(commandArgs())

}


dir_tree <- directoryTree(path=file.path('tests', 'output'))

# rif is BRD-K01507359-001-19-5
concensus_data <- concensusDataSetFromFile(data_filename=file_inputs$data,
                                           annotation_filename=file_inputs$annotation,
                                           output_path=dir_tree$top,
                                           controls=controls,
                                           threshold=0,
                                           checkpoint=TRUE)

# ConcensusGLM model
# split on strain
concensus_data1 <- scatter(concensus_data, 'strain')

# get rough dispersions
concensus_data1 <- getRoughDispersions(concensus_data1)

# get batch effects
concensus_data1 <- getBatchEffects(concensus_data1)

# get final dispersions
concensus_data1 <- getFinalDispersions(concensus_data1)

# resample untreated
#concensus_data <- resampleNegative(concensus_data)

# final models
concensus_data1 <- getFinalModel(concensus_data1)

concensus_data1 <- execute(concensus_data1, locality='local', parallel=FALSE)

# gather back together
concensus_data1 <- gather(concensus_data1)

stopifnot('model_parameters' %in% names(concensus_data1$pipelines[[1]]$data))

# save
write_concensusDataSet(concensus_data1, paste0(concensus_data1$pipelines[[1]]$data$output_prefix, '-concensus-data.rds'))
write_concensusDataSet(concensus_data1, paste0(concensus_data1$pipelines[[1]]$data$output_prefix, '-concensus-data.rds'),
                       output_matrix=TRUE)

# test analyze
concensus_data2 <- analyze(concensus_data1)

stopifnot('model_parameters' %in% names(concensus_data2$pipelines[[1]]$data))

