library(concensusGLM)

file_inputs <- list(data=file.path('tests/input/count-data.csv'),
                      annotation=file.path('tests/input/annotations.csv'))
controls    <- list(positive='BRD-K01507359',
                      negative='untreated')

dir_tree <- directoryTree(path=file.path('tests', 'output', 'sge'))

# rif is BRD-K01507359-001-19-5
concensus_data <- concensusDataSetFromFile(data_filename=file_inputs$data,
                                           annotation_filename=file_inputs$annotation,
                                           output_path=dir_tree$top,
                                           controls=controls,
                                           threshold=0,
                                           checkpoint=TRUE)

# ConcensusGLM model
# split on strain
concensus_data <- scatter(concensus_data, 'strain')

# get rough dispersions
concensus_data <- getRoughDispersions(concensus_data)

# get batch effects
concensus_data <- getBatchEffects(concensus_data)

# get final dispersions
concensus_data <- getFinalDispersions(concensus_data)

# resample untreated
#concensus_data <- resampleNegative(concensus_data)

# final models
concensus_data <- getFinalModel(concensus_data)

concensus_data <- execute(concensus_data, locality='sge',
                          submit_script=file.path('tests', 'input', 'sge-template.sh'))

# gather back together
concensus_data <- gather(concensus_data)

# save
write_concensusDataSet(concensus_data, paste0(concensus_data$pipelines[[1]]$data$output_prefix, '-concensus-data.rds'))

stopifnot('model_parameters' %in% names(concensus_data$pipelines[[1]]$data))
