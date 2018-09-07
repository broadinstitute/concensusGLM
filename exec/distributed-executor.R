
library(workflows)

checkpoint_file <- commandArgs(trailingOnly=TRUE)[1]
clobber         <- as.logical(commandArgs(trailingOnly=TRUE)[2])
libraries        <- commandArgs(trailingOnly=TRUE)[3:length(commandArgs(trailingOnly=TRUE))]

for ( lib in libraries ) library(lib, character.only=TRUE)

this_pipeline   <- readRDS(checkpoint_file)
this_pipeline   <- execute(this_pipeline, clobber=clobber)

checkpoint(this_pipeline)

this_success_filename <- file.path(this_pipeline$working_directory, '_sge-logs',
                                   paste0(basename(this_pipeline$checkpoint_filename), '.success'))

writeLines(this_success_filename, this_success_filename)
