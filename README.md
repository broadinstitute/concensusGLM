# concensusGLM

R package to facilitate microbial fitness inference (using Generalized Linear Models) from barcode counts based on the Negative Binomial (NB) distribution. 

Designed for ultra-high throughput screening with (more than) tens of millions of observations. This means:

- Estimation of batch effects based on known experimental metadata. This is used in final estimation of effect size, and also in the estimation of the NB dispersion parameter to guard against overly conservative p-values
- Computational short-cuts which allow very large datasets to be analyzed without large `glm` objects causing memory explosions
- Support for parallelization across cores or GridEngine cluster jobs
- Support for checkpointing in case of intermediate failure

## Installation

Not yet on CRAN or BioconductoR, so you need to install using devtools.

In an R session:

`devtools::install_github('eachanjohnson/concensusGLM')`

## Usage

This package provides the objects:

- `concensusDataSet`, a container for ConcensusGLM analysis
- `concenusWorkflow`, a wrapper to the `workflow` object from [workflows](https://github.com/eachanjohnson/workflows), allowing scale-up to multi-core on a laptop or a GridEngine style cluster

The `concensusDataSet` object has methods to carry out ConcensusGLM analysis.

## Citation

[New inhibitors of Mycobacterium tuberculosis identified using systems chemical biology](https://doi.org/10.1101/396440), bioRxiv, Aug 2018

Eachan O Johnson, Emily LaVerriere, Mary Stanley, Emma Office, Elisabeth Meyer, Tomohiko Kawate, James Gomez, Rebecca E Audette, Nirmalya Bandyopadhyay, Natalia Betancourt, Kayla Delano, Israel Da Silva, Joshua Davis, Christina Gallo, Michelle Gardner,  Aaron Golas,  Kristine M Guinn, Rebecca Korn, Jennifer A McConnell, Caitlin E Moss,  Kenan C Murphy, Ray Nietupski,  Kadamba G Papavinasasundaram, Jessica T Pinkham, Paula A Pino,  Megan K Proulx, Nadine Ruecker, Naomi Song, Matthew Thompson, Carolina Trujillo, Shoko Wakabayashi, Joshua B Wallach, Christopher Watson,  Thomas R Ioerger, Eric S Lander, Brian K Hubbard,  Michael H Serrano-Wu,  Sabine Ehrt, Michael Fitzgerald, Eric J Rubin,  Christopher M Sassetti,  Dirk Schnappinger,  Deborah T Hung

## See also

Dependencies:

- R
- [broom](https://github.com/hadley/broom)
- [dplyr](https://github.com/hadley/dplyr)
- [errR](https://github.com/eachanjohnson/errR)
- [future](https://github.com/HenrikBengtsson/future)
- [magrittr](https://github.com/hadley/magrittr)
- MASS
- parallel
- [readr](https://github.com/hadley/readr) for fast loading of large CSV files
- [workflows](https://github.com/eachanjohnson/workflows)

Similar approaches include but are not limited to:

- [DESeq2](https://github.com/mikelove/DESeq2)
- [edgeR]()
- [limma](https://github.com/mikelove/DESeq2)
