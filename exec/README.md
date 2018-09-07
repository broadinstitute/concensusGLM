# concensusGLM scripts

This directory contains scripts for default analysis.

## Fitness inference

The script `concensusGLM.R` will check input files, 

```
Usage:
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
```

