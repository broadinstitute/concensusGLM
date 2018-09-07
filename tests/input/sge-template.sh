#$ -cwd
#$ -V
#$ -q broad
#$ -N __name__
#$ -j y
#$ -o __logfile__
#$ -l h_vmem=__virtualmemory__m,h_rt=__runtime__

source /broad/software/scripts/useuse
reuse R-3.2

Rscript tests/input/distributed-executor.R __checkpointfile__ __clobber__ __libraries__
