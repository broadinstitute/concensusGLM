#$ -cwd
#$ -V
#$ -q broad
#$ -N __name__
#$ -j y
#$ -o __logfile__
#$ -l h_vmem=__virtualmemory__g,h_rt=2:00:00

source /broad/software/scripts/useuse
reuse R-3.2

Rscript tests/input/distributed-executor.R __checkpointfile__ __clobber__ __libraries__
