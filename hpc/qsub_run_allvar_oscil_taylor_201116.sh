#!/bin/bash
#------------------------------------------
#$ -q all.q # Select the serial queue
#$ -l h_vmem=4G # Max. virt. Mem. = 4G per slot/core
#$ -cwd # Change to current working directory
#$ -V # Export environment variables into script
#$ -N myprogram # A name for the job
#$ -o myprogram.out # Name of the SGE-Output File, default: $JOB_NAME.o$JOB_ID
#$ -e myprogram.err # Name of the SGE-Erro File, default: $JOB_NAME.e$JOB_ID
#$ -m e # E-Mail message when Job has




#$ -q all.q
#$ -l h_vmem=4G
#$ -N pts_all
#$ -m e

#$ -l h_rt=4:00:00 


#$ -j y
#$ -pe smp 1
#$ -t 1-10

date
module load R

cd /cl_tmp/clarka/
./run_allvar_oscil_taylor_200429.R $SGE_TASK_ID
date
