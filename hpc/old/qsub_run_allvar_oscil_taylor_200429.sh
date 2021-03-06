#!/bin/bash
#------------------------------------------
#$ -N pts_all
#$ -o iman.out$JOB_ID
#$ -j y
#$ -S /bin/bash
#$ -l h_rt=2:00:00 
#$ -l h_vmem=6G
#$ -pe smp 1
#$ -t 1001-2000

date
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load R/3.5.1

cd /home/clarka/pts_r_package/hpc/
./run_allvar_oscil_taylor_200429.R $SGE_TASK_ID
date
