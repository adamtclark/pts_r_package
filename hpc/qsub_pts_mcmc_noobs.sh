#!/bin/bash
#------------------------------------------
#$ -N pts_all
#$ -o iman.out$JOB_ID
#$ -j y
#$ -S /bin/bash
#$ -l h_rt=1:00:00 
#$ -l h_vmem=6G
#$ -pe smp 1
#$ -t 1-1000

date
module load R
cd /home/clarka/pts_r_package/hpc/
./run_ptstab_mcmc_grad_noobs.R $SGE_TASK_ID
date