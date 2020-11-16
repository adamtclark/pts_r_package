#!/bin/bash
#------------------------------------------
#$ -q all.q
#$ -l h_vmem=4G
#$ -N pts_all
#$ -m e
#$ -cwd
#$ -l h_rt=4:00:00 
#$ -pe smp 1
#$ -t 1-10

date
module load R

cd /cl_tmp/clarka/pts_r_package/hpc
./run_allvar_oscil_taylor_200429.R $SGE_TASK_ID
date
