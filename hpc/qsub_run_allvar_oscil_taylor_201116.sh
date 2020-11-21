#!/bin/bash
#------------------------------------------
#$ -S /bin/bash
#$ -V
#$ -j yes
#$ -q all.q
#$ -l h_vmem=4G
#$ -N pts_all
#$ -m e
#$ -M adam.clark@uni-graz.at
#$ -cwd
#$ -l h_rt=14:00:00 
#$ -pe smp 1
#$ -t 747-778

date
module load R

cd /cl_tmp/clarka/pts_r_package/hpc
./run_allvar_oscil_taylor_201116.R $SGE_TASK_ID
date
