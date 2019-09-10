#!/bin/bash
#------------------------------------------
#$ -N pts_example
#$ -o iman.out$JOB_ID
#$ -j y
#$ -S /bin/bash
#$ -l h_rt=6:00:00 
#$ -l h_vmem=6G
#$ -pe smp 1

date
module load R
cd /home/clarka/pts_r_package/pttstability/
../pttstability_man_example.R
date
