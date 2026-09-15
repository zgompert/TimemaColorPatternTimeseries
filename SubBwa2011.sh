#!/bin/bash 
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --mem=120000
#SBATCH --account=gompert
#SBATCH --qos=notchpeak
#SBATCH --partition=notchpeak
#SBATCH --job-name=bwa
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

ml bwa
## version Version: 0.7.19-r1273

cd /scratch/general/nfs1/u6000989/tcr_fha_timeseries/dat_fha_2011

perl /uufs/chpc.utah.edu/common/home/gompert-group4/projects/timema_color_pattern_complexity/gbs_time_series/data_scripts/bwa_aln_fork.pl *astq
