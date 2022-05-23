#!/bin/bash -l

#SBATCH --mem=200GB
#SBATCH --time=4320
#SBATCH --output=outinfo/isimip.%A_%a.out
#SBATCH --error=outinfo/isimip.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-1
#SBATCH --qos=long

##########################################################################
# need to make sure number in array matches number of elements in mod_arr
#           SLURM_ARRAY_TASK_ID=1  # for testing
# (qos normal, time 360) (--time=4320 --qos=long)
##########################################################################

echo This is task $SLURM_ARRAY_TASK_ID

module load scitools

declare -a mod_arr=("ISIMIP3a_GSWP2-W5E5")
# get length of an array
arraylength=${#mod_arr[@]}

# Run this task.
echo This is SLURM task $SLURM_ARRAY_TASK_ID, mipName is ${mod_arr[$SLURM_ARRAY_TASK_ID-1]}
python -u process_jules.py ${mod_arr[$SLURM_ARRAY_TASK_ID-1]} --l_backfill_missing_files

date
