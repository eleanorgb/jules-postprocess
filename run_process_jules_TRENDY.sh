#!/bin/bash -l

#SBATCH --mem=100GB
#SBATCH --time=4320
#SBATCH --output=outinfo/trendy.%A_%a.out
#SBATCH --error=outinfo/trendy.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-1
#SBATCH --qos=long

##########################################################################
# need to make sure number in --array matches number of elements in mod_arr
# (qos normal, time 360)
##########################################################################

echo This is task $SLURM_ARRAY_TASK_ID

module load scitools/experimental-current

declare -a mod_arr=("TRENDY_JULES-ES_S2")
# get length of an array
arraylength=${#mod_arr[@]}

# Run this task.
echo This is SLURM task $SLURM_ARRAY_TASK_ID, mipName is ${mod_arr[$SLURM_ARRAY_TASK_ID-1]}
python -u process_jules.py  ${mod_arr[$SLURM_ARRAY_TASK_ID-1]} --l_backfill_missing_files

date
