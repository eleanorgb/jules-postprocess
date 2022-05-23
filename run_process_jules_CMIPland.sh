#!/bin/bash -l

#SBATCH --mem=100GB
#SBATCH --time=360
#SBATCH --output=outinfo/cmipland.%A_%a.out
#SBATCH --error=outinfo/cmipland.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-8
#SBATCH --qos=normal

##########################################################################
# need to make sure number in array matches number of elements in mod_arr
#           SLURM_ARRAY_TASK_ID=1  # for testing
# (qos normal, time 360), --time=4320 --qos=long
##########################################################################

echo This is task $SLURM_ARRAY_TASK_ID

module load scitools

#"CMIPland_UKESM_land-hist-altStartYear" 

declare -a mod_arr=("CMIPland_UKESM_land-hist" 
"CMIPland_UKESM_land-noLu" "CMIPland_UKESM_land-cClim" 
"CMIPland_UKESM_land-cCO2" "CMIPland_UKESM_land-hist-cruNcep" 
"CMIPland_HadGEM_land-hist" "CMIPland_HadGEM_land-noLu"
"CMIPland_HadGEM_land-hist-cruNcep")

# get length of an array
arraylength=${#mod_arr[@]}

# Run this task.
echo This is SLURM task $SLURM_ARRAY_TASK_ID, mipName is ${mod_arr[$SLURM_ARRAY_TASK_ID-1]}
python -u process_jules.py ${mod_arr[$SLURM_ARRAY_TASK_ID-1]} --l_backfill_missing_files

date
