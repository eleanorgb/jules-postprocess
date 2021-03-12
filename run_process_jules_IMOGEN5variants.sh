#!/bin/bash -l

#SBATCH --mem=100GB
#SBATCH --time=4320
#SBATCH --output=outinfo/imogen5.%A_%a.out
#SBATCH --error=outinfo/imogen5.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-2
#SBATCH --qos=long

##########################################################################
# need to make sure number in array matches number of elements in mod_arr
#           SLURM_ARRAY_TASK_ID=1  # for testing
# (qos normal, time 360) (--time=4320 --qos=long)
##########################################################################

echo This is task $SLURM_ARRAY_TASK_ID

module load scitools/experimental-current

#declare -a mod_arr=("IMOGEN5variant_C_co2only_85")
declare -a mod_arr=("IMOGEN5variant_CN_PIndep_26" "IMOGEN5variant_CN_PIndep_45"
                    "IMOGEN5variant_CN_PIndep_85" "IMOGEN5variant_C_co2only_26"
                    "IMOGEN5variant_C_co2only_45" "IMOGEN5variant_C_co2only_85"
                    "IMOGEN5variant_CN_co2only_26" "IMOGEN5variant_CN_co2only_45"
                    "IMOGEN5variant_CN_co2only_85" "IMOGEN5variant_C_climonly_26"
                    "IMOGEN5variant_C_climonly_45" "IMOGEN5variant_C_climonly_85"
                    "IMOGEN5variant_CN_climonly_26" "IMOGEN5variant_CN_climonly_45"
                    "IMOGEN5variant_CN_climonly_85" "IMOGEN5variant_CN_26"
                    "IMOGEN5variant_CN_45" "IMOGEN5variant_CN_85"
                    "IMOGEN5variant_C_26" "IMOGEN5variant_C_45"
                    "IMOGEN5variant_C_85")
declare -a mod_arr=("IMOGEN5variant_CN_85" "IMOGEN5variant_C_85")
# get length of an array
arraylength=${#mod_arr[@]}

# Run this task.
echo This is SLURM task $SLURM_ARRAY_TASK_ID, mipName is ${mod_arr[$SLURM_ARRAY_TASK_ID-1]}
python -u process_jules.py ${mod_arr[$SLURM_ARRAY_TASK_ID-1]} --l_backfill_missing_files

date
