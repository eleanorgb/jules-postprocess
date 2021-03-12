#!/bin/bash -l

#SBATCH --mem=100GB
#SBATCH --time=360
#SBATCH --output=outinfo/imogen5.%A_%a.out
#SBATCH --error=outinfo/imogen5.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-21
#SBATCH --qos=normal

##########################################################################
# need to make sure number in array matches number of elements in mod_arr
#           SLURM_ARRAY_TASK_ID=1  # for testing
# (qos normal, time 360) (--time=4320 --qos=long)
##########################################################################

echo This is task $SLURM_ARRAY_TASK_ID

module load scitools/experimental-current

declare -a mod_arr=("IMOGEN5_aj257_26" "IMOGEN5_aj257_45" "IMOGEN5_aj257_85"
"IMOGEN5_aj257_26_ninorg_cemissions_excludeoldc" "IMOGEN5_aj257_45_ninorg_cemissions_excludeoldc" "IMOGEN5_aj257_85_ninorg_cemissions_excludeoldc"
"IMOGEN5_aj257_26_ninorg_cemissions" "IMOGEN5_aj257_45_ninorg_cemissions" "IMOGEN5_aj257_85_ninorg_cemissions"
"IMOGEN5_aj257_26_excludeoldc" "IMOGEN5_aj257_45_excludeoldc" "IMOGEN5_aj257_85_excludeoldc"
"IMOGEN5_aj257_26_isunfrozen" "IMOGEN5_aj257_45_isunfrozen" "IMOGEN5_aj257_85_isunfrozen"
"IMOGEN5_ah181_26_excludeoldc" "IMOGEN5_ah181_45_excludeoldc" "IMOGEN5_ah181_85_excludeoldc"
"IMOGEN5_ah181_26" "IMOGEN5_ah181_45" "IMOGEN5_ah181_85")
# get length of an array
arraylength=${#mod_arr[@]}

# Run this task.
echo This is SLURM task $SLURM_ARRAY_TASK_ID, mipName is ${mod_arr[$SLURM_ARRAY_TASK_ID-1]}
python -u process_jules.py ${mod_arr[$SLURM_ARRAY_TASK_ID-1]} --l_backfill_missing_files

date
