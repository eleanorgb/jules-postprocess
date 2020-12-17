#!/bin/bash -l

#SBATCH --mem=100GB
#SBATCH --time=4320
#SBATCH --output=isimip2b_outinfo/isimip2b.%A_%a.out
#SBATCH --error=isimip2b_outinfo/isimip2b.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-12
#SBATCH --qos=long

##########################################################################
# need to make sure number in array matches number of elements in mod_arr
#           SLURM_ARRAY_TASK_ID=1  # for testing
# (qos normal, time 360) (--time=4320 --qos=long)
##########################################################################

echo This is task $SLURM_ARRAY_TASK_ID

module load scitools/experimental-current

declare -a mod_arr=("ISIMIP2b_GFDL-ESM2M_C20C" "ISIMIP2b_HADGEM2-ES_C20C" 
                   "ISIMIP2b_IPSL-CM5A-LR_C20C" "ISIMIP2b_MIROC5_C20C"
                   "ISIMIP2b_GFDL-ESM2M_RCP2P6" "ISIMIP2b_HADGEM2-ES_RCP2P6" 
                   "ISIMIP2b_IPSL-CM5A-LR_RCP2P6" "ISIMIP2b_MIROC5_RCP2P6"
                   "ISIMIP2b_GFDL-ESM2M_RCP6P0" "ISIMIP2b_HADGEM2-ES_RCP6P0" 
                   "ISIMIP2b_IPSL-CM5A-LR_RCP6P0" "ISIMIP2b_MIROC5_RCP6P0")
# get length of an array
arraylength=${#mod_arr[@]}

# Run this task.
echo This is SLURM task $SLURM_ARRAY_TASK_ID, mipName is ${mod_arr[$SLURM_ARRAY_TASK_ID-1]}
python -u process_jules.py ${mod_arr[$SLURM_ARRAY_TASK_ID-1]}

date
