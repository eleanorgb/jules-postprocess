"""
cmip specific functions for postprocessing
dont need to look at if not processing cmip
"""

import sys
import subprocess
from iris.coords import DimCoord
import numpy as np
from read_input import parse_args
from read_input import config_parse_args
from read_input import read_mip_info_no_rose

MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE = parse_args()

if L_JULES_ROSE:
    CONFIG_ARGS = config_parse_args(MIPNAME)
elif not L_JULES_ROSE:
    MIP_INFO = read_mip_info_no_rose(MIPNAME)


# #############################################################################
# #############################################################################
def make_global_grid_n96e():
    """
    N96E degrees global grid
    """
    nlat = 144
    latitude = DimCoord(
        np.linspace(-89.375, 89.375, nlat), standard_name="latitude", units="degrees"
    )
    nlon = 192
    longitude = DimCoord(
        np.linspace(0.9375, 359.0625, nlon), standard_name="longitude", units="degrees"
    )
    return latitude, longitude, nlat, nlon


# #############################################################################
# #############################################################################
def make_outfilename_cmip(out_dir, outprofile, var, syr, eyr):
    """
    make outfilename as required by mip convert
    """
    if outprofile is not None:
        if "day" in outprofile:
            out_dir_time = out_dir + "/lnd"
        elif "yr" in outprofile:
            out_dir_time = out_dir + "/lna"
        else:
            out_dir_time = out_dir + "/lnm"
        retcode = subprocess.call("mkdir -p " + out_dir_time, shell=True)
        if retcode != 0:
            print("mkdir -p " + out_dir_time)
            sys.exit("mkdir broken")

        if not L_JULES_ROSE:
            outfilename = (
                f"{out_dir_time}/"
                f"{var}_{outprofile}_"
                f"{MIP_INFO['model'][MIPNAME]}_"
                f"{MIP_INFO['out_scenario'][MIPNAME]}_r1i1p1f1_"
                f"{syr}01-{eyr}12.nc"
            )
        else:
            outfilename = (
                f"{out_dir_time}/"
                f"{var}_{outprofile}_"
                f"{CONFIG_ARGS['OUT_INFO']['model_out_id']}-{CONFIG_ARGS['MODEL_INFO']['configname']}_"
                f"{CONFIG_ARGS['MODEL_INFO']['climate_scenario']}_r1i1p1f1_"
                f"{syr}01-{eyr}12.nc"
            )
    else:
        print(f"ERROR: no outprofile defined for {var}")
        print(f"ERROR: no outfilename defined for {var}")
        outfilename = ""
    return outfilename


# #############################################################################
# #############################################################################
def define_years_daily_cmip(syrin, eyrin):
    """
    break into 10 year batches
    """
    syrall = np.arange(syrin[0], eyrin[0], 10)
    eyrall = syrall + 9
    eyrall[eyrall > eyrin] = eyrin
    print(MIPNAME, ": ", syrall)
    print(MIPNAME, ": ", eyrall)
    return syrall, eyrall
