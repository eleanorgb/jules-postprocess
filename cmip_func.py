"""
cmip specific functions for postprocessing
dont need to look at if not processing cmip
"""
import sys
import subprocess
from iris.coords import DimCoord
import numpy as np
from read_input import parse_args
from read_input import read_mip_info_no_rose

MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE = parse_args()

if not L_JULES_ROSE:
    MIP_INFO = read_mip_info_no_rose(MIPNAME)

# #############################################################################
def make_global_grid_n96e():
    """
    N96E degrees global grid
    """
    nlat = 144
    latitude = DimCoord(np.linspace(-89.375, 89.375, nlat),
                        standard_name="latitude", units="degrees")
    nlon = 192
    longitude = DimCoord(np.linspace(0.9375, 359.0625, nlon),
                         standard_name="longitude", units="degrees")
    return latitude, longitude, nlat, nlon
# #############################################################################


# #############################################################################
def make_outfilename_cmip(out_dir, outprofile, var, syr, eyr):
    """
    make outfilename as required by mip convert
    """
    if "day" in outprofile:
        out_dir_time = out_dir + "/lnd"
    elif "yr" in outprofile:
        out_dir_time = out_dir + "/lna"
    else:
        out_dir_time = out_dir + "/lnm"
    retcode = subprocess.call("mkdir -p "+ out_dir_time, shell=True)
    if retcode != 0:
        print("mkdir -p "+ out_dir_time)
        sys.exit("mkdir broken")
    try:
        if not L_JULES_ROSE:
            outfilename = out_dir_time+"/"+var+"_"+outprofile+"_"+\
                  MIP_INFO["model"][MIPNAME]+"_"+\
                  MIP_INFO["out_scenario"][MIPNAME]+"_r1i1p1f1_"+\
                  str(syr)+"01-"+str(eyr)+"12.nc"
    except:
        print("broken "+var+": check make_outfilename_cmip")
        raise

    return outfilename
# #############################################################################


# #############################################################################
def define_years_daily_cmip():
    """
    break into 10 year batches
    """
    syrall = np.arange(1850, 2020, 10)
    eyrall = syrall + 9
    eyrall[eyrall == 2019] = 2014
    print(MIPNAME, ": ", syrall)
    print(MIPNAME, ": ", eyrall)
    return syrall, eyrall
# #############################################################################
