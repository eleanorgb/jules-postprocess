"""
imogen specific functions for postprocessing
dont need to look at if not processing imogen
"""

import sys
import os
import glob
import iris
import iris.coords as icoords
import numpy as np
from read_input import parse_args
from read_input import config_parse_args
from read_input import read_mip_info_no_rose
import jules_xarray

MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE, JSONMAPFILE = parse_args()

if L_JULES_ROSE:
    CONFIG_ARGS = config_parse_args(MIPNAME)
elif not L_JULES_ROSE:
    MIP_INFO = read_mip_info_no_rose(MIPNAME)

# #############################################################################
imodel_dict_cmip6 = {
    "ACCESS-CM2": 1,
    "ACCESS-ESM1-5": 2,
    "CAS-ESM2-0": 3,
    "CMCC-ESM2": 4,
    "CNRM-CM6-1-HR": 5,
    "CNRM-CM6-1": 6,
    "CNRM-ESM2-1": 7,
    "CanESM5": 8,
    "E3SM-1-0": 9,
    "EC-Earth3-Veg": 10,
    "EC-Earth3": 11,
    "FGOALS-g3": 12,
    "GFDL-ESM4": 13,
    "GISS-E2-1-G": 14,
    "HadGEM3-GC31-LL": 15,
    "HadGEM3-GC31-MM": 16,
    "INM-CM4-8": 17,
    "INM-CM5-0": 18,
    "IPSL-CM6A-LR": 19,
    "MIROC-ES2L": 20,
    "MIROC6": 21,
    "MPI-ESM1-2-HR": 22,
    "MPI-ESM1-2-LR": 23,
    "MRI-ESM2-0": 24,
    "NorESM2-MM": 25,
    "UKESM1-0-LL": 26,
}


def add_ensemble(imodel_dict, drive_model, cube):
    """
    Add our own realization coordinate if it doesn't already exist
    Also add time coordinate year and transfer time to middle of period
    """
    # #####################################################################
    # BEGIN add our own realization coordinate if it doesn't already exist.
    if not cube.coords("realization"):
        realization = imodel_dict[drive_model]
        ensemble_coord = icoords.AuxCoord(
            realization, standard_name="realization", var_name="realization"
        )
        cube.add_aux_coord(ensemble_coord)
    for key in list(cube.attributes.keys()):
        del cube.attributes[key]
    return cube


# #############################################################################
def read_ensemble(files_in, variable_cons, time_cons, diag_in, drive_model):
    """
    Read imogen ensembles defined by imodel_dict
    """

    # #########################################################################
    def model_ensemble_callback(cube, field, filename):
        """
        Add our own realization coordinate if it doesn't already exist
        Also add time coordinate year and transfer time to middle of period
        """
        # #####################################################################
        # BEGIN add our own realization coordinate if it doesn't already exist.
        if not cube.coords("realization"):
            realization = imodel_dict[filename.split("/")[-2]]
            ensemble_coord = icoords.AuxCoord(
                realization, standard_name="realization", var_name="realization"
            )
            cube.add_aux_coord(ensemble_coord)
        for key in list(cube.attributes.keys()):
            del cube.attributes[key]
        # from add_time_middle function in process jules
        cube.coord("time").bounds = 3600 * (
            np.rint(cube.coord("time").bounds.astype("int64") / 3600.0)
        )
        cube.coord("time").points = (
            cube.coord("time").bounds[:, 0] + cube.coord("time").bounds[:, 1]
        ) / 2.0
        # END add our own realization coordinate if it doesn't already exist.
        # #####################################################################

    keyall = []
    cubeall = iris.cube.CubeList([])
    try:
        if "imogen" in MIP_INFO["model"][MIPNAME].lower():
            imodel_dict = imodel_dict_cmip6
    except:
        if "imogen" in CONFIG_ARGS["MODEL_INFO"]["mipname"].lower():
            imodel_dict = imodel_dict_cmip6
    keys = imodel_dict.keys()

    for key in keys:
        print("INFO: ensembles - reading files made using patterns from:", key)
        files_read = list()
        if len(drive_model) == 0:
            files_tmp = [f.replace("*", key) if "*" in f else f for f in files_in]
        else:
            files_tmp = [f.replace(drive_model, key) for f in files_in]
        files_tmp = [glob.glob(f) for f in files_tmp]
        files_read = [f for f in files_tmp if f]
        if len(files_read) > 0:
            if isinstance(files_read[0], list):
                files_read = [f for sublist in files_read for f in sublist]
            files_read = [f for f in files_read if key in f]

        if len(files_read) > 0:
            cube = jules_xarray.load(files_read, variable_cons & time_cons)
            coord_names = [coord.name() for coord in cube.coords()]
            if "scpool" in coord_names:
                cube = cube.collapsed("scpool", iris.analysis.SUM)
            if diag_in in [
                "imogen_radf",
                "dtemp_g",
                "c_emiss_out",
                "d_ocean_atmos",
                "d_land_atmos_co2",
            ]:
                cube = cube.collapsed(["latitude", "longitude"], iris.analysis.MEAN)
                cube.remove_coord("latitude")
                cube.remove_coord("longitude")

            cube = add_ensemble(imodel_dict, key, cube)
            keyall.append(key)
            cubeall.append(cube)
    sel_dict = {
        "model: " + key: imodel_dict[key] for key in keyall if key in imodel_dict
    }
    cubeall = cubeall.merge_cube()
    cubeall.attributes = sel_dict
    coord_names = [coord.name() for coord in cubeall.coords()]
    if "latitude" in coord_names:
        cubeall.coord("latitude").guess_bounds()
        cubeall.coord("latitude").long_name = "latitude"
    if "longitude" in coord_names:
        cubeall.coord("longitude").guess_bounds()
        cubeall.coord("longitude").long_name = "longitude"
    cubeall.coord("time").points = cube.coord("time").points.astype(np.int32)

    return cubeall


# #############################################################################


# #############################################################################
def make_outfilename_imogen(out_dir, outprofile, var, syr, eyr, diag_dic):
    """
    sort out filename for outputfile
    """
    errorcode = 0
    outfilename = ""

    if outprofile not in ["monthly"]:
        key = f"{var}_{outprofile}"
        outprofile = outprofile + "_"
    else:
        key = var

    # remove timestep of data from output filename
    outprofile = ""

    if not key in diag_dic.keys():
        print(f"ERROR: {key} not in diag_dic")
        errorcode = 1
        return outfilename, errorcode

    if not L_JULES_ROSE:
        if diag_dic[key]["cmip_varname"] is not None:
            outfilename = (
                out_dir
                + "/"
                + MIP_INFO["model"][MIPNAME]
                + "_"
                + MIP_INFO["out_scenario"][MIPNAME]
                + "_"
                + diag_dic[key]["cmip_varname"]
                + "_"
                + outprofile
                + str(syr)
                + "_"
                + str(eyr)
                + ".nc"
            )
    else:
        if diag_dic[key]["cmip_varname"] is not None:
            outfilename = (
                out_dir
                + "/"
                + CONFIG_ARGS["OUT_INFO"]["model_out_id"]
                + "_"
                + CONFIG_ARGS["MODEL_INFO"]["configname"]
                + "_"
                + CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]
                + "_"
                + diag_dic[key]["cmip_varname"]
                + "_"
                + outprofile
                + str(syr)
                + "_"
                + str(eyr)
                + ".nc"
            )
        else:
            print(f"ERROR: cmip_varname for {var}, {outprofile} is none in diag_dic")
            errorcode = 1

    return outfilename, errorcode


# #############################################################################


# #############################################################################
def make_infilename_imogen(src_dir, jules_profname, years):
    """
    sort out filename for input file
    """

    files_in = []
    errorcode = 0

    if not L_JULES_ROSE:
        if isinstance(jules_profname, list):
            for profname in jules_profname:
                files_in.extend(
                    [
                        [
                            f"{src_dir}/*_{MIP_INFO['in_scenario'][MIPNAME]}{MIP_INFO['run_name'][MIPNAME]}{profname}.{year}.nc"
                            for year in years
                        ]
                    ]
                )
        else:
            files_in = [
                src_dir
                + "/*_"
                + MIP_INFO["in_scenario"][MIPNAME]
                + MIP_INFO["run_name"][MIPNAME]
                + jules_profname
                + "."
                + year
                + ".nc"
                for year in years
            ]
    else:
        if isinstance(jules_profname, list):
            for profname in jules_profname:
                files_in.extend(
                    [
                        src_dir
                        + "/*/*"
                        + CONFIG_ARGS["MODEL_INFO"]["driving_model"]
                        + "_"
                        + CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]
                        + "_"
                        + CONFIG_ARGS["MODEL_INFO"]["configname"]
                        + "."
                        + profname
                        + "."
                        + year
                        + ".nc"
                        for year in years
                    ]
                )
        else:
            files_in = [
                src_dir
                + "/*/*"
                + CONFIG_ARGS["MODEL_INFO"]["driving_model"]
                + "_"
                + CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]
                + "_"
                + CONFIG_ARGS["MODEL_INFO"]["configname"]
                + "."
                + jules_profname
                + "."
                + year
                + ".nc"
                for year in years
            ]

    existing_files = [glob.glob(file) for file in files_in]
    if len(existing_files) > 0:
        if isinstance(existing_files[0], list):
            existing_files = [file for sublist in existing_files for file in sublist]

    if len(existing_files) == 0:
        print("ERROR: No input files found of form:")
        if not L_JULES_ROSE:
            print(
                [
                    f"{src_dir}*/*_{MIP_INFO['in_scenario'][MIPNAME]}{MIP_INFO['run_name'][MIPNAME]}{profname}.*.nc"
                    for profname in jules_profname
                ]
            )
        else:
            print(
                f"{src_dir}/{CONFIG_ARGS['MODEL_INFO']['driving_model']}/{CONFIG_ARGS['MODEL_INFO']['driving_model']}_{CONFIG_ARGS['MODEL_INFO']['climate_scenario']}_{CONFIG_ARGS['MODEL_INFO']['configname']}.{jules_profname}.*.nc"
            )
        errorcode = 1

    return existing_files, errorcode
