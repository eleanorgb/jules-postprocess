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
from process_jules import USE_JULES_PY

sys.path.append("/home/h03/hadea/bin")
import jules

MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE = parse_args()

if L_JULES_ROSE:
    CONFIG_ARGS = config_parse_args(MIPNAME)
elif not L_JULES_ROSE:
    MIP_INFO = read_mip_info_no_rose(MIPNAME)


imodel_dict_cmip5 = {
    "cen_bcc_mod_bcc-csm1-1": 1,
    "cen_ipsl_mod_ipsl-cm5a-lr": 2,
    "cen_mri_mod_mri-cgcm3": 3,
    "cen_noaa-gfdl_mod_gfdl-esm2m": 4,
    "cen_bcc_mod_bcc-csm1-1-m": 5,
    "cen_ipsl_mod_ipsl-cm5a-mr": 6,
    "cen_nasa-giss_mod_giss-e2-h": 7,
    "cen_nsf-doe-ncar_mod_cesm1-bgc": 8,
    "cen_bnu_mod_bnu-esm": 9,
    "cen_ipsl_mod_ipsl-cm5b-lr": 10,
    "cen_nasa-giss_mod_giss-e2-h-cc": 11,
    "cen_nsf-doe-ncar_mod_cesm1-cam5": 12,
    "cen_cccma_mod_canesm2": 13,
    "cen_miroc_mod_miroc-esm": 14,
    "cen_nasa-giss_mod_giss-e2-r": 15,
    "cen_nsf-doe-ncar_mod_cesm1-waccm": 16,
    "cen_cmcc_mod_cmcc-cms": 17,
    "cen_miroc_mod_miroc-esm-chem": 18,
    "cen_nasa-giss_mod_giss-e2-r-cc": 19,
    "cen_cnrm-cerfacs_mod_cnrm-cm5": 20,
    "cen_miroc_mod_miroc5": 21,
    "cen_ncar_mod_ccsm4": 22,
    "cen_csiro-bom_mod_access1-0": 23,
    "cen_mohc_mod_hadgem2-cc": 24,
    "cen_ncc_mod_noresm1-m": 25,
    "cen_csiro-bom_mod_access1-3": 26,
    "cen_mohc_mod_hadgem2-es": 27,
    "cen_ncc_mod_noresm1-me": 28,
    "cen_csiro-qccce_mod_csiro-mk3-6-0": 29,
    "cen_mpi-m_mod_mpi-esm-lr": 30,
    "cen_noaa-gfdl_mod_gfdl-cm3": 31,
    "cen_inm_mod_inmcm4": 32,
    "cen_mpi-m_mod_mpi-esm-mr": 33,
    "cen_noaa-gfdl_mod_gfdl-esm2g": 34,
}

# imodel_dict_cmip5 = {'cen_cnrm-cerfacs_mod_cnrm-cm5':1,
#              'cen_csiro-bom_mod_access1-3':2,
#              'cen_ipsl_mod_ipsl-cm5a-lr':3}

imodel_dict_cmip6 = {
    "ACCESS-CM2": 1,
    "ACCESS-ESM1-5": 2,
    "CNRM-CM6-1-HR": 3,
    "CNRM-CM6-1": 4,
    "CNRM-ESM2-1": 5,
    "CanESM5": 6,
    "EC-Earth3-Veg": 7,
    "FGOALS-g3": 8,
    "HadGEM3-GC31-LL": 9,
    "HadGEM3-GC31-MM": 10,
    "INM-CM4-8": 11,
    "IPSL-CM6A-LR": 12,
    "MIROC-ES2L": 13,
    "MIROC6": 14,
    "MPI-ESM1-2-HR": 15,
    "MPI-ESM1-2-LR": 16,
    "MRI-ESM2-0": 17,
    "UKESM1-0-LL": 18,
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
    imodel_dict = imodel_dict_cmip5
    try:
        if "imogen6" in MIP_INFO["model"][MIPNAME].lower():
            imodel_dict = imodel_dict_cmip6
    except:
        if "imogen6" in CONFIG_ARGS["MODEL_INFO"]["mipname"].lower():
            imodel_dict = imodel_dict_cmip6
    keys = imodel_dict.keys()
    for key in keys:
        files_read = list()
        if drive_model is None:
            files_tmp = [f.replace("*", key) for f in files_in]
        else:
            files_tmp = [f.replace(drive_model, key) for f in files_in]
        files_tmp = [glob.glob(f) for f in files_tmp]
        files_read = [f for f in files_tmp if f]
        files_read = [f for sublist in files_read for f in sublist]
        if len(files_read) > 0:
            if USE_JULES_PY:
                cube = jules.load(
                    files_read,
                    variable_cons & time_cons,
                    missingdata=np.ma.masked,
                    callback=model_ensemble_callback,
                )
                for ijk, icube in enumerate(cube):
                    if ijk > 0:
                        cube[ijk].coord("time").convert_units(cube[0].coord("time").units)
                cube = cube.concatenate_cube()
            else:
                cube = jules_xarray.load(files_read,variable_cons & time_cons)
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
    sel_dict = {"model: " + key: imodel_dict[key] for key in keyall if key in imodel_dict}
    cubeall = cubeall.merge_cube()
    cubeall.attributes = sel_dict
    coord_names = [coord.name() for coord in cubeall.coords()]
    if "latitude" in coord_names:
        cubeall.coord("latitude").guess_bounds()
        cubeall.coord("latitude").long_name = "latitude"
    if "longitude" in coord_names:
        cubeall.coord("longitude").guess_bounds()
        cubeall.coord("longitude").long_name = "longitude"

    return cubeall


# #############################################################################


# #############################################################################
def make_outfilename_imogen(out_dir, outprofile, var, syr, eyr, diag_dic):
    """
    sort out filename for outputfile
    """
    errorcode = 0
    outfilename = ""

    if outprofile != "monthly":
        key = f"{var}_{outprofile}"
    else:
        key = var

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
                + "_"
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
                + "_"
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
                        + "/"
                        + CONFIG_ARGS["MODEL_INFO"]["driving_model"]
                        + "/"
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
                + "/"
                + CONFIG_ARGS["MODEL_INFO"]["driving_model"]
                + "/"
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

    existing_files = [file for file in files_in if os.path.isfile(file)]

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
