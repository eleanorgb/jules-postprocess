"""
script to postprocess jules output 
suitable for ISIMIP, TRENDY, ILAMB, IMOGEN, CMIP
BE CAREFUL with nitrogen on and off and npp and nlim
sometimes needs lots of memory - run on spice
"""

USE_JULES_PY = False  # false - trying a quicker way of reading data
READ_JSON = False  # true - trying to read from json file

import time
import os
import sys
import glob
import warnings
import getpass
import subprocess
import numpy.ma as ma
import numpy as np
import iris
import iris.coord_categorisation
from iris.time import PartialDateTime
from iris.cube import Cube
from cf_units import Unit
from cftime import num2date
from landcover_types import select_vegfrac
from landcover_types import add_tile_info
from sel_diags_test import sel_diags_test
from read_input import parse_args
from read_input import config_parse_args
from read_input import read_mip_info_no_rose
from sort_varout_outprofile_name import sort_varout_outprofile_name
import isimip_func  # isimip runs only
import imogen_func  # imogen runs only
import cmip_func  # cmip runs only
import json

if not USE_JULES_PY:
    import jules_xarray

from cubelist_functions import (
    burntarea_func,
    sum_func,
    nee_func,
    div_func,
    rhums_func,
    nbp_func,
    annmax_func,
    minus_func,
    fracweight_func,
    rflow_func,
    mult_func,
    top10cm_func,
    sth_func,
    layered_soilbgc_func,
    conv_360days_to_sec,
    neefire_func,
    select_toplayer_soilmois,
)

warnings.filterwarnings("ignore")

# #############################################################################
# #########################################################################
# read command line arguments and finsh imports
global MIPNAME, L_TESTING, L_JULES_ROSE
MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE = parse_args()
# #########################################################################

if "imogen" in MIPNAME.lower():
    READ_JSON = True
if "isimip" in MIPNAME.lower():
    READ_JSON = True
if "cmip" in MIPNAME.lower():
    READ_JSON = True

if not L_JULES_ROSE:
    if not READ_JSON:
        import ISIMIP2_variables
        import ISIMIP3_variables
        import CMIP_variables
        import TRENDY_variables

        # import IMOGEN_variables
        # import IMOGENvariant_variables
        # import IMOGEN6_variables

if L_JULES_ROSE:
    if not READ_JSON:
        import suite_postprocessed_variables

        diag_dic = suite_postprocessed_variables.get_var_dict()

if not L_JULES_ROSE:
    CONTACT_EMAIL = "eleanor.burke@metoffice.gov.uk"  # dont use my email!
    OUT_BASEDIR = "/scratch/" + getpass.getuser() + "/jules_postprocess/"  # out dir
    INSTITUTION = "Met Office"
    COMMENT = ""

if USE_JULES_PY:
    sys.path.append("/home/h03/hadea/bin")
    import jules  # https://code.metoffice.gov.uk/svn/utils/smstress_jpeg/trunk/jules.py


# #############################################################################
# #############################################################################
def define_inout_paths():
    if L_JULES_ROSE:
        global CONFIG_ARGS
        CONFIG_ARGS = config_parse_args(MIPNAME)  # MIPNAME.ini
        src_dir = CONFIG_ARGS["MODEL_INFO"]["model_input_dir"]
        out_dir = CONFIG_ARGS["MODEL_INFO"]["output_dir"]
        mip = CONFIG_ARGS["MODEL_INFO"]["mipname"]
        if "CMIP" in MIPNAME.upper():
            out_dir = f"{out_dir}/{CONFIG_ARGS['MODEL_INFO']['climate_scenario']}"
    else:  # not L_JULES_ROSE:
        global MIP_INFO
        mip = MIPNAME.split("_")[0]
        MIP_INFO = read_mip_info_no_rose(MIPNAME)
        src_dir = MIP_INFO["src_dir"][MIPNAME] + MIP_INFO["suite_id"][MIPNAME] + "/"
        out_dir = OUT_BASEDIR + MIP_INFO["suite_id"][MIPNAME]
        if "CMIP" in MIPNAME.upper():
            out_dir = f"{OUT_BASEDIR}/{MIP_INFO['out_scenario'][MIPNAME]}"

    # make directory
    try:
        retcode = subprocess.call("mkdir -p " + out_dir, shell=True)
        if retcode != 0:
            raise OSError(f"ERROR: Failed to create directory {out_dir}")
    except OSError as e:
        print(e)
        raise
    return src_dir, out_dir, mip


# #########################################################################
# #########################################################################
def get_diag_for_output(mip):
    # get diagnostics for processing
    if READ_JSON:  # newer version
        if "ISIMIP3" in MIPNAME.upper():
            json_file = "ISIMIP3_variables.json"
        elif "imogen" in MIPNAME.lower():
            json_file = "imogen6_variables.json"
        else:
            json_file = mip + "_variables.json"
        try:
            with open(json_file, "r") as json_file:
                diag_dic = json.load(json_file)
        except FileNotFoundError:
            raise FileNotFoundError(f"ERROR: File {json_file} not found.")
        except json.JSONDecodeError:
            raise json.JSONDecodeError(
                f"ERROR: problem decoding JSON from file {json_file}."
            )
    else:  # older version
        if L_JULES_ROSE:
            diag_dic = suite_postprocessed_variables.get_var_dict()
        elif not L_JULES_ROSE:
            if "ISIMIP2" in mip:
                diag_dic = ISIMIP2_variables.get_var_dict()
            if "ISIMIP3" in mip or "ISIMIP3".lower() in mip:
                diag_dic = ISIMIP3_variables.get_var_dict()
            elif "CMIP" in mip:
                diag_dic = CMIP_variables.get_var_dict()
            # elif "IMOGEN5variant" in mip:
            #     diag_dic = IMOGENvariant_variables.get_var_dict()
            # elif "IMOGEN5" in mip:
            #     diag_dic = IMOGEN_variables.get_var_dict()
            # elif "IMOGEN6" in mip:
            #     diag_dic = IMOGEN6_variables.get_var_dict()
            elif "TRENDY" in mip:
                diag_dic = TRENDY_variables.get_var_dict()
            else:
                raise ValueError("ERROR: need to define dictionary for mip")

    # define all prognistics
    diags = diag_dic.keys()  ## all diagnostics

    if READ_JSON:
        # select required prognostics
        # REMOVE ALL DAILY DATA POSTPROCESSING
        diags = [s for s in diags if "daily" not in s]
        # REMOVE ALL DAILY DATA POSTPROCESSING

    # select specific ones here for testing (L_TESTING=True)
    if L_TESTING:
        diags = sel_diags_test()
    print("####################################################")
    print("Diagnostics:", diags)
    print("####################################################")

    return diags, diag_dic


# #########################################################################
# #########################################################################
def print_run_information(src_dir, out_dir):
    """print information for the model run"""
    print("")
    print("####################################################")
    print("INFORMATION FOR MODEL RUNS TO BE PROCESSED")
    print("MIPNAME", MIPNAME)
    print("L_TESTING", L_TESTING)
    print("L_BACKFILL_MISSING_FILES", L_BACKFILL_MISSING_FILES)
    print("src_dir", src_dir)
    print("out_dir", out_dir)
    print("####################################################")
    print("")


# #############################################################################
# #############################################################################
def main():
    """
    this is the main code reads information and loops through all diagnostics
    """
    # #########################################################################
    # read jules info for relevant mip and define input and output directories
    # #########################################################################
    src_dir, out_dir, mip = define_inout_paths()

    # #########################################################################
    # get diagnstics for processing
    diags, diag_dic = get_diag_for_output(mip)

    # #########################################################################
    # write out information
    print_run_information(src_dir, out_dir)

    # #########################################################################
    # loop through variables
    # #########################################################################
    var_error = "\n"
    outfile_error = "\n"
    for var in diags:
        errorcode = 0
        print("# ##### " + var + " ##########")
        start_time = time.time()
        # ######################################################################
        # define start and end years
        syrall, eyrall, l_alltimes_in_one_file = define_syr_eyr_output(var, MIPNAME)

        for i, syr in enumerate(syrall):  # loop through all years
            eyr = eyrall[i]
            print(
                "INFO: start year: "
                + str(syr)
                + "  end year: "
                + str(eyr)
                + " for: "
                + var
            )

            if not l_alltimes_in_one_file:
                time_cons = iris.Constraint(
                    time=lambda cell: PartialDateTime(syr) <= cell.point
                    and PartialDateTime(eyr) >= cell.point
                )
            else:
                time_cons = None

            # first call to define filenames - check file exists or not
            outfilename, varout, errorcode = write_out_final_cube(
                diag_dic, None, var, out_dir, syr, eyr, l_onlymakefname=True
            )

            if "isimip" in MIPNAME.lower() or "crujra" in MIPNAME.lower():
                outfilename = isimip_func.sort_outfilename_isimip(
                    outfilename, var, varout
                )

            # if file exists do nothing unless L_BACKFILL_MISSING_FILES is false.
            if (
                os.path.exists(outfilename)
                and bool(L_BACKFILL_MISSING_FILES)
                and errorcode == 0
            ):
                print("INFO: file exists: " + os.path.basename(outfilename))
                continue
            elif len(outfilename) == 0:
                print("ERROR: outfilename is empty look for previous error")
                continue
            else:
                print("INFO: outfilename: " + os.path.basename(outfilename))

            cube, errorcode = make_gridded_files(
                src_dir, diag_dic, time_cons, var, syr, eyr
            )
            if errorcode == 1:
                var_error = f"{var_error} ERROR: {var}\n"
                outfile_error = outfile_error + outfilename + "\n"
                continue

            if L_JULES_ROSE:
                cube.attributes["contact"] = CONFIG_ARGS["OUT_INFO"]["contact_email"]
                cube.attributes["institution"] = CONFIG_ARGS["OUT_INFO"]["institution"]
                comment = CONFIG_ARGS["OUT_INFO"]["comment"]
                cube.attributes["comment"] = (
                    comment + "rose-suite is " + CONFIG_ARGS["MODEL_INFO"]["suite_id"]
                )
            elif not L_JULES_ROSE:
                cube.attributes["contact"] = CONTACT_EMAIL
                cube.attributes["institution"] = INSTITUTION
                comment = COMMENT
                cube.attributes["comment"] = (
                    COMMENT + "rose-suite is " + MIP_INFO["suite_id"][MIPNAME]
                )

            outfilename, varout, errorcode = write_out_final_cube(
                diag_dic, cube, var, out_dir, syr, eyr, l_onlymakefname=False
            )
            print(f"SUCCESS: {var} written to {outfilename}")

        end_time = time.time()
        elapsed_time = end_time - start_time
        print("INFO: elapsed time:", elapsed_time, "seconds")
        print("# ##################################")
        print("")

    # show list of broken variables at the end and write out missing files
    print("ERROR: these variables broke: " + var_error)
    retcode = subprocess.call("mkdir -p outinfo", shell=True)
    if retcode != 0:
        raise OSError(f"Failed to create directory outinfo")
    with open("outinfo/errors_warnings_" + MIPNAME + ".outinfo", "w") as text_file:
        text_file.write("These variables broke (ERROR): " + var_error)
        text_file.write("These files are missing (WARNING):\n" + outfile_error)


# #############################################################################
# #############################################################################
def make_gridded_files(src_dir, diag_dic, time_cons, var, syr, eyr):
    """
    make the gridded jules data
    """

    errorcode = 0

    # get input diagnostics
    try:
        if READ_JSON:
            inputdiags = diag_dic[var]["var_name"]
        else:
            inputdiags = diag_dic[var][0]
    except:
        print("ERROR: check variable translations (either diag_dic or json file)")
        print("ERROR:" + var + " is not linked to a specified jules output variable")
        errorcode = 1  # problem with variable name
        return None, errorcode

    if READ_JSON:
        jules_profname = diag_dic[var]["jules_profile"]
        # jules_profname should be a list here
        if len(jules_profname) == 1:
            jules_profname = jules_profname[0]
    else:
        jules_profname = diag_dic[var][4]

    profname_print = (
        ", ".join(jules_profname)
        if isinstance(jules_profname, list)
        else jules_profname
    )
    print(
        f"INFO: {var} -- jules variable: {str(inputdiags)} -- jules_profile: {profname_print}"
    )

    if L_JULES_ROSE:
        print(
            f"INFO: {var} -- scenario: {CONFIG_ARGS['MODEL_INFO']['climate_scenario']}"
        )
        drive_model = CONFIG_ARGS["MODEL_INFO"]["driving_model"]
    else:
        print(f"INFO: {var} -- scenario: {MIP_INFO['in_scenario'][MIPNAME]}")
        drive_model = ""

    # get input filenames and check they exist
    files_in, errorcode = make_infilename(src_dir, jules_profname, syr, eyr)
    if errorcode == 1:  # not all files exist
        return None, errorcode

    func = None
    if READ_JSON:  # get pre-processing function
        process_func = diag_dic.get(var, {}).get("process_func", {})
        func_name = process_func.get("func")
        func = globals()[func_name] if func_name is not None else None
    else:
        if isinstance(diag_dic[var], str):
            print("INFO: " + diag_dic[var])
        else:
            print("INFO: " + ", ".join(map(str, diag_dic[var])))
        if diag_dic[var][5] is not None:
            func = globals()[diag_dic[var][5]]

    if func is not None:
        print(f"INFO: function for pre-processing {func}")
    else:
        if inputdiags.__class__.__name__ == "list":
            print("ERROR: missing function for pre-processing")
            errorcode = 1  # problem with function name
            return None, errorcode

    # read jules data
    if inputdiags.__class__.__name__ == "list":
        # more than one variable read and combined together with a defined func
        cubelist = iris.cube.CubeList([])
        for diag_in in inputdiags:
            cube, errorcode = get_jules_cube(
                diag_in, files_in, MIPNAME, drive_model, time_cons=time_cons
            )
            if errorcode == 1:
                print(f"ERROR: can't find {diag_in} in {files_in[0]}")
                return None, errorcode
            cubelist.append(cube)
        cube, errorcode = func(cubelist, var)
        if errorcode == 1:
            print(f"ERROR: problem with {inputdiags} in {files_in[0]}")
    else:
        # one variable read from the files
        cube, errorcode = get_jules_cube(
            inputdiags, files_in, MIPNAME, drive_model, time_cons=time_cons
        )
        if errorcode == 1:
            print(f"ERROR: cant find {inputdiags} in {files_in[0]}")
            return None, errorcode
        if func is not None:
            cube, errorcode = func(cube, var)
            if errorcode == 1:
                print(f"ERROR: problem wth {func} for preprocessing")
                return None, errorcode

    # sort out units for cube
    if READ_JSON:
        units_for_cube = Unit(diag_dic[var]["units"])
    else:
        units_for_cube = Unit(diag_dic[var][3])

    if cube.units != units_for_cube:
        cube, _ = conv_360days_to_sec(cube, cube.var_name)
        try:
            cube.convert_units(units_for_cube)
        except:
            print(
                f"ERROR: {var}: units are: {str(cube.units)}, required: {units_for_cube}"
            )
            return None, errorcode

    # change variable long name
    if READ_JSON:
        longname_for_cube = diag_dic[var]["long_name"]
    else:
        longname_for_cube = diag_dic[var][1]

    if longname_for_cube is not None:
        cube.long_name = longname_for_cube

    # sort out vegetation fractions when need to separate them
    if READ_JSON:
        fracname = diag_dic[var]["var_name"]
    else:
        fracname = diag_dic[var][0]

    if fracname == "frac" and var not in [
        "landCoverFrac",
        "pft-pft",
        "pft-pft_annual",
        "frac_annual",
    ]:
        cube, errorcode = select_vegfrac(cube, var)
        if errorcode == 1:
            print(f"ERROR: problem with select_vegfrac function")
            return None, errorcode

    # expand grids to global (not coded for imogen)
    cube = expand_to_global(cube, MIPNAME)

    return cube, errorcode


# #############################################################################
# #############################################################################
def define_syr_eyr_output(var, mip_name):
    """
    define start and end year array for output data
    determine whether we want all data for variable in one file
    currently anything _daily in varname will be in multiple files
    """
    if L_JULES_ROSE:
        syrall = [int(CONFIG_ARGS["MODEL_INFO"]["start_year"])]
        eyrall = [int(CONFIG_ARGS["MODEL_INFO"]["end_year"])]
    else:  # not L_JULES_ROSE:
        syrall = [int(MIP_INFO["start_year"][mip_name])]
        eyrall = [int(MIP_INFO["end_year"][mip_name])]
    l_alltimes_in_one_file = True

    if L_TESTING:
        eyrall = [syrall[0] + 11]

    mip_name = mip_name.upper()
    if "daily" in var:
        if "CMIP" in mip_name:
            syrall, eyrall = cmip_func.define_years_daily_cmip(syrall, eyrall)
        elif "ISIMIP2" in mip_name:
            syrall, eyrall = isimip_func.define_years_daily_isimip2b()
        elif "CRUJRA" in mip_name:
            syrall, eyrall = isimip_func.define_years_daily_wrpmip()
        else:
            syrall = [-1]
            eyrall = [-1]
            raise ValueError("need to define how to split daily data")
        l_alltimes_in_one_file = False  # split output files for daily data

    return syrall, eyrall, l_alltimes_in_one_file


# #############################################################################
# #############################################################################
def add_time_middle(cube, field, filename):
    """
    Round time bounds to nearest day, change time coord points to
    be midway between the bounds
    """
    cube.coord("time").bounds = 3600 * (
        np.rint(cube.coord("time").bounds.astype("int64") / 3600.0)
    )
    cube.coord("time").points = (
        cube.coord("time").bounds[:, 0] + cube.coord("time").bounds[:, 1]
    ) / 2.0


# #############################################################################
# #############################################################################
def make_infilename(src_dir, jules_profname, syr, eyr):
    """
    make filename of all the required infiles
    """
    # only files between start_year and end_year
    # CONFIG_ARGS = config_parse_args(MIPNAME)

    errorcode = 0
    files_in = []

    years = [str(year) for year in np.arange(syr, eyr + 1)]
    if "imogen" in MIPNAME.lower():
        files_in, errorcode = imogen_func.make_infilename_imogen(
            src_dir, jules_profname, years
        )
        if errorcode == 1:
            return files_in, errorcode
    else:
        if not L_JULES_ROSE:
            if isinstance(jules_profname, list):
                for profname in jules_profname:
                    files_in.extend(
                        [
                            f"{src_dir}{MIP_INFO['run_name'][MIPNAME]}{MIP_INFO['in_scenario'][MIPNAME]}.{profname}.{year}.nc"
                            for year in years
                        ]
                    )
            else:
                files_in = [
                    f"{src_dir}{MIP_INFO['run_name'][MIPNAME]}{MIP_INFO['in_scenario'][MIPNAME]}.{jules_profname}.{year}.nc"
                    for year in years
                ]
        elif L_JULES_ROSE:
            if (
                "isimip" in CONFIG_ARGS["MODEL_INFO"]["mipname"].lower()
                or "crujra" in CONFIG_ARGS["MODEL_INFO"]["mipname"].lower()
            ):
                if isinstance(jules_profname, list):
                    for profname in jules_profname:
                        files_in.extend(
                            [
                                (
                                    f"{src_dir}"
                                    f"{CONFIG_ARGS['MODEL_INFO']['mipname'].lower()}_"
                                    f"{CONFIG_ARGS['MODEL_INFO']['configname']}_"
                                    f"{CONFIG_ARGS['MODEL_INFO']['driving_model']}_"
                                    f"{CONFIG_ARGS['MODEL_INFO']['climate_scenario']}."
                                    f"{profname}.{year}.nc"
                                )
                                for year in years
                            ]
                        )
                else:
                    files_in = [
                        (
                            f"{src_dir}"
                            f"{CONFIG_ARGS['MODEL_INFO']['mipname'].lower()}_"
                            f"{CONFIG_ARGS['MODEL_INFO']['configname']}_"
                            f"{CONFIG_ARGS['MODEL_INFO']['driving_model']}_"
                            f"{CONFIG_ARGS['MODEL_INFO']['climate_scenario']}."
                            f"{jules_profname}.{year}.nc"
                        )
                        for year in years
                    ]
            else:
                if isinstance(jules_profname, list):
                    for profname in jules_profname:
                        files_in.extend(
                            [
                                (
                                    f"{src_dir}"
                                    f"{CONFIG_ARGS['MODEL_INFO']['driving_model']}_"
                                    f"{CONFIG_ARGS['MODEL_INFO']['climate_scenario']}."
                                    f"{profname}.{year}.nc"
                                )
                                for year in years
                            ]
                        )
                else:
                    files_in = [
                        (
                            f"{src_dir}"
                            f"{CONFIG_ARGS['MODEL_INFO']['driving_model']}_"
                            f"{CONFIG_ARGS['MODEL_INFO']['climate_scenario']}."
                            f"{jules_profname}.{year}.nc"
                        )
                        for year in years
                    ]
            print(f"INFO: First input file: {files_in[0]}")

    # check files exist
    existing_files = [glob.glob(file) for file in files_in]
    if len(existing_files) > 0:
        if isinstance(existing_files[0], list):
            existing_files = [file for sublist in existing_files for file in sublist]

    if len(existing_files) == 0:
        print("ERROR: No input files found")
        errorcode = 1

    return existing_files, errorcode


# #############################################################################
# #############################################################################
def make_outfilename(out_dir, outprofile, var, syr, eyr, diag_dic):
    """
    make filename of outfilename
    """
    errorcode = 0
    if "CMIP" in MIPNAME.upper() and "IMOGEN" not in MIPNAME.upper():
        outfilename = cmip_func.make_outfilename_cmip(
            out_dir, outprofile, var, syr, eyr
        )
    elif "TRENDY" in MIPNAME.upper():
        if not L_JULES_ROSE:
            outfilename = (
                out_dir
                + "/"
                + MIP_INFO["model"][MIPNAME]
                + "_"
                + MIP_INFO["out_scenario"][MIPNAME]
                + "_"
                + var
                + "_"
                + outprofile
                + ".nc"
            )
        if L_JULES_ROSE:
            outfilename = (
                out_dir
                + "/"
                + CONFIG_ARGS["MODEL_INFO"]["model"]
                + "_"
                + CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]
                + "_"
                + var
                + "_"
                + outprofile
                + ".nc"
            )
    elif "isimip" in MIPNAME.lower() or "crujra" in MIPNAME.lower():
        outfilename = isimip_func.make_outfilename_isimip(
            out_dir, outprofile, var, syr, eyr
        )
    elif "imogen" in MIPNAME.lower():
        outfilename, errorcode = imogen_func.make_outfilename_imogen(
            out_dir, outprofile, var, syr, eyr, diag_dic
        )
    else:
        sys.exit("need outfilename for " + MIPNAME)
    return outfilename, errorcode


# #############################################################################
# #############################################################################
def write_out_final_cube(diag_dic, cube, var, out_dir, syr, eyr, l_onlymakefname=True):
    """
    write out cube
    """
    errorcode = 0
    # print stats if interested
    if L_TESTING and not l_onlymakefname:
        print(cube)
        print("INFO: lazy data still? ", cube.has_lazy_data())

    # sort out varout and outprofile
    varout, outprofile = sort_varout_outprofile_name(var)
    if "cmip" in MIPNAME.lower() and not "imogen" in MIPNAME.lower():
        if READ_JSON:
            outprofile = diag_dic[var]["cmip_profile"]
        else:
            outprofile = diag_dic[var][6]

    if not l_onlymakefname:
        iris.coord_categorisation.add_year(cube, "time")
        syr_data = int(np.min(cube.coord("year").points))
        if (cube.coord("year").points[1] - cube.coord("year").points[0]) / 2.0 > 0.9:
            syr_data = int(
                syr_data
                - (cube.coord("year").points[1] - cube.coord("year").points[0]) / 2.0
            )
        if syr_data != syr:
            print("INFO: start years (data, filename)", syr_data, syr)
            print("ERROR: start years in data and filenames are incompatible")
            raise ValueError(
                "ERROR: start years in data and filenames are incompatible"
            )
        eyr_data = int(np.max(cube.coord("year").points))
        if (
            abs(cube.coord("year").points[-1] - cube.coord("year").points[-2]) / 2.0
        ) > 0.9:
            eyr_data = int(
                eyr_data
                - (cube.coord("year").points[-1] - cube.coord("year").points[-2]) / 2.0
            )
        if eyr_data != eyr:
            print("INFO: end years  (data, filename)", eyr_data, eyr)
            print("ERROR: end years in data and filenames are incompatible")
            print("INFO: changing end years in filenames and continuing")
            eyr = eyr_data
        cube.remove_coord("year")
        fill_value = np.float32(
            1.0e20
        )  # this might need to change with different mips?

        # sort out formatting of ISIMIP output
        if "isimip" in MIPNAME.lower() or "crujra" in MIPNAME.lower():
            cube = isimip_func.sort_isimip_cube(cube, outprofile)

    if "cmip" in MIPNAME.lower():
        varout = diag_dic[var]["cmip_varname"]
    outfilename, errorcode = make_outfilename(
        out_dir, outprofile, varout, syr, eyr, diag_dic
    )
    if errorcode == 1:
        return outfilename, varout, errorcode

    if not l_onlymakefname:
        if "imogen" in MIPNAME.lower():
            if READ_JSON:
                cube.var_name = diag_dic[var]["cmip_varname"]
        else:
            cube.var_name = varout
        # might have to do some work ISIMIP2 output files
        # https://www.isimip.org/protocol/preparing-simulation-files/
        # quality-check-of-your-simulation-data
        if "pft" in var and (
            "ISIMIP" in MIPNAME.upper() or "crujra" in MIPNAME.lower()
        ):
            # need to separate cube by pft
            for ipft in cube.coord("vegtype").points:
                isimip_func.sort_and_write_pft_cube(
                    varout, cube, outfilename, ipft, fill_value
                )
        elif "-pool" in var and (
            "ISIMIP" in MIPNAME.upper() or "crujra" in MIPNAME.lower()
        ):
            for ipool in cube.coord("scpool").points:
                isimip_func.sort_and_write_pool_cube(
                    varout, cube, outfilename, ipool, fill_value
                )
        else:
            divide_files = False
            chunksizes = define_chunksizes(cube)
            print("INFO: cube should still have lazy data", cube.has_lazy_data())

            coord_names = [coord.name() for coord in cube.coords()]

            if "depth" in coord_names and len(cube.coord("depth").points) > 10:
                print("INFO: >10 soil levels which means files are very big")
                divide_files = True  # divide into separate files

            if "sclayer" in coord_names and len(cube.coord("sclayer").points) > 10:
                print("INFO: >10 soil bgc levels which means files can be very big")
                if READ_JSON:  # get pre-processing function
                    process_func = diag_dic.get(var, {}).get("process_func", {})
                    func_name = process_func.get("func")
                else:
                    if diag_dic[var][5] is not None:
                        func_name = diag_dic[var][5]
                if func_name == "layered_soilbgc_func":
                    divide_files = True  # divide into separate files

            if "frac_name" in coord_names:
                cube.remove_coord("frac_name")

            if not L_TESTING:
                if divide_files:
                    subtimes = 10  # change this in anger
                    cube_count = -1
                    iris.coord_categorisation.add_year(cube, "time")
                    for i in range(0, len(cube.coord("time").points), subtimes):
                        if "imogen" in MIPNAME.lower():
                            times = cube.coord("time")[i : i + subtimes]
                            calendar = cube.coord("time").units.calendar
                            tmptimes = [
                                num2date(
                                    time,
                                    units=cube.coord("time").units.name,
                                    calendar=calendar,
                                )
                                for time in times.points
                            ]
                            years = [date.year for date in tmptimes]
                            time_cons = iris.Constraint(
                                year=lambda cell: years[0] <= cell.point
                                and years[-1] >= cell.point
                            )
                            cubeout = cube.extract(time_cons)
                        else:
                            cubeout = cube[i : i + subtimes].copy()
                        cube_count = cube_count + 1
                        print("INFO: lazy data still? ", cubeout.has_lazy_data())
                        if "latitude" in coord_names or "lat" in coord_names:
                            cubeout.data.mask[np.isnan(cubeout.data)] = (
                                True  # breaks lazy data
                            )
                        outfilenametmp = (
                            outfilename[:-3] + "_separate" + str(cube_count) + ".nc"
                        )
                        netcdf_format = "NETCDF4_CLASSIC"
                        if "cmip" in MIPNAME.lower():
                            netcdf_format = "NETCDF4"
                        iris.save(
                            cubeout,
                            outfilenametmp,
                            fill_value=fill_value,
                            zlib=True,
                            netcdf_format=netcdf_format,
                            chunksizes=chunksizes,
                            contiguous=False,
                            complevel=9,
                        )
                else:
                    if "latitude" in coord_names or "lat" in coord_names:
                        cube.data.mask[np.isnan(cube.data)] = True  # breaks lazy data
                    netcdf_format = "NETCDF4_CLASSIC"
                    if "cmip" in MIPNAME.lower():
                        netcdf_format = "NETCDF4"
                    iris.save(
                        cube,
                        outfilename,
                        fill_value=fill_value,
                        zlib=True,
                        netcdf_format=netcdf_format,
                        chunksizes=chunksizes,
                        contiguous=False,
                        complevel=9,
                    )

                if "ISIMIP2" in MIPNAME.upper():
                    retcode = subprocess.call(
                        "mv " + outfilename + " " + outfilename + "4", shell=True
                    )
                    if retcode != 0:
                        print("ERROR: mv " + outfilename + " " + outfilename + "4")
                        raise RuntimeError("mv command failed")

                if "ISIMIP3" in MIPNAME.upper():
                    isimip_func.rename_cfcompliant_to_isimip(outfilename, cube)

    return outfilename, varout, errorcode


# #############################################################################
# #############################################################################
def expand_to_global(cube, mip_name):
    """
    Define new cube with global grid for output
    """
    mip_name = mip_name.upper()
    if "IMOGEN" not in mip_name:
        if "ISIMIP" in mip_name or "CRUJRA" in mip_name:
            latitude, longitude, nlat, nlon = isimip_func.make_global_grid_0p5()
        elif "CMIP" in mip_name or "TRENDY" in mip_name:
            latitude, longitude, nlat, nlon = cmip_func.make_global_grid_n96e()
        else:
            raise ValueError(
                "add definition of the full global grid in expand_to_global"
            )

        cube_fullgrid = Cube(
            np.zeros((nlat, nlon), np.float32),
            dim_coords_and_dims=[(latitude, 0), (longitude, 1)],
        )
        cube = cube.regrid(
            cube_fullgrid, iris.analysis.Nearest(extrapolation_mode="mask")
        )
    return cube


# #############################################################################
# #############################################################################
def add_soil_info(cube):
    """
    add soil depth information
    """
    soil_coord = cube.coord("soil")
    if len(soil_coord.points) == 4:  # standard jules
        dzsoil_io = [0.1, 0.25, 0.65, 2.0]
    elif len(soil_coord.points) == 14:
        dzsoil_io = (
            0.1,
            0.2,
            0.2,
            0.2,
            0.3,
            0.3,
            0.3,
            0.4,
            0.4,
            0.4,
            1.0,
            1.0,
            3.0,
            3.0,
        )
    elif len(soil_coord.points) == 20:
        dzsoil_io = (
            0.05000000,
            0.09012505,
            0.12721053,
            0.16245048,
            0.19637876,
            0.22929710,
            0.26139865,
            0.29281714,
            0.32365039,
            0.35397289,
            0.38384314,
            0.41330824,
            0.44240686,
            0.47117131,
            0.49962893,
            0.52780316,
            0.55571432,
            0.58338013,
            0.61081623,
            0.63803647,
        )
    elif len(soil_coord.points) == 10:
        dzsoil_io = [0.1, 0.15, 0.22, 0.33, 0.5, 0.74, 1.1, 1.64, 2.45, 3.66]
    else:
        raise ValueError("need depth coordinates to be set properly")

    bottom_of_layer = np.cumsum(dzsoil_io)
    top_of_layer = np.append([0.0], bottom_of_layer[:-1])
    bounds = np.column_stack((top_of_layer, bottom_of_layer))
    soil_coord.points = calc_mid_depth(bottom_of_layer)
    soil_coord.bounds = bounds

    soil_coord.units = "m"
    soil_coord.rename("depth")
    cube.coord("depth").long_name = "depth of middle of soil layer"

    return cube


# #############################################################################
# #############################################################################
def get_jules_cube(diag_in, files_in, mip_name, drive_model, time_cons=None):
    """
    read and pre-process jules cube
    """
    errorcode = 0

    load_method = "jules load" if USE_JULES_PY else "jules_xarray load"
    print(f"INFO: load cube using {load_method} for {diag_in}")

    variable_cons = iris.Constraint(cube_func=(lambda c: c.var_name == diag_in))

    if "imogen" in mip_name.lower():  # and "ens" in mip_name.lower():
        # read ensembles for imogen if available
        try:
            cube = imogen_func.read_ensemble(
                files_in, variable_cons, time_cons, diag_in, drive_model
            )
        except:
            errorcode = 1
            return None, errorcode
            # raise ValueError(f"{diag_in} not available in {files_in[0]}")
    else:
        # any mips which are not imogen ensembles
        try:
            if USE_JULES_PY:
                cube = jules.load(
                    files_in,
                    variable_cons & time_cons,
                    missingdata=ma.masked,
                    callback=add_time_middle,
                )
            else:
                cube = jules_xarray.load(files_in, variable_cons & time_cons)
        except:
            errorcode = 1
            return None, errorcode
        # raise ValueError(diag_in + " not working for " + files_in[0])

        if isinstance(cube, iris.cube.CubeList):
            for ijk in np.arange(0, int(len(cube))):
                for key in list(cube[ijk].attributes.keys()):
                    del cube[ijk].attributes[key]
                if ijk > 0:
                    cube[ijk].coord("time").convert_units(cube[0].coord("time").units)
            cube = cube.concatenate_cube()

    # sort out coordinates
    for coord_name in ["type", "tile", "pft", "soil"]:
        if coord_name in [coord.name() for coord in cube.coords()]:
            cube = (
                add_tile_info(cube, coord_name)
                if coord_name != "soil"
                else add_soil_info(cube)
            )

    # might need these for TRENDY
    # if var == "tsl":
    #    cube.coord("soil").long_name = "stlayer"
    # if var == "msl":
    #    cube.coord("soil").long_name = "smlayer"

    return cube, errorcode


# #############################################################################
# #############################################################################
def stat(cube):
    """
    print out cube statistics
    """
    print(" ~~~ stat start ~~~ ")
    print(cube)
    print(" ~~~ stat ~~~ ")
    print(
        "min, mean, max data: "
        + str(np.min(cube.data))
        + ","
        + str(np.mean(cube.data))
        + ","
        + str(np.max(cube.data))
    )
    print(
        "NAN, min, mean, max data: "
        + str(np.nanmin(cube.data))
        + ","
        + str(np.nanmean(cube.data))
        + ","
        + str(np.nanmax(cube.data))
    )
    print("type of cube data: " + cube.data.__class__.__name__)
    if cube.data.__class__.__name__ == "MaskedArray":
        print("number of masked points: " + str(np.sum(cube.data.mask)))
        print("number of NOT masked points: " + str(np.sum(~cube.data.mask)))
    print("shape of data: " + str(np.shape(cube.data)))
    print(
        "NANS?: "
        + str(np.sum(np.isnan(cube.data)))
        + "   "
        + "not NANS?: "
        + str(np.sum(~np.isnan(cube.data)))
    )
    print(" ~~~ stat end ~~~ ")


# #############################################################################
# #############################################################################
def calc_mid_depth(bottom_of_layer):
    """
    given bottom_of_layer of soil layers calculate the middle of each layer
    """
    mid_depth = np.zeros_like(bottom_of_layer)
    mid_depth[0] = bottom_of_layer[0] / 2.0
    mid_depth[1:] = bottom_of_layer[:-1] + np.diff(bottom_of_layer) / 2.0

    return mid_depth


# #############################################################################
# #############################################################################
def define_chunksizes(cube):
    """
    define the chunksizes for the output file
    """
    chunksizes = None
    all_coord_names = [coord.name() for coord in cube.coords()]
    ndims = len(cube.shape)
    ndims_orig = ndims

    has_realization = "realization" in all_coord_names and ndims_orig == len(
        all_coord_names
    )
    if has_realization:
        ndims -= 1

    if 2 < ndims < 6:
        chunksizes = [1] * ndims
        chunksizes[-3:] = cube.shape[-3:]
    else:
        print("ERROR: chunksizes are undefined set to shape of data")
        chunksizes = list(cube.shape)

    if has_realization and chunksizes:
        chunksizes = [1] + chunksizes

    if "time" == [coord.name() for coord in cube.coords()][0]:
        chunksizes[0] = 1

    return chunksizes


# #############################################################################
# #############################################################################
if __name__ == "__main__":
    main()
