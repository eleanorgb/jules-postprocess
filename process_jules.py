'''
script to postprocess jules output 
suitable for ISIMIP, TRENDY, ILAMB, IMOGEN, CMIP
BE CAREFUL with nitrogen on and off and npp and nlim
sometimes needs lots of memory - run on spice
'''
import os
import sys
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
from landcover_types import select_vegfrac
from landcover_types import add_tile_info
from sel_diags_test import sel_diags_test
from read_input import parse_args
from read_input import config_parse_args
from read_input import read_mip_info_no_rose
from sort_varout_outprofile_name import sort_varout_outprofile_name
import isimip_func   #isimip runs only
import imogen_func   #imogen runs only
import cmip_func     # cmip runs only

warnings.filterwarnings("ignore")

# #########################################################################
# read command line arguments
global MIPNAME, L_TESTING, L_JULES_ROSE
MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE = parse_args()
# #########################################################################

if not L_JULES_ROSE:
    import ISIMIP2_variables
    import ISIMIP3_variables
    import CMIP_variables
    import TRENDY_variables
    import IMOGEN_variables
    import IMOGENvariant_variables
    import IMOGEN6_variables
    sys.path.append("/home/h03/hadea/bin")
    import jules #https://code.metoffice.gov.uk/svn/utils/smstress_jpeg/trunk/jules.py
elif L_JULES_ROSE:
    import suite_postprocessed_variables
    diag_dic = suite_postprocessed_variables.get_var_dict()
    import jules
    
# #############################################################################
if not L_JULES_ROSE:
    CONTACT_EMAIL = "eleanor.burke@metoffice.gov.uk"  # dont use my email!
    OUT_BASEDIR = "/scratch/"+getpass.getuser()+"/jules_postprocess/" # out dir
    INSTITUTION = "Met Office"
    COMMENT = ""

# #############################################################################
def main():
    """
    this is the main code reads information and loops through all diagnostics
    """
    # #########################################################################
    # read jules info for relevant mip and define input and output directories
    # #########################################################################
    if L_JULES_ROSE:
        global CONFIG_ARGS
        CONFIG_ARGS = config_parse_args(MIPNAME)  # MIPNAME.ini
        src_dir = CONFIG_ARGS['MODEL_INFO']['model_input_dir']
        out_dir =  CONFIG_ARGS['MODEL_INFO']['output_dir']
        mip = CONFIG_ARGS['MODEL_INFO']['mipname']
        if "CMIP" in MIPNAME.upper():
            out_dir = out_dir+CONFIG_ARGS["MODEL_INFO"]["suite_id"]+"_"+\
                   CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]

    # #########################################################################
    elif not L_JULES_ROSE:
        global MIP_INFO
        mip = MIPNAME.split('_')[0]
        MIP_INFO = read_mip_info_no_rose(MIPNAME)
        src_dir = MIP_INFO["src_dir"][MIPNAME]+MIP_INFO["suite_id"][MIPNAME]+"/"
        out_dir = OUT_BASEDIR+MIP_INFO["suite_id"][MIPNAME]
        if "CMIP" in MIPNAME.upper():
            out_dir = OUT_BASEDIR+MIP_INFO["suite_id"][MIPNAME]+"_"+\
                   MIP_INFO["out_scenario"][MIPNAME]

    # #########################################################################
    retcode = subprocess.call("mkdir -p "+ out_dir, shell=True)
    if retcode != 0:
        print("mkdir -p "+ out_dir)
        sys.exit("mkdir broken")
    error = "\n"
    outfile_error = ""

    # #########################################################################
    # get diagnstics for processing
    if L_JULES_ROSE:
        diag_dic = suite_postprocessed_variables.get_var_dict()
    elif not L_JULES_ROSE:
        if "ISIMIP2" in mip:
            diag_dic = ISIMIP2_variables.get_var_dict()
        if "ISIMIP3" in mip:
            diag_dic = ISIMIP3_variables.get_var_dict()
        elif "CMIP" in mip:
            diag_dic = CMIP_variables.get_var_dict()
        elif "IMOGEN5variant" in mip:
            diag_dic = IMOGENvariant_variables.get_var_dict()
        elif "IMOGEN5" in mip:
            diag_dic = IMOGEN_variables.get_var_dict()
        elif "IMOGEN6" in mip:
            diag_dic = IMOGEN6_variables.get_var_dict()
        elif "TRENDY" in mip:
            diag_dic = TRENDY_variables.get_var_dict()
        else:
            print("need to define dictionary for mip")
            sys.exit()

    diags = diag_dic.keys() ## all diagnostics

    # select specific ones here for testing (L_TESTING=True)
    if L_TESTING:
        diags = sel_diags_test()
    print("Diagnostics:", diags)
    # #########################################################################

    # #########################################################################
    # write out information
    print("####################################################")
    print("INFORMATION FOR MODEL RUNS TO BE PROCESSED")
    print("MIPNAME", MIPNAME)
    print("L_TESTING", L_TESTING)
    print("L_BACKFILL_MISSING_FILES", L_BACKFILL_MISSING_FILES)
    print("src_dir", src_dir)
    print("out_dir", out_dir)
    print("####################################################")

    # #########################################################################
    # loop through variables
    for var in diags:
        # ######################################################################
        # define start and end years
        syrall, eyrall, l_alltimes_in_one_file = define_syr_eyr_output(var)

        for i, syr in enumerate(syrall): # loop through all years
            eyr = eyrall[i]
            print("start year: "+str(syr)+"  end year: "+str(eyr)+" for: "+var)

            if not l_alltimes_in_one_file:
                time_cons = iris.Constraint(time=lambda cell: \
                                            PartialDateTime(syr) \
                                            <= cell.point and \
                                            PartialDateTime(eyr) >= cell.point)
            else:
                time_cons = None

            # first call to define filenames - check file exists or not
            outfilename, varout = write_out_final_cube(diag_dic, None, var,
                                                       out_dir, syr, eyr,
                                                       l_onlymakefname=True)
            if "ISIMIP" in MIPNAME.upper():
                outfilename = isimip_func.sort_outfilename_isimip(outfilename,
                                                                  var, varout)
            print("outfilename: "+outfilename)

            # if file exists do nothing unless BACKFILL is false.
            if os.path.exists(outfilename) and bool(L_BACKFILL_MISSING_FILES):
                print("File exists: "+outfilename)
                continue

            try:
                cube = make_gridded_files(src_dir, diag_dic,
                                                       time_cons, var, syr, eyr)
                if L_JULES_ROSE:
                    cube.attributes["contact"] = \
                                        CONFIG_ARGS['OUT_INFO']['contact_email']
                    cube.attributes["institution"] = \
                                          CONFIG_ARGS['OUT_INFO']['institution']
                    comment = CONFIG_ARGS['OUT_INFO']['comment']
                    cube.attributes["comment"] = comment+"rose-suite is "+\
                                           CONFIG_ARGS["MODEL_INFO"]["suite_id"]
                elif not L_JULES_ROSE:
                    cube.attributes["contact"] = \
                                        CONTACT_EMAIL
                    cube.attributes["institution"] = \
                                          INSTITUTION
                    comment = COMMENT
                    cube.attributes["comment"] = COMMENT+"rose-suite is "+\
                                           MIP_INFO["suite_id"][MIPNAME]

                outfilename, varout = write_out_final_cube(diag_dic, cube,
                                                          var, out_dir,
                                                          syr, eyr,
                                                          l_onlymakefname=False)
            except:
                print("ERROR "+MIPNAME+": "+var)
                error = error+" "+var+"\n"
                outfile_error = outfile_error+outfilename+"\n"
                if not L_TESTING:
                    pass
                else:
                    raise

    # show list of broken variables at the end and write out missing files
    print("THESE variables broke: "+error)
    retcode = subprocess.call("mkdir -p outinfo", shell=True)
    with open("outinfo/"+MIPNAME+".outinfo", "w") as text_file:
        text_file.write("These variables broke (ERROR): "+error)
        text_file.write("These files are missing (WARNING):\n"+outfile_error)
# #############################################################################


# #############################################################################
def make_gridded_files(src_dir, diag_dic, time_cons, var, syr, eyr):
    """
    make the gridded jules data
    """

    # read input diagnostics
    try:
        inputdiags = diag_dic[var][0]
    except:
        print("check dictionary of variable translations (see diag_dic)")
        print(var + " is not linked to a specified jules output variable")
        sys.exit()
    #print(inputdiags.__class__.__name__)
    jules_profname = diag_dic[var][4]
    print(var+": variable: "+str(inputdiags)+" jules_profile: "+jules_profname)
    files_in = make_infilename(src_dir, jules_profname, syr, eyr)
    if inputdiags.__class__.__name__ == "list":
        # more than one variable read and combined together with a defined func
        cubelist = iris.cube.CubeList([])
        for diag_in in inputdiags:
            cubelist.append(get_jules_cube(diag_in, files_in,
                                           time_cons=time_cons))
        if diag_dic[var][5] is not None:
            print("function for pre-processing "+diag_dic[var][5])
        else:
            sys.exit("MISSING function for pre-processing")
        func = globals()[diag_dic[var][5]]
        cube = func(cubelist, var)  # return a cube and then work with it!
    else:
        # one variable read from the files
        cube = get_jules_cube(inputdiags, files_in, time_cons=time_cons)
        if diag_dic[var][5] is not None:
            print("function for pre-processing: "+diag_dic[var][5])
            func = globals()[diag_dic[var][5]]
            cube = func(cube, var)  # return a cube and then work with it!

    # sort out units and names for jules cube
    if cube.units != Unit(diag_dic[var][3]):
        try:
            cube.convert_units(diag_dic[var][3])
        except:
            if diag_dic[var][2] is not None:
                print("multiply cube data by "+str(diag_dic[var][2]))
                print(var+" units required = "+diag_dic[var][3])
            if diag_dic[var][2] is not None and diag_dic[var][3] is not None:
                cube.data = cube.core_data() * diag_dic[var][2]
                cube.units = Unit(diag_dic[var][3])
                pass
            else:
                if diag_dic[var][3] is not None:
                    print("units are : "+str(cube.units)+\
                          " required: "+diag_dic[var][3])
                else:
                    print("units are : "+str(cube.units)+" diag_dic is broken")
                print("fix units/scaling in get_var_dict for "+var)
                raise
    if diag_dic[var][1] is not None: # change variable long name
        cube.long_name = diag_dic[var][1]

   # sort out vegetation fractions when need to separate them
    if diag_dic[var][0] == "frac" and \
                 var not in ["landCoverFrac", "pft-pft",
                             "pft-pft_annual", "frac_annual"]:
        cube = select_vegfrac(cube, var)

    # need to expand grids to global (not done for imogen, maybe wont bother)
    if "IMOGEN" not in MIPNAME.upper():
        cube = expand_to_global(cube)

    return cube
# #############################################################################


# #############################################################################
def define_syr_eyr_output(var):
    """
    define start and end year array for output data
    determine whether we want all data for variable in one file
    currently anything _daily in varname will be in multiple files
    """
    if L_JULES_ROSE:
        syrall = [int(CONFIG_ARGS["MODEL_INFO"]["start_year"])]
        eyrall = [int(CONFIG_ARGS["MODEL_INFO"]["end_year"])]
    elif not L_JULES_ROSE:
        syrall = [int(MIP_INFO["start_year"][MIPNAME])]
        eyrall = [int(MIP_INFO["end_year"][MIPNAME])]
    l_alltimes_in_one_file = True

    if L_TESTING:
        eyrall = [syrall[0] + 11]

    if "daily" in var:
        if "CMIP" in MIPNAME.upper():
            syrall, eyrall = cmip_func.define_years_daily_cmip()
        elif "ISIMIP2B" in MIPNAME.upper():
            syrall, eyrall = isimip_func.define_years_daily_isimip2b()
        else:
            syrall = [-1]
            eyrall = [-1]
            print("need to define how to split daily data")
            sys.exit("need to define how to split daily data")
        l_alltimes_in_one_file = False    # split output files for daily data

    return syrall, eyrall, l_alltimes_in_one_file

# #############################################################################
def add_time_middle(cube, field, filename):
    """
    Round time bounds to nearest day, change time coord points to
    be midway between the bounds
    """
    cube.coord("time").bounds = \
           3600*(np.rint(cube.coord("time").bounds.astype("int64")/3600.0))
    cube.coord("time").points = \
           (cube.coord("time").bounds[:, 0] + \
            cube.coord("time").bounds[:, 1])/2.0
# #############################################################################


# #############################################################################
def make_infilename(src_dir, jules_profname, syr, eyr):
    """
    make filename of all the required infiles
    """
    # only files between start_year and end_year
    # CONFIG_ARGS = config_parse_args(MIPNAME)
    years = [str(year) for year in np.arange(syr, eyr+1)]
    if "IMOGEN" in MIPNAME.upper():
        files_in = imogen_func.make_infilename_imogen(src_dir, jules_profname,
                                                      years)
    else:
        if not L_JULES_ROSE:
            files_in = [src_dir+MIP_INFO["run_name"][MIPNAME]+\
                    MIP_INFO["in_scenario"][MIPNAME]+\
                    "."+jules_profname+"."+year+".nc" for year in years]
        elif L_JULES_ROSE:
            if "isimip" in CONFIG_ARGS["MODEL_INFO"]["mipname"]:
                files_in = [src_dir+\
                            CONFIG_ARGS["MODEL_INFO"]["mipname"]+"_"+\
                            CONFIG_ARGS["MODEL_INFO"]["configname"]+"_"+\
                            CONFIG_ARGS["MODEL_INFO"]["driving_model"]+"_"+\
                            CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]+\
                            "."+jules_profname+"."+year+".nc" for year in years]
            else:
                files_in = [src_dir+CONFIG_ARGS["MODEL_INFO"]["driving_model"]+"_"+\
                            CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]+\
                            "."+jules_profname+"."+year+".nc" for year in years]
    print("First input file:")
    print(files_in[0])

    return files_in
# #############################################################################


# #############################################################################
def make_outfilename(out_dir, outprofile, var, syr, eyr):
    """
    make filename of outfilename
    """
    if "CMIP" in MIPNAME.upper():
        outfilename = cmip_func.make_outfilename_cmip(out_dir, outprofile,
                                                      var, syr, eyr)
    elif "TRENDY" in MIPNAME.upper():
        if not L_JULES_ROSE:
            outfilename = out_dir+"/"+MIP_INFO["model"][MIPNAME]+"_"+\
                  MIP_INFO["out_scenario"][MIPNAME]+"_"+var+"_"+\
                  outprofile+".nc"
        if L_JULES_ROSE:
            outfilename = out_dir+"/"+CONFIG_ARGS["MODEL_INFO"]["model"]+"_"+\
                  CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]+"_"+var+"_"+\
                  outprofile+".nc"
    elif "ISIMIP" in MIPNAME.upper():
        outfilename = isimip_func.make_outfilename_isimip(out_dir, outprofile,
                                                          var, syr, eyr)
    elif "IMOGEN" in MIPNAME.upper():
        outfilename = imogen_func.make_outfilename_imogen(out_dir, outprofile,
                                                          var, syr, eyr)

    else:
        sys.exit("need outfilename for "+MIPNAME)
    return outfilename
# #############################################################################


# #############################################################################
def write_out_final_cube(diag_dic, cube, var, out_dir, syr,
                         eyr, l_onlymakefname=True):
    """
    write out cube
    """
    #print stats if interested
    if L_TESTING:
        if not l_onlymakefname:
            #this realises lazy data
            try:
                print(cube)
            except:
                print("error in print cube statement")
            print("lazy data still? ",cube.has_lazy_data())

    # sort out varout and outprofile
    varout, outprofile = sort_varout_outprofile_name(var)
    if "CMIP" in MIPNAME.upper():
        outprofile = diag_dic[var][6]

    if not l_onlymakefname:
        iris.coord_categorisation.add_year(cube, "time")
        if np.min(cube.coord("year").points) != syr:
            print("start years", np.min(cube.coord("year").points), syr)
            sys.exit("start years in data and filenames are incompatible")
        if np.max(cube.coord("year").points) != eyr:
            print("end years", np.max(cube.coord("year").points), eyr)
            sys.exit("end years in data and filenames are incompatible")
        cube.remove_coord("year")
        fill_value = 1.e+20  # this might need to change with different mips?
        # sort out formatting of ISIMIP output
        if "ISIMIP" in MIPNAME.upper():
            cube = isimip_func.sort_isimip_cube(cube, outprofile)

    outfilename = make_outfilename(out_dir, outprofile, varout, syr, eyr)

    if not l_onlymakefname:
        cube.var_name = varout
        # might have to do some work ISIMIP2 output files
        # https://www.isimip.org/protocol/preparing-simulation-files/
        # quality-check-of-your-simulation-data
        if "pft" in var and "ISIMIP" in MIPNAME.upper():
            # need to separate cube by pft
            for ipft in cube.coord("vegtype").points:
                isimip_func.sort_and_write_pft_cube(varout, cube, outfilename,
                                                    ipft, fill_value)
        else:
            chunksizes = define_chunksizes(cube)
            print("cube should still have lazy data ",cube.has_lazy_data())
            coord_names = [coord.name() for coord in cube.coords()]
            if "latitude" in coord_names or "lat" in coord_names:
                cube.data.mask[np.isnan(cube.data)] = True    # breaks lazy data
            if not L_TESTING:
                iris.save(cube, outfilename, fill_value=fill_value, zlib=True,
                      netcdf_format='NETCDF4_CLASSIC', chunksizes=chunksizes,
                      contiguous=False, complevel=9)
                if "ISIMIP2" in MIPNAME.upper():
                    retcode = subprocess.call("mv "+outfilename+" "+\
                                                   outfilename+"4", shell=True)
                    if retcode != 0:
                        print("mv "+outfilename+" "+outfilename+"4")
                        sys.exit("mv broken")

    return outfilename, varout

# #############################################################################
def expand_to_global(cube):
    """
    define new cube with global grid for output
    """
    if "ISIMIP" in MIPNAME.upper():
        latitude, longitude, nlat, nlon = isimip_func.make_global_grid_0p5()
    elif "CMIP" in MIPNAME.upper() or "TRENDY" in MIPNAME.upper():
        latitude, longitude, nlat, nlon = cmip_func.make_global_grid_n96e()
    else:
        sys.exit("add definition of the full global grid in expand_to_global")
    cube_fullgrid = Cube(np.zeros((nlat, nlon), np.float32),
                         dim_coords_and_dims=[(latitude, 0),
                                              (longitude, 1)])
    cube = cube.regrid(cube_fullgrid,
                       iris.analysis.Nearest(extrapolation_mode="mask"))
    return cube
# #############################################################################


# #############################################################################
def add_soil_info(cube):
    """
    add soil depth information
    """
    if len(cube.coord("soil").points) == 4:  # standard jules
        cube.coord("soil").points = [0.05, 0.225, 0.675, 2.0]
        cube.coord("soil").bounds = [[0.0, 0.1], [0.1, 0.35],
                                     [0.35, 1.0], [1.0, 3.0]]
    elif len(cube.coord("soil").points) == 14:  # ejb modified jules
        dzsoil_io = (0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4,
                     1.0, 1.0, 3.0, 3.0)
        bottom_of_layer = np.cumsum(dzsoil_io)
        top_of_layer = np.array([0.0])
        top_of_layer = np.append(top_of_layer, np.cumsum(dzsoil_io)[:-1])
        bounds = np.vstack((top_of_layer, bottom_of_layer))
        bounds = np.transpose(bounds)
        cube.coord("soil").points = calc_mid_depth(bottom_of_layer)
        cube.coord("soil").bounds = bounds
    else:
        print("look in add_soil_info - check the relevant information is coded")
        sys.exit("need depth coordinates to be set properly")
    cube.coord("soil").units = "m"
    cube.coord("soil").rename("depth")
    cube.coord("depth").long_name = "depth of middle of soil layer"

    return cube
# #############################################################################


# #############################################################################
def get_jules_cube(diag_in, files_in, time_cons=None):
    """
    read and pre-process jules cube
    """
    print("load cube using jules load for "+diag_in)
    variable_cons = iris.Constraint(cube_func=(lambda c: c.var_name == diag_in))

    if "IMOGEN" in MIPNAME.upper():
        # read ensembles for imogen
        try:
            cube = imogen_func.read_ensemble(files_in, variable_cons,
                                             time_cons, diag_in)
        except:
            print(diag_in+" not available in "+files_in[0])
            raise
    else:
        # any mips which are not imogen ensembles
        try:
            cube = jules.load(files_in, variable_cons & time_cons,
                          missingdata=ma.masked, callback=add_time_middle)
        except:
            print(diag_in+" not available in "+files_in[0])
            raise
        if isinstance(cube, iris.cube.CubeList):
            for ijk in np.arange(0, int(len(cube))):
                for key in list(cube[ijk].attributes.keys()):
                    del cube[ijk].attributes[key]
                if ijk > 0:
                    cube[ijk].coord("time").convert_units \
                                (cube[0].coord("time").units)
            cube = cube.concatenate_cube()

    # sort out coordinates
    all_coord_names = [coord.name() for coord in cube.coords()]
    if "type" in all_coord_names:
        cube = add_tile_info(cube, "type")
    if "tile" in all_coord_names:
        cube = add_tile_info(cube, "tile")
    if "pft" in all_coord_names:
        cube = add_tile_info(cube, "pft")
    if "soil" in all_coord_names:
        cube = add_soil_info(cube)
    #might need these for TRENDY
    #if var == "tsl":
    #    cube.coord("soil").long_name = "stlayer"
    #if var == "msl":
    #    cube.coord("soil").long_name = "smlayer"

    return cube
# #############################################################################

# #############################################################################
def twsa_func(cubelist, var):
    """
    variation in total water mass
    (-1 * depth to water table) * water density +
    canopy water + snow + soil moisture
    should check that the input cubes are the variables expected.
    """
    #print([cube.units for cube in cubelist])
    # which element of cubelist is depth to water table? - could select by units
    watertabledepthname = "zw"
    idx = [i for i, cube in enumerate(cubelist)
           if cube.var_name == watertabledepthname]
    try:
        zw_cube = cubelist[idx[0]]
    except:
        print("name of watertabledepth is not recognised")
        raise
    water_density = 1000.0 #kg m-3
    zw_cube = zw_cube * water_density
    zw_cube.units = Unit("kg m-2")

    cubelist_minuszw = iris.cube.CubeList([])
    for cube in cubelist:
        if cube.var_name != watertabledepthname:
            cubelist_minuszw.append(cube)
    cube = sum_func(cubelist_minuszw, var)

    return cube - zw_cube
# #############################################################################

# #############################################################################
def burntarea_func(cube, var):
    """
    converts units from "fraction of land per second"
    to "% of land per month"
    """
    cube.data = cube.core_data()*30.0*86400.0*100.0
    cube.units = Unit("%")
    return cube
# #############################################################################


# #############################################################################
def fracweight_func(cubelist, var):
    """
    weight by fractional cover
    """
    # find out how many pfts
    npft = np.array([])
    for cube in cubelist:
        npft = np.append(npft, [len(coord.points) for coord in \
                                cube.coords() if coord.name() == "vegtype"])
    npft = np.min(npft)

    # which element of cubelist is landcover fraction?
    fracname = "frac"
    idx = [i for i, cube in enumerate(cubelist) if cube.var_name == fracname]
    try:
        weights = cubelist[idx[0]]
    except:
        print("name of landcover fraction is not recognised")
        raise
    weights = weights.extract(iris.Constraint(
              vegtype=lambda cell: cell.point < npft-0.5))

    cubelist_minusfrac = iris.cube.CubeList([])
    for cube in cubelist:
        if cube.var_name != fracname:
            cubelist_minusfrac.append(cube)
    cube = sum_func(cubelist_minusfrac, var)
    cube = cube.collapsed("vegtype", iris.analysis.SUM, weights=weights.data)
    return cube
# #############################################################################


# #############################################################################
def top10cm_func(cube, var):
    """"
    define mean value in top 10cm of soil
    """
    # need to put lambda functions here
    if cube.coord("depth").bounds[0, 1] == 0.1 and \
                  cube.coord("depth").points[0] == 0.05:
        cube = cube.extract(iris.Constraint(depth=lambda cell: cell == 0.05))
    else:
        sys.exit("need to sort top 10cm soil variables")
    return cube
# #############################################################################


# #############################################################################
def minus_func(cubelist, var):
    """
    subtract from first element in cubelist
    """
    out_cube = cubelist[0]
    for cube in cubelist[1:]:
        out_cube = out_cube - cube
    return out_cube
# #############################################################################


# #############################################################################
def sth_func(cubelist, var):
    """
    get water content in a layer in kg m-2 from fraction of saturation
    and saturated watercontent
    """
    soil_thick = cubelist[0].coord("depth").bounds[:, 1] - \
                      cubelist[0].coord("depth").bounds[:, 0]
    soil_thick = np.array(soil_thick)
    # print(len(soil_thick), cubelist[0].core_data().shape)
    soil_thick = np.broadcast_to(soil_thick[0], cubelist[0].core_data().shape)
    out_cube = cubelist[0] * cubelist[1] * soil_thick * 1000.0
    out_cube.units = "kg m-2"
    return out_cube
# #############################################################################


# #############################################################################
def nbp_func(cubelist, var):
    """
    should check that the input cubes are the variables expected.
    assumes first cube is npp and all others are loss terms
    """
    for i, cube in enumerate(cubelist):
        if cube.var_name in ["resp_s_to_atmos_gb", "WP_fast_out",
                             "WP_med_out", "WP_slow_out",
                             "npp_n_gb", "harvest_gb", "npp_gb",
                             "veg_c_fire_emission_gb", "burnt_carbon_dpm",
                             "burnt_carbon_rpm"]:
            if cube.units != Unit("kg m-2 s-1"):
                if "360" in str(cube.units):
                    print("nbp_func: units from "+
                          str(cube.units)+" "+cube.var_name)
                    cube.data = cube.core_data() * 1.0/(86400.0*360.0)
                    cube.units = "kg m-2 s-1"
                else:
                    sys.exit("problem with units in nbp_func")
            if cubelist[0].var_name != "npp_n_gb":
                if cubelist[0].var_name != "npp_gb":
                    sys.exit("check nbp function - cubes in wrong order")
            cubelist[i] = cube

    if len(cubelist) != 9:
        sys.exit("check nbp function - wrong number of cubes - have we added fire or not?")

    out_cube = minus_func(cubelist, var)
    out_cube.units="kg m-2 s-1"
    return out_cube
# #############################################################################


# #############################################################################
def sum_func(cubelist, var):
    """
    add cubes in cubelist
    """
    out_cube = cubelist[0]
    for cube in cubelist[1:]:
        out_cube = out_cube + cube
    return out_cube
# #############################################################################


# #############################################################################
def mult_func(cubelist, var):
    """
    multiply cubes in cubelist
    """
    out_cube = cubelist[0]
    for cube in cubelist[1:]:
        out_cube = out_cube * cube
    return out_cube
# #############################################################################


# #############################################################################
def rflow_func(cube, var):
    """
    used to convert river flow (kg/m2/s) to discharge (m3/s)
    """
    cube_area = cube.copy()
    if cube_area.coord("latitude").bounds is None:
        cube_area.coord('latitude').guess_bounds()
        cube_area.coord('longitude').guess_bounds()
    # area of grid cells
    area_weights = iris.analysis.cartography.area_weights(cube_area)
    # 1000 kg/s = 1 m3/sec
    cube = cube * 1000. * area_weights
    cube.units = "m3 s-1"
    return cube

# #############################################################################


# #############################################################################
def annmax_func(cube, var):
    """
    used e.g for annual maximum thaw depth
    """
    iris.coord_categorisation.add_year(cube, "time")
    cube = cube.aggregated_by("year", iris.analysis.MAX)
    cube.remove_coord("year")
    return cube

# #############################################################################


# #############################################################################
def cs_func(cube, var):
    """
    used for outputting cs/ns/rh into layers
    """
    all_coord_names = [ coord.name() for coord in cube.coords() ]
    if "sclayer" in all_coord_names and "scpool" in all_coord_names:
            cube = cube.collapsed("scpool", iris.analysis.SUM)
    cube.coord("sclayer").rename("soil")
    return cube

# #############################################################################


# #############################################################################
def conv360_func(cubelist, var):
    """
    used to add N fluxes with different units
    """
    cubelist_sameunits = iris.cube.CubeList([])
    for cube in cubelist:
        print(cube.units.__str__())
        if "360" in cube.units.__str__():
            cube = cube/(86400.0*360.0)
            cube.units="kg/m2/s"
        cubelist_sameunits.append(cube)
    cube = sum_func(cubelist_sameunits, var)
    return cube

# #############################################################################


# #############################################################################
def stat(cube):
    """
    print out cube statistics
    """
    print(" ~~~ stat start ~~~ ")
    print(cube)
    print(" ~~~ stat ~~~ ")
    print("min, mean, max data: "+str(np.min(cube.data))+","+
          str(np.mean(cube.data))+","+str(np.max(cube.data)))
    print("NAN, min, mean, max data: "+str(np.nanmin(cube.data))+","+
          str(np.nanmean(cube.data))+","+str(np.nanmax(cube.data)))
    print("type of cube data: "+cube.data.__class__.__name__)
    if cube.data.__class__.__name__ == "MaskedArray":
        print("number of masked points: "+str(np.sum(cube.data.mask)))
        print("number of NOT masked points: "+str(np.sum(~cube.data.mask)))
    print("shape of data: "+str(np.shape(cube.data)))
    print("NANS?: "+str(np.sum(np.isnan(cube.data)))+"   "+\
          "not NANS?: "+str(np.sum(~np.isnan(cube.data))))
    print(" ~~~ stat end ~~~ ")
# #############################################################################


# #############################################################################
def calc_mid_depth(bottom_of_layer):
    """
    given bottom_of_layer of soil layers calculate the middle of each layer
    """
    mid_depth = np.zeros(len(bottom_of_layer))
    for idim, bot_depth in enumerate(bottom_of_layer):
        if idim == 0:
            mid_depth[idim] = bot_depth/2.0
        if idim > 0:
            mid_depth[idim] = bottom_of_layer[idim-1] + \
                (bot_depth - bottom_of_layer[idim-1]) / 2.0

    return mid_depth
# #############################################################################

 #############################################################################
def define_chunksizes(cube):
    """"
    define the chunksizes for the output file
    """
    data_shape = False
    chunksizes = None
    all_coord_names = [ coord.name() for coord in cube.coords() ]
    ndims = len(cube.shape)
    ndims_orig = ndims
    if "realization" in all_coord_names and ndims_orig == len(all_coord_names):
        ndims = ndims - 1
    if ndims == 3:
        chunksizes = [1, cube.shape[-2], cube.shape[-1]]
    elif ndims == 4:
        chunksizes = [1, cube.shape[-3], cube.shape[-2], cube.shape[-1]]
    elif ndims == 5:  # hmm makes a big file
        chunksizes = [1, cube.shape[-4], cube.shape[-3], cube.shape[-2],
                      cube.shape[-1]]
    else:
        data_shape = True
        print("chunksizes are undefined set to shape of data")
    if "realization" in all_coord_names and ndims_orig == len(all_coord_names):
        if chunksizes is not None:
            chunksizes = np.append([1], chunksizes)
    if data_shape:
        if len(cube.shape) == 2:
            chunksizes = [cube.shape[0], cube.shape[1]]
    return chunksizes
# #############################################################################


# #############################################################################
if __name__ == "__main__":
    main()
