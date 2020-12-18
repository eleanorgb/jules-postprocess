# needs lots of memory - run on spice
"""
script to postprocess jules output
BE CAREFUL with nitrogen on and off and npp and nlim
"""
import os
import sys
import warnings
import getpass
import subprocess
from re import sub
import numpy.ma as ma
import numpy as np
import pandas
# import matplotlib.pyplot as plt
# import iris.quickplot  as qplt
import iris
import iris.coord_categorisation
from iris.time import PartialDateTime
from iris.cube import Cube
from cf_units import Unit
#from landcover_types import UKESM_TYPES
#from landcover_types import HADGEM_TYPES
#from landcover_types import PERMAFROST10_TYPES
from landcover_types import add_tile_info
from read_mip_name import read_mip_name
# this version of JULES is currently in karinas directory
# sys.path.append("/home/h03/kwilliam/other_fcm/jules_py/trunk/jules")
sys.path.append("/home/h03/hadea/bin")
import jules


warnings.filterwarnings("ignore")

# #############################################################################
# #############################################################################
# #############################################################################
# this will change depending on the run we want to process
# see first column in mip_info.csv for currently available simulations
MIPNAME = read_mip_name()
L_TESTING = True  # change this to TRUE for testing specific diagnostics
L_BACKFILL_MISSING_FILES = True # change this to TRUE to just add missing files
# redefine diag in main() to the specific diagnistics required for testing
CONTACT_EMAIL = "eleanor.burke@metoffice.gov.uk"  # dont use my email!
# #############################################################################
# #############################################################################
# #############################################################################
OUT_BASEDIR = "/scratch/"+getpass.getuser()+"/jules_postprocess/" #output dir
INSTITUTION = "Met Office"
COMMENT = ""

print(MIPNAME)
if "ISIMIP" in MIPNAME:
    from ISIMIP_variables import get_var_dict
    import isimip_func
elif "CMIP" in MIPNAME:
    from CMIP_variables import get_var_dict
    import cmip_func
elif "TRENDY" in MIPNAME:
    from TRENDY_variables import get_var_dict
elif "IMOGEN5" in MIPNAME:
    import imogen_func
    if "variant" in MIPNAME:
        from IMOGEN_variant_variables import get_var_dict
    else:
        from IMOGEN_variables import get_var_dict
else:
    sys.exit("variable list not available for "+MIPNAME)


# #############################################################################


# #############################################################################
def main():
    """
    this is the main code reads information and loops through all diagnostics
    """
    # import dictionaries specific to mips - should be the same for all mips now
    diag_dic = get_var_dict()
    diags = get_var_dict().keys() ## all diagnostics
    print(diags)
    # ######################################################
    # select specific ones here for testing (L_TESTING=True)
    #if L_TESTING:
    #     diags=["tas","tsl"] #,"tsl_daily","snc","snc_daily"]
    # select specific ones here for testing (L_TESTING=True)
    # ######################################################

    # read dictionary of jules info for relevant mip
    # may need to add information to this file
    mip_info = pandas.read_csv("mip_info.csv", skiprows=16)
    mip_info = mip_info.set_index("MIPNAME")
    print("####################################################")
    print("INFORMATION FOR MODEL RUNS TO BE PROCESSED")
    print(mip_info.loc[MIPNAME])
    print("outputdir = "+OUT_BASEDIR)
    print("####################################################")
    mip_info = mip_info.to_dict()

    # input and output directories
    try:
        src_dir = mip_info["src_dir"][MIPNAME]
    except:
        print("ADD INFORMATION for mip to file: mip_info.csv")
        raise
    # src_dir = "/scratch/hadea/tmp/" # temporary for testing
    src_dir = src_dir+mip_info["suite_id"][MIPNAME]+"/"
    out_dir = OUT_BASEDIR+mip_info["suite_id"][MIPNAME]
    if "CMIP" in MIPNAME:
        out_dir = OUT_BASEDIR+mip_info["suite_id"][MIPNAME]+"_"+\
                   mip_info["out_scenario"][MIPNAME]
    retcode = subprocess.call("mkdir -p "+ out_dir, shell=True)
    if retcode != 0:
        print("mkdir -p "+ out_dir)
        sys.exit("mkdir broken")

    error = "\n"
    outfile_error = ""
    # loop through variables
    for var in diags:
        # determine whether we want all data for variable in one file or not
        # currently anything _daily in varname will be in multiple files
        syrall = [mip_info["start_year"][MIPNAME]]
        eyrall = [mip_info["end_year"][MIPNAME]]
        time_cons = None
        l_alltimes_in_one_file = True
        if "daily" in var:
            l_alltimes_in_one_file = False # split output files for daily data
            syrall, eyrall = define_syr_eyr_daily_output()

        for i, syr in enumerate(syrall): # loop through all years
            eyr = eyrall[i]
            print("start year: "+str(syr)+"  end year: "+str(eyr)+" "+var)
            if not l_alltimes_in_one_file:
                time_cons = iris.Constraint(time=lambda cell: \
                                            PartialDateTime(syr) \
                                            <= cell.point and \
                                            PartialDateTime(eyr) >= cell.point)

            # first call to define filenames - check file exists or not
            l_onlymakefname = True
            outfilename, varout = write_out_final_cube(mip_info, diag_dic, None,
                                                       var, out_dir, syr, eyr,
                                                       l_onlymakefname)
            if "ISIMIP" in MIPNAME:
                outfilename = isimip_func.sort_outfilename_isimip(outfilename,
                                                                  var, varout)

            # if file exists do nothing unless BACKFILL is false.
            if os.path.exists(outfilename) and L_BACKFILL_MISSING_FILES == True:
                print("File exists: "+outfilename)
                continue

            try:
                cube = make_gridded_files(src_dir, mip_info, diag_dic,
                                          time_cons, var, syr, eyr)
                l_onlymakefname = False
                outfilename, varout = write_out_final_cube(mip_info,
                                                           diag_dic, cube,
                                                           var, out_dir,
                                                           syr, eyr,
                                                           l_onlymakefname)
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
    with open("outinfo/"+MIPNAME+".txt", "w") as text_file:
        text_file.write("These files are missing:\n"+outfile_error)
# #############################################################################


# #############################################################################
def make_gridded_files(src_dir, mip_info, diag_dic, time_cons, var, syr, eyr):
    """
    make the gridded jules data
    """

    # read input diagnostics
    inputdiags = diag_dic[var][0]
    #print(inputdiags.__class__.__name__)
    jules_profname = diag_dic[var][4]
    print(var+": variable: "+str(inputdiags)+" jules_profile: "+jules_profname)
    files_in = make_infilename(mip_info, src_dir, jules_profname, syr, eyr)
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
        cube = func(cubelist)  # return a cube and then work with it!
    else:
        # one variable read from the files
        cube = get_jules_cube(inputdiags, files_in, time_cons=time_cons)
        if diag_dic[var][5] is not None:
            print("function for pre-processing "+diag_dic[var][5])
            func = globals()[diag_dic[var][5]]
            cube = func(cube)  # return a cube and then work with it!

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

    # any more attributes to add?
    cube.attributes["contact"] = CONTACT_EMAIL
    cube.attributes["institution"] = INSTITUTION
    cube.attributes["comment"] = COMMENT+"rose-suite is "+\
                                                mip_info["suite_id"][MIPNAME]

   # sort out vegetation fractions when need to separate them
    if diag_dic[var][0] == "frac" and \
                 var not in ["landCoverFrac", "pft-pft", "pft-pft_annual"]:
        cube = select_vegfrac(cube, var)

    # need to expand grids to global (not done for imogen, maybe wont bother)
    if "IMOGEN" not in MIPNAME:
        cube = expand_to_global(cube)

    return cube
# #############################################################################


# #############################################################################
def define_syr_eyr_daily_output():
    """
    define start and end year array for daily output
    """
    if "CMIP" in MIPNAME:
        syrall, syrall =   cmip_func.define_years_daily_cmip()
    elif "ISIMIP2" in MIPNAME:  #ejb change to 2b
        syrall, syrall =   isimip_func.define_years_daily_isimip2b()
    else:
        syrall = [-1]
        eyrall = [-1]
        print("need to define how to split daily data")
        sys.exit("need to define how to split daily data")
    return syrall, eyrall

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
def make_infilename(mip_info, src_dir, jules_profname, syr, eyr):
    """
    make filename of all the required infiles
    """
    # only files between start_year and end_year
    print(MIPNAME, syr, eyr)
    years = [str(year) for year in np.arange(syr, eyr+1)]
    if "IMOGEN" in MIPNAME:
        files_in = imogen_func.make_infilename_imogen(MIPNAME, mip_info,
                                                 src_dir, jules_profname, years)
    else:
        files_in = [src_dir+mip_info["run_name"][MIPNAME]+\
                    mip_info["in_scenario"][MIPNAME]+\
                    "."+jules_profname+"."+year+".nc" for year in years]
    print("First input file:")
    print(files_in[0])

    return files_in
# #############################################################################


# #############################################################################
def make_outfilename(mip_info, out_dir, outprofile, var, syr, eyr):
    """
    make filename of outfilename
    """
    print(outprofile)
    if "CMIP" in MIPNAME:
        outfilename = cmip_func.make_outfilename_cmip(mip_info,
                                             out_dir, outprofile, var, syr, eyr)
    elif "TRENDY" in MIPNAME:
        outfilename = out_dir+"/"+mip_info["model"][MIPNAME]+"_"+\
                  mip_info["out_scenario"][MIPNAME]+"_"+var+"_"+\
                  outprofile+".nc"
    elif "ISIMIP" in MIPNAME:
        outfilename = isimip_func.make_outfilename_isimip(mip_info,
                                             out_dir, outprofile, var, syr, eyr)
    elif "IMOGEN" in MIPNAME:
        outfilename = imogen_func.make_outfilename_imogen(mip_info,
                                             out_dir, outprofile, var, syr, eyr)

    else:
        sys.exit("need outfilename for "+MIPNAME)
    print("outfilename: "+outfilename)
    return outfilename
# #############################################################################


# #############################################################################
def write_out_final_cube(mip_info, diag_dic, cube, var, out_dir, syr,
                         eyr, l_onlymakefname):
    """
    write out cube
    """
    #print stats if interested
    if L_TESTING:
        if not l_onlymakefname:
            #this realises lazy data
            #stat(cube)
            print("after stat lazy data is realised ",cube.has_lazy_data())

    varout = var
    if "npp" in var and "_nlim" in var:
        varout = "npp"  #npp is different in Nitrogen/non-Nitrogen cases
        if "pft" in var:
            varout = varout+"-pft"
    elif "npp" in var:
        varout = "npp_noNlimitation"
        if "pft" in var:
            varout = varout+"-pft"
    if "ecoatmflux" in var and "_nlim" in var and "ISIMIP" in MIPNAME:
        varout = "ecoatmflux"  #npp is different in Nitrogen/non-Nitrogen cases
    elif "ecoatmflux" in var and "ISIMIP" in MIPNAME:
        varout = "ecoatmflux_noNlimitation"  #npp different in N/non-N cases
    if "nbp" in var and "_nlim" in var and "CMIP" in MIPNAME:
        varout = "nbp"  # npp different in N/non-N cases
    elif "nbp" in var and "CMIP" in MIPNAME:
        varout = "nbp_noNlimitation"   # npp different in N/non-N cases
    outprofile = "monthly"
    if "daily" in var:
        varout = sub("_daily$", "", var)
        outprofile = "daily"
    if "annual" in var:
        outprofile = "annual"
        if not "npp" in varout:
            varout = sub("_annual$", "", var)
        if "npp_nlim" in var:
            varout = "npp"
            if "pft" in var:
                varout = varout+"-pft"
        elif "npp" in var:
            varout = "npp_noNlimitation"
            if "pft" in var:
                varout = varout+"-pft"

    if "CMIP" in MIPNAME:
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
        if "ISIMIP" in MIPNAME:
            cube = isimip_func.sort_isimip_cube(cube, outprofile)

    outfilename = make_outfilename(mip_info, out_dir, outprofile,
                                   varout, syr, eyr)

    if not l_onlymakefname:
        cube.var_name = varout
        # might have to do some work ISIMIP2 output files
        # https://www.isimip.org/protocol/preparing-simulation-files/
        # quality-check-of-your-simulation-data
        if "pft" in var and "ISIMIP" in MIPNAME:
            # need to separate cube by pft
            for ipft in cube.coord("vegtype").points:
                isimip_func.sort_and_write_pft_cube(varout, cube, outfilename,
                                                    ipft, fill_value)
        else:
            chunksizes = define_chunksizes(cube)
            print("cube should still have lazy data ",cube.has_lazy_data())
            cube.data.mask[np.isnan(cube.data)] = True    # breaks lazy data
            iris.save(cube, outfilename, fill_value=fill_value, zlib=True,
                      netcdf_format='NETCDF4_CLASSIC', chunksizes=chunksizes,
                      contiguous=False, complevel=9)
            if "ISIMIP" in MIPNAME:
                retcode = subprocess.call("mv "+outfilename+" "+outfilename+\
                                          "4", shell=True)
                if retcode != 0:
                    print("mv "+outfilename+" "+outfilename+"4")
                    sys.exit("mv broken")

    return outfilename, varout

# #############################################################################
def expand_to_global(cube):
    """
    define new cube with global grid for output
    """
    if "ISIMIP" in MIPNAME:
        latitude, longitude, nlat, nlon = isimip_func.make_global_grid_0p5()
    elif "CMIP" in MIPNAME or "TRENDY" in MIPNAME:
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
def select_vegfrac(cube, var):
    """
    make mapping of vegtypes that are defined by the output requirements
    """
    lengthoftype = len(cube.coord("vegtype").points)
    print(lengthoftype)
    if lengthoftype == 17:  # JULES-ES type
        vegtype_mapping = {"BdlDcd": "BdlDcd",
                           "treeFracBdlDcd": "BdlDcd",
                           "treeFracBdlEvg": ["BdlEvgTemp", "BdlEvgTrop"],
                           "treeFracNdlDcd": "NdlDcd",
                           "treeFracNdlEvg": "NdlEvg",
                           "NdlDcd": "NdlDcd",
                           "NdlEvg": "NdlEvg",
                           "BdlEvg": ["BdlEvgTemp", "BdlEvgTrop"],
                           "grassFracC3": "c3grass",
                           "cropFracC3": "c3crop",
                           "pastureFracC3": "c3pasture",
                           "grassFracC4": "c4grass",
                           "cropFracC4": "c4crop",
                           "pastureFracC4": "c4pasture",
                           "treeFrac": ["BdlDcd", "BdlEvgTemp", "BdlEvgTrop",
                                        "NdlDcd", "NdlEvg"],
                           "c3PftFrac": ["c3grass", "c3pasture", "c3crop"],
                           "c4PftFrac": ["c4grass", "c4pasture", "c4crop"],
                           "shrubFrac": ["shrubDcd", "shrubEvg"],
                           "baresoilFrac": "soil",
                           "residualFrac": ["urban", "lake", "ice"]
                           }
    elif lengthoftype == 9:  # JULES_GL7 type
        vegtype_mapping = {"treeFrac": ["evgTree", "dcdTree"],
                           "c3PftFrac": "c3",
                           "c4PftFrac": "c4",
                           "shrubFrac": "shrub",
                           "baresoilFrac": "soil",
                           "residualFrac": ["urban", "lake", "ice"]
                           }
    else:
        sys.exit("output not setup for this number of pfts")
    try:
        print(cube.coord("frac_name"))
        print(var)
        print(vegtype_mapping[var])
        cube = cube.extract(iris.Constraint(frac_name=vegtype_mapping[var]))
        if len(cube.coord("vegtype").points) > 1:
            cube = cube.collapsed("vegtype", iris.analysis.SUM)
    except:
        print("vegetation type not recognised")
        print("need to add it to vegtype_mapping dictionary")
        raise
    return cube
# #############################################################################



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
    # print(files_in)
    if "IMOGEN" in MIPNAME:
        # read ensembles for imogen
        try:
            cube = imogen_func.read_ensemble(files_in, variable_cons, time_cons)
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
def twsa_func(cubelist):
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
    zw_cube.units = Unit("kg/m2")

    cubelist_minuszw = iris.cube.CubeList([])
    for cube in cubelist:
        if cube.var_name != watertabledepthname:
            cubelist_minuszw.append(cube)
    cube = sum_func(cubelist_minuszw)

    return cube - zw_cube
# #############################################################################

# #############################################################################
def burntarea_func(cubelist):
    """
    converts units from "fraction of land per second"
    to "% of land per year"
    """
    sys.exit("check burntarea function")
    cubelist[0].data = cubelist[0].data*365.0*86400.0*100.0
    cubelist[0].units = Unit("%")
    return cubelist[0]
# #############################################################################


# #############################################################################
def fracweight_func(cubelist):
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
    cube = sum_func(cubelist_minusfrac)
    cube = cube.collapsed("vegtype", iris.analysis.SUM, weights=weights.data)
    return cube
# #############################################################################


# #############################################################################
def top10cm_func(cube):
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
def minus_func(cubelist):
    """
    subtract from first element in cubelist
    """
    out_cube = cubelist[0]
    for cube in cubelist[1:]:
        out_cube = out_cube - cube
    return out_cube
# #############################################################################

# #############################################################################
def sth_func(cubelist):
    """
    get water content in a layer in kg/m2 from fraction of saturation
    and saturated watercontent
    """
    soil_thick = cubelist[0].coord("depth").bounds[:, 1] - \
                      cubelist[0].coord("depth").bounds[:, 0]
    soil_thick = np.array(soil_thick)
    print(len(soil_thick), cubelist[0].data.shape)
    soil_thick = np.broadcast_to(soil_thick[0], cubelist[0].data.shape)
    out_cube = cubelist[0] * cubelist[1] * soil_thick * 1000.0
    out_cube.units = "kg/m2"
    return out_cube
# #############################################################################


# #############################################################################
def nbp_func(cubelist):
    """
    should check that the input cubes are the variables expected.
    assumes first cube is npp and all others are loss terms
    """
    for i, cube in enumerate(cubelist):
        if cube.var_name in ["resp_s_to_atmos_gb", "WP_fast_out",
                             "WP_med_out", "WP_slow_out",
                             "npp_n_gb", "harvest_gb", "npp_gb"]:
            if cube.units != Unit("kg/m2/s"):
                if "360" in str(cube.units):
                    print("nbp_func: units from "+
                          str(cube.units)+" "+cube.var_name)
                    cube.data = cube.core_data() * 1.0/(86400.0*360.0)
                    cube.units = "kg/m2/s"
                else:
                    sys.exit("problem with units in nbp_func")
            if cubelist[0].var_name != "npp_n_gb":
                if cubelist[0].var_name != "npp_gb":
                    sys.exit("check nbp function - cubes in wrong order")
            cubelist[i] = cube

    if len(cubelist) != 6:
        sys.exit("check nbp function - wrong number of cubes")

    out_cube = minus_func(cubelist)
    return out_cube
# #############################################################################


# #############################################################################
def sum_func(cubelist):
    """
    add cubes in cubelist
    """
    out_cube = cubelist[0]
    for cube in cubelist[1:]:
        out_cube = out_cube + cube
    return out_cube
# #############################################################################


# #############################################################################
def annmax_func(cube):
    """
    used e.g for annual maximum thaw depth
    """
    iris.coord_categorisation.add_year(cube, "time")
    cube = cube.aggregated_by("year", iris.analysis.MAX)
    cube.remove_coord("year")
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
    all_coord_names = [ coord.name() for coord in cube.coords() ]
    ndims = len(cube.shape)
    if "realization" in all_coord_names:
        ndims = ndims - 1
    if ndims == 3:
        chunksizes = [1, cube.shape[-2], cube.shape[-1]]
    elif ndims == 4:
        chunksizes = [1, cube.shape[-3], cube.shape[-2], cube.shape[-1]]
    elif ndims == 5:  # hmm makes a big file
        chunksizes = [1, cube.shape[-4], cube.shape[-3], cube.shape[-2],
                      cube.shape[-1]]
    else:
        sys.exit("chunksizes are undefined")
    if "realization" in all_coord_names:
        chunksizes = np.append([1], chunksizes)
    return chunksizes
# #############################################################################


# #############################################################################
if __name__ == "__main__":
    main()
