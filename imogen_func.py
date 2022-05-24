"""
imogen specific functions for postprocessing
dont need to look at if not processing imogen
"""
import sys
import glob
import iris
import iris.coords as icoords
import numpy as np
from read_input import parse_args
from read_input import config_parse_args
from read_input import read_mip_info_no_rose
sys.path.append("/home/h03/hadea/bin")
import jules

MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE = parse_args()

if L_JULES_ROSE:
    CONFIG_ARGS = config_parse_args(MIPNAME)
elif not L_JULES_ROSE:
    MIP_INFO = read_mip_info_no_rose(MIPNAME)


imodel_dict_cmip5 = {'cen_bcc_mod_bcc-csm1-1':1,
              'cen_ipsl_mod_ipsl-cm5a-lr':2,
              'cen_mri_mod_mri-cgcm3':3,
              'cen_noaa-gfdl_mod_gfdl-esm2m':4,
              'cen_bcc_mod_bcc-csm1-1-m':5,
              'cen_ipsl_mod_ipsl-cm5a-mr':6,
              'cen_nasa-giss_mod_giss-e2-h':7,
              'cen_nsf-doe-ncar_mod_cesm1-bgc':8,
              'cen_bnu_mod_bnu-esm':9,
              'cen_ipsl_mod_ipsl-cm5b-lr':10,
              'cen_nasa-giss_mod_giss-e2-h-cc':11,
              'cen_nsf-doe-ncar_mod_cesm1-cam5':12,
              'cen_cccma_mod_canesm2':13,
              'cen_miroc_mod_miroc-esm':14,
              'cen_nasa-giss_mod_giss-e2-r':15,
              'cen_nsf-doe-ncar_mod_cesm1-waccm':16,
              'cen_cmcc_mod_cmcc-cms':17,
              'cen_miroc_mod_miroc-esm-chem':18,
              'cen_nasa-giss_mod_giss-e2-r-cc':19,
              'cen_cnrm-cerfacs_mod_cnrm-cm5':20,
              'cen_miroc_mod_miroc5':21,
              'cen_ncar_mod_ccsm4':22,
              'cen_csiro-bom_mod_access1-0':23,
              'cen_mohc_mod_hadgem2-cc':24,
              'cen_ncc_mod_noresm1-m':25,
              'cen_csiro-bom_mod_access1-3':26,
              'cen_mohc_mod_hadgem2-es':27,
              'cen_ncc_mod_noresm1-me':28,
              'cen_csiro-qccce_mod_csiro-mk3-6-0':29,
              'cen_mpi-m_mod_mpi-esm-lr':30,
              'cen_noaa-gfdl_mod_gfdl-cm3':31,
              'cen_inm_mod_inmcm4':32,
              'cen_mpi-m_mod_mpi-esm-mr':33,
              'cen_noaa-gfdl_mod_gfdl-esm2g':34}

#imodel_dict_cmip5 = {'cen_cnrm-cerfacs_mod_cnrm-cm5':1,
#              'cen_csiro-bom_mod_access1-3':2,
#              'cen_ipsl_mod_ipsl-cm5a-lr':3}

imodel_dict_cmip6 = {
    "ACCESS-CM2":1,
    "ACCESS-ESM1-5":2,
    "CNRM-CM6-1-HR":3,
    "CNRM-CM6-1":4,
    "CNRM-ESM2-1":5,
    "CanESM5":6,
    "EC-Earth3-Veg":7,
    "FGOALS-g3":8,
    "HadGEM3-GC31-LL":9,
    "HadGEM3-GC31-MM":10,
    "INM-CM4-8": 11,
    "IPSL-CM6A-LR":12,
    "MIROC-ES2L":13,
    "MIROC6": 14,
    "MPI-ESM1-2-HR": 15,
    "MPI-ESM1-2-LR": 16,
    "MRI-ESM2-0": 17,
    "UKESM1-0-LL": 18 }

# #############################################################################
def read_ensemble(files_in, variable_cons, time_cons):
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
        if not cube.coords('realization'):
            realization = imodel_dict[filename.split("/")[-2]]
            ensemble_coord = icoords.AuxCoord(realization,
                                              standard_name='realization',
                                              var_name='realization')
            cube.add_aux_coord(ensemble_coord)
        for key in list(cube.attributes.keys()):
            del cube.attributes[key]
        # from add_time_middle function in process jules
        cube.coord("time").bounds = \
                3600*(np.rint(cube.coord("time").bounds.astype("int64")/3600.0))
        cube.coord("time").points = \
                (cube.coord("time").bounds[:, 0] + \
                cube.coord("time").bounds[:, 1])/2.0
        # END add our own realization coordinate if it doesn't already exist.
        # #####################################################################


    cubeall = iris.cube.CubeList([])
    imodel_dict = imodel_dict_cmip5
    try:
      if "imogen6" in MIP_INFO["model"][MIPNAME].lower() :
        imodel_dict = imodel_dict_cmip6
    except:
      if "imogen6" in CONFIG_ARGS["MODEL_INFO"]["mipname"].lower() :
        imodel_dict = imodel_dict_cmip6
    print(imodel_dict)
    keys = imodel_dict.keys()
    for key in keys:
        files_read=list()
        files_tmp = [f.replace("*", key) for f in files_in]
        files_tmp = [glob.glob(f) for f in files_tmp]
        files_read = [f for f in files_tmp if f]
        files_read = [f for sublist in files_read for f in sublist]
        #print(files_read)
        try:
            cube = jules.load(files_read, variable_cons & time_cons,
                        missingdata=np.ma.masked,
                        callback=model_ensemble_callback)
            for ijk, icube in enumerate(cube):
                if ijk >0:
                    cube[ijk].coord('time').convert_units \
                            (cube[0].coord('time').units)
            cube=cube.concatenate_cube()
            coord_names = [coord.name() for coord in cube.coords()]
            if "scpool" in coord_names:
                cube = cube.collapsed("scpool", iris.analysis.SUM)
            #print(cube)
            cubeall.append(cube)
        except:
            continue
    cubeall = cubeall.merge_cube()
    cubeall.coord('latitude').guess_bounds()
    cubeall.coord('longitude').guess_bounds()
    cubeall.coord("latitude").long_name ="latitude"
    cubeall.coord("longitude").long_name ="longitude"
    #iris.coord_categorisation.add_year(cubeall, "time")
    
    return cubeall
# #############################################################################


# #############################################################################
def make_outfilename_imogen(out_dir, outprofile, var, syr, eyr):
    """
    sort out filename for outputfile
    """
    if not L_JULES_ROSE:
        outfilename = out_dir+"/"+MIP_INFO["model"][MIPNAME]+"_"+\
                      MIP_INFO["out_scenario"][MIPNAME]+"_"+var+"_"+\
                      outprofile+"_"+str(syr)+"_"+str(eyr)+".nc"
    else:
        outfilename = out_dir+CONFIG_ARGS["OUT_INFO"]["model_out_id"]+"_"+\
                       CONFIG_ARGS["MODEL_INFO"]["configname"]+"_"+\
                       CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]+"_"+\
                       var+"_"+outprofile+"_"+str(syr)+"_"+str(eyr)+".nc"

    return outfilename
# #############################################################################


# #############################################################################
def make_infilename_imogen(src_dir, jules_profname, years):
    """
    sort out filename for input file
    """
    # print(MIP_INFO["run_name"][MIPNAME])
    if not L_JULES_ROSE:
        files_in = [src_dir+"*/*_"+MIP_INFO["in_scenario"][MIPNAME]+\
                    MIP_INFO["run_name"][MIPNAME]+\
                    jules_profname+"."+year+".nc" for year in years]
    else:
        files_in = [src_dir+"*/*_"+CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]+"_"+\
                    CONFIG_ARGS["MODEL_INFO"]["configname"]+"."+\
                    jules_profname+"."+year+".nc" for year in years]

    return files_in
