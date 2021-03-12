"""
imogen specific functions for postprocessing
dont need to look at if not processing imogen
"""
import sys
import iris
import iris.coords as icoords
import numpy as np
from read_input import parse_args
from read_input import read_mip_info_no_rose
sys.path.append("/home/h03/hadea/bin")
import jules

MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE = parse_args()

if not L_JULES_ROSE:
    MIP_INFO = read_mip_info_no_rose(MIPNAME)


imodel_dict = {'cen_bcc_mod_bcc-csm1-1'.lower():1,
              'cen_ipsl_mod_ipsl-cm5a-lr'.lower():2,
              'cen_mri_mod_mri-cgcm3'.lower():3,
              'cen_noaa-gfdl_mod_gfdl-esm2m'.lower():4,
              'cen_bcc_mod_bcc-csm1-1-m'.lower():5,
              'cen_ipsl_mod_ipsl-cm5a-mr'.lower():6,
              'cen_nasa-giss_mod_giss-e2-h'.lower():7,
              'cen_nsf-doe-ncar_mod_cesm1-bgc'.lower():8,
              'cen_bnu_mod_bnu-esm'.lower():9,
              'cen_ipsl_mod_ipsl-cm5b-lr'.lower():10,
              'cen_nasa-giss_mod_giss-e2-h-cc'.lower():11,
              'cen_nsf-doe-ncar_mod_cesm1-cam5'.lower():12,
              'cen_cccma_mod_canesm2'.lower():13,
              'cen_miroc_mod_miroc-esm'.lower():14,
              'cen_nasa-giss_mod_giss-e2-r'.lower():15,
              'cen_nsf-doe-ncar_mod_cesm1-waccm'.lower():16,
              'cen_cmcc_mod_cmcc-cms'.lower():17,
              'cen_miroc_mod_miroc-esm-chem'.lower():18,
              'cen_nasa-giss_mod_giss-e2-r-cc'.lower():19,
              'cen_cnrm-cerfacs_mod_cnrm-cm5'.lower():20,
              'cen_miroc_mod_miroc5'.lower():21,
              'cen_ncar_mod_ccsm4'.lower():22,
              'cen_csiro-bom_mod_access1-0'.lower():23,
              'cen_mohc_mod_hadgem2-cc'.lower():24,
              'cen_ncc_mod_noresm1-m'.lower():25,
              'cen_csiro-bom_mod_access1-3'.lower():26,
              'cen_mohc_mod_hadgem2-es'.lower():27,
              'cen_ncc_mod_noresm1-me'.lower():28,
              'cen_csiro-qccce_mod_csiro-mk3-6-0'.lower():29,
              'cen_mpi-m_mod_mpi-esm-lr'.lower():30,
              'cen_noaa-gfdl_mod_gfdl-cm3'.lower():31,
              'cen_inm_mod_inmcm4'.lower():32,
              'cen_mpi-m_mod_mpi-esm-mr'.lower():33,
              'cen_noaa-gfdl_mod_gfdl-esm2g'.lower():34}

#imodel_dict = {'cen_cnrm-cerfacs_mod_cnrm-cm5':1,
#              'cen_csiro-bom_mod_access1-3':2,
#              'cen_ipsl_mod_ipsl-cm5a-lr':3}

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
    keys = imodel_dict.keys()
    for key in keys:
        files_read = [f.replace("*", key) for f in files_in]
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
    return files_in
