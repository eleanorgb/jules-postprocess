"""
isimip specific functions for postprocessing
dont need to look at if not processing isimip
"""

import sys
import subprocess
import glob
import numpy as np
import iris
import netCDF4
import os
import dask.array as da
from iris.coords import DimCoord
from cf_units import Unit
from landcover_types import UKESM_TYPES
from landcover_types import SOIL_POOLS
from read_input import parse_args
from read_input import config_parse_args
from read_input import read_mip_info_no_rose

MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE, JSONMAPFILE = parse_args()

if L_JULES_ROSE:
    CONFIG_ARGS = config_parse_args(MIPNAME)
elif not L_JULES_ROSE:
    MIP_INFO = read_mip_info_no_rose(MIPNAME)


# #############################################################################
# #############################################################################
def make_global_grid_0p5():
    """
    0p5 degrees global grid
    """
    nlat = 360
    latitude = DimCoord(
        np.linspace(-89.75, 89.75, nlat), standard_name="latitude", units="degrees"
    )
    nlon = 720
    longitude = DimCoord(
        np.linspace(-179.75, 179.75, nlon), standard_name="longitude", units="degrees"
    )
    return latitude, longitude, nlat, nlon


# #############################################################################
# #############################################################################
def rename_cfcompliant_to_isimip(outfilename, cube, divide_files=False):
    """
    annoyingly iris save changes some variable and dimension names
    -total is automatically changed to _total
    """
    if divide_files:
        try:
            filesall = glob.glob(outfilename[:-3] + "_separate*.nc")
        # subprocess.run("conda activate cdo; cdo cat "+' '.join(filesall)+" "+outfilename+"; conda deactivate", shell=True)
        except:
            pass

    # changes from _total to -total
    if "-total" in outfilename:
        var_old = cube.var_name.split("-")
        if len(var_old) == 2:
            var_old = "".join([var_old[0], "_", var_old[1]])
            dataset = netCDF4.Dataset(outfilename, "r+")
            try:
                dataset.renameVariable(var_old, cube.var_name)
                dataset[var_old].long_name = cube.var_name
            except:
                print("varname already cf compliant when looking at _/-")
            dataset.close()
    return


# #############################################################################
# #############################################################################
def define_years_daily_isimip2b():
    """
    break into 10 year batches
    """
    if "hist" in MIPNAME:
        syrall = np.arange(1861, 2011, 10)
        eyrall = syrall + 9
        eyrall[eyrall == 2010] = 2005
    if "rcp" in MIPNAME:
        syrall = np.arange(2001, 2101, 10)
        eyrall = syrall + 9
        syrall[syrall == 2001] = 2006
    print(MIPNAME, ": ", syrall)
    print(MIPNAME, ": ", eyrall)
    return syrall, eyrall


# #############################################################################
# #############################################################################
def define_years_daily_wrpmip():
    """
    break into 1 year batches
    """
    syrall = np.arange(2000, 2023)
    eyrall = syrall
    return syrall, eyrall


# #############################################################################
# #############################################################################
def make_outfilename_isimip(out_dir, outprofile, var, syr, eyr):
    """
    define outfilename for isimip2b:
    <modelname>_<gcm>_<bias-correction>_<climate-scenario>_
    <soc-scenario>_<co2sens-scenarios>_<variable>_<region>_
    <timestep>_<start-year>_<end-year>.nc4

    define outfilename for isimip3a: appendix nc not nc4
    <model>_<climate-forcing>_<climate-scenario>_
    <soc-scenario>_<sens-scenario>_<variable>(-<crop>-<irrigation>|-<pft>)
    _<region>_<time-step>_<start-year>_<end-year>.nc

    define outfilename for isimip3b: appendix nc not nc4
    <model>_<climate-forcing>_<bias-adjustment>_<climate-scenario>
    _<soc-scenario>_<sens-scenario>_<variable>(-<crop>-<irrigation>|-<pft>)
    _<region>_<time-step>_<start-year>_<end-year>.nc
    """
    if not L_JULES_ROSE:
        # print("scenario: ",MIP_INFO["in_scenario"][MIPNAME])
        if MIP_INFO["in_scenario"][MIPNAME] in "obsclim":
            soc = "histsoc_co2"
        elif MIP_INFO["in_scenario"][MIPNAME] in "c20c":
            soc = "histsoc_co2"
        elif MIP_INFO["in_scenario"][MIPNAME] in "historical":
            soc = "default"
        elif MIP_INFO["in_scenario"][MIPNAME] in "rcp2p6":
            soc = "rcp26soc_co2"
        elif MIP_INFO["in_scenario"][MIPNAME] in "rcp6p0":
            soc = "rcp60soc_co2"
        elif MIP_INFO["in_scenario"][MIPNAME] in "rcp8p5":
            soc = "nosoc_co2"
        # these below are only for 2015 land use with transient historical land use
        # plus changing co2
        elif MIP_INFO["in_scenario"][MIPNAME] in "ssp126":
            soc = "2015soc-from-histsoc_default"
        elif MIP_INFO["in_scenario"][MIPNAME] in "ssp370":
            soc = "2015soc-from-histsoc_default"
        elif MIP_INFO["in_scenario"][MIPNAME] in "ssp585":
            soc = "2015soc-from-histsoc_default"
        # these above are only for 2015 land use with transient historical land use
        # plus changing co2
        else:
            sys.exit("no soc")
        if "isimip2" in MIPNAME:
            add_drive_info = "ewembi_"
        elif "isimip3b" in MIPNAME:
            add_drive_info = "w5e5_"
        else:
            add_drive_info = ""
        outfilename = (
            out_dir
            + "/"
            + MIP_INFO["model"][MIPNAME].lower()
            + "_"
            + MIP_INFO["run_name"][MIPNAME].lower()
            + add_drive_info
            + MIP_INFO["out_scenario"][MIPNAME]
            + "_"
            + soc
            + "_"
            + var
            + "_global_"
            + outprofile
            + "_"
            + str(syr)
            + "_"
            + str(eyr)
            + ".nc"
        )

    if L_JULES_ROSE:
        outfilename = (
            out_dir
            + CONFIG_ARGS["OUT_INFO"]["model_out_id"].lower()
            + "_"
            + CONFIG_ARGS["MODEL_INFO"]["driving_model"].lower()
            + CONFIG_ARGS["MODEL_INFO"]["bias_correction"].lower()
            + CONFIG_ARGS["MODEL_INFO"]["climate_scenario"]
            + "_"
            + CONFIG_ARGS["MODEL_INFO"]["soc_scenario"]
            + "_"
            + CONFIG_ARGS["MODEL_INFO"]["co2_scenario"]
            + "_"
            + var
            + "_global_"
            + outprofile
            + "_"
            + str(syr)
            + "_"
            + str(eyr)
            + ".nc"
        )
    return outfilename


# #############################################################################
# #############################################################################
def sort_outfilename_isimip(outfilename, var, varout):
    """
    JULES-ES assumed here for vegetation types
    do we need to return varout?
    """
    if "isimip2" in MIPNAME:
        outfilename = outfilename + "4"
    if "-pool" in var:  # only checking filename containing first soil pool
        var_minus_pool = varout.replace("-pool", "-")
        print(varout, ": only checking file exists for first soil pool")
        outfilename = outfilename.replace(
            varout, var_minus_pool + SOIL_POOLS[0].lower()
        )

    if "pft" in var:  # only checking filename containing first pft
        var_minus_pft = varout.replace("-pft", "-")
        print(varout, ": only checking file exists for first pft")
        outfilename = outfilename.replace(
            varout, var_minus_pft + UKESM_TYPES[0].lower()
        )
    if "npp" in var and "_nlim" in var:
        outfilename = outfilename.replace(varout, "npp")
        if "total" in var:
            outfilename = outfilename.replace("npp", "npp-total")
        # npp is different in Nitrogen/non-Nitrogen cases
    elif "npp" in var:
        outfilename = outfilename.replace(varout, "npp_nonlimitation")
        if "total" in var:
            outfilename = outfilename.replace(
                "npp_nonlimitation", "npp_nonlimitation-total"
            )
        # npp is different in Nitrogen/non-Nitrogen cases
    return outfilename


# #############################################################################
# #############################################################################
def sort_isimip_cube(cube, outprofile):
    """
    make isimip data past validity tests
    isimip reference dates:
    ISIMIP3a 	1901-01-01, 00:00:00, 'proleptic_gregorian'
    ISIMIP3b 	1601-01-01, 00:00:00, 'proleptic_gregorian'
    """
    tcoord = cube.coord("time")
    tcoord.units = Unit(tcoord.units.origin, calendar="gregorian")

    ref_year = 1901
    if "isimip2b" in MIPNAME.lower():
        ref_year = 1661
    elif "isimip3a" in MIPNAME.lower():
        ref_year = 1901
    elif "isimip3b" in MIPNAME.lower():
        ref_year = 1601
    tcoord.convert_units("days since " + str(ref_year) + "-01-01 00:00:00")
    tcoord.units = Unit(tcoord.units.origin, calendar="proleptic_gregorian")
    cube.remove_coord("time")
    cube.add_dim_coord(tcoord, 0)  # might need to find this dimension
    cube.coord("time").long_name = "Time"
    cube.coord("time").bounds = None
    cube.coord("latitude").long_name = "latitude"
    cube.coord("longitude").long_name = "longitude"
    cube.coord("latitude").var_name = "lat"
    cube.coord("longitude").var_name = "lon"
    cube.coord("latitude").points = np.flip(cube.coord("latitude").points, axis=0)
    # which coordinate is latitude

    for coord in cube.coords():
        if coord.name() == "vegtype":
            coord.bounds = None
    for ilat, coord in enumerate(cube.coords()):
        if coord.name() == "latitude":
            lat_coord_idx = ilat
    cube.data = np.flip(cube.core_data(), axis=lat_coord_idx)
    cube.data = cube.core_data().astype("float32")
    cube.attributes["missing_value"] = np.float32(1.0e20)
    for key in list(cube.attributes.keys()):
        if key == "coordinates":
            del cube.attributes[key]
    cube.cell_methods = None
    cube = iris.util.squeeze(cube)

    # depth coordinate
    for coord in cube.coords():
        if coord.name() == "depth":
            coord.long_name = "Depth of vertical layer center below land"
        if coord.name() == "vegtype":
            if len(cube.coord("vegtype").points) == 1:
                cube.remove_coord("vegtype")
    if cube.units == "kg/m2":
        cube.units = "kg m-2"
    if cube.units == "kg/m2/s":
        cube.units = "kg m-2 s-1"
    if cube.units == "m-2.kg.s-1":
        cube.units = "kg m-2 s-1"

    return cube


# #############################################################################
# #############################################################################
def sort_and_write_pft_cube(varout, cube, outfilename, ipft, fill_value):
    """
    for each pft sort out the cube
    """
    cubeout = cube.extract(iris.Constraint(vegtype=ipft))
    var_minus_pft = varout.replace("-pft", "-")
    cubeout.var_name = var_minus_pft + cubeout.coord("frac_name").points[0]
    cubeout.var_name = cubeout.var_name.lower()
    outfilename = outfilename.replace(varout, cubeout.var_name.lower())
    # print("outfilename: " + os.path.basename(outfilename))
    wrong_name = varout.replace("-pft", "_") + cubeout.coord("frac_name").points[0]
    wrong_name = wrong_name.lower()

    # remove some stuff
    for coord in cubeout.coords():
        if coord.name() == "frac_name":
            cubeout.remove_coord("frac_name")
        if coord.name() == "vegtype":
            if len(cubeout.coord("vegtype").points) == 1:
                cubeout.remove_coord("vegtype")
    if "vegtype" in list(cubeout.attributes.keys()):
        del cubeout.attributes["vegtype"]
    # end remove some stuff

    # print("cube should still have lazy data ", cubeout.has_lazy_data())
    cubeout.data.mask[np.isnan(cubeout.data)] = True  # breaks lazy data
    # cubeout.data = da.where(da.isnan(cubeout.data), True, cubeout.data)  # think fixes above issue
    chunksizes = [1, cubeout.shape[1], cubeout.shape[2]]
    if not L_TESTING:
        iris.save(
            cubeout,
            outfilename,
            fill_value=fill_value,
            zlib=True,
            netcdf_format="NETCDF4_CLASSIC",
            chunksizes=chunksizes,
            contiguous=False,
            complevel=9,
        )
        # dont really understand why this below happens
        retcode = subprocess.call(
            "ncrename -h -v " + wrong_name + "," + cubeout.var_name + " " + outfilename,
            shell=True,
        )
        if retcode != 0:
            raise ValueError("ERROR: ncrename variable broken " + outfilename)

        if "ISIMIP2B" in MIPNAME.upper():
            retcode = subprocess.call(
                "mv " + outfilename + " " + outfilename + "4", shell=True
            )
            if retcode != 0:
                raise ValueError("ERROR: mv " + outfilename + " " + outfilename + "4")


# #############################################################################
# #############################################################################
def sort_and_write_pool_cube(varout, cube, outfilename, ipool, fill_value):
    """
    for each pool sort out the cube
    """
    cubeout = cube.extract(iris.Constraint(scpool=ipool))
    var_minus_pool = varout.replace("-pool", "-")
    cubeout.var_name = var_minus_pool + SOIL_POOLS[ipool - 1]
    cubeout.var_name = cubeout.var_name.lower()
    outfilename = outfilename.replace(varout, cubeout.var_name.lower())
    # print("outfilename: " + os.path.basename(outfilename))
    wrong_name = varout.replace("-pool", "_") + SOIL_POOLS[ipool - 1]
    wrong_name = wrong_name.lower()

    # remove some stuff
    for coord in cubeout.coords():
        if coord.name() == "scpool":
            cubeout.remove_coord("scpool")
        if coord.name() == "sclayer":
            if len(cubeout.coord("sclayer").points) > 1:
                cubeout = cubeout.collapsed("sclayer", iris.analysis.SUM)
            cubeout.remove_coord("sclayer")
    # end remove some stuff

    print("cube should still have lazy data ", cubeout.has_lazy_data())
    cubeout.data.mask[np.isnan(cubeout.data)] = True  # breaks lazy data
    chunksizes = [1, cubeout.shape[1], cubeout.shape[2]]
    if not L_TESTING:
        iris.save(
            cubeout,
            outfilename,
            fill_value=fill_value,
            zlib=True,
            netcdf_format="NETCDF4_CLASSIC",
            chunksizes=chunksizes,
            contiguous=False,
            complevel=9,
        )
        # dont really understand why this below happens
        retcode = subprocess.call(
            "ncrename -h -v " + wrong_name + "," + cubeout.var_name + " " + outfilename,
            shell=True,
        )
        if retcode != 0:
            raise ValueError("ERROR: ncrename variable broken " + outfilename)
        if "isimip2b" in MIPNAME:
            retcode = subprocess.call(
                "mv " + outfilename + " " + outfilename + "4", shell=True
            )
            if retcode != 0:
                raise ValueError("ERROR: mv " + outfilename + " " + outfilename + "4")
