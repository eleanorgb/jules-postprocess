import iris
import numpy as np
from cf_units import Unit


# #############################################################################
# #####################################################################
def conv_360days_to_sec(cube, var):
    errorcode = 0
    def conv360(cube):
        cube.data = cube.core_data() / (86400.0 * 360.0)
        cube.units = "kg m-2 s-1"
        return cube

    if isinstance(cube, iris.cube.Cube):
        if "360" in cube.units.__str__():
            cube = conv360(cube)
    elif isinstance(cube, iris.cube.CubeList):
        for i, c in enumerate(cube):
            if "360" in c.units.__str__():
                cube[i] = conv360(c)
    return cube, errorcode


# #############################################################################
# #############################################################################
def annmax_func(cube, var):
    """
    annual maximum function
    used e.g for annual maximum thaw depth
    """
    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    if "year" not in [coord.name() for coord in cube.coords()]:
        iris.coord_categorisation.add_year(cube, "time")
    cube = cube.aggregated_by("year", iris.analysis.MAX)
    cube.remove_coord("year")
    return cube, errorcode


# #############################################################################
# #############################################################################
def burntarea_func(cube, var):
    """
    converts units from "fraction of land per second"
    to "% of land per month"
    """
    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    cube.data = cube.core_data() * 30.0 * 86400.0 * 100.0
    cube.units = Unit("%")
    cube.long_name = var
    return cube, errorcode


# #############################################################################
# #############################################################################
def burntarea_pftfunc(cube, var):
    """
    converts units from "fraction of pft per second"
    to "% of land per month"
    """
    errorcode = 0
    raise ValueError("this should be redundant")
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube.")

    cube.data = cube.core_data() * 30.0 * 86400.0 * 100.0
    cube.units = Unit("%")
    cube.long_name = var
    return cube, errorcode


# #############################################################################
# #############################################################################
def fracweight_func(cubelist, var):
    """
    weight by fractional cover
    """
    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    # find out how many pfts
    npft = np.array([])
    for cube in cubelist:
        npft = np.append(
            npft,
            [len(coord.points) for coord in cube.coords() if coord.name() == "vegtype"],
        )
    npft = int(np.min(npft))

    # which element of cubelist is landcover fraction?
    fracname = "frac"
    idx = [i for i, cube in enumerate(cubelist) if cube.var_name == fracname]
    try:
        weights = cubelist[idx[0]]
    except:
        raise ValueError("name of landcover fraction is not recognised")
    weights = weights.extract(
        iris.Constraint(vegtype=lambda cell: cell.point < npft + 0.5)
    )

    cubelist_withoutfrac = iris.cube.CubeList([])
    for cube in cubelist:
        if cube.var_name != fracname:
            cubelist_withoutfrac.append(cube)
    if len(cubelist_withoutfrac) != 0:
        cube, errorcode = sum_func(cubelist_withoutfrac, var)
    else:
        errorcode = 1
        print("ERROR: frac weight function does not have enough variabes")
        print(cubelist)
        return None, errorcode

    cube = cube.collapsed("vegtype", iris.analysis.SUM, weights=weights.core_data())
    return cube, errorcode


# #############################################################################
# #############################################################################
def minus_func(cubelist, var):
    """
    subtract from first element in cubelist
    """

    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    cubelist = conv_360days_to_sec(cubelist, var)

    out_cube = cubelist[0]
    for cube in cubelist[1:]:
        out_cube = out_cube - cube
    out_cube.long_name = var
    return out_cube, errorcode


# #############################################################################
# #############################################################################
def mult_func(cubelist, var):
    """
    multiply all cubes in cubelist together
    """

    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    cubelist = conv_360days_to_sec(cubelist, var)

    out_cube = cubelist[0]
    for cube in cubelist[1:]:
        out_cube = out_cube * cube
    return out_cube, errorcode


# #############################################################################


# #############################################################################
def nbp_func(cubelist, var, fire=True):
    """
    should check that the input cubes are the variables expected.
    first cube is npp and all others are loss terms
    """

    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    npp_found = False
    cubelist_used = iris.cube.CubeList([])

    # make sure npp is first
    for i, cube in enumerate(cubelist):
        if cube.var_name in ["npp_n_gb", "npp_gb"]:
            if npp_found == False:
                cubelist_used.append(cube)
                npp_found = True
            else:
                raise ValueError("check nbp function - 2 variables with npp in cube")

    for i, cube in enumerate(cubelist):
        if "sclayer" in [coord.name() for coord in cube.coords()]:
            cube = cube.collapsed("sclayer", iris.analysis.SUM)

        # check right variables
        if cube.var_name in [
            "resp_s_to_atmos_gb",
            "WP_fast_out",
            "WP_med_out",
            "WP_slow_out",
            "harvest_gb",
            "veg_c_fire_emission_gb",
            "burnt_carbon_dpm",
            "burnt_carbon_rpm",
        ]:
            cubelist_used.append(cube)

    if len(cubelist_used) != 9 and fire:
        raise ValueError(
            "check nbp function - wrong number of cubes - have we added fire or not?"
        )

    if cubelist_used[0].var_name != "npp_n_gb":
        if cubelist_used[0].var_name != "npp_gb":
            raise ValueError("check nbp function - cubes in wrong order")

    out_cube, errorcode = minus_func(cubelist_used, var)
    out_cube.units = "kg m-2 s-1"
    return out_cube, errorcode


# #############################################################################
# #############################################################################
def nee_func(cubelist, var):
    """
    should check that the input cubes are the variables expected.
    assumes first cube is npp and all others are loss terms
    """

    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    if cubelist[0].var_name != "npp_n_gb":
        if cubelist[0].var_name != "npp_gb":
            raise ValueError("check nee function - cubes in wrong order")

    cubelist_used = iris.cube.CubeList([])
    for i, cube in enumerate(cubelist):
        if "sclayer" in [coord.name() for coord in cube.coords()]:
            cube = cube.collapsed("sclayer", iris.analysis.SUM)

        if cube.var_name in ["resp_s_to_atmos_gb", "npp_n_gb"]:
            cubelist_used.append(cube)

    if len(cubelist_used) != 2:
        raise ValueError("check nee function - wrong number of cubes")

    out_cube, errorcode = minus_func(cubelist, var)

    return out_cube, errorcode


# #############################################################################


# #############################################################################
def sth_func(cubelist, var):
    """
    get water content in a layer in kg m-2 from fraction of saturation
    and saturated watercontent
    """
    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    if len(cubelist) != 2:
        raise ValueError("cubelist must contain exactly two elements.")

    soil_thick = (
        cubelist[0].coord("depth").bounds[:, 1]
        - cubelist[0].coord("depth").bounds[:, 0]
    )
    soil_thick = np.array(soil_thick)
    # print(len(soil_thick), cubelist[0].core_data().shape)
    # print(cubelist[0])
    # print(cubelist[1])
    # raise
    soil_thick = np.broadcast_to(soil_thick[0], cubelist[0].core_data().shape)
    out_cube = cubelist[0] * cubelist[1] * soil_thick * 1000.0
    out_cube.units = "kg m-2"
    return out_cube, errorcode


# #############################################################################


# #############################################################################
def div_func(cubelist, var):
    """
    divide cubes in cubelist, first check there are only 2 cubes present
    """

    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    cubelist = conv_360days_to_sec(cubelist, var)

    if len(cubelist) != 2:
        raise ValueError("cubelist must contain exactly two elements.")
    else:
        out_cube = cubelist[0] / cubelist[1]

    return out_cube, errorcode


# #############################################################################


# #############################################################################
def top10cm_func(cube, var):
    """ "
    define mean value in top 10cm of soil only one special case here
    """
    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    # need to put lambda functions here
    try:
        if (
            cube.coord("depth").bounds[0, 1] == 0.1
            and cube.coord("depth").points[0] == 0.05
        ):
            cube = cube.extract(iris.Constraint(depth=lambda cell: cell == 0.05))
        else:
            raise ValueError("Need to sort top 10cm soil variables")
    except iris.exceptions.CoordinateNotFoundError:
        raise ValueError("The 'depth' coordinate was not found in the cube.")
    return cube, errorcode


# #############################################################################


# #############################################################################
def rflow_func(cube, var):
    """
    used to convert river flow (kg/m2/s) to discharge (m3/s)
    for a specific point need to multiply by area of grid cell
    """
    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    cube_area = cube.copy()
    if cube_area.coord("latitude").bounds is None:
        cube_area.coord("latitude").guess_bounds()
        cube_area.coord("longitude").guess_bounds()
    # area of grid cells
    area_weights = iris.analysis.cartography.area_weights(cube_area)
    # 1000 kg/s = 1 m3/sec
    cube = cube / 1000.0 * area_weights
    cube.units = "m3 s-1"
    return cube, errorcode


# #############################################################################


# #############################################################################
def layered_soilbgc_func(cube, var):
    """
    used for outputting layered cs/ns/rh but adding pools together
    """
    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    all_coord_names = [coord.name() for coord in cube.coords()]
    if "sclayer" in all_coord_names and "scpool" in all_coord_names:
        cube = cube.collapsed("scpool", iris.analysis.SUM)
    # raise ValueError("Need to sort out layered soilbgc function")
    # cube.coord("sclayer").rename("depth")
    # add depth coordinate
    return cube, errorcode


# #############################################################################


# #############################################################################
def soilbgc_pool_func(cube, var):
    """
    used for outputting cs/ns as pools but removing layers
    """
    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    if "sclayer" in [coord.name() for coord in cube.coords()]:
        cube = cube.collapsed("sclayer", iris.analysis.SUM)
    return cube, errorcode


# #############################################################################
# #############################################################################
def rhums_func(cubelist, var):
    """
    calculate relative humidity from specific humidity
    inputs are 1.5m q, 1.5m T and p*
    http://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    """

    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    if len(cubelist) != 3:
        raise ValueError("cubelist must contain exactly three elements.")

    q1p5m_constraint = iris.Constraint(cube_func=lambda cube: "q1p5m" in cube.var_name)
    q1p5m = cubelist.extract_cube(q1p5m_constraint)
    if not q1p5m:
        raise ValueError("Cube 'q1p5m' not found in cubelist")
    # q1p5m = cubelist[0]  # specific humidity
    # if "q1p5m" not in q1p5m.var_name:
    #     raise ValueError(
    #         f"name of specific humidity is not recognised {q1p5m.var_name}"
    #     )
    t1p5m_constraint = iris.Constraint(cube_func=lambda cube: "t1p5m" in cube.var_name)
    t1p5m = cubelist.extract_cube(t1p5m_constraint)
    if not t1p5m:
        raise ValueError("Cube 't1p5m' not found in cubelist")
    # t1p5m = cubelist[1]  # air temperature
    # if "t1p5m" not in t1p5m.var_name:
    #     raise ValueError(f"name of air temperature is not recognised {t1p5m.var_name}")
    t1p5m.convert_units("celsius")

    pstar_constraint = iris.Constraint(cube_func=lambda cube: "pstar" in cube.var_name)
    pstar = cubelist.extract_cube(pstar_constraint)
    if not pstar:
        raise ValueError("Cube 'pstar' not found in cubelist")
    #    pstar = cubelist[2]  # surface pressure
    #    if "pstar" not in pstar.var_name:
    #        raise ValueError(f"name of air temperature is not recognised {pstar.var_name}")
    pstar.convert_units("mbar")

    # convert q1p5m to vapour pressure in millibars
    vp = q1p5m * pstar / (0.378 * q1p5m + 0.622)
    # saturated vapor pressure = svp
    svp = iris.analysis.maths.exp((17.26938818 * t1p5m) / (237.3 + t1p5m)) * 6.1078
    # calculate relative humidity
    hurs = vp / svp * 100.0
    hurs.core_data()[hurs.core_data() > 100.0] = 100.0
    hurs.core_data()[hurs.core_data() < 0.0] = 0.0
    hurs.units = Unit("%")
    hurs.var_name = var

    return hurs, errorcode


# #############################################################################
# #############################################################################
def get_vegtype_frac(cube, var):

    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    if "vegtype" not in [coord.name() for coord in cube.coords()]:
        cube = add_tile_info(cube, "type")
    cube = select_vegfrac(cube, var)
    return cube, errorcode


# #############################################################################
# #############################################################################
def sum_func(cubelist, var, collapse_sclayer=True):
    """
    add cubes in cubelist
    """
    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    if collapse_sclayer:
        for i, cube in enumerate(cubelist):
            if "sclayer" in [coord.name() for coord in cube.coords()]:
                cubelist[i] = cube.collapsed("sclayer", iris.analysis.SUM)

    cubelist = conv_360days_to_sec(cubelist, var)

    out_cube = cubelist[0]
    for cube in cubelist[1:]:
        out_cube = out_cube + cube
    out_cube.long_name = var
    return out_cube, errorcode


# #############################################################################


# ############################################################################
# #############################################################################
def netatmco2flux_func(cubelist, var):
    """
    should check that the input cubes are the variables expected.
    assumes first cube is npp and all others are loss terms
    """
    errorcode = 0
    if not isinstance(cubelist, iris.cube.CubeList):
        raise ValueError("cubelist must be an instance of iris.cube.CubeList")

    for i, cube in enumerate(cubelist):
        if cubelist[0].var_name != "npp_n_gb":
            if cubelist[0].var_name != "npp_gb":
                sys.exit("check nbp function - cubes in wrong order")
        cubelist[i] = cube

    if len(cubelist) != 9:
        sys.exit(
            "check nbp function - wrong number of cubes - have we added fire or not?"
        )

    out_cube = minus_func(cubelist, var)
    out_cube.units = "kg m-2 s-1"
    out_cube.long_name = var
    return out_cube, errorcode


# #############################################################################
# #############################################################################
def pf_annmaxthaw_func(cube, var):
    """
    used to caluclate permafrost extent from annual maximum thawdepth
    """
    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    if "year" not in [coord.name() for coord in cube.coords()]:
        iris.coord_categorisation.add_year(cube, "time")
    # qplt.pcolormesh(cube[0])
    # plt.show()
    cube = cube.aggregated_by("year", iris.analysis.MAX)
    cube.data.mask[cube.data > MAX_THAW_DEPTH] = True
    mask = cube.data.mask.copy()
    cube.data[~cube.data.mask] = 1.0
    cube.data[cube.data.mask] = 0.0
    cube.data.mask = mask.copy()
    cube = iris.util.squeeze(cube)
    # print(cube)
    # qplt.pcolormesh(cube)
    # plt.show()
    cube.remove_coord("year")
    cube.long_name = var
    cube.units = Unit("1e6 km2")
    return cube, errorcode


# #############################################################################
# #############################################################################
def pf_deepsoilt_func(cube, var):
    """
    used to caluclate permafrost extent from depsoil temperture
    """
    errorcode = 0
    if not isinstance(cube, iris.cube.Cube):
        raise ValueError("cube must be an instance of iris.cube.Cube")

    if "year" not in [coord.name() for coord in cube.coords()]:
        iris.coord_categorisation.add_year(cube, "time")
    cube = cube.extract(
        iris.Constraint(soil=lambda s: s == max(cube.coord("soil").points))
    )
    cube = cube.aggregated_by("year", iris.analysis.MAX)
    cube = iris.util.squeeze(cube)
    cube.data.mask[cube.data > 273.2] = True  # remove magic number
    # print(cube)
    # qplt.pcolormesh(cube)
    # plt.show()
    mask = cube.data.mask.copy()
    cube.data[~cube.data.mask] = 1.0
    cube.data[cube.data.mask] = 0.0
    cube.data.mask = mask.copy()
    cube.remove_coord("year")
    cube.long_name = var
    cube.units = Unit("1e6 km2")
    return cube, errorcode
