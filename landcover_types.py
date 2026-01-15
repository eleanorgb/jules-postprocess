"""
define land cover types
"""

import sys
import iris

# land cover types - are there standard names here?
UKESM_TYPES = [
    "BdlDcd",
    "BdlEvgTrop",
    "BdlEvgTemp",
    "NdlDcd",
    "NdlEvg",
    "c3grass",
    "c3crop",
    "c3pasture",
    "c4grass",
    "c4crop",
    "c4pasture",
    "shrubDcd",
    "shrubEvg",
    "urban",
    "lake",
    "soil",
    "ice",
]
PERMAFROST10_TYPES = [
    "BdlDcd",
    "BdlEvgTrop",
    "BdlEvgTemp",
    "NdlDcd",
    "NdlEvg",
    "c3grass",
    "c3arctic",
    "c4grass",
    "shrubDcd",
    "shrubEvg",
    "urban",
    "lake",
    "soil",
    "ice",
]
PERMAFROST14_TYPES = [
    "BdlDcd",
    "BdlEvgTrop",
    "BdlEvgTemp",
    "NdlDcd",
    "NdlEvg",
    "c3grass",
    "c3arctic",
    "c3crop",
    "c3pasture",
    "c4grass",
    "c4crop",
    "c4pasture",
    "shrubDcd",
    "shrubEvg",
    "urban",
    "lake",
    "soil",
    "ice",
]
HADGEM_TYPES = [
    "evgTree",
    "dcdTree",
    "c3",
    "c4",
    "shrub",
    "urban",
    "lake",
    "soil",
    "ice",
]
SOIL_POOLS = ["dpm", "rpm", "bio", "hum"]


# #############################################################################
def add_tile_info(cube, typename):
    """
    sort out tile information depending on pfts available
    """
    errorcode = 0
    
    cube.coord(typename).long_name = "vegtype"
    # cube.coord(typename).var_name = "vegtype"
    lengthoftype = len(cube.coord("vegtype").points)
    if lengthoftype in (13, 17):  # JULES-ES type
        # below broken for use with cdo
        tilecoord = iris.coords.AuxCoord(
            UKESM_TYPES[0:lengthoftype], long_name="frac_name"
        )
        cube.attributes["vegtype"] = "; ".join(
            [(str(i) + "." + x) for i, x in enumerate(UKESM_TYPES[0:lengthoftype])]
        )
    elif lengthoftype in (5, 9):  # JULES-GL7 type
        tilecoord = iris.coords.AuxCoord(
            HADGEM_TYPES[0:lengthoftype], long_name="frac_name"
        )
        cube.attributes["vegtype"] = "; ".join(
            [(str(i) + "." + x) for i, x in enumerate(HADGEM_TYPES[0:lengthoftype])]
        )
    elif lengthoftype in (10, 14):  # PERMAFROST10 NO LANDUSE TILES
        tilecoord = iris.coords.AuxCoord(
            PERMAFROST10_TYPES[0:lengthoftype], long_name="frac_name"
        )
        cube.attributes["vegtype"] = "; ".join(
            [
                (str(i) + "." + x)
                for i, x in enumerate(PERMAFROST10_TYPES[0:lengthoftype])
            ]
        )
    else:
        print("[ERROR]: look in add_tile_info - make sure the relevant info is coded")
        print("[ERROR]: coordinate " + typename + " not defined for this configuration")
        errorcode = 1
        return cube, errorcode

    idxtile = [i for i, t in enumerate(cube.core_data().shape) if t == lengthoftype]
    cube.add_aux_coord(tilecoord, idxtile)

    return cube, errorcode


# #############################################################################
# #############################################################################
def select_vegfrac(cube, var):
    """
    make mapping of vegtypes that are defined by the output requirements
    """
    errorcode = 0

    lengthoftype = len(cube.coord("vegtype").points)
    print(f"Number of vegetation types {lengthoftype}")
    if lengthoftype == 17:  # JULES-ES type
        vegtype_mapping = {
            "BdlDcd": "BdlDcd",
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
            "grassFrac": [
                "c3grass",
                "c3pasture",
                "c3crop",
                "c4grass",
                "c4pasture",
                "c4crop",
            ],
            "treeFrac": ["BdlDcd", "BdlEvgTemp", "BdlEvgTrop", "NdlDcd", "NdlEvg"],
            "c3PftFrac": ["c3grass", "c3pasture", "c3crop"],
            "c4PftFrac": ["c4grass", "c4pasture", "c4crop"],
            "shrubFrac": ["shrubDcd", "shrubEvg"],
            "baresoilFrac": "soil",
            "residualFrac": ["urban", "lake", "ice"],
        }
    elif lengthoftype == 9:  # JULES_GL7 type
        vegtype_mapping = {
            "treeFrac": ["evgTree", "dcdTree"],
            "c3PftFrac": "c3",
            "c4PftFrac": "c4",
            "shrubFrac": "shrub",
            "baresoilFrac": "soil",
            "residualFrac": ["urban", "lake", "ice"],
        }
    elif lengthoftype == 14:  # ADDED c3arctic removed pasture and crop
        vegtype_mapping = {
            "BdlDcd": "BdlDcd",
            "treeFracBdlDcd": "BdlDcd",
            "treeFracBdlEvg": ["BdlEvgTemp", "BdlEvgTrop"],
            "treeFracNdlDcd": "NdlDcd",
            "treeFracNdlEvg": "NdlEvg",
            "NdlDcd": "NdlDcd",
            "NdlEvg": "NdlEvg",
            "BdlEvg": ["BdlEvgTemp", "BdlEvgTrop"],
            "grassFracC3": "c3grass",
            "arcticFracC3": "c3arctic",
            "grassFracC4": "c4grass",
            "treeFrac": ["BdlDcd", "BdlEvgTemp", "BdlEvgTrop", "NdlDcd", "NdlEvg"],
            "c3PftFrac": ["c3grass", "c3arctic"],
            "c4PftFrac": ["c4grass"],
            "shrubFrac": ["shrubDcd", "shrubEvg"],
            "baresoilFrac": "soil",
            "residualFrac": ["urban", "lake", "ice"],
        }
    else:
        print("ERROR: output not setup for this number of pfts")
        errorcode = 1
    try:
        cube = cube.extract(iris.Constraint(frac_name=vegtype_mapping[var]))
        if len(cube.coord("vegtype").points) > 1:
            cube = cube.collapsed("vegtype", iris.analysis.SUM)
    except:
        print(f"ERROR: vegetation type not in mapping: {var}")
        print("ERROR: add it to vegtype_mapping dictionary if you think it is needed")
        errorcode = 1
    return cube, errorcode
