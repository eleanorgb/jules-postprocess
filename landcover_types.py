import iris 

# land cover types - are there standard names here?
UKESM_TYPES = ["BdlDcd", "BdlEvgTrop", "BdlEvgTemp", "NdlDcd", "NdlEvg",
               "c3grass", "c3crop", "c3pasture", "c4grass", "c4crop",
               "c4pasture", "shrubDcd", "shrubEvg", "urban", "lake",
               "soil", "ice"]
PERMAFROST10_TYPES = ["BdlDcd", "BdlEvgTrop", "BdlEvgTemp", "NdlDcd", "NdlEvg",
               "c3grass", "c3arctic", "c4grass", "shrubDcd", "shrubEvg", 
               "urban", "lake", "soil", "ice"]
PERMAFROST14_TYPES = ["BdlDcd", "BdlEvgTrop", "BdlEvgTemp", "NdlDcd", "NdlEvg",
               "c3grass", "c3arctic", "c3crop", "c3pasture", "c4grass",
                "c4crop", "c4pasture", "shrubDcd", "shrubEvg", "urban", "lake",
               "soil", "ice"]
HADGEM_TYPES = ["evgTree", "dcdTree", "c3", "c4", "shrub",
                "urban", "lake", "soil", "ice"]

# #############################################################################
def add_tile_info(cube, typename):
    """
    sort out tile information depending on pfts available
    """
    cube.coord(typename).long_name = "vegtype"
    #cube.coord(typename).var_name = "vegtype"
    lengthoftype = len(cube.coord("vegtype").points)
    if lengthoftype in (13, 17): # JULES-ES type
        # below broken for use with cdo
        tilecoord = iris.coords.AuxCoord(UKESM_TYPES[0:lengthoftype],
                                         long_name="frac_name")
        cube.attributes["vegtype"] = "; ".join([(str(i)+"."+x) for i, x in
                                       enumerate(UKESM_TYPES[0:lengthoftype])])
    elif lengthoftype in (5, 9): # JULES-GL7 type
        tilecoord = iris.coords.AuxCoord(HADGEM_TYPES[0:lengthoftype],
                                         long_name="frac_name")
        cube.attributes["vegtype"] = "; ".join([(str(i)+"."+x) for i, x in
                                       enumerate(HADGEM_TYPES[0:lengthoftype])])
    elif lengthoftype in (10, 14): # PERMAFROST10 NO LANDUSE TILES
        tilecoord = iris.coords.AuxCoord(PERMAFROST10_TYPES[0:lengthoftype],
                                         long_name="frac_name")
        cube.attributes["vegtype"] = "; ".join([(str(i)+"."+x) for i, x in
                                 enumerate(PERMAFROST10_TYPES[0:lengthoftype])])
    else:
        print("look in add_tile_info - make sure the relevant info is coded")
        sys.exit("coordinate "+typename+" not defined for this configuration")

    idxtile = [i for i, t in enumerate(cube.core_data().shape)
                                                           if t == lengthoftype]
    cube.add_aux_coord(tilecoord, idxtile)

    return cube
