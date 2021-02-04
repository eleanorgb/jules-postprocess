# #############################################################################
# no daily data here
def get_var_dict():
#               "outVar":[0.inputVars, 1.longName, 2.scale, 3.units, 4.outprofile, 5.func, 6.CMIPprofileName]
    var_dict = {
        # ilamb - MONTHLY
        "cSoil":         ["cs_gb", "carbon mass in model soil pool", None, "kg/m2", "ilamb", None, None],
        "cVeg":          ["cv", "carbon mass in vegetation", None, "kg/m2", "ilamb", None, None],
        "hfss":          ["ftl_gb", "sensible heat flux", None, "W/m2", "ilamb", None, None],
        "gpp":           ["gpp_gb", "gross primary production", None, "kg/m2/s", "ilamb", None, None],
        "fHarvest":      ["harvest_gb", "Carbon Mass Flux into Atmosphere Due to Crop Harvesting", 1.0/(86400.0*360.0), 'kg m-2 s-1', "ilamb", None, None],
        "lai":           ["lai_gb", "leaf area index", None, "1", "ilamb", None, None],
        "hfls":          ["latent_heat", "latent heat flux", None, "W/m2", "ilamb", None, None],
        "rlds":          ["lw_down", "downward longwave radiation", None, "W/m2", "ilamb", None, None],
        "rlus":          ["lw_up", "upward longwave radiation", None, "W/m2", "ilamb", None, None],
        "npp_nlim":      ["npp_n_gb", "Net Primary Production on Land as Carbon Mass Flux", 1.0/(86400.0*360.0), "kg m-2 s-1", "ilamb", None, None],
        "npp":           ["npp_gb", "Net Primary Production on Land as Carbon Mass Flux", None, "kg m-2 s-1", "ilamb", None, None],
        "pr":            ["precip", "precipitation rate", None, "kg/m2/s", "ilamb", None, None],
        "ps":            ["pstar", "surface air pressure", None, "Pa", "ilamb", None, None],
        "ra":            ["resp_p_gb", "autotrophic respiration", None, "kg/m2/s", "ilamb", None, None],
        "rh":            ["resp_s_to_atmos_gb", "Total Heterotrophic Respiration on Land as Carbon Mass Flux", 1.0/(86400.0*360.0), "kg m-2 s-1", "ilamb", None, None],# is this right?
        "mrso":          ["smc_tot", "total soil moisture", None, "kg/m2", "ilamb", None, None],
        "snw":           ["snow_mass_gb", "snow water equivalent", None, "kg/m2", "ilamb", None, None],
        "rsds":          ["sw_down", "downward shortwave radiation", None, "W/m2", "ilamb", None, None],
        "tas":           ["t1p5m_gb", "1.5m air temeprature", None, "K", "ilamb", None, None],
        "tsl":           ["t_soil", "Temperature of soil - layer", None, "K", "ilamb", None, None],
        "mrro":          ["runoff", "total runoff", None, "kg/m2/s", "ilamb", None, None],
        "wtd":           ["zw", "water table depth", None, "m", "ilamb", None, None],
        "mrtws":         [["zw", "canopy_gb", "snow_mass_gb", "smc_tot"], "Total water storage", None, "kg/m2", "ilamb", "twsa_func", None],
        "nbp_nlim":      [["npp_n_gb", "resp_s_to_atmos_gb", "harvest_gb", "WP_fast_out", "WP_med_out", "WP_slow_out"], 'Carbon Mass Flux out of Atmosphere Due to Net Biospheric Production on Land', None, "kg m-2 s-1", "ilamb", "nbp_func", None],# downwards positive
        "nbp":           [["npp_gb", "resp_s_to_atmos_gb", "harvest_gb", "WP_fast_out", "WP_med_out", "WP_slow_out"], 'Carbon Mass Flux out of Atmosphere Due to Net Biospheric Production on Land', None, "kg m-2 s-1", "ilamb", "nbp_func", None],# downwards positive
        "rsus":          [["sw_down", "sw_net"], "upward shortwave radiation", None, "W/m2", "ilamb", "minus_func", None],
        "fLuc":          [["WP_fast_out", "WP_med_out", "WP_slow_out"], 'Net Carbon Mass Flux into Atmosphere Due to Land-Use Change', 1.0/(86400.0*360.0), "kg m-2 s-1", "ilamb", "sum_func", None],
        # gen_mon_gb - MONTHLY
        "evspsblveg":    ["ecan_gb", "evaporation from the canopy", None, 'kg m-2 s-1', "gen_mon_gb", None, None],
        "trans":         ["et_stom_gb", "vegetation transpiration", None, "kg/m2/s", "gen_mon_gb", None, None],
        "snd":           ["snow_depth_gb", "depth of snow layer", None, "m", "gen_mon_gb", None, None],
        "snc":           ["snow_frac", "Snow fraction", None, "1", "gen_mon_gb", None, None],
        # not available "prsn":          ["snowfall", "snowfall rate", None, "kg/m2/s", "gen_mon_gb", None, None],
        "mrros":         ["surf_roff", "surface runoff", None, "kg/m2/s", "gen_mon_gb", None, None],
        "mrrob":         ["sub_surf_roff", "sub-surface runoff", None, "kg/m2/s", "gen_mon_gb", None, None],
        "ts":            ["tstar_gb", "surface temperature", None, "K", "gen_mon_gb", None, None],
        "evspsblsoi":    [["esoil_gb", "et_stom_gb"], "water evaporation from the soil", None, 'kg m-2 s-1', "gen_mon_gb", "minus_func", None],
        "evspsbl":       [["esoil_gb", "ei_gb"], "Evaporation Including Sublimation and Transpiration", None, 'kg m-2 s-1', "gen_mon_gb", "sum_func", None],
        # gen_mon_pft - MONTHLY
        "treeFrac":      ["frac", None, None, "%", "gen_mon_pft", None, None], 
        "c3PftFrac":     ["frac", "Percentage Cover by C3 Plant Functional Type", 100, "%", "gen_mon_pft", None, None],
        "c4PftFrac":     ["frac", "Percentage Cover by C3 Plant Functional Type", 100, "%", "gen_mon_pft", None, None],
        "shrubFrac":     ["frac", "Percentage Cover by Shrub", 100, "%", "gen_mon_pft", None, None],
        "baresoilFrac":  ["frac", "Bare Soil Percentage Area Coverage", 100, "%", "gen_mon_pft", None, None],
        "residualFrac":  ["frac", "Percentage of Grid Cell That Is Land but neither Vegetation Covered nor Bare Soil", 100, "%", "gen_mon_pft", None, None, ],
        # gen_mon_layer - MONTHLY
        "mrsol":         ["smcl", "average layer soil moisture", None, "kg/m2", "gen_mon_layer", None, None],
        "mrsos":         ["smcl", "moisture in top soil (10cm) layer", None, "kg/m2", "gen_mon_layer", "top10cm_func", None],
        # UKESM below
        "landCoverFrac": ["frac", "Plant Functional Type Grid Fraction", 100, "%", "gen_mon_pft", None, None],
        "treeFracBdlDcd":   ["frac", "Broadleaf Deciduous Tree Area Percentage", 100, "%", "gen_mon_pft", None, None],
        "treeFracBdlEvg":   ["frac", "Broadleaf Evergreen Tree Area Percentage", 100, "%", "gen_mon_pft", None, None],
        "treeFracNdlDcd":   ["frac", "Needleleaf Deciduous Tree Area Percentage", 100, "%", "gen_mon_pft", None, None],
        "treeFracNdlEvg":   ["frac", "Needleleaf Evergreen Tree Area Percentage", 100, "%", "gen_mon_pft", None, None],
        "grassFracC3":      ["frac", "C3 Natural Grass Area Percentage", 100, "%", "gen_mon_pft", None, None],
        "cropFracC3":       ["frac", "Percentage Cover by C3 Crops", 100, "%", "gen_mon_pft", None, None],
        "pastureFracC3":    ["frac", "C3 Pasture Area Percentage", 100, "%", "gen_mon_pft", None, None],
        "grassFracC4":      ["frac", "C4 Natural Grass Area Percentage", 100, "%", "gen_mon_pft", None, None],
        "cropFracC4":       ["frac", "Percentage Cover by C4 Crops", 100, "%", "gen_mon_pft", None, None],
        "pastureFracC4":    ["frac", "C4 Pasture Area Percentage", 100, "%", "gen_mon_pft", None, None],
        "cProduct":         [["wood_prod_fast", "wood_prod_med", "wood_prod_slow"], 'Carbon Mass in Products of Land-Use Change', None, "kg/m2", "cn_mon_gb", "sum_func", None],
        #cn_mon_gb - MONTHLY
        "fBNF":          ["n_fix_gb", "Biological Nitrogen Fixation", 1.0/(86400.0*360.0), 'kg m-2 s-1', "cn_mon_gb", None, None],
        "fNup":          ["n_uptake_gb", "Total Plant Nitrogen Uptake", 1.0/(86400.0*360.0), 'kg m-2 s-1', "cn_mon_gb", None, None],
        "fVegSoil":      ["lit_c_mean", "total carbon mass from vegetation directly into the soil", 1.0/(86400.0*360.0), 'kg m-2 s-1', "cn_mon_gb", None, None],
        #n_ann_gb - ANNUAL
        "fNdep_annual":         ["deposition_n", "Dry and Wet Deposition of Reactive Nitrogen onto Land", None, 'kg m-2 s-1', "n_ann_gb", None, None],
           }
    return var_dict
# #############################################################################
