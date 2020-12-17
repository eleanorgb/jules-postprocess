# #############################################################################
def get_var_dict():
#               "outVar":[0.inputVars,1.longName,2.scale,3.units,4.outprofile, 5.func, 6.CMIPprofileName]
    var_dict = {
        # ilamb - MONTHLY
        "cSoil":         ["cs_gb","carbon mass in model soil pool",None,"kg/m2","ilamb",None,"Emon"],
        "cVeg":          ["cv","carbon mass in vegetation",None,"kg/m2","ilamb",None,"Lmon"], 
        "hfss":          ["ftl_gb","sensible heat flux",None,"W/m2","ilamb",None,"Amon"],
        "gpp":           ["gpp_gb","Carbon mass flux out of the atmosphere due to gross primary production",None,"kg/m2/s","ilamb",None,"Lmon"], 
        "fHarvest":      ["harvest_gb","Carbon Mass Flux into Atmosphere Due to Crop Harvesting",1.0/(86400.0*360.0),'kg m-2 s-1',"ilamb",None,"Lmon"],
        "lai":           ["lai_gb","leaf area index",None,"1","ilamb",None,"Lmon"],
        "hfls":          ["latent_heat","latent heat flux",None,"W/m2","ilamb",None,"Amon"],
        "rlds":          ["lw_down","downward longwave radiation",None,"W/m2","ilamb",None,"Amon"],
        "rlus":          ["lw_up","upward longwave radiation",None,"W/m2","ilamb",None,"Amon"],
        "npp_nlim":      ["npp_n_gb","Net Primary Production on Land as Carbon Mass Flux",1.0/(86400.0*360.0),"kg m-2 s-1","ilamb",None,"Lmon"], 
        "npp":           ["npp_gb","Net Primary Production on Land as Carbon Mass Flux",None,"kg m-2 s-1","ilamb",None,"Lmon"], 
        "pr":            ["precip","precipitation rate",None,"kg/m2/s","ilamb",None,"Amon"],
        "ps":            ["pstar","surface air pressure",None,"Pa","ilamb",None,"Emon"],
        "ra":            ["resp_p_gb","Carbon mass flux out of the atmosphere due to autotrophic respiration",None,"kg/m2/s","ilamb",None,"Lmon"], 
        "rh":            ["resp_s_to_atmos_gb","Total Heterotrophic Respiration on Land as Carbon Mass Flux",1.0/(86400.0*360.0),"kg m-2 s-1","ilamb",None,"Lmon"], # is this right?
        "mrso":          ["smc_tot","total soil moisture",None,"kg/m2","ilamb",None,"Lmon"],
        "snw":           ["snow_mass_gb","snow water equivalent",None,"kg/m2","ilamb",None,"LImon"], 
        "rsds":          ["sw_down","downward shortwave radiation",None,"W/m2","ilamb",None,"Amon"],
        "tas":           ["t1p5m_gb","1.5m air Temperature",None,"K","ilamb",None,"Amon"],
        "tsl":           ["t_soil","Temperature of soil - layer",None,"K","ilamb",None,"Lmon"], 
        "mrro":          ["runoff","total runoff",None,"kg/m2/s","ilamb",None,"Lmon"],
        "wtd":           ["zw","water table depth",None,"m","ilamb",None,"Emon"],
        "mrtws":         [["zw","canopy_gb","snow_mass_gb","smc_tot"],"Total water storage",None,"kg/m2","ilamb","twsa_func","Emon"],
        "nbp_nlim":      [["npp_n_gb","resp_s_to_atmos_gb","harvest_gb","WP_fast_out","WP_med_out","WP_slow_out"],'Carbon Mass Flux out of Atmosphere Due to Net Biospheric Production on Land',None,"kg m-2 s-1","ilamb","nbp_func","Lmon"], # downwards positive
        "nbp":           [["npp_gb","resp_s_to_atmos_gb","harvest_gb","WP_fast_out","WP_med_out","WP_slow_out"],'Carbon Mass Flux out of Atmosphere Due to Net Biospheric Production on Land',None,"kg m-2 s-1","ilamb","nbp_func","Lmon"], # downwards positive
        "rsus":          [["sw_down","sw_net"],"upward shortwave radiation",None,"W/m2","ilamb","minus_func","Amon"],
        "fLuc":          [["WP_fast_out","WP_med_out","WP_slow_out"],'Net Carbon Mass Flux into Atmosphere Due to Land-Use Change',1.0/(86400.0*360.0),"kg m-2 s-1","ilamb","sum_func", "Emon"], 
        # gen_mon_gb - MONTHLY
        "evspsblveg":    ["ecan_gb","evaporation from the canopy",None,'kg m-2 s-1',"gen_mon_gb",None,"Lmon"],
        "tran":          ["et_stom_gb","vegetation transpiration",None,"kg/m2/s","gen_mon_gb",None,"Lmon"],
        "snd":           ["snow_depth_gb","depth of snow layer",None,"m","gen_mon_gb",None,"LImon"], 
        "snc":           ["snow_frac","Snow Area Percentage",100,"%","gen_mon_gb",None,"LImon"],
        # not available "prsn":          ["snowfall","snowfall rate",None,"kg/m2/s","gen_mon_gb",None,"Amon"], 
        "mrros":         ["surf_roff","surface runoff",None,"kg/m2/s","gen_mon_gb",None,"Lmon"],
        "mrrob":         ["sub_surf_roff","sub-surface runoff",None,"kg/m2/s","gen_mon_gb",None,"Lmon"],
        "ts":            ["tstar_gb","surface temperature",None,"K","gen_mon_gb",None,"Amon"],
        "evspsblsoi":    [["esoil_gb","et_stom_gb"],"water evaporation from the soil",None,'kg m-2 s-1',"gen_mon_gb","minus_func","Lmon"],
        "evspsbl":       [["esoil_gb","ei_gb"],"Evaporation Including Sublimation and Transpiration",None,'kg m-2 s-1',"gen_mon_gb","sum_func","Amon"],
        # gen_mon_pft - MONTHLY
        "treeFrac":      ["frac","Tree cover percentage",100,"%","gen_mon_pft",None,"Lmon"],   
        "c3PftFrac":     ["frac","Percentage Cover by C3 Plant Functional Type",100,"%","gen_mon_pft",None,"Lmon"],
        "c4PftFrac":     ["frac","Percentage Cover by C4 Plant Functional Type",100,"%","gen_mon_pft",None,"Lmon"],
        "shrubFrac":     ["frac", "Percentage Cover by Shrub",100,"%","gen_mon_pft",None,"Lmon"],
        "baresoilFrac":  ["frac","Bare Soil Percentage Area Coverage",100,"%","gen_mon_pft",None,"Lmon"],
        "residualFrac":  ["frac","Percentage of Grid Cell That Is Land but neither Vegetation Covered nor Bare Soil",100,"%","gen_mon_pft",None,"Lmon",], 
        # gen_mon_layer - MONTHLY
        "mrsol":         ["smcl","average layer soil moisture",None,"kg/m2","gen_mon_layer",None,"Emon"], 
        "mrsos":         ["smcl","moisture in top soil (10cm) layer",None,"kg/m2","gen_mon_layer","top10cm_func","Lmon"],
        # UKESM below
        "landCoverFrac": ["frac","Plant Functional Type Grid Fraction",100,"%","gen_mon_pft",None,"Lmon"],
        "treeFracBdlDcd":   ["frac","Broadleaf Deciduous Tree Area Percentage",100,"%","gen_mon_pft",None,"Lmon"],
        "treeFracBdlEvg":   ["frac","Broadleaf Evergreen Tree Area Percentage",100,"%","gen_mon_pft",None,"Lmon"], 
        "treeFracNdlDcd":   ["frac","Needleleaf Deciduous Tree Area Percentage",100,"%","gen_mon_pft",None,"Lmon"],
        "treeFracNdlEvg":   ["frac","Needleleaf Evergreen Tree Area Percentage",100,"%","gen_mon_pft",None,"Lmon"],
        "grassFracC3":      ["frac","C3 Natural Grass Area Percentage",100,"%","gen_mon_pft",None,"Emon"],
        "cropFracC3":       ["frac","Percentage Cover by C3 Crops",100,"%","gen_mon_pft",None,"Emon"],
        "pastureFracC3":    ["frac","C3 Pasture Area Percentage",100,"%","gen_mon_pft",None,"Emon"],
        "grassFracC4":      ["frac","C4 Natural Grass Area Percentage",100,"%","gen_mon_pft",None,"Emon"],
        "cropFracC4":       ["frac","Percentage Cover by C4 Crops",100,"%","gen_mon_pft",None,"Emon"],
        "pastureFracC4":    ["frac","C4 Pasture Area Percentage",100,"%","gen_mon_pft",None,"Emon"],
        "cProduct":         [["wood_prod_fast","wood_prod_med","wood_prod_slow"],'Carbon Mass in Products of Land-Use Change',None,"kg/m2","cn_mon_gb","sum_func","Lmon"], 
        # drive_day_gb - DAILY
        "mrsos_daily":         ["smcl","moisture in top soil (10cm) layer",None,"kg/m2","gen_day_layer","top10cm_func","day"],
        "rlds_daily":          ["lw_down","downward longwave radiation",None,"W/m2","drive_day_gb",None,"day"],
        "rlus_daily":          ["lw_up","upward longwave radiation",None,"W/m2","drive_day_gb",None,"day"],
        "rsds_daily":          ["sw_down","downward shortwave radiation",None,"W/m2","drive_day_gb",None,"day"],
        "tasmax_daily":        ["t1p5m_gb_dailymax","Daily Maximum Near-Surface Air Temperature",None,"K","drive_day_gb",None,"day"],
        "tasmin_daily":        ["t1p5m_gb_dailymin","Daily Minimum Near-Surface Air Temperature",None,"K","drive_day_gb",None,"day"], 
        "tas_daily":           ["t1p5m_gb","Daily Near-Surface Air Temperature",None,"K","drive_day_gb",None,"day"], 
        "rsus_daily":          [["sw_down","sw_net"],"upward shortwave radiation",None,"W/m2","drive_day_gb","minus_func","day"],
        "tsl_daily":           ["t_soil","Temperature of soil - layer",None,"K","gen_day_layer",None,"Eday"], 
        # gen_day_gb - DAILY
        "hfss_daily":          ["ftl_gb","sensible heat flux",None,"W/m2","gen_day_gb",None,"day"], 
        "hfls_daily":          ["fqw_gb","latent heat flux",2.5e6,"W/m2","gen_day_gb",None,"day"],  
        "mrso_daily":          ["smc_tot","total soil moisture",None,"kg/m2","gen_day_gb",None,"day"],
        "snc_daily":           ["snow_frac","Snow Area Percentage",100,"%","gen_day_gb",None,"day"],
        "snd_daily":           ["snow_depth_gb","depth of snow layer",None,"m","gen_day_gb",None,"Eday"],
        "snw_daily":           ["snow_mass_gb"," Surface Snow Amount",None,"kg/m2","gen_day_gb",None,"day"],
        "mrro_daily":          [["surf_roff","sub_surf_roff"],"Runoff (surface + subsurface)",None,"kg/m2/s","gen_day_gb","sum_func","day"],
        #cn_mon_gb - MONTHLY
        "fBNF":          ["n_fix_gb","Biological Nitrogen Fixation",1.0/(86400.0*360.0),'kg m-2 s-1',"cn_mon_gb",None,"Emon"],
        "fNup":          ["n_uptake_gb","Total Plant Nitrogen Uptake",1.0/(86400.0*360.0),'kg m-2 s-1',"cn_mon_gb",None,"Emon"],
        "fVegSoil":      ["lit_c_mean","total carbon mass from vegetation directly into the soil",1.0/(86400.0*360.0),'kg m-2 s-1',"cn_mon_gb",None,"Lmon"],
        #n_ann_gb - ANNUAL
        "fNdep_annual":         ["deposition_n","Dry and Wet Deposition of Reactive Nitrogen onto Land",None,'kg m-2 s-1',"n_ann_gb",None,"Eyr"],
           }
    return var_dict
# #############################################################################
