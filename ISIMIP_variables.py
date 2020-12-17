def get_var_dict():
#############"outVar":[0.inputVars,1.longName,2.scale,3.units,4.outprofile,5.func]
    var_dict = {
        ##################### ilamb - MONTHLY #####################
        #"burntarea":     ["burnt_area_gb","Burnt area fraction",None,"%","ilamb","burntarea_func"], # check fire currently switched off
        "thawdepth":      ["depth_unfrozen","Thaw Depth",None,"m","ilamb",None],
        "evap":           ["fqw_gb","Evapotranspiration",None,"kg/m2/s","ilamb",None],
        "gpp":            ["gpp_gb","Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land",None,"kg/m2/s","ilamb",None],
        "lai":            ["lai_gb","Leaf Area Index",None,"1","ilamb",None],
        "ra":             ["resp_p_gb","Carbon Mass Flux into Atmosphere due to Autotrophic (plant) Respiration on Land",None,"kg/m2/s","ilamb",None],
        "rh":             ["resp_s_to_atmos_gb","Carbon Mass Flux into Atmosphere due to Heterotrophic Respiration on Land",1.0/(86400.0*360.0),"kg/m2/s","ilamb",None],
        "qtot":           ["runoff","Total Runoff",None,"kg/m2/s","ilamb",None],
        "swe":            ["snow_mass_gb","Snow Water Equivalent",None,"kg/m2","ilamb",None],
        "tas":            ["t1p5m_gb","Near-Surface Air Temperature",None,"K","ilamb",None],
        "tsl":            ["t_soil","Temperature of soil",None,"K","ilamb",None],
        "ecoatmflux":     [["npp_gb","resp_s_to_atmos_gb","harvest_gb","WP_fast_out","WP_med_out","WP_slow_out"],"Carbon Mass Flux out of Atmosphere due to Net biome Production on Land (NBP)",None,"kg/m2/s","ilamb","nbp_func"],
        "ecoatmflux_nlim":[["npp_n_gb","resp_s_to_atmos_gb","harvest_gb","WP_fast_out","WP_med_out","WP_slow_out"],"Carbon Mass Flux out of Atmosphere due to Net biome Production on Land (NBP) Nitrogen limitied",None,"kg/m2/s","ilamb","nbp_func"],
        "tws":           [["zw","canopy_gb","snow_mass_gb","smc_tot"],"Total Water Storage",None,"kg/m2","ilamb","twsa_func"],
        "thawdepth_annual":  ["depth_unfrozen","Annual Maximum Thaw Depth",None,"m","ilamb","annmax_func"],
        "npp_nlim":      ["npp_n_gb","Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land",1.0/(86400.0*360.0),"kg m-2 s-1","ilamb",None,"Lmon"], 
        "npp":           ["npp_gb","Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land",None,"kg m-2 s-1","ilamb",None,"Lmon"], 
        ##################### gen_mon_gb - MONTHLY #####################
        "albedo":        ["albedo_land","Albedo",None,"1","gen_mon_gb",None],
        "intercep":     ["ecan_gb","Evaporation from Canopy (interception)",None,"kg/m2/s","gen_mon_gb",None],
        "sublim":        ["ei_gb","Gridbox sublimation from lying snow or sea-ice",None,"kg/m2/s","gen_mon_gb",None],
        "trans":         ["et_stom_gb","Transpiration",None,"kg/m2/s","gen_mon_gb",None],
        "dis":            ["rflow","Discharge",None,"kg/m2/s","gen_mon_gb",None],
        "snd":           ["snow_depth_gb","Snow Depth",None,"m","gen_mon_gb",None],
        "qs":            ["surf_roff","Surface Runoff",None,"kg/m2/s","gen_mon_gb",None],
        "qsb":           ["sub_surf_roff","Subsurface Runoff",None,"kg/m2/s","gen_mon_gb",None],
        "esoil":         [["esoil_gb","et_stom_gb"],"Water Evaporation from Soil",None,"kg/m2/s","gen_mon_gb","minus_func"], # check esoil-et_stom???????? incldue sublimation
        ##################### gen_mon_pft - MONTHLY #####################
        "intercep-pft": ["ecan","PFT Evaporation from Canopy (interception)",None,"kg/m2/s","gen_mon_pft",None],
        "sublim-pft":   ["ei","Tile sublimation from lying snow or sea-ice",None,"kg/m2/s","gen_mon_pft",None],
        "trans-pft":     ["et_stom","PFT transpiration",None,"kg/m2/s","gen_mon_pft",None],
        "pft-pft":       ["frac","Plant Functional Type Grid Fraction",None,"1","gen_mon_pft",None],
        "esoil-pft":     [["esoil","et_stom"],"PFT Evaporation from soil",None,"kg/m2/s","gen_mon_pft","minus_func"], 
        ##################### gen_mon_layer - MONTHLY #####################
        "soilmoist":     ["smcl","Soil moisture for each layer",None,"kg/m2","gen_mon_layer",None],
        "soilmoistfroz": [["sthf","sm_sat"],"Frozen soil moisture for each layer",None,"kg/m2","gen_mon_layer","sth_func"],
        ##################### drive_day_gb - DAILY #####################
        "rlds_daily":           ["lw_down","Downward longwave radiation",None,"W/m2","drive_day_gb",None,"day"],
        "rlus_daily":           ["lw_up","Upward longwave radiation",None,"W/m2","drive_day_gb",None,"day"],
        "precip_daily":         ["precip","Precipitation rate",None,"kg/m2/s","drive_day_gb",None,"day"],
        "rsds_daily":           ["sw_down","Downward shortwave radiation",None,"W/m2","drive_day_gb",None,"day"],
        "tas_daily":            ["t1p5m_gb","Daily Near-Surface Air Temperature",None,"K","drive_day_gb",None],
        "tasmin_daily":         ["t1p5m_gb_dailymin","Daily Minimum Near-Surface Air Temperature",None,"K","drive_day_gb",None],
        "tasmax_daily":         ["t1p5m_gb_dailymax","Daily Maximum Near-Surface Air Temperature",None,"K","drive_day_gb",None],
        "rsus_daily":           [["sw_down","sw_net"],"Upward shortwave radiation",None,"W/m2","drive_day_gb","minus_func","day"],
        ##################### gen_day_gb - DAILY #####################
        "intercep_daily":      ["ecan_gb","Evaporation from Canopy (interception)",None,"kg/m2/s","gen_day_gb",None],
        "sublim_daily":         ["ei_gb","Gridbox sublimation from lying snow or sea-ice",None,"kg/m2/s","gen_day_gb",None],
        "trans_daily":          ["et_stom_gb","Transpiration",None,"kg/m2/s","gen_day_gb",None], 
        "hfss_daily":           ["ftl_gb","sensible heat flux",None,"W/m2","gen_day_gb",None,"day"],  #also Eday in tables
        "evap_daily":           ["fqw_gb","Evapotranspiration",None,"kg/m2/s","gen_day_gb",None],
        "gpp_daily":            ["gpp_gb","Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land",None,"kg/m2/s","gen_day_gb",None],
        "lai_daily":            ["lai_gb","TOTAL Leaf Area Index",None,"1","gen_day_gb",None],
        "ra_daily":             ["resp_p_gb","Carbon Mass Flux into atmosphere due to Autotrophic (Plant) Respiration on Land",None,"kg/m2/s","gen_day_gb",None],
        "dis_daily":            ["rflow","Discharge",None,"kg/m2/s","gen_day_gb",None],
        "mrso_daily":           ["smc_tot","total soil moisture",None,"kg/m2","gen_day_gb",None,"day"],
        "snd_daily":            ["snow_depth_gb","Snow depth",None,"m","gen_day_gb",None],
        "snc_daily":            ["snow_frac","Snow fraction",None,"1","gen_day_gb",None,"day"],
        "swe_daily":            ["snow_mass_gb","Snow water equivalent",None,"kg/m2","gen_day_gb",None],
        "qsb_daily":            ["sub_surf_roff","Subsurface runoff",None,"kg/m2/s","gen_day_gb",None],
        "qs_daily":             ["surf_roff","Surface runoff",None,"kg/m2/s","gen_day_gb",None],
        "qtot_daily":           [["surf_roff","sub_surf_roff"],"Runoff (surface + subsurface)",None,"kg/m2/s","gen_day_gb","sum_func"],
        ##################### gen_day_pft - DAILY #####################
        "intercep-pft_daily":  ["ecan","PFT Evaporation from Canopy (interception)",None,"kg/m2/s","gen_day_pft",None],
        "sublim-pft_daily":     ["ei","Tile sublimation from lying snow for land tiles",None,"kg/m2/s","gen_day_pft",None],
        "trans-pft_daily":      ["et_stom","PFT transpiration",None,"kg/m2/s","gen_day_pft",None], # et_stom 
        "evap-pft_daily":       ["fqw","Evapotranspiration",None,"kg/m2/s","gen_day_pft",None],
        #gen_day_layer - DAILY
        "tsl_daily":            ["t_soil","Temperature of soil - layer",None,"K","gen_day_layer",None],
        "soilmoist_daily":      ["smcl","Soil moisture for each layer",None,"kg/m2","gen_day_layer",None],
        # cn_day_nbp - DAILY
        "npp_daily":            ["npp_gb","TOTAL Carbon Mass flux out of Atmosphere due to Net Primary Production on Land",None,"kg/m2/s","cn_day_nbp",None],
        "npp_daily_nlim":       ["npp_n_gb","TOTAL Carbon Mass flux out of Atmosphere due to Net Primary Production on Land Nitrogen limited",1.0/(86400.0*360.0),"kg/m2/s","cn_day_nbp",None], 
        "rh_daily":             ["resp_s_to_atmos_gb","Carbon Mass Flux into atmosphere due to Heterotrophic Respiration on Land",1.0/(86400.0*360.0),"kg/m2/s","cn_day_nbp",None],
        "ecoatmflux_daily_nlim":[["npp_n_gb","resp_s_to_atmos_gb","harvest_gb","WP_fast_out","WP_med_out","WP_slow_out"],"Carbon Mass Flux out of Atmosphere due to Net biome Production on Land (NBP)  n_limited",None,"kg/m2/s","cn_day_nbp","nbp_func"],
        "ecoatmflux_daily":     [["npp_gb","resp_s_to_atmos_gb","harvest_gb","WP_fast_out","WP_med_out","WP_slow_out"],"Carbon Mass Flux out of Atmosphere due to Net biome Production on Land (NBP)",None,"kg/m2/s","cn_day_nbp","nbp_func"],
        #gen_ann_pftlayer - ANNUAL
        "pft-pft_annual":    ["frac","Plant Functional Type Grid Fraction",None,"1","gen_ann_pftlayer",None],
        "gpp-pft_annual":    ["gpp","PFT Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land",None,"kg/m2/s","gen_ann_pftlayer",None],
        "lai-pft_annual":    ["lai","PFT Leaf Area Index",None,"1","gen_ann_pftlayer",None],
        "npp-pft_annual":    ["npp","PFT Carbon Mass flux out of Atmosphere due to Net Primary Production on Land",1.0/(86400.0*360.0),"kg/m2/s","gen_ann_pftlayer",None],
        #c_ann_gb - ANNUAL
        "csoil_annual":      ["cs_gb","Carbon Mass in soil",None,"kg/m2","c_ann_gb",None],
        "cveg_annual":       ["cv","Carbon Mass in Vegetation",None,"kg/m2","c_ann_gb",None],
        #c_ann_pftlayer - ANNUAL
        "cvegag-pft_annual": [["woodC","leafC"],"PFT Carbon Mass in above ground Vegetation",None,"kg/m2","c_ann_pftlayer","sum_func"],
        "cveg-pft_annual":   ["c_veg","PFT Carbon Mass in Vegetation",None,"kg/m2","c_ann_pftlayer",None],
        "cvegag_annual":     [["frac","woodC","leafC"],"Carbon Mass in Above Ground Vegetation Biomass",None,"kg/m2","c_ann_pftlayer","fracweight_func"],
        "cvegbg-pft_annual": ["rootC","PFT Carbon Mass in below ground Vegetation",None,"kg/m2","c_ann_pftlayer",None],
        "cvegbg_annual":     [["frac","rootC"],"TOTAL Carbon Mass in below ground Vegetation",None,"kg/m2","c_ann_pftlayer","fracweight_func"],
        #"csoil_annual":      ["cs","Carbon Mass in each soil pool",None,"kg/m2","c_ann_pftlayer",None],
        #n_ann_pftlayer - ANNUAL
        "npp_nlim_annual":      ["npp_n_gb","Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land",1.0/(86400.0*360.0),"kg m-2 s-1","c_ann_gb",None], 
        "npp_annual":           ["npp_gb","Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land",None,"kg m-2 s-1","gen_ann_gb",None], 
        "gpp_annual":            ["gpp_gb","Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land",None,"kg/m2/s","gen_ann_gb",None],
        "npp-pft_nlim_annual":    ["npp_n","PFT Carbon Mass flux out of Atmosphere due to Net Primary Production on Land nitrogen limited",1.0/(86400.0*360.0),"kg/m2/s","n_ann_pftlayer",None]
            }
    return var_dict
