def get_var_dict():
#############"outVar":[0.inputVars, 1.longName, 2.scale, 3.units, 4.outprofile, 5.func]
    var_dict = {
        ##################### ilamb - MONTHLY #####################
        # "burntarea-total":        ["burnt_area_gb", "Burnt area fraction", None, "%/year", "ilamb", "burntarea_func"], # need specific function to convert units
        # "thawdepth":        ["depth_unfrozen", "Thaw Depth", None, "m", "ilamb", None],
        # "evap-total":       ["fqw_gb", "Evapotranspiration", None, "kg m-2 s-1", "ilamb", None],
        # "gpp-total":        ["gpp_gb", "Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land", None, "kg m-2 s-1", "ilamb", None],
        # "lai":              ["lai_gb", "Leaf Area Index", None, "1", "ilamb", None],
        # "ra-total":         ["resp_p_gb", "Carbon Mass Flux into Atmosphere due to Autotrophic (plant) Respiration on Land", None, "kg m-2 s-1", "ilamb", None],
        # "rh-total":         ["resp_s_to_atmos_gb", "Carbon Mass Flux into Atmosphere due to Heterotrophic Respiration on Land", 1.0/(86400.0*360.0), "kg m-2 s-1", "ilamb", None],
        # "qtot":             ["runoff", "Total Runoff", None, "kg m-2 s-1", "ilamb", None],
        # "swe":              ["snow_mass_gb", "Snow Water Equivalent", None, "kg m-2", "ilamb", None],
        # "tas":              ["t1p5m_gb", "Near-Surface Air Temperature", None, "K", "ilamb", None],
        # "tsl":              ["t_soil", "Temperature of Soil", None, "K", "ilamb", None],
        # "nbp-total":        [["npp_gb", "resp_s_to_atmos_gb", "harvest_gb", "WP_fast_out", "WP_med_out", "WP_slow_out"], "Carbon Mass Flux out of Atmosphere due to Net Biospheric Production on Land", None, "kg m-2 s-1", "ilamb", "nbp_func"],
        # "nbp-total_nlim":   [["npp_n_gb", "resp_s_to_atmos_gb", "harvest_gb", "WP_fast_out", "WP_med_out", "WP_slow_out"], "Carbon Mass Flux out of Atmosphere due to Net Biospheric Production on Land Nitrogen limited", None, "kg m-2 s-1", "ilamb", "nbp_func"],
        # "tws":              [["zw", "canopy_gb", "snow_mass_gb", "smc_tot"], "Total Water Storage", None, "kg m-2", "ilamb", "twsa_func"],
        # "thawdepth_annual": ["depth_unfrozen", "Annual Maximum Thaw Depth", None, "m", "ilamb", "annmax_func"], 
        # "npp-total_nlim":   ["npp_n_gb", "Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land", 1.0/(86400.0*360.0), "kg m-2 s-1", "ilamb", None, "Lmon"],
        # "npp-total":        ["npp_gb", "Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land", None, "kg m-2 s-1", "ilamb", None, "Lmon"],
        # ##################### gen_mon_gb - MONTHLY #####################
        # "landalbedo":       ["albedo_land", "Surface Albedo of Land", None, "1", "gen_mon_gb", None],
        # "intercep-total":   ["ecan_gb", "Evaporation from Canopy (interception)", None, "kg m-2 s-1", "gen_mon_gb", None],
        # "sublim":           ["ei_gb", "Gridbox sublimation from lying snow or sea-ice", None, "kg m-2 s-1", "gen_mon_gb", None],
        # "trans-total":      ["et_stom_gb", "Transpiration", None, "kg m-2 s-1", "gen_mon_gb", None],
        # "dis":              ["rflow", "Discharge", None, "m3 s-1", "gen_mon_gb", "rflow_func"],
        # "snd":              ["snow_depth_gb", "Snow depth", None, "m", "gen_mon_gb", None],
        # "qs":               ["surf_roff", "Surface Runoff", None, "kg m-2 s-1", "gen_mon_gb", None],
        # "qsb":              ["sub_surf_roff", "Subsurface Runoff", None, "kg m-2 s-1", "gen_mon_gb", None],
        # "esoil-total":      [["esoil_gb", "et_stom_gb"], "Water Evaporation from Soil", None, "kg m-2 s-1", "gen_mon_gb", "minus_func"], # check esoil-et_stom???????? incldue sublimation
        # ##################### gen_mon_pft - MONTHLY #####################
        # "intercep-pft":     ["ecan", "Evaporation from Canopy (interception)", None, "kg m-2 s-1", "gen_mon_pft", None],
        # "sublim-pft":       ["ei", "Tile sublimation from lying snow or sea-ice", None, "kg m-2 s-1", "gen_mon_pft", None],
        # "trans-pft":        ["et_stom", "Transpiration", None, "kg m-2 s-1", "gen_mon_pft", None],
        # "pft-pft":          ["frac", "Plant Functional Type Grid Fraction", 100., "%", "gen_mon_pft", None],
        # "esoil-pft":        [["esoil", "et_stom"], "Water Evaporation from Soil", None, "kg m-2 s-1", "gen_mon_pft", "minus_func"],
        # ##################### gen_mon_layer - MONTHLY #####################
        # "soilmoist":        ["smcl", "Total Soil Moisture Content", None, "kg m-2", "gen_mon_layer", None],
        # "soilmoistfroz":    [["sthf", "sm_sat"], "Frozen Soil Moisture Content", None, "kg m-2", "gen_mon_layer", "sth_func"],
        # ##################### drive_day_gb - DAILY #####################
        # "rlds_daily":           ["lw_down", "Downward longwave radiation", None, "W m-2", "drive_day_gb", None, "day"],
        # "rlus_daily":           ["lw_up", "Upward longwave radiation", None, "W m-2", "drive_day_gb", None, "day"],
        # "precip_daily":         ["precip", "Precipitation rate", None, "kg m-2 s-1", "drive_day_gb", None, "day"],
        # "rsds_daily":           ["sw_down", "Downward shortwave radiation", None, "W m-2", "drive_day_gb", None, "day"],
        # "tas_daily":            ["t1p5m_gb", "Daily Near-Surface Air Temperature", None, "K", "drive_day_gb", None],
        # "tasmin_daily":         ["t1p5m_gb_dailymin", "Daily Minimum Near-Surface Air Temperature", None, "K", "drive_day_gb", None],
        # "tasmax_daily":         ["t1p5m_gb_dailymax", "Daily Maximum Near-Surface Air Temperature", None, "K", "drive_day_gb", None],
        # "rsus_daily":           [["sw_down", "sw_net"], "Upward shortwave radiation", None, "W m-2", "drive_day_gb", "minus_func", "day"],
        # ##################### gen_day_gb - DAILY #####################
        # "intercep_daily":      ["ecan_gb", "Evaporation from Canopy (interception)", None, "kg m-2 s-1", "gen_day_gb", None],
        # "sublim_daily":         ["ei_gb", "Gridbox sublimation from lying snow or sea-ice", None, "kg m-2 s-1", "gen_day_gb", None],
        # "trans_daily":          ["et_stom_gb", "Transpiration", None, "kg m-2 s-1", "gen_day_gb", None],
        # "hfss_daily":           ["ftl_gb", "sensible heat flux", None, "W m-2", "gen_day_gb", None, "day"], # also Eday in tables
        # "evap_daily":           ["fqw_gb", "Evapotranspiration", None, "kg m-2 s-1", "gen_day_gb", None],
        # "gpp_daily":            ["gpp_gb", "Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land", None, "kg m-2 s-1", "gen_day_gb", None],
        # "lai_daily":            ["lai_gb", "TOTAL Leaf Area Index", None, "1", "gen_day_gb", None],
        # "ra_daily":             ["resp_p_gb", "Carbon Mass Flux into atmosphere due to Autotrophic (Plant) Respiration on Land", None, "kg m-2 s-1", "gen_day_gb", None],
        # "dis_daily":            ["rflow", "Discharge", None, "m3 s-1", "gen_day_gb", "rflow_func"], 
        # "mrso_daily":           ["smc_tot", "total soil moisture", None, "kg m-2", "gen_day_gb", None, "day"],
        # "snd_daily":            ["snow_depth_gb", "Snow depth", None, "m", "gen_day_gb", None],
        # "snc_daily":            ["snow_frac", "Snow fraction", None, "1", "gen_day_gb", None, "day"], 
        # "swe_daily":            ["snow_mass_gb", "Snow water equivalent", None, "kg m-2", "gen_day_gb", None],
        # "qsb_daily":            ["sub_surf_roff", "Subsurface runoff", None, "kg m-2 s-1", "gen_day_gb", None],
        # "qs_daily":             ["surf_roff", "Surface runoff", None, "kg m-2 s-1", "gen_day_gb", None],
        # "qtot_daily":           [["surf_roff", "sub_surf_roff"], "Total Runoff", None, "kg m-2 s-1", "gen_day_gb", "sum_func"],
        # ##################### gen_day_pft - DAILY #####################
        # "intercep-pft_daily":  ["ecan", "Evaporation from Canopy (interception)", None, "kg m-2 s-1", "gen_day_pft", None],
        # "sublim-pft_daily":     ["ei", "Tile sublimation from lying snow for land tiles", None, "kg m-2 s-1", "gen_day_pft", None],
        # "trans-pft_daily":      ["et_stom", "Transpiration", None, "kg m-2 s-1", "gen_day_pft", None], # et_stom 
        # "evap-pft_daily":       ["fqw", "Evapotranspiration", None, "kg m-2 s-1", "gen_day_pft", None],
        # #gen_day_layer - DAILY
        # "tsl_daily":            ["t_soil", "Temperature of soil - layer", None, "K", "gen_day_layer", None],
        # "soilmoist_daily":      ["smcl", "Total Soil Moisture Content", None, "kg m-2", "gen_day_layer", None],
        # # cn_day_nbp - DAILY
        # "npp_daily":            ["npp_gb", "TOTAL Carbon Mass flux out of Atmosphere due to Net Primary Production on Land", None, "kg m-2 s-1", "cn_day_nbp", None],
        # "npp_daily_nlim":       ["npp_n_gb", "TOTAL Carbon Mass flux out of Atmosphere due to Net Primary Production on Land Nitrogen limited", 1.0/(86400.0*360.0), "kg m-2 s-1", "cn_day_nbp", None],
        # "rh_daily":             ["resp_s_to_atmos_gb", "Carbon Mass Flux into atmosphere due to Heterotrophic Respiration on Land", 1.0/(86400.0*360.0), "kg m-2 s-1", "cn_day_nbp", None],
        # "nbp-total_daily_nlim":[["npp_n_gb", "resp_s_to_atmos_gb", "harvest_gb", "WP_fast_out", "WP_med_out", "WP_slow_out"], "Carbon Mass Flux out of Atmosphere due to Net biome Production on Land (NBP)  n_limited", None, "kg m-2 s-1", "cn_day_nbp", "nbp_func"],
        # "nbp-total_daily":     [["npp_gb", "resp_s_to_atmos_gb", "harvest_gb", "WP_fast_out", "WP_med_out", "WP_slow_out"], "Carbon Mass Flux out of Atmosphere due to Net biome Production on Land (NBP)", None, "kg m-2 s-1", "cn_day_nbp", "nbp_func"],
        #gen_ann_pftlayer - ANNUAL
        #"pft-pft_annual":    ["frac", "Plant Functional Type Grid Fraction", 100., "%", "gen_ann_pftlayer", None],
        #"gpp-pft_annual":    ["gpp", "Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land", None, "kg m-2 s-1", "gen_ann_pftlayer", None],
        #"lai-pft_annual":    ["lai", "Leaf Area Index", None, "1", "gen_ann_pftlayer", None],
        #"npp-pft_annual":    ["npp", "PFT Carbon Mass flux out of Atmosphere due to Net Primary Production on Land", 1.0/(86400.0*360.0), "kg m-2 s-1", "gen_ann_pftlayer", None],
        #c_ann_gb - ANNUAL
        #"csoil-total_annual":      ["cs_gb", "Carbon Mass in Soil Pool", None, "kg m-2", "c_ann_gb", None],
        #"cveg-total_annual":       ["cv", "Carbon Mass in Vegetation", None, "kg m-2", "c_ann_gb", None],
        #c_ann_pftlayer - ANNUAL
        #"cvegag-pft_annual": [["woodC", "leafC"], "Carbon Mass in Above Ground Vegetation Biomass", None, "kg m-2", "c_ann_pftlayer", "sum_func"],
        #"cveg-pft_annual":   ["c_veg", "Carbon Mass in Vegetation", None, "kg m-2", "c_ann_pftlayer", None],
        #"cvegag-total_annual":     [["frac", "woodC", "leafC"], "Carbon Mass in Above Ground Vegetation Biomass", None, "kg m-2", "c_ann_pftlayer", "fracweight_func"],
        #"cvegbg-pft_annual": ["rootC", "Carbon Mass in Below Ground Vegetation Biomass", None, "kg m-2", "c_ann_pftlayer", None],
        #"cvegbg-total_annual":     [["frac", "rootC"], "Carbon Mass in Below Ground Vegetation Biomass", None, "kg m-2", "c_ann_pftlayer", "fracweight_func"],
        #n_ann_pftlayer - ANNUAL
        #"npp_nlim_annual":      ["npp_n_gb", "Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land", 1.0/(86400.0*360.0), "kg m-2 s-1", "c_ann_gb", None],
        #"npp-total_annual":           ["npp_gb", "Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land", None, "kg m-2 s-1", "gen_ann_gb", None],
         #"gpp-total_annual":            ["gpp_gb", "Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land", None, "kg m-2 s-1", "gen_ann_gb", None],
         #"npp-pft_nlim_annual":    ["npp_n", "PFT Carbon Mass flux out of Atmosphere due to Net Primary Production on Land nitrogen limited", 1.0/(86400.0*360.0), "kg m-2 s-1", "n_ann_pftlayer", None]
        #NEW ADDITIONS
        #"thawdepth_daily":        ["t_soil", "Thaw Depth", None, "m", "gen_day_layer", "thawdepth_func"],
         "csoillayer-total_annual":["cs", "Carbon Mass in Soil Pool in each layer", None, "kg m-2", "c_ann_pftlayer", "cs_func"],
        "nsoil-total_annual":     ["ns_gb", "Nitrogen Mass in Soil Pool", None, "kg m-2", "n_ann_gb", None],
        "nsoillayer-total_annual":["ns", "Nitrogen Mass in Soil Pool in each layer", None, "kg m-2", "n_ann_pftlayer", "cs_func"],
        "nmineral_annual":        ["n_inorg", "Mineral Nitrogen in the soil", None, "kg m-2", "n_ann_pftlayer", None],
        "wetlandfrac":            ["fwetl", "Wetland Fraction", 100.0, "%", "gen_mon_gb", None],
        "cproduct-total_annual":  [["wood_prod_fast", "wood_prod_med", "wood_prod_slow"], "Carbon in Products of Land Use Change", None, "kg/m2", "c_ann_gb", "sum_func"],
        "cleaf-total_annual":     [["frac","leafC"], "Carbon Mass in Leaves", None, "kg m-2", "c_ann_pftlayer", "fracweight_func"],
        "cwood-total_annual":     [["frac","woodC"], "Carbon Mass in Wood", None, "kg m-2", "c_ann_pftlayer", "fracweight_func"],	
        "croot-total_annual":     [["frac","rootC"], "Carbon Mass in Roots", None, "kg m-2", "c_ann_pftlayer", "fracweight_func"],
        "nveg_annual":            ["n_veg_gb", "Nitrogen Mass in Vegetation", None, "kg m-2", "cn_mon_gb", None],
        "ra-pft":                 ["resp_p", "Carbon Mass Flux into Atmosphere due to Autotrophic (plant) Respiration on Land", None, "kg m-2 s-1", "gen_mon_pft", None],
        "rhlayer-total":          ["resp_s", "Carbon Mass Flux into Atmosphere due to Heterotrophic Respiration on Land in each layer", None, "kg m-2 s-1", "c_mon_pftlayer", "cs_func"],
        "wetlandch4":             ["fch4_wetl_cs", "Grid averaged methane emissions from wetland (NPP)", None, "kg m-2 s-1", "gen_mon_gb", None],
        "ch4":                    [["fwetl", "fch4_wetl_cs"], "Total Surface Carbon Mass Flux into CH4 emissions", None, "kg m-2 s-1", "gen_mon_gb", None], # do we add fire here? what do we do about microbial model?
        "ffire-total":            ["fire_em_CO2_gb", "Carbon Mass Flux into Atmosphere due to CO₂ Emission from Fire", None, "kg m-2 s-1", "gen_mon_gb", None],
        "ffire-total_annual":     ["fire_em_CO2_gb", "Carbon Mass Flux into Atmosphere due to CO₂ Emission from Fire", None, "kg m-2 s-1", "gen_ann_gb", None],
        "ffire-pft_annual":       ["fire_em_CO2", "Carbon Mass Flux into Atmosphere due to CO₂ Emission from Fire", None, "kg m-2 s-1", "gen_ann_pftlayer", None],
        "fluc-total":             [["harvest_gb", "WP_fast_out", "WP_med_out", "WP_slow_out"], "CO2 Flux to Atmosphere from Land Use Change", 1.0/(86400.0*360.0), "kg m-2 s-1", "ilamb", "conv360_func"],
        "fvegsoil-total":         ["lit_c_mean", "Total Carbon Flux from Vegetation Directly to Soil", 1.0/(86400.0*360.0), "kg m-2 s-1", "cn_mon_gb", None],
        "fvegsoil-pft":           ["lit_c", "Total Carbon Flux from Vegetation Directly to Soil", 1.0/(86400.0*360.0), "kg m-2 s-1", "c_mon_pftlayer", None],
        "fcropharvest-total":     ["harvest_gb", "CO2 Flux to Atmosphere from Crop Harvesting", 1.0/(86400.0*360.0), "kg m-2 s-1", "ilamb", None],
        "fngasfire":              ["fire_em_NOx_gb", "Total Nitrogen lost to the atmosphere (including NHx, N2O, N2) from fire", None, "kg m-2 s-1", "gen_mon_gb", None],
        "fngasfire_annual":       ["fire_em_NOx_gb", "Total Nitrogen lost to the atmosphere (including NHx, N2O, N2) from fire", None, "kg m-2 s-1", "gen_ann_gb", None],
        "fngas":                  ["n_gas_gb", "Total nitrogen lost to the atmosphere (sum of NHx, NOx, N2O, N2)", 1.0/(86400.0*360.0), "kg m-2 s-1", "cn_mon_gb", None], #should we add fire here?
        "fngas_annual":           ["n_gas_gb", "Total nitrogen lost to the atmosphere (sum of NHx, NOx, N2O, N2)", 1.0/(86400.0*360.0), "kg m-2 s-1", "n_ann_gb", None],  #should we add fire here?
        "fnleach_annual":         ["n_leach", "Total nitrogen loss due to leaching or runoff (sum of ammonium, nitrite and nitrate)", None, "kg m-2 s-1", "n_ann_gb", None],
        "fnleach":                ["n_leach", "Total nitrogen loss due to leaching or runoff (sum of ammonium, nitrite and nitrate)", None, "kg m-2 s-1", "cn_mon_gb", None],
        "fnproduct":              [["frac", "harvest_n"], "Nitrogen in deforested or harvested biomass as a result of anthropogenic land-use or change", 1.0/(86400.0*360.0), "kg m-2 s-1", "cn_mon_gb", "fracweight_func"],	
        "fnproduct_annual":       ["harvest_n_gb", "Nitrogen in deforested or harvested biomass as a result of anthropogenic land-use or change", 1.0/(86400.0*360.0), "kg m-2 s-1", "n_ann_gb", None],
        "fbnf":                   ["n_fix_gb", "Biological nitrogen fixation", 1.0/(86400.0*360.0), "kg m-2 s-1", "cn_mon_gb", None],
        "fbnf_annual":            ["n_fix_gb", "Biological nitrogen fixation", 1.0/(86400.0*360.0), "kg m-2 s-1", "n_ann_gb", None],
        "fndep_annual":           ["deposition_n", "Dry and wet deposition of reactive nitrogen onto land", None, "kg m-2 s-1", "n_ann_gb", None],
        "fnetmin":                [["minl_n_gb", "immob_n_gb"], "Net nitrogen release from soil and litter as the outcome of nitrogen immobilisation and gross mineralisation", 1.0/(86400.0*360.0), "kg m-2 s-1", "cn_mon_gb", "minus_func"],
        "fnetmin_annual":         [["minl_n_gb", "immob_n_gb"], "Net nitrogen release from soil and litter as the outcome of nitrogen immobilisation and gross mineralisation", 1.0/(86400.0*360.0), "kg m-2 s-1", "n_ann_gb", "minus_func"],
        "nloss":                  [["n_leach", "n_loss", "n_gas_gb"], "Total nitrogen lost (including NHx, NOx, N2O, N2 and leaching)", None, "kg m-2 s-1", "cn_mon_gb", "conv360_func"],
        "nloss_annual":           [["n_leach", "n_loss", "n_gas_gb"], "Total nitrogen lost (including NHx, NOx, N2O, N2 and leaching)", None, "kg m-2 s-1", "n_ann_gb", "conv360_func"],
        "fapar-pft":              ["fapar", "Fraction of Absorbed Photosynthetically Active Radiation", 100, "%", "ilamb", None],
        "fapar":                  [["frac", "fapar"], "Fraction of Absorbed Photosynthetically Active Radiation", 100, "%", "ilamb", "fracweight_func"],
        "burntarea-pft":          ["burnt_area", "Burnt area fraction", None, "%/month", "c_mon_pftlayer", "burntarea_func"], # need specific function to convert units
        "evap-pft":               ["fqw", "Evapotranspiration", None, "kg m-2 s-1", "gen_mon_pft", None]
            }
    return var_dict
