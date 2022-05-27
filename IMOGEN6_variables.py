# #############################################################################
def get_var_dict():
#               "outVar":[0.inputVars, 1.longName, 2.scale, 3.units, 4.outprofile, 5.func, 6.CMIPprofileName]
    var_dict = {
         # "tas":      ["t1p5m_gb", "1.5m air Temperature", None, "K", "driving", None, None],
         # "pstar":    ["pstar", "Gridbox surface pressure", None, "Pa", "driving", None, None],
         # "wind":     ["wind", "Gridbox wind speed", None, "m s-1", "driving", None, None],
         # "rainfall": ["rainfall", "Gridbox rainfall rate", None, "kg m-2 s-1", "driving", None, None],
         # "snowfall": ["snowfall", "Gridbox snowfall rate", None, "kg m-2 s-1", "driving", None, None],
         # "precip":   ["precip", "Gridbox precipitation rate", None, "kg m-2 s-1", "driving", None, None],
         # "lwdown":   ["lw_down", "Gridbox surface downward LW radiation", None, "W m-2", "driving", None, None],
         # "swdown":   ["sw_down", "Gridbox surface downward SW radiation", None, "W m-2", "driving", None, None],
         "dtempg_annual":   ["dtemp_g", "Gloabl temperature change in IMOGEN", None, "K", "radf", None, None],
         "radf_annual":   ["imogen_radf", "Radiative forcing in IMOGEN", None, "W m-2", "radf", None, None]
         # "gpp_annual":      ["gpp_gb", "Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land", None, "kg m-2 s-1", "Annual", None],
         # "npp_nlim_annual": ["npp_n_gb", "N limited net Primary Production on Land as Carbon Mass Flux", 1.0/(86400.0*360.0), "kg m-2 s-1", "Annual", None, None],
         # "rh_annual":       ["resp_s_to_atmos_gb", "Total Heterotrophic Respiration on Land as Carbon Mass Flux", 1.0/(86400.0*360.0), "kg m-2 s-1", "Annual", None, None],
         # "exudates_annual": ["exudates_gb", "exudates_gb", 1.0/(86400.0*360.0), "kg m-2 s-1", "Annual", None, None],
         # "mrro_annual":     ["runoff", "total runoff", None, "kg m-2 s-1", "Annual", None, None],
         # "cs_gb_annual":    ["cs_gb", "carbon mass in model soil pool", None, "kg/m2", "Annual", None, None],
         # "ns_gb_annual":    ["ns_gb", "nitrogen mass in model soil pool", None, "kg/m2", "Annual", None, None],
         # "cv_annual":       ["cv", "carbon mass in vegetation", None, "kg m-2", "Annual", None, None],
         # "n_inorg_gb_annual":["n_inorg_gb", "gb inorganic nitrogen mass", None, "kg/m2", "Annual", None, None]
        }
    return var_dict
# #############################################################################
