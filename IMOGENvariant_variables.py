# #############################################################################
def get_var_dict():
#               "outVar":[0.inputVars, 1.longName, 2.scale, 3.units, 4.outprofile, 5.func, 6.CMIPprofileName]
    var_dict = {
        #"cVeg_annual":             ["cv", "carbon mass in vegetation", None, "kg m-2", "Annual", None, None],
        "cVeg_tiles_annual":        ["c_veg", "carbon mass in vegetation", None, "kg m-2", "Annual", None, None],
        "cSoil_annual":            ["cs_gb", "carbon mass in model soil pool", None, "kg m-2", "Annual", None, None],
        "ns_annual":               ["ns_gb", "nitrogen mass in model soil pool", None, "kg m-2", "Annual", None, None],
        "n_demand_gb_annual":       ["n_demand_gb", "inorganic nitrogen demand", 1.0/(86400.0*360.0), "kg m-2 s-1", "Annual", None, None],
        "n_inorg_gb_annual":       ["n_inorg_gb", "gb inorganic nitrogen mass", None, "kg m-2", "Annual", None, None],
        "n_inorg_annual":          ["n_inorg", "pft inorganic nitrogen mass", None, "kg m-2", "Annual", None, None],
        "n_inorg_avail_pft_annual":["n_inorg_avail_pft", "available inorganic nitrogen mass", None, "kg m-2", "Annual", None, None],
        "tas_annual":              ["t1p5m_gb", "1.5m air Temperature", None, "K", "Annual", None, "Amon"],
        "gpp_annual":              ["gpp_gb", "Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land", None, "kg m-2 s-1", "Annual", None, None],
        "npp_nlim_tiles_annual":    ["npp_n", "Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land", 1.0/(86400.0*360.0), "kg m-2 s-1", "Annual", None, None],
        "gpp_tiles_annual":         ["gpp", "Carbon Mass Flux out of Atmosphere due to Gross Primary Production on Land", None, "kg m-2 s-1", "Annual", None, None],
        "npp_nlim_annual":         ["npp_n_gb", "Net Primary Production on Land as Carbon Mass Flux", 1.0/(86400.0*360.0), "kg m-2 s-1", "Annual", None, None],
        "rh":                      ["resp_s_to_atmos_gb", "Total Heterotrophic Respiration on Land as Carbon Mass Flux", 1.0/(86400.0*360.0), "kg m-2 s-1", "Annual", None, None],
        "npp_annual":              ["npp_gb", "Carbon Mass Flux out of Atmosphere due to Net Primary Production on Land", None, "kg m-2 s-1", "Annual", None, None],
        "exudates_annual":          ["exudates_gb", "exudates_gb", 1.0/(86400.0*360.0), "kg m-2 s-1", "Annual", None]
        "frac_annual":              ["frac", "Plant Functional Type Grid Fraction", 100., "%", "Annual", None]
        }
    return var_dict
# #############################################################################
