# #############################################################################
def get_var_dict():
#               "outVar":[0.inputVars, 1.longName, 2.scale, 3.units, 4.outprofile, 5.func, 6.CMIPprofileName]
    var_dict = {
          "tas":       ["t1p5m_gb", "1.5m air Temperature", None, "K", "driving", None, None],
         # "pstar":    ["pstar", "Gridbox surface pressure", None, "Pa", "driving", None, None],
         # "wind":     ["wind", "Gridbox wind speed", None, "m s-1", "driving", None, None],
         # "rainfall": ["rainfall", "Gridbox rainfall rate", None, "kg m-2 s-1", "driving", None, None],
         # "snowfall": ["snowfall", "Gridbox snowfall rate", None, "kg m-2 s-1", "driving", None, None],
         # "precip":   ["precip", "Gridbox precipitation rate", None, "kg m-2 s-1", "driving", None, None],
         # "lwdown":   ["lw_down", "Gridbox surface downward LW radiation", None, "W m-2", "driving", None, None],
         # "swdown":   ["sw_down", "Gridbox surface downward SW radiation", None, "W m-2", "driving", None, None],
         "dtempg":     ["dtemp_g", "Gloabl temperature change in IMOGEN", None, "K", "radf", None, None],
        }
    return var_dict
# #############################################################################
