"""
get name for output
"""

from re import sub
from read_input import parse_args

MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE, JSONMAPFILE = parse_args()


def sort_varout_outprofile_name(var):
    """
    get name for output
    """
    varout = var
    if "isimip" in MIPNAME.lower():
        if "npp" in var and "_nlim" in var:
            varout = "npp"  # npp is different in Nitrogen/non-Nitrogen cases
            if "tiles" in var:
                varout = varout + "_tiles"
            elif "pft" in var:
                varout = varout + "-pft"
            elif "total" in var:
                varout = varout + "-total"
        elif "npp" in var:
            varout = "npp_nonlimitation"
            if "pft" in var:
                varout = varout + "-pft"
            elif "total" in var:
                varout = varout + "-total"
        if "ecoatmflux" in var and "_nlim" in var:
            varout = "ecoatmflux"  # npp is different in Nitrogen/non-Nitrogen cases
        elif "ecoatmflux" in var:
            varout = "ecoatmflux_nonlimitation"  # npp different in N/non-N cases
        if "nbp" in var and "_nlim" in var:
            varout = "nbp"  # npp is different in Nitrogen/non-Nitrogen cases
            if "pft" in var:
                varout = varout + "-pft"
            elif "total" in var:
                varout = varout + "-total"
        elif "nbp" in var:
            varout = "nbp_nonlimitation"  # npp different in N/non-N cases
            if "pft" in var:
                varout = varout + "-pft"
            elif "total" in var:
                varout = varout + "-total"
    elif "cmip" in MIPNAME.lower() and "imogen" not in MIPNAME.lower():
        if "nbp" in var and "_nlim" in var:
            varout = "nbp"  # npp different in N/non-N cases
            if "pft" in var:
                varout = varout + "-pft"
            elif "total" in var:
                varout = varout + "-total"
        elif "nbp" in var:
            varout = "nbp_nonlimitation"  # npp different in N/non-N cases
            if "pft" in var:
                varout = varout + "-pft"
            elif "total" in var:
                varout = varout + "-total"
    outprofile = "monthly"
    if "daily" in var:
        varout = sub("_daily$", "", var)
        outprofile = "daily"
    if "annual" in var:  # still not right for all npps
        outprofile = "annual"
        if "imogen" not in MIPNAME.lower():
            if not "npp" in varout:
                varout = sub("_annual$", "", var)
            if "npp_nlim" in var:
                varout = "npp"
                if "pft" in var:
                    varout = varout + "-pft"
                elif "tiles" in var:
                    varout = varout + "_tiles"
                elif "total" in var:
                    varout = varout + "-total"
            elif "npp" in var:
                varout = "npp_nonlimitation"
                if "pft" in var:
                    varout = varout + "-pft"
                elif "tiles" in var:
                    varout = varout + "_tiles"
                elif "total" in var:
                    varout = varout + "-total"
        else:
            varout = sub("_annual$", "", var)
    return varout, outprofile
