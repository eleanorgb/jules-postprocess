import sys

def read_mip_name():
    MIPNAME = None
    #MIPNAME = "IMOGEN5_aj257_26" 
    #"CMIPland_HadGEM_land-hist"
    #MIPNAME = "IMOGEN5_variant_C_26"
    L_TESTING = False  # change this to TRUE for testing specific diagnostics
    L_BACKFILL_MISSING_FILES = True # change this to TRUE to just add missing files
    MIPNAME = sys.argv[1]

    return MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES
