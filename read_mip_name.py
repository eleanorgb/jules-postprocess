import sys

def read_mip_name():
    MIPNAME = None
    #MIPNAME = "IMOGEN5_aj257_26" 
    #"CMIPland_HadGEM_land-hist"
    #MIPNAME = "IMOGEN5_variant_C_26"
    MIPNAME = sys.argv[1]

    return MIPNAME
