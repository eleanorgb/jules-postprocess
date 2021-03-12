'''
read relevant model information
'''
import argparse
import configparser
import pandas

# #############################################################################
def parse_args():
    """
    command line parser
    """
    description = "Postprocessing jules files"
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=formatter_class)
    parser.add_argument("MIPNAME", default=None,
                        help="MIPNAME for model definitions")
    parser.add_argument("--l_testing", default=False, action="store_true",
                        help="If true: time is shortened and no output written")
    parser.add_argument("--l_backfill_missing_files", default=False,
                        action="store_true",
                        help="If true: only process files which dont exit")
    parser.add_argument("--l_jules_rose", default=False, action="store_true",
                        help="If true: code reads an MIPNAME.ini"+ \
                             "file generated automatically by a rose suite")
    MIPNAME = parser.parse_args().MIPNAME
    L_TESTING = parser.parse_args().l_testing
    L_BACKFILL_MISSING_FILES = parser.parse_args().l_backfill_missing_files
    L_JULES_ROSE = parser.parse_args().l_jules_rose

    return MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE
# #############################################################################


# #############################################################################
def config_parse_args(MIPNAME):
    """
    parser from configuration files
    """
    config = configparser.ConfigParser()
    fileini = MIPNAME+".ini"
    #print('using: '+fileini)
    config.read(fileini)
    #for key in config['OUT_INFO']: print(key)
    return config
# #############################################################################


# #############################################################################
def read_mip_info_no_rose(MIPNAME):
    """
    get mip details from mip_info.csv
    """
    MIP_INFO = pandas.read_csv("mip_info.csv", skiprows=16) # details all mips
    try:
        MIP_INFO = MIP_INFO[MIP_INFO["MIPNAME"] == MIPNAME]
    except:
        print("ADD INFORMATION for mip to file: mip_info.csv")
        raise
    MIP_INFO = MIP_INFO.set_index("MIPNAME")
    MIP_INFO = MIP_INFO.to_dict()
    return MIP_INFO
