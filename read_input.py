"""
read relevant model information
"""

import argparse
import configparser
import pandas


# #############################################################################
def parse_args():
    '''
    Parses command line arguments for the JULES postprocessing script.

    Returns
    -------
    tuple
        A tuple containing:
        - MIPNAME (str): The MIPNAME for model definitions (required positional argument).
        - L_TESTING (bool): If True, time is shortened and no output is written (set with --l_testing).
        - L_BACKFILL_MISSING_FILES (bool): If True, only process files which do not exist (set with --l_backfill_missing_files).
        - L_JULES_ROSE (bool): If True, code reads an MIPNAME.ini file generated automatically by a rose suite (set with --l_jules_rose).
        - JSONMAPFILE (str or None): Optional input JSON mapping file if provided, else None.

    Command Line Usage
    ------------------
    python read_input.py MIPNAME [--l_testing] [--l_backfill_missing_files] [--l_jules_rose] [--jsonmapfile JSONMAPFILE]
    Examples
    --------
    python read_input.py my_mipname
    python read_input.py my_mipname --l_testing
    python read_input.py my_mipname --l_backfill_missing_files --jsonmapfile input.json
    ''' 

    description = "Postprocessing jules files"
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(
        description=description, formatter_class=formatter_class
    )
    parser.add_argument("MIPNAME", default=None, help="MIPNAME for model definitions")
    parser.add_argument(
        "--l_testing",
        default=False,
        action="store_true",
        help="If true: time is shortened and no output written",
    )
    parser.add_argument(
        "--l_backfill_missing_files",
        default=False,
        action="store_true",
        help="If true: only process files which dont exit",
    )
    parser.add_argument(
        "--l_jules_rose",
        default=False,
        action="store_true",
        help="If true: code reads an MIPNAME.ini"
        + "file generated automatically by a rose suite",
    )
    parser.add_argument(
        "--jsonmapfile",
        default=None,
        help="Optional input json mapping file. If not provided,"
        + "the default mapping file will be used. (variable_mapping.json)",
    )

    args = parser.parse_args()
    MIPNAME = args.MIPNAME
    L_TESTING = args.l_testing
    L_BACKFILL_MISSING_FILES = args.l_backfill_missing_files
    L_JULES_ROSE = args.l_jules_rose
    JSONMAPFILE = args.jsonmapfile

    return MIPNAME, L_TESTING, L_BACKFILL_MISSING_FILES, L_JULES_ROSE, JSONMAPFILE


# #############################################################################


# #############################################################################
def config_parse_args(MIPNAME):
    """
    parser from configuration files
    """
    config = configparser.ConfigParser()
    fileini = MIPNAME + ".ini"
    # print('using: '+fileini)
    config.read(fileini)
    # for key in config['OUT_INFO']: print(key)
    return config


# #############################################################################


# #############################################################################
def read_mip_info_no_rose(MIPNAME):
    """
    get mip details from mip_info.csv
    """
    MIP_INFO = pandas.read_csv("mip_info.csv", skiprows=16)  # details all mips
    try:
        MIP_INFO = MIP_INFO[MIP_INFO["MIPNAME"] == MIPNAME]
    except:
        print("ADD INFORMATION for mip to file: mip_info.csv")
        raise
    MIP_INFO = MIP_INFO.set_index("MIPNAME")
    MIP_INFO = MIP_INFO.to_dict()
    return MIP_INFO
