#############################################################################
# RUNNING postprocessing
###########################################################################
# add the following row to file mip_info.csv
# MIPNAME - identifier used in code to select the relevant information from mip_info.csv
#          - include 'CMIP' or 'TRENDY' or 'ISIMIP' in this name depending on what the MIP is.
#          - you might need to add for another mip, but try and use one of these if possible
#          - the differences between them are mainly in the syntax of output file 
#          - plus small formatting details and sometimes variable names
# model - model name for output file (e.g. JULES-ES, UKESM1-0-LL, HadGEM3-GC31-LL)
# suite_id - usually rose suite name (src_dir+suite_name = path to data - this is not always rose suite id)
# src_dir - base path to input data (suite_name will be added to find actual data)
# run_name - start of filename for input data (as defined in JULES output profiles)
# out_scenario - scenario for output file name (e.g. land-hist, historical, rcp2p6, c20c)
# in_scenario - scenario for input file name (e.g. S2, S3)
# start_year_in - start year for output (if not all available files are used)
#                 set to less than 0 to read all available files
# end_year_in - end year for output (if not all available files are used)
#                set to less than 0 to read all available files

On the Met Office linux system
module load scitools
python process_jules.py [-h] [--l_testing] [--l_backfill_missing_files]
                        [--l_jules_rose]
                        MIPNAME

l_testing switch: test runs for the first 10 years of data for the variables
                  in sel_diags_test.
                  
l_backfill_missing_files: this determines whether the file for that variable exists already 
                          if it exists it doesnt make a new file
                          
l_jules_rose: this is used when need the postprocessing to run within the jules rose suite.

MIPNAME: is used to look up information in mip_info which gives details of the
         model to be postprocessed. This is used outside rose

The variables to be postprocessed are defined in ***_variables.py where *** is the 
first term (before the underscore) of the MIPNAME. This is only used outside rose.

#############################################################################
# ISSUES
###################################################################
CMIP6:
   fNdep: should be monthly (Emon) and not annual 
      - needs to be changed in the output profiles and in the variable namelists
      
###################################################################
ISIMIP:
  need to get final year of future simulation
   - this can be done by adding another year to the driving by
      making links to the last year.



#############################################################################
# RUNNING MIP_CONVERT to make CMIP data for submission
###########################################################################
Getting hold of the code;
fcm co https://code.metoffice.gov.uk/svn/cdds/main/branches/dev/matthewmizielinski/r8439_1898_JULES_land_only

Activate code;
cd to r8439_1898_JULES_land_only
run “source ./setup_env_for_devel”

working_directory/
               input/
                              <suite_id>/
                                             lnd/       (daily data files)
                                             lnm/       (monthly data files)
               output/
               mip_convert.cfg

Get hold of config file; download from 
https://code.metoffice.gov.uk/trac/cdds/ticket/1898
adapt as needed; change experiment_id, variant_label, suite_id, run_bounds, etc.

Run mip_convert for monthlies
   ../r8439_1898_JULES_land_only/mip_convert/bin/mip_convert mip_convert.cfg -s lnm
check for CRITICAL messages in the output files (mip_convert_<datestamp>.log)

Run mip_convert for dailies
               mip_convert mip_convert.cfg -s lnd
check for CRITICAL messages in the output files (mip_convert_<datestamp>.log)

Check the files look OK. (last chance to review)
