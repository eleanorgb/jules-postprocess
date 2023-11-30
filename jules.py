# -*- coding: iso-8859-1 -*-
'''

A collection of functions to facilitate the use of IRIS with JULES-related NETCDF files.

Although some testing has been carried out, use at your own risk!

Crown Copyright 2014, Met Office and licensed under LGPL v3.0 or later: 
https://www.gnu.org/licenses/lgpl.html 

The unittest suite for these functions is only available internally in the 
Met Office (contact karina.williams@metoffice.gov.uk for more information).

Richard Gilham and Karina Williams, Met Office

'''

from __future__ import absolute_import, division, print_function

import os
import subprocess
import tempfile
import datetime
import copy

import iris
if int(iris.__version__[0]) == 1: # made default in Iris 2.0.0
    iris.FUTURE.netcdf_promote = True


import numpy as np
import numpy.ma as ma

# version number refers to revision number on the internal Met Office repository,
# not the mirror on the external rep.
__version__ = "last known fcm revision: " + ''.join(filter(str.isdigit, "$Rev: 31131 $"))

# 9.8.2021: updated to JULES6.1     
JULES_DIM_NAMES_DEFAULT = {
                   'land_dim_name':'land',  
                   'pft_dim_name':'pft',
                   'cpft_dim_name':'cpft',
                   'nvg_dim_name':'nvg',
                   'type_dim_name':'type',
                   'tile_dim_name':'tile',
                   'soil_dim_name':'soil',
                   'snow_dim_name':'snow',
                   'scpool_dim_name':'scpool',   
                   'scalar_dim_name':"scalar",   
                   'sclayer_dim_name':"sclayer",                                              
                   'nolevs_dim_name':"olevs",      
                   'nfarray_dim_name':'nfarray',    
                   'seed_dim_name':'seed',             
                   'bedrock_dim_name':'bedrock',
                   'p_rivers_dim_name':'p_rivers',
                   'bl_level_dim_name':'bllevel',  
                   'soilt_dim_name':'soilt',  
                   'soil_n_pool_dim_name':'snpool',  
                   'tracer_dim_name':'tracer',  
                   'ch4layer_dim_name':'ch4layer',  
                   'ch4subgrid_dim_name':'ch4subgrid',  
                  }

# 9.8.2021: updated to JULES6.1     
VAR_NO_LATLON_DIMS = (
                      'co2_ppmv', 'co2_change_ppmv', 'dtemp_o', 'fa_ocean', 'seed_rain', 
                      'frac_irr_all_tiles', 'irrtiles', 'nirrtile', 'co2_mmr',
                      'ch4_ppbv', 
                      #'rivers_lat_rp', 'rivers_lon_rp', 'rivers_sto_rp', 
                      'rfm_surfstore_rp', 
                      'rfm_substore_rp', 'rfm_flowin_rp', 'rfm_bflowin_rp'
                      'set_irrfrac_on_irrtiles'
                      )               

SOME_POSSIBLE_FRAC_NAMES = ("frac", "fractions", "field1391")

STANDARD_TILE_LIST = ['bleaftree', 'nleaftree', 'c3grass', 'c4grass', 'shrub', 'urban', 'lake', 'soil', 'ice'] 

# if two lats or lons are within this amount of each other, they are considered to be the same
# This needs to be quite large to allow for rotated pole use cases where there seems to be 
#quite large rounding errors.
TOLERANCE_LATLON = 0.5E-2 # in degrees
TOLERANCE_LATLONIND = 1.0E-2


def load(filename, constraint=None, cubename=None, conv_to_grid=True, latlon_str=None, 
         gridfile=None, missingdata=np.nan, nsmax=None, ntiles=None,
         jules_dim_names=None, latlon_flatten_function=None, coord_system=None,
         **kwargs):
    """
    Loads an iris compatible file using iris.load and returns a cubelist instance. 
    By default assumes a JULES land point format and converts to grid format. This conversion needs 
    a regular grid.
    
    Args:
    
    * filename:
        filenames or list of filenames to read from. Can contain wildcards.
    
    Kwargs:
    
    * constraint:
        iris.Constraint object
    * cubename:
        load only cubes with this name i.e. same as giving constraint=iris.Constraint(name=cubename)
    * conv_to_grid:
        conv_to_grid = True (note this is the default) assumes that the file being read in is
        land points only and puts it on a grid.
    * latlon_str:
        A tuple of strings (lat_str,lon_str) containing additional names for longitude and latitude 
        (if latlon_str=None, the standard names e.g 'longitude', 'lon' are searched for).
    * gridfile:
        The file to read longitudes and latitudes from. Ususally only useful if the input file 
        does not contain longitudes and latitudes for some reason.
    * missingdata:
        what the sea points should be set to, if points-only data is converted to a grid.
        e.g. 0.0, np.nan, ma.masked. (If missingdata=ma.masked, then the data in the resulting 
        cube is a np.MaskedArray, otherwise it is a np.ndarray)
    * jules_dim_names:
        a dictionary containing dimension names e.g. jules_dim_names = {'pft_dim_name':'my_pft'}.
        The default names in the JULES User Guide are used for any dimensions not specified here.
    * latlon_flatten_function:
        function that can be applied to flatten 2D arrays of latitudes and longitudes. 
        e.g. latlon_flatten_function=np.hstack
    * coord_system:
        an iris.coord_system object to facilitate pole rotation. Used when the longitudes and latitudes in the cube are 
        unrotated and the cube needs to be rotated to a regular grid.
        
    Any kwarg not described above is passed straight to iris.load.
    
    """
    
    if nsmax is not None:
        print('Warning: nsmax is an obsolete option')        
    
    if ntiles is not None:
        print('Warning: ntiles is an obsolete option')
    
    comb_constraint = _combine_constraints(cubename, constraint)
    
    kwargs['callback'] = _combine_callbacks(kwargs, jules_dim_names=jules_dim_names)
       
    cubelist = iris.load(filename, comb_constraint, **kwargs)
    
    cubelist = _tidy_cubelist_after_loading(
                   cubelist,   
                   filename=filename,
                   conv_to_grid=conv_to_grid,
                   latlon_str=latlon_str, 
                   gridfile=gridfile, 
                   missingdata=missingdata, 
                   jules_dim_names=jules_dim_names, 
                   latlon_flatten_function=latlon_flatten_function,
                   coord_system=coord_system)
        
    return cubelist

def load_cube(filename, constraint=None, cubename=None, conv_to_grid=True, latlon_str=None, 
         gridfile=None, missingdata=np.nan, nsmax=None, ntiles=None,
         jules_dim_names=None, latlon_flatten_function=None, coord_system=None,
          **kwargs):
    """     
    Loads an iris compatible file using iris.load_cube and returns a cube instance. 
    By default assumes a JULES land point format and converts to grid format. This conversion needs 
    a regular grid.
    
    Args:
    
    * filename:
        filenames or list of filenames to read from. Can contain wildcards.
    
    Kwargs:
    
    * constraint:
        iris.Constraint object
    * cubename:
        load only cubes with this name i.e. same as giving constraint=iris.Constraint(name=cubename)
    * conv_to_grid:
        conv_to_grid = True (note this is the default) assumes that the file being read in is
        land points only and puts it on a grid.
    * latlon_str:
        A tuple of strings (lat_str,lon_str) containing additional names for longitude and latitude 
        (if latlon_str=None, the standard names e.g 'longitude', 'lon' are searched for).
    * gridfile:
        The file to read longitudes and latitudes from. Ususally only useful if the input file 
        does not contain longitudes and latitudes for some reason.
    * missingdata:
        what the sea points should be set to, if points-only data is converted to a grid.
        e.g. 0.0, np.nan, ma.masked. (If missingdata=ma.masked, then the data in the resulting 
        cube is a np.MaskedArray, otherwise it is a np.ndarray)
    * jules_dim_names:
        a dictionary containing dimension names e.g. jules_dim_names = {'pft_dim_name':'my_pft'}.
        The default names in the JULES User Guide are used for any dimensions not specified here.
    * latlon_flatten_function:
        function that can be applied to flatten 2D arrays of latitudes and longitudes. 
        e.g. latlon_flatten_function=np.hstack
    * coord_system:
        an iris.coord_system object to facilitate pole rotation. Used when the longitudes and latitudes in the cube are 
        unrotated and the cube needs to be rotated to a regular grid.
        
    Any kwarg not described above is passed straight to iris.load_cube.
    
    """
    
    if nsmax is not None:
        print('Warning: nsmax is an obsolete option')       
    
    if ntiles is not None:
        print('Warning: ntiles is an obsolete option')
        
    comb_constraint = _combine_constraints(cubename, constraint)
    
    kwargs['callback'] = _combine_callbacks(kwargs, jules_dim_names=jules_dim_names)
    
    cube = iris.load_cube(filename, comb_constraint, **kwargs)
    
    cubelist = _tidy_cubelist_after_loading(
                   iris.cube.CubeList([cube]),  
                   filename=filename,
                   conv_to_grid=conv_to_grid,
                   latlon_str=latlon_str, 
                   gridfile=gridfile, 
                   missingdata=missingdata, 
                   jules_dim_names=jules_dim_names, 
                   latlon_flatten_function=latlon_flatten_function,
                   coord_system=coord_system)
   
    if len(cubelist) != 1:
        raise UserWarning('expecting a cubelist of length 1, not length '+str(len(cubelist)))          
  
    return cubelist[0]

def _tidy_cubelist_after_loading(cubelist, filename=None, 
         conv_to_grid=True, latlon_str=None, 
         gridfile=None, missingdata=np.nan,
         jules_dim_names=None, latlon_flatten_function=None, coord_system=None):    
    """
    Takes the cubelist as read in by Iris and tidies it up e.g. by converting from list of points to grid.
    
    Args:
    
    * cubelist:
        cubelist to tidy.
    
    Kwargs:
    
    * filename:
        filenames or list of filenames to read from. Can contain wildcards.
    * constraint:
        iris.Constraint object
    * cubename:
        load only cubes with this name i.e. same as giving constraint=iris.Constraint(name=cubename)
    * conv_to_grid:
        conv_to_grid = True (note this is the default) assumes that the file being read in is
        land points only and puts it on a grid.
    * latlon_str:
        A tuple of strings (lat_str,lon_str) containing additional names for longitude and latitude 
        (if latlon_str=None, the standard names e.g 'longitude', 'lon' are searched for).
    * gridfile:
        The file to read longitudes and latitudes from. Ususally only useful if the input file 
        does not contain longitudes and latitudes for some reason.
    * missingdata:
        what the sea points should be set to, if points-only data is converted to a grid.
        e.g. 0.0, np.nan, ma.masked. (If missingdata=ma.masked, then the data in the resulting 
        cube is a np.MaskedArray, otherwise it is a np.ndarray)
    * jules_dim_names:
        a dictionary containing dimension names e.g. jules_dim_names = {'pft_dim_name':'my_pft'}.
        The default names in the JULES User Guide are used for any dimensions not specified here.
    * latlon_flatten_function:
        function that can be applied to flatten 2D arrays of latitudes and longitudes. 
        e.g. latlon_flatten_function=np.hstack
    * coord_system:
        an iris.coord_system object to facilitate pole rotation. Used when the longitudes and latitudes in the cube are 
        unrotated and the cube needs to be rotated to a regular grid.
    
    """
    
    if not cubelist:
        raise UserWarning('cubelist has no cubes')      

    if conv_to_grid:
        # have a list of points to be converted to a grid
    
        gridcubelist = iris.cube.CubeList()

        while cubelist:

            temp_cube = cubelist.pop()
            
            if temp_cube.core_data().dtype.kind in ['S', 'a'] or temp_cube.var_name in VAR_NO_LATLON_DIMS:
                    gridcubelist.append(temp_cube)
            else:
                try:    
                
                    gridcubelist.append(
                        points_to_grid(temp_cube, latlon_str=latlon_str, gridfile=gridfile,
                        missingdata=missingdata, jules_dim_names=jules_dim_names,
                        latlon_flatten_function=latlon_flatten_function, coord_system=coord_system))
                        
                except _LonlatException: 
                
                    gridfile = filename
                        
                    gridcubelist.append(
                        points_to_grid(temp_cube, latlon_str=latlon_str, gridfile=gridfile,
                        missingdata=missingdata, jules_dim_names=jules_dim_names,
                        latlon_flatten_function=latlon_flatten_function, coord_system=coord_system))      
                                      
        return gridcubelist

    else:
        # no points-to-grid conversion necessary, either because want to keep as points (e.g.
        # to save computer memory) or because data is already on a lat,lon grid
        
        for cube in cubelist:
        
            try:
                # check whether there are lat, lon coords in cube
            
                all_coord_names = [ coord.name() for coord in cube.coords() ]
                
                (_lat_str, _lon_str) = _get_latlon_str(all_coord_names, latlon_str=latlon_str)
                
                dim_coord_names = [ coord.name() for coord in cube.coords(dim_coords=True) ]
                                                                        
                if ( len(cube.coord(_lat_str).points) > 1 and
                     len(cube.coord(_lon_str).points) > 1 and
                     not ( _lat_str in dim_coord_names )  and
                     not ( _lon_str in dim_coord_names )  
                   ):
                     
                    # this special case is needed for e.g. JULES 3.4.1 output files with 
                    # land_only=F, where the lat and lon is in the files, but given 
                    # as a 2D array. Iris reads both the lat and lon as a 2-dim aux coords, 
                    # so want to reduce these to 1-dim and add as dimcoords.                   
                   
                    try:
                        (latcoord, loncoord, latpts, lonpts) = _parse_grid_file(
                                   latlon_str=latlon_str, gridfile=filename, latlon_flatten_function=np.hstack,
                                   coord_system=coord_system)
                                                
                        # don't want to get rid of the auxcoord lat and lon
                        # until we know whether the new dimcoords are ok
                        # so give the new dimcoords a temporary name 
                        
                        latcoord.rename('temporary_latitude')
                        loncoord.rename('temporary_longitude')
                        
                        cube.add_dim_coord(latcoord,  cube.ndim-2) # note this assumes the position of lat and lon
                        cube.add_dim_coord(loncoord,  cube.ndim-1)
                                                                                                
                        cube.remove_coord(_lat_str)
                        cube.remove_coord(_lon_str)
                                                                        
                        cube.coord('temporary_latitude').rename('latitude')
                        cube.coord('temporary_longitude').rename('longitude')
                        
                    except (ValueError, IndexError):
                    
                        pass                          
                
            except _LonlatException:
                
                
                try:
                    # this one is handy for old-style JULES land_only=F files, because they're
                    # not set up so that lat and lon are recognised as coords and instead are
                    # treated as separate variables
                    
                    (latcoord, loncoord, latpts, lonpts) = _parse_grid_file(
                            latlon_str=latlon_str, gridfile=filename, latlon_flatten_function=np.hstack,
                            coord_system=coord_system)
                    
                    cube.add_dim_coord(latcoord,  cube.ndim-2) # note this assumes the position of lat and lon
                    cube.add_dim_coord(loncoord,  cube.ndim-1)
                    
                except (ValueError, IndexError, UserWarning) :    
                
                    pass                    
        
        return cubelist


def load_dump(*args, **kwargs):
    """Loads a JULES dump file with jules.load and returns an iris cubelist instance, by calling jules.load 
    and then setting the long_name of each variable to be the same as it's var_name.
    Takes same arguments as jules.load.
    
    Args:
    
    * filename:
        filenames or list of filenames to read from. Can contain wildcards.
    
    Kwargs:
    
    * constraint:
        iris.Constraint object
    * cubename:
        load only cubes with this name i.e. same as giving constraint=iris.Constraint(name=cubename)
    * conv_to_grid:
        conv_to_grid = True (note this is the default) assumes that the file being read in is
        land points only and puts it on a grid.
    * latlon_str:
        A tuple of strings (lat_str,lon_str) containing additional names for longitude and latitude 
        (if latlon_str=None, the standard names e.g 'longitude', 'lon' are searched for).
    * gridfile:
        The file to read longitudes and latitudes from. Usually only useful if the input file 
        does not contain longitudes and latitudes for some reason.
    * missingdata:
        what the sea points should be set to, if points-only data is converted to a grid.
        e.g. 0.0, np.nan, ma.masked. (If missingdata=ma.masked, then the data in the resulting 
        cube is a np.MaskedArray, otherwise it is a np.ndarray)
    * jules_dim_names:
        a dictionary containing dimension names e.g. jules_dim_names = {'pft_dim_name':'my_pft'}.
        The default names in the JULES User Guide are used for any dimensions not specified here.
    * latlon_flatten_function:
        function that can be applied to flatten 2D arrays of latitudes and longitudes. 
        e.g. latlon_flatten_function=np.hstack
    * coord_system:
        an iris.coord_system object to facilitate pole rotation. Used when the longitudes and latitudes in the cube are 
        unrotated and the cube needs to be rotated to a regular grid.
        
    Any kwarg not described above is passed straight to iris.load.
    
    """
        
    cubelist = load(*args, **kwargs)
                          
    for cube in cubelist:
        cube.long_name = getattr(cube, 'var_name', None)
    
    return cubelist
         
         
def _combine_constraints(cubename, constraint):
    """combines the cubename constraint and any other specified constraint"""
    
    if cubename is None:
        combined_constraint  = constraint
    elif constraint is None:
        combined_constraint  = cubename
    else:
        combined_constraint  = iris.Constraint(name=cubename) & constraint
        
    return combined_constraint   


def _combine_callbacks(load_kwargs, jules_dim_names=None):
    """combines the callback that adds recognised dim coords and any other specified callback"""    
    
    jules_recognised_dim_names = {'Psuedo':'tile'} # this is needed for some of the older files. Note spelling.
    
    for key,val in JULES_DIM_NAMES_DEFAULT.items():
       jules_recognised_dim_names[val] = val
       if jules_dim_names is not None:      
           if key in jules_dim_names:
               jules_recognised_dim_names[val] = jules_dim_names[key]
               jules_recognised_dim_names[jules_dim_names[key]] = jules_dim_names[key]
               
    if jules_dim_names is not None:      
       for key,val in jules_dim_names.items(): 
           if key not in jules_recognised_dim_names:
               jules_recognised_dim_names[val] = val

    def _add_recognised_dim_coords(cube, field, filename):
        '''
        Adds a dim coord to each anonymous dimension with a recognised name in the netCDF file.
        Values of the dim coord are integers starting from 0.
        '''
        
        for i,dim_length in enumerate(cube.shape): 
            if not cube.coords(dim_coords=True, contains_dimension=i):
                if field.dimensions[i] in jules_recognised_dim_names:
                    new_dim_coord = iris.coords.DimCoord(list(range(dim_length)), long_name=jules_recognised_dim_names[field.dimensions[i]])
                    cube.add_dim_coord(new_dim_coord, (i,))
    
    if 'callback' in load_kwargs.keys():
        user_callback = load_kwargs['callback']  
    else:
        user_callback = None
        
    def comb_callback(cube, field, filename):
        _callback_discard_empty(cube, field, filename)
        _add_recognised_dim_coords(cube, field, filename)
        
        if user_callback is not None: 
            user_callback(cube, field, filename)
            
    return comb_callback


def _callback_discard_empty(cube, field, filename):
    '''
    Callback to discard any empty cubes i.e. cubes where one or more dimensions have zero length
    '''
    if 0 in cube.shape:
         raise iris.exceptions.IgnoreCubeException
         
         
def load_pftfrac_cube(*args, **kwargs):
    """Deprecated. Please use load_frac_cube instead."""

    raise UserWarning("Deprecated function. Please use load_frac_cube instead.")
    #load_frac_cube(*args, **kwargs)


def load_frac_cube(filename, frac_names=None, tilelist=None,
                   **kwargs):
    """Loads a JULES tile fractions mask file with jules.load_cube and returns it as an iris cube. 
    Take same arguments as jules.load_cube except gridfile is required, not optional) 
    and also there are two additional optional arguments: frac_names and tilelist.
    
    Args:
    
    * filename:
        filenames or list of filenames to read from. Can contain wildcards.
    
    Kwargs:
    
    * frac_names: 
        list of possible variable names for the variable storing the tile fractions
        within the netcdf file
    * tilelist:    
        list of tile names to be added to the cube as an aux coord.
    * constraint:
        iris.Constraint object
    * cubename:
        load only cubes with this name i.e. same as giving constraint=iris.Constraint(name=cubename)
    * conv_to_grid:
        conv_to_grid = True (note this is the default) assumes that the file being read in is
        land points only and puts it on a grid.
    * latlon_str:
        A tuple of strings (lat_str,lon_str) containing additional names for longitude and latitude 
        (if latlon_str=None, the standard names e.g 'longitude', 'lon' are searched for).
    * gridfile:
        The file to read longitudes and latitudes from. Ususally only useful if the input file 
        does not contain longitudes and latitudes for some reason.
    * missingdata:
        what the sea points should be set to, if points-only data is converted to a grid.
        e.g. 0.0, np.nan, ma.masked. (If missingdata=ma.masked, then the data in the resulting 
        cube is a np.MaskedArray, otherwise it is a np.ndarray)
    * jules_dim_names:
        a dictionary containing dimension names e.g. jules_dim_names = {'pft_dim_name':'my_pft'}.
        The default names in the JULES User Guide are used for any dimensions not specified here.
    * latlon_flatten_function:
        function that can be applied to flatten 2D arrays of latitudes and longitudes. 
        e.g. latlon_flatten_function=np.hstack
    * coord_system:
        an iris.coord_system object to facilitate pole rotation. Used when the longitudes and latitudes in the cube are 
        unrotated and the cube needs to be rotated to a regular grid.
        
    Any kwarg not described above is passed straight to iris.load_cube.
    
    """    

    if not "gridfile" in kwargs:
        raise UserWarning('Error: Grid file needed to load a frac file.')         
        
    if frac_names is None:  
        frac_names = SOME_POSSIBLE_FRAC_NAMES
          
    # put the possible names for the variable containing the fraction into lowercase      
    # don't need _nice_lower here because none of x should be None anyway
    frac_names = [x.lower() for x in frac_names]
       
    ##ensure that only the variable containing the tile fractions are read in to the cube.
    frac_var = iris.Constraint(cube_func=lambda cube: _nice_lower(cube.var_name) in frac_names)
    
    cube = load_cube(filename, constraint=frac_var, **kwargs)
    
    if tilelist is not None:
        tilecoord = iris.coords.AuxCoord(tilelist, long_name="frac_name")
        cube.add_aux_coord(tilecoord, 0)    
        
    return cube


def apply_pft_mask(*args, **kwargs):
    """Deprecated. Please use apply_frac_mask instead"""
    
    raise UserWarning("Deprecated function. Please use apply_frac_mask instead.")
    #apply_frac_mask(*args, **kwargs)
    

def apply_frac_mask(cube, fraccube, targetfrac, threshold, missingdata=np.nan, pos_mask=True):
    """Returns a modified IRIS cube where some data is multiplied by a 'missingdata' value,
    according to whether the corresponding element in a tile fraction cube is above or below a threshold
    value.
    
    Args:
    
    * cube:
        Cube of data to be masked
    * fraccube:
        Cube of tile fractions. The fraction dimension should have the name 'frac_name'.
    * targetfrac:
        The name of the fraction to use for masking.
    * threshold:
        threshold value
    
    Kwargs:
    
    * missingdata:
        what to multiply the unwanted data values by. Recommended values are np.nan (default) or 0.0.
    * pos_mask:
        If pos_mask=True, data less than the threshold is multiplied by the 'missingdata' value. The rest is multiplied by 1.
        If pos_mask=False, data greater than the threshold is multiplied by the 'missingdata' value. The rest is multiplied by 1.
    
    """
    
    maskcube = fraccube.extract(iris.Constraint(frac_name=targetfrac))
    
    if pos_mask:
        maskcube.data[maskcube.data <  threshold] = missingdata
        maskcube.data[maskcube.data >= threshold] = 1.0
    else:
        maskcube.data[maskcube.data >  threshold] = missingdata
        maskcube.data[maskcube.data <= threshold] = 1.0
        
    #Deal with the possibility of gridded and point cubes being masked        
    if len(maskcube.data.shape) == 2:
        cube.data[..., :, :] *= maskcube.data[:, :]
    elif len(maskcube.data.shape) == 1:
        cube.data[...,  :]   *= maskcube.data[:]
    else:
        raise ValueError("Unexpected maskcube dimensionality") 
    
    return cube


def points_to_grid(pcube, gridfile=None, missingdata=np.nan, latlon_str=None, latlon_flatten_function=None, coord_system=None, **kwargs):
    """
    Takes an IRIS cube with land points only data and returns the data as a gridded cube. 
    Designed to work with JULES output. The output grid should be a regular grid.
    Some functionality for rotated poles has been added.
    
    Args:
    
    * pcube:
        Iris cube before being put on a grid i.e. containing a list of land point data.
    
    Kwargs:
    
    * constraint:
        iris.Constraint object
    * cubename:
        load only cubes with this name i.e. same as giving constraint=iris.Constraint(name=cubename)
    * conv_to_grid:
        conv_to_grid = True (note this is the default) assumes that the file being read in is
        land points only and puts it on a grid.
    * latlon_str:
        A tuple of strings (lat_str,lon_str) containing additional names for longitude and latitude 
        (if latlon_str=None, the standard names e.g 'longitude', 'lon' are searched for).
    * gridfile:
        The file to read longitudes and latitudes from. Ususally only useful if the input file 
        does not contain longitudes and latitudes for some reason.
    * missingdata:
        what the sea points should be set to, if points-only data is converted to a grid.
        e.g. 0.0, np.nan, ma.masked. (If missingdata=ma.masked, then the data in the resulting 
        cube is a np.MaskedArray, otherwise it is a np.ndarray)
    * jules_dim_names:
        a dictionary containing dimension names e.g. jules_dim_names = {'pft_dim_name':'my_pft'}.
        The default names in the JULES User Guide are used for any dimensions not specified here.
    * latlon_flatten_function:
        function that can be applied to flatten 2D arrays of latitudes and longitudes. 
        e.g. latlon_flatten_function=np.hstack
    * coord_system:
        an iris.coord_system object to facilitate pole rotation. Used when the longitudes and latitudes in the cube are 
        unrotated and the cube needs to be rotated to a regular grid.
    
    """

    if 'nsmax' in kwargs:
        print('Warning: nsmax is an obsolete option')        
    
    if 'ntiles' in kwargs:
        print('Warning: ntiles is an obsolete option')

    # assume last element is the land-points dimension
    # start list of pcube.data dimensions that should not be copied straight to new cube
    not_to_be_copied_over = [len(pcube.core_data().shape)-1]

    # list of locations of dim coords which have names (i.e. are not anonymous) in pcube.data
    named_coord_loc = []
    for coord in pcube.coords(dim_coords=True):
        named_coord_loc.extend(pcube.coord_dims(coord))
    
    # list of locations of anonymous dim coords in pcube.data
    unnamed_coord_loc = []
    for i in range(len(pcube.data.shape)):
        if not i in named_coord_loc:
            unnamed_coord_loc.append(i)
                           
    # list of the names of all named coord (both dim and aux)                       
    all_coord_names = [ coord.name() for coord in pcube.coords() ]

    # in a special case, want to get rid of one of the anonymous dimensions right away:
    if (
       ( len(pcube.core_data().shape) >= 2) and 
       ( pcube.core_data().shape[-2] == 1 ) and 
       ( len(pcube.core_data().shape)-2 in unnamed_coord_loc ) and 
       ( 'latitude' in [coord.name() for coord in pcube.aux_coords]) and 
       ( len(pcube.core_data().shape)-2 in pcube.coord_dims(pcube.coords('latitude')[0]) )
       ):
        data = np.mean(pcube.core_data(), axis=-2) 
        not_to_be_copied_over.append(len(pcube.core_data().shape)-2)
    elif ( # as above, except dimension we want to get rid of is a dim_coord called 'y'
       ( len(pcube.core_data().shape) >= 2) and 
       ( pcube.core_data().shape[-2] == 1 ) and 
       ( 'y' in [coord.name() for coord in pcube.coords(dim_coords=True)] ) and 
       ( pcube.coord_dims(pcube.coord('y')) == (len(pcube.core_data().shape)-2 ,) ) and 
       ( 'latitude' in [coord.name() for coord in pcube.aux_coords]) and 
       ( len(pcube.core_data().shape)-2 in pcube.coord_dims(pcube.coords('latitude')[0]) )
       ):
        data = np.mean(pcube.core_data(), axis=-2) 
        not_to_be_copied_over.append(len(pcube.core_data().shape)-2)
    else:
        data = pcube.core_data().copy()  

    # list of dim coords that are named 
    reduced_dim_coords = [ coord for coord in pcube.coords(dim_coords=True) 
                           if not set(pcube.coord_dims(coord)).intersection(set(not_to_be_copied_over)) ] 
    
    #Set up an empty list to hold iris.coord instances for each named dim coord and its location in pcube.data
    cubedimlist = [] 
    
    #Set up each named 'dim coord and dim' instance we want to copy over and append it onto cubedimlist
    for coord in reduced_dim_coords:
        cubedimlist.append((coord, pcube.coord_dims(coord)[0]))
 
    # append the unnamed ones we want too, and give them names
    idim = 1
    n_unknown_dim = len([i for i in unnamed_coord_loc if not i in not_to_be_copied_over])
    for i in unnamed_coord_loc:
        if not i in not_to_be_copied_over:
            ntiles        = pcube.core_data().shape[i]
            long_name = 'tiles_or_layers_dim_'+ str(idim)
                
            tilecoord     = iris.coords.DimCoord(list(range(ntiles)), long_name=long_name) # units='no_unit' is for strings only
            cubedimlist.append((tilecoord, i))
            idim         += 1
           
    # now make the new longitude and latitude coords
    if gridfile is None:
        (_lat_str, _lon_str) = _get_latlon_str(all_coord_names, latlon_str=latlon_str)
        
        for coord_str in [_lat_str, _lon_str]: 
            if np.ma.is_masked(pcube.coord(coord_str).points):
                print(pcube.coord(coord_str))
                raise UserWarning('This function does not work when coordinate is masked. Unmask beforehand.')
        
        #flatten for the point case (treated as 1x1)
        latpts   = pcube.coord(_lat_str).points.flatten()
        lonpts   = pcube.coord(_lon_str).points.flatten()

        #Determine whether lat and lon points need modifiying, eg for a rotated pole
        if isinstance(type(coord_system), type(iris.coord_systems.RotatedGeogCS)):
            lonpts, latpts = iris.analysis.cartography.rotate_pole(lonpts, latpts, 
                                                                   coord_system.grid_north_pole_longitude, 
                                                                   coord_system.grid_north_pole_latitude)
            
        # case where lat,lon coords are 2D and neither dim is 1
        if ( ( len(pcube.coord(_lat_str).points.shape) == 2 ) & ( not 1 in pcube.coord(_lat_str).points.shape )
           & ( len(pcube.coord(_lon_str).points.shape) == 2 ) & ( not 1 in pcube.coord(_lon_str).points.shape ) ):

            # check first whether data looks already gridded
            test_latcoord = _make_spatial_coord(latpts, 'latitude', coord_system=coord_system)
            test_loncoord = _make_spatial_coord(lonpts, 'longitude', coord_system=coord_system)

            if test_latcoord.points.shape == np.unique(latpts).shape:
                if test_loncoord.points.shape == np.unique(lonpts).shape:
                    if np.allclose(test_latcoord.points, np.unique(latpts)):
                        if np.allclose(test_loncoord.points, np.unique(lonpts)):
                            raise _MaybeAlreadyGriddedException

            #Flatten to a point pseudocube
            raise UserWarning("UNTESTED CODE- Multipoint pseudo-gridded data- attempting to flattening and regridding.")
            pcube    = _flatten_pseudogrid(pcube)
            latpts   = pcube.coord(_lat_str).points.flatten()
            lonpts   = pcube.coord(_lon_str).points.flatten()
           
        latcoord = _make_spatial_coord(latpts, 'latitude',  coord_system=coord_system)
        loncoord = _make_spatial_coord(lonpts, 'longitude', coord_system=coord_system)
        
    else:
        (latcoord, loncoord, latpts, lonpts) = _parse_grid_file(
                latlon_str=latlon_str, gridfile=gridfile, latlon_flatten_function=latlon_flatten_function,
                coord_system=coord_system)
    
    #Get the lat/lon value for each index.
    #This is where we need to work in rotated pole world.
    if len(latcoord.points) == 1:
        latind_int = np.zeros(len(latpts), dtype=int)
    else:
        latind_float = (latpts - latcoord.points[0]) / (latcoord.points[1] - latcoord.points[0])
        latind_int = np.rint(latind_float).astype(int)
        
        if not np.allclose( latind_int.astype('float'), latind_float, atol=TOLERANCE_LATLONIND ):
            raise UserWarning('lats not interpreted correctly')
        
    if len(loncoord.points) == 1:
        lonind_int = np.zeros(len(lonpts), dtype=int)
    else:    
        lonind_float = (lonpts - loncoord.points[0]) / (loncoord.points[1] - loncoord.points[0])
        lonind_int = np.rint(lonind_float).astype(int)
        
        if not np.allclose( lonind_int.astype('float'), lonind_float, atol=TOLERANCE_LATLONIND ):
            raise UserWarning('lons not interpreted correctly')
    
    cubedimlist.append((latcoord, len(cubedimlist)))
    cubedimlist.append((loncoord, len(cubedimlist))) 

    griddata_shape = [dim[0].shape[0] for dim in cubedimlist]
    
    if missingdata is ma.masked: # n.b. need the 'is' here rather than ==
        griddata = ma.masked_all(griddata_shape)
        data_for_copying = data
    else:    
        griddata    = np.zeros(griddata_shape)
        griddata[:] = missingdata
        if np.ma.is_masked(data):
            data_for_copying = data.filled(missingdata)
        else:
            data_for_copying = data

    if griddata.shape[:-2] != data.shape[:-1]:
        print(cubedimlist)
        print(griddata.shape)
        print(data.shape)
        raise UserWarning("either griddata or data has an unexpected shape")

    if data.shape[-1] - 1 > len(latind_int):
        raise _MaybeAlreadyGriddedException('Latitudes are not consistent with the data.')
    
    for gridpt in range(data.shape[-1]):
        griddata[..., latind_int[gridpt], lonind_int[gridpt]] = data_for_copying[..., gridpt]

    # now to copy over any aux coords we want
    # n.b. do not want to copy over any aux coords which depend on one we're getting rid of
    aux_coords_and_dims = []         
    for coord in pcube.aux_coords:
        dim_tuple = pcube.coord_dims(coord)
        if set(dim_tuple).intersection(set(not_to_be_copied_over)) == set([]):
            aux_coords_and_dims.append((coord, dim_tuple))     

    #Make the cube and go home
    cube = iris.cube.Cube(griddata, 
        dim_coords_and_dims = cubedimlist,
        aux_coords_and_dims = aux_coords_and_dims,
        standard_name       = pcube.standard_name, 
        long_name           = pcube.long_name, 
        var_name            = pcube.var_name, 
        units               = pcube.units, 
        cell_methods        = pcube.cell_methods,
        attributes          = pcube.attributes
        )
        #aux_factories       = pcube.aux_factories, # need to deal with this in a better way
                           
    for name in ['latitude', 'longitude']:
        if [coord.standard_name for coord in cube.coords(dim_coords=True)].count(name) > 1:
            raise _MaybeAlreadyGriddedException('a cube has ended up with more than one dim coord with '
                              'the standard_name "' + name + '"')
    
                              
    return cube


def _parse_grid_file(latlon_str=None, gridfile='dummy', latlon_flatten_function=None, coord_system=None):
    """Loads a JULES land point-only grid definition file and returns a tuple of iris coordinates for lat and lon"""

    cubelist = iris.load(gridfile)
    
    try:
        var_name_list = [cube.var_name for cube in cubelist]
        (_lat_str, _lon_str) = _get_latlon_str(var_name_list, latlon_str=latlon_str)   
        
        latdim = var_name_list.index(_lat_str)
        londim = var_name_list.index(_lon_str)
        
        if len(cubelist[latdim].data.shape) == 1 :
            latdata = cubelist[latdim].data
        elif latlon_flatten_function != None:
            try:
                latdata = latlon_flatten_function(cubelist[latdim].data)
            except (ValueError, IndexError):
                print('this latlon_flatten_function did not work')
                raise _LonlatException 
        else:    
            print('warning: latitude coord found but not list of points. Maybe think about setting a latlon_flatten_function?')
            raise _LonlatException
            
        if len(cubelist[londim].data.shape) == 1 :
            londata = cubelist[londim].data
        elif latlon_flatten_function != None:
            try:
                londata = latlon_flatten_function(cubelist[londim].data)
            except (ValueError, IndexError):
                print('this latlon_flatten_function did not work')
                raise _LonlatException
        else:
            print('warning: longitude coord found but not list of points. Maybe think about setting a latlon_flatten_function?')
            raise _LonlatException        
        
        latpts = latdata
        lonpts = londata
        
    except _LonlatException:
        cube0 = cubelist[0]    #fixme: generalise so that it checks each cube in cubelist 
        coord_name_list = [ coord.name() for coord in cube0.coords() ]
        
        (_lat_str, _lon_str) = _get_latlon_str(coord_name_list, latlon_str=latlon_str)
        
        # flatten() is needed because might have a (1,npoints) array or lat, lon as auxcoords
        # that are in the process of being reduced to dim coords
            
        latpts = cube0.coord(_lat_str).points.flatten()
        lonpts = cube0.coord(_lon_str).points.flatten()
    
    #Determine whether lat and lon points need modifiying, eg for a rotated pole
    if isinstance(type(coord_system), type(iris.coord_systems.RotatedGeogCS)):   
        lonpts, latpts = iris.analysis.cartography.rotate_pole(lonpts, latpts, 
                                                                coord_system.grid_north_pole_longitude, 
                                                                coord_system.grid_north_pole_latitude)
                                                                   
    latcoord = _make_spatial_coord(latpts, 'latitude',  coord_system=coord_system)
    loncoord = _make_spatial_coord(lonpts, 'longitude', coord_system=coord_system)
    
    return (latcoord, loncoord, latpts, lonpts)


def _make_spatial_coord(coord, name, units='degrees', tolerance = TOLERANCE_LATLON, coord_system=None):
    """Takes an unordered list of grid box lat/lons and generates an iris coords instance assuming a regular grid
    and assuming that two points are next to each other in the grid
    """
    coordset = np.sort(np.unique(coord))
    
    if len(coordset) == 1: 
        step      = 0.0
    elif len(coordset) == 2:
        step      = abs(coordset[1] - coordset[0]) 
    else:
        steparr   = np.diff( coordset )
        steparr   = steparr[np.nonzero( steparr > tolerance )]
        stepmin   = np.min( steparr ) 
        step      = np.mean( steparr[np.nonzero( steparr < stepmin + tolerance )] )
        
        normalised_step = steparr / step
        if not np.allclose( np.rint(normalised_step).astype('float'), normalised_step, atol=tolerance):
            raise Exception('problem converting this coord - check grid is regular and that two points are next to each other')
        
    if step < tolerance:
        coordlist = coordset[0] 
    else:
        # old method using arange: 
        #coordlist = np.arange(min(coordset), max(coordset) + step, step)
        # new method using linspace:
        coordlist = np.linspace(min(coordset), max(coordset), np.rint((max(coordset) - min(coordset)) / step).astype(int) + 1)
        
    iriscoord = iris.coords.DimCoord(coordlist, standard_name=name, units=units, coord_system=coord_system)
    return iriscoord

def save(source, target, lsmask=None, missingdata=np.nan, lsmask_missingdata_str='nan', 
         landpointsonly=True, latlon_str=None, data2D = False, user_def_start_corner=None, **kwargs):
    """
    Saves an iris cube/cubelist/cube sequence. Requires CDO. Has the option of just outputting
    landpoints, in a format suitable for JULES i.e. data stored in a 1 x npoints array where npoints
    is the number of land points.

    Args:
    
     * source:
         A iris cube/cubelist/cube sequence
     * target:
         Output filename.
     
    Kwargs:
    
     * lsmask:
         The land-sea mask as a cube
     * landpointsonly:  
         landpointsonly = True means gridded data will be converted to a list of land points before outputting.
     * missingdata:
         What the missing data in the output array will be set to.
     * lsmask_missingdata_str:
         Specifies how sea points are labelled in lsmask. Can be either 'nan' or 'zeros'. For an array
               of Trues and Falses, where False labels sea points, pick 'zeros'.
     * latlon_str:
         A tuple (lat_str,lon_str) containing the variable names for longitude and latitude in lsmask.
     * data2D:
         If data2D=True, there will be two spatial dimensions with sizes (1,npoints), rather than just one with size (npoints).  
         This is a similar format to the JULES output files.
     * user_def_start_corner:
         When converting to a list of points, a particular order of points is chosen as default. This keyword allows 
         this choice to be overridden. See grid_to_points for more information.
        
    """
    print('jules.save: just entering function at '+str(datetime.datetime.now()))

    # force netcdf output
    saver           = iris.fileformats.netcdf.save
    kwargs          = kwargs.copy()
    kwargs['saver'] = saver
    
    if landpointsonly:
        if lsmask is None:
            raise UserWarning('if landpointsonly=True, then lsmask should be given')
            
        source_for_saving = iris.cube.CubeList()
        
        # Single cube?
        if isinstance(source, iris.cube.Cube):
            try:
                source_for_saving.append(grid_to_points(source, lsmask,
                    lsmask_missingdata_str=lsmask_missingdata_str, latlon_str=latlon_str,
                    data2D =data2D, user_def_start_corner=user_def_start_corner))
            except _LonlatException: 
                source_for_saving.append(source)
       
        # CubeList or sequence of cubes?
        elif isinstance(source, iris.cube.CubeList) or \
          (isinstance(source, (list, tuple)) and all([type(i)==iris.cube.Cube for i in source])):
           
            for cube in source:
                try:
                    source_for_saving.append(grid_to_points(cube, lsmask,
                        lsmask_missingdata_str=lsmask_missingdata_str, latlon_str=latlon_str,
                        data2D =data2D, user_def_start_corner=user_def_start_corner))
                except _LonlatException: 
                    source_for_saving.append(cube)
        else:
            raise ValueError("Cannot save; non Cube found in source")
       
    else:
        source_for_saving = source
    
    tempfilename = _generate_temp_file_name() # we need a temporary file for the next bit
    
    print('jules.save: about to enter try/except at '+str(datetime.datetime.now()))
    try:   
        iris.save(source_for_saving, tempfilename, **kwargs)
        print('jules.save: finished iris.save at '+str(datetime.datetime.now()))
        _run_command('mv %s %s', tempfilename, target)
        print('jules.save: finished first mv at '+str(datetime.datetime.now()))
        
        if np.isfinite(missingdata): #stops it being used for NaNs
            _run_command('cdo setmissval,%s %s %s', str(missingdata), target, tempfilename)
            _run_command('mv %s %s', tempfilename, target)

# may need this back in:            
#        time_str_arr=[]
#        for cube in source_for_saving:
#            element_coord_names = [ coord.name() for coord in cube.coords(dim_coords=True) ]  
#            time_str_arr.append(_get_time_str(element_coord_names))
#        time_dim_always_present = all([a != None for a in time_str_arr])
#        
#        if not time_dim_always_present:
#            _run_command('nccopy -u %s %s', target, tempfilename) 
#            _run_command('mv %s %s', tempfilename, target)

        
    finally:
        if os.path.isfile(tempfilename):
            os.remove(tempfilename)
    print('jules.save: about to leave function at '+str(datetime.datetime.now()))
    
     
def grid_to_points(cube, lsmask, latlon_str=None, lsmask_missingdata_str='nan', data2D = False, user_def_start_corner=None):  
    """
    Takes a cube containing data on a regular grid and returns the data as a cube of land points only.

    Args:
    
     * cube:
         Cube with data on a regular lat lon grid.
     * lsmask:
         The land-sea mask as a cube
     
    Kwargs:
    
     * missingdata:
         What the missing data in the output array will be set to.
     * lsmask_missingdata_str:
         Specifies how sea points are labelled in lsmask. Can be either 'nan' or 'zeros'. For an array
               of Trues and Falses, where False labels sea points, pick 'zeros'.
     * latlon_str:
         A tuple (lat_str,lon_str) containing the variable names for longitude and latitude in lsmask.
     * data2D:
         If data2D=True, there will be two spatial dimensions with sizes (1,npoints), rather than just one with size (npoints).  
         This is a similar format to the JULES output files.
     * user_def_start_corner:
         When converting to a list of points, a particular order of points is chosen as default. This keyword allows 
         this choice to be overridden. Allowed values are 'topleft' (starting in the top left and going downwards first) and 
         'bottomleft' (starting in the bottom left and going towards the right first).
        
    """
    
    # Need a test to make sure that we're starting with a gridded cube
        
    # tidy mask up    
        
    possible_lsmask_missingdata_str = ['nan', 'zeros']
    possible_start_corner = ['topleft', 'bottomleft']
    
    if not (lsmask_missingdata_str in possible_lsmask_missingdata_str):
        raise UserWarning("allowed values for lsmask_missingdata_str are " + 
                  str(possible_lsmask_missingdata_str)+
                  '. The lsmask_missingdata_str given by the user is not one of these allowed values.')
        
    if lsmask_missingdata_str == 'nan':  
        mask = np.isfinite(lsmask.data)   
       
    if lsmask_missingdata_str == 'zeros': 
        mask = np.ma.make_mask(lsmask.data)   
    
    if len(mask.shape) != 2:
        print(mask.shape)
        raise UserWarning('mask should be 2-dimensional (lat and lon)')        
    
    # get the lat and lon strings for the cube and mask
    
    mask_coord_names = [coord.name() for coord in lsmask.coords(dim_coords=True)]
    (_mask_lat_str, _mask_lon_str) = _get_latlon_str(mask_coord_names, latlon_str=latlon_str)
    
    cube_coord_names = [coord.name() for coord in cube.coords(dim_coords=True)]
    (_cube_lat_str, _cube_lon_str) = _get_latlon_str(cube_coord_names, latlon_str=latlon_str)
    
    if len(cube.coord(_cube_lat_str).points) != len(lsmask.coord(_mask_lat_str).points) or \
     not np.allclose(cube.coord(_cube_lat_str).points, lsmask.coord(_mask_lat_str).points, atol=1e-05):
        raise UserWarning('Latitude DimCoord of cube does not match mask')
        
    if len(cube.coord(_cube_lon_str).points) != len(lsmask.coord(_mask_lon_str).points) or \
     not np.allclose(cube.coord(_cube_lon_str).points, lsmask.coord(_mask_lon_str).points, atol=1e-05):
        raise UserWarning('Longitude DimCoord of cube does not match mask')
    
    if cube.coord_dims(cube.coord(_cube_lat_str)) != (len(cube.data.shape)-2, ):
        raise UserWarning('Latitude DimCoord of cube should be last-but-one element of data array')
        
    if cube.coord_dims(cube.coord(_cube_lon_str)) != (len(cube.data.shape)-1, ):
        raise UserWarning('Longitude DimCoord of cube should be last element of data array')
    
    if lsmask.coord_dims(lsmask.coord(_mask_lat_str)) != (len(lsmask.data.shape)-2, ):
        raise UserWarning('Latitude DimCoord of mask should be last-but-one element of data array')
        
    if lsmask.coord_dims(lsmask.coord(_mask_lon_str)) != (len(lsmask.data.shape)-1, ):
        raise UserWarning('Longitude DimCoord of mask should be last element of data array')
    
    if user_def_start_corner is None:
        start_corner = 'topleft'
    else:
        if user_def_start_corner in possible_start_corner:
            start_corner = user_def_start_corner
        else:
            raise UserWarning('value given for user_def_start_corner is one of possible_start_corner')
        
    
    # the lat and lon dimension coordinates will not be copied over.  
    not_to_be_copied_over = (len(cube.data.shape)-2, len(cube.data.shape)-1)
        
    # now to copy over any dim coords we do want
    # n.b. do not want to copy over the lat,lon dim coords
    dim_coords_and_dims = []         
    for coord in cube.coords(dim_coords=True):
        dim_tuple = cube.coord_dims(coord)
        if set(dim_tuple).intersection(set(not_to_be_copied_over)) == set([]):
            dim_coords_and_dims.append((coord, dim_tuple))    
    
    # now to copy over any aux coords we want
    # n.b. do not want to copy over any aux coords which depend on one we're getting rid of
    aux_coords_and_dims = []         
    for coord in cube.aux_coords:
        dim_tuple = cube.coord_dims(coord)
        if set(dim_tuple).intersection(set(not_to_be_copied_over)) == set([]):
            aux_coords_and_dims.append((coord, dim_tuple)) 
    
    # set some stuff that will be useful when making the new lat,lon aux coords
     
    nlat = len(lsmask.coord(_mask_lat_str).points)
    nlon = len(lsmask.coord(_mask_lon_str).points)
    
    npoints = np.sum(mask)
    
    if data2D:
        sizeofdata = [1, npoints]
        ycoord = iris.coords.DimCoord(np.arange(1), long_name='y')
        xcoord = iris.coords.DimCoord(np.arange(npoints), long_name='x')
    else:    
        sizeofdata = [npoints]
        pointscoord = iris.coords.DimCoord(np.arange(npoints), long_name='land')
        
    latdata = np.empty(sizeofdata, dtype=float)
    londata = np.empty(sizeofdata, dtype=float)
    
    # initialise the new data array for the pseudo cube: points data
    
    coord_len = cube.data.shape[:-2]
    coord_len += tuple(sizeofdata)

    # Not the best way to deal with a masked input cube, which has the possibility
    # that some land points are masked, but at least it's safe for now.
    # Means that pcube.data is always an ndarray. 
    if np.ma.is_masked(cube.data):
        data_for_copying = cube.data.filled(np.nan)
    else:
        data_for_copying = cube.data
    pointsdata = np.empty(coord_len, dtype=float)         
    
    ipoint = 0
    if start_corner == 'topleft':
    
        # want to fill the points array from the top left corner e.g. for a 3x3 grid the order is:  
        #    lat
        #     ^     1 4 7     
        #     |     2 5 8   
        #     |     3 6 9      
        #     ------> lon
        
        lat_generator = range(nlat-1, 0-1, -1)
        lon_generator = range(nlon)
        
        if data2D:
            for ilon in lon_generator:
                for ilat in lat_generator: 
                    if mask[ilat, ilon]: 
                        pointsdata[..., 0, ipoint] = data_for_copying[..., ilat, ilon] 
                        latdata[..., ipoint]      = cube.coord(_cube_lat_str).points[ilat]
                        londata[..., ipoint]      = cube.coord(_cube_lon_str).points[ilon]                  
                        ipoint += 1    
        else:
            for ilon in lon_generator:
                for ilat in lat_generator: 
                    if mask[ilat, ilon]: 
                        pointsdata[..., ipoint] = data_for_copying[..., ilat, ilon] 
                        latdata[..., ipoint]    = cube.coord(_cube_lat_str).points[ilat]
                        londata[..., ipoint]    = cube.coord(_cube_lon_str).points[ilon]                  
                        ipoint += 1
        
    elif start_corner == 'bottomleft':
    
        # want to fill the points array from the bottom left corner e.g. for a 3x3 grid the order is:  
        #    lat
        #     ^     7 8 9     
        #     |     4 5 6   
        #     |     1 2 3      
        #     ------> lon
    
        lat_generator = range(nlat)
        lon_generator = range(nlon)
        
        if data2D:
            for ilat in lat_generator: 
                for ilon in lon_generator:
                    if mask[ilat, ilon]: 
                        pointsdata[..., 0, ipoint] = data_for_copying[..., ilat, ilon] 
                        latdata[..., ipoint]      = cube.coord(_cube_lat_str).points[ilat]
                        londata[..., ipoint]      = cube.coord(_cube_lon_str).points[ilon]                  
                        ipoint += 1    
        else:
            for ilat in lat_generator: 
                for ilon in lon_generator:
                    if mask[ilat, ilon]: 
                        pointsdata[..., ipoint] = data_for_copying[..., ilat, ilon] 
                        latdata[..., ipoint]    = cube.coord(_cube_lat_str).points[ilat]
                        londata[..., ipoint]    = cube.coord(_cube_lon_str).points[ilon]                  
                        ipoint += 1       
                        
    else:
        raise UserWarning('start_corner not recognised')
                            
    # create the pseudo cube
    
    pcube = iris.cube.Cube(pointsdata, 
        dim_coords_and_dims = dim_coords_and_dims,
        aux_coords_and_dims = aux_coords_and_dims,
        standard_name       = cube.standard_name, 
        long_name           = cube.long_name, 
        var_name            = cube.var_name, 
        units               = cube.units, 
        cell_methods        = cube.cell_methods, 
        attributes          = cube.attributes
        )
        #aux_factories       = cube.aux_factories,  # need to deal with this in a better way
    
    
    if data2D:
        pcube.add_dim_coord(ycoord, data_dim = (len(pcube.data.shape)-2,) )
        pcube.add_dim_coord(xcoord, data_dim = (len(pcube.data.shape)-1,) )
        pass
    else:    
        pcube.add_dim_coord(pointscoord, data_dim = (len(pcube.data.shape)-1, ) )
    
    # now add the lat and lon aux coords
    
    lat_auxcoord = iris.coords.AuxCoord(latdata, standard_name='latitude',  long_name='latitude',  units='degrees')
    lon_auxcoord = iris.coords.AuxCoord(londata, standard_name='longitude', long_name='longitude', units='degrees')
     
    if data2D:
        pcube.add_aux_coord(lat_auxcoord, data_dims = (len(pcube.data.shape)-2, len(pcube.data.shape)-1))
        pcube.add_aux_coord(lon_auxcoord, data_dims = (len(pcube.data.shape)-2, len(pcube.data.shape)-1))
    else:    
        pcube.add_aux_coord(lat_auxcoord, data_dims = (len(pcube.data.shape)-1, ) )
        pcube.add_aux_coord(lon_auxcoord, data_dims = (len(pcube.data.shape)-1, ) )
    
    return pcube
    
def _nice_lower(str_in):
    '''
    Use when you don't want to raise an Exception when str_in is None
    '''

    if str_in is None:
        lower_str = 'unknown'
    else:
        lower_str = str_in.lower()
    
    return lower_str

def _get_time_str(strarr):
    """Finds the name for time (if present) by checking some common names (possible_time_str)
    against the strings in the array strarr.
    """
    possible_time_str = ['time', 't']     
    
    try: 
        time_str = (next(x for x in strarr if (_nice_lower(x) in possible_time_str) ) )
    except StopIteration:
        time_str = None
        
    return time_str

def _get_latlon_str(strarr, latlon_str=None):
    """Finds the name for latiitude and longitude by checking latlon_str and some common names (possible_lon_str
    and possible_lat_str) against the strings in the array strarr.
    """
    possible_lat_str = ['latitude', 'lat', 'lats']   
    possible_lon_str = ['longitude', 'lon', 'lons']
    
    if latlon_str is None or len(latlon_str) != 2: # len(None) raises an exception but this case is checked for beforehand
        user_lat_str = None
        user_lon_str = None
    else:    
        user_lat_str = latlon_str[0]
        user_lon_str = latlon_str[1]
    
    if isinstance(user_lat_str, str) and user_lat_str in strarr:
        found_lat_str = user_lat_str
    else:        
        try: 
            found_lat_str = (next(x for x in strarr if (_nice_lower(x) in possible_lat_str) ) )
        except StopIteration:
            raise _LonlatException("_get_latlon_str: none of the possible variable names for latitude match")
            
    if isinstance(user_lon_str, str) and user_lon_str in strarr:
        found_lon_str = user_lon_str
    else:      
        try: 
            found_lon_str = (next(x for x in strarr if (_nice_lower(x) in possible_lon_str) ) )

        except StopIteration:
            raise _LonlatException("_get_latlon_str: none of the possible variable names for longitude match")
            
    return (found_lat_str, found_lon_str)

class _LonlatException(UserWarning):
    pass
    
class _MaybeAlreadyGriddedException(UserWarning):
    def __init__(self, value=None):
        self.value = value
    def __str__(self):
        msg_string = 'This data may already be gridded. (Are you calling this function ' \
                     'via jules.load, jules.load_cube, jules.load_dump or ' \
                     'jules.load_frac_cube? If so, try setting conv_to_grid=False )' 
        if self.value is not None:
            msg_string = repr(self.value) + ' ' + msg_string
        return msg_string

def _generate_temp_file_name():
    ''' generates a name that can be used for temporary files
    '''
    
    fd, filename = tempfile.mkstemp(suffix='.nc')
    os.close(fd)
        
    return filename


def _run_command(cmd, *args):
    '''wrap the subprocess.call bit, so that program is easier to read'''

    cmdline = cmd % args
    
    retcode = subprocess.call(cmdline, shell=True)
    
    if retcode != 0:
        raise UserWarning(cmdline+'\n did not run successfully, return code = '+str(retcode))
 
def _flatten_pseudogrid(cube, latlon_str=None):
    """Takes a 'gridded' JULES pseudocube and returns a vectorised version
    We'll take the approach that we use elsewhere and build the bits for a new cube
    and then put it together.
    """
    
    cube_coord_names              = [coord.name() for coord in cube.coords()]
    (cube_lat_str, cube_lon_str)   = _get_latlon_str(cube_coord_names, latlon_str=latlon_str)
   
    # the lat and lon dimension coordinates will not be copied over.  
    not_to_be_copied_over = (len(cube.data.shape)-2, len(cube.data.shape)-1)
        
    # now to copy over any dim coords we do want
    # n.b. do not want to copy over the lat,lon dim coords
    dim_coords_and_dims = []         
    for coord in cube.coords(dim_coords=True):
        dim_tuple = cube.coord_dims(coord)
        if set(dim_tuple).intersection(set(not_to_be_copied_over)) == set([]):
            dim_coords_and_dims.append((coord, dim_tuple))    
    
    # now to copy over any aux coords we want
    # n.b. do not want to copy over any aux coords which depend on one we're getting rid of
    aux_coords_and_dims = []         
    for coord in cube.aux_coords:
        dim_tuple = cube.coord_dims(coord)
        if set(dim_tuple).intersection(set(not_to_be_copied_over)) == set([]):
            aux_coords_and_dims.append((coord, dim_tuple)) 
    
    # set some stuff that will be useful when making the new lat,lon aux coords
    nlat       = cube.coord(cube_lat_str).points.shape[0]
    nlon       = cube.coord(cube_lat_str).points.shape[1]
    sizeofdata = [nlat * nlon]
    latdata    = np.empty(sizeofdata, dtype=float)
    londata    = np.empty(sizeofdata, dtype=float)
    
    # initialise the new data array for the pseudo cube: points data
    coord_len  = cube.data.shape[:-2]
    coord_len += tuple(sizeofdata)
    pointsdata = np.empty(coord_len, dtype=float) 
    
    ipoint = 0
    for ilat in range(nlat):
        for ilon in range(nlon): 
            pointsdata[..., ipoint] = cube.data[..., ilat, ilon] 
            latdata[..., ipoint]    = cube.coord(cube_lat_str).points[ilat, ilon]
            londata[..., ipoint]    = cube.coord(cube_lon_str).points[ilat, ilon]                  
            ipoint                += 1    
    
    #create the pcube
    pcube = iris.cube.Cube(pointsdata, 
    dim_coords_and_dims = dim_coords_and_dims,
    aux_coords_and_dims = aux_coords_and_dims,
    standard_name       = cube.standard_name, 
    long_name           = cube.long_name, 
    var_name            = cube.var_name, 
    units               = cube.units, 
    cell_methods        = cube.cell_methods, 
    aux_factories       = cube.aux_factories, 
    attributes          = cube.attributes
    )

    # now add the lat and lon aux coords
    lat_auxcoord = iris.coords.AuxCoord(latdata, standard_name='latitude',  long_name='latitude',  units='degrees')
    lon_auxcoord = iris.coords.AuxCoord(londata, standard_name='longitude', long_name='longitude', units='degrees')
    
    pcube.add_aux_coord(lat_auxcoord, data_dims = (len(pcube.data.shape)-1))
    pcube.add_aux_coord(lon_auxcoord, data_dims = (len(pcube.data.shape)-1))
    
    return pcube
    
def pick_1point(cubelist_full, desired_lonlat, latlon_str=None, tol=3.0):
    ''' 
    Takes a cubelist as argument, returns a cubelist, with just one lon,lat point 
    (specified by desired_lonlat).
    
    Args:
    
    * cubelist_full:
        cubelist of data to pick a point from
    * desired_lonlat:
        A tuple (latitude, longitude) specifying which point is to be picked
    
    Kwargs:
    
    * latlon_str:
        A tuple (lat_str,lon_str) containing the variable names for longitude and latitude in lsmask.
    * tol:
        The cube method nearest_neighbour_index is used to extract the point. 
        If the resulting latitude and longitude is further from the desired_lonlat than the value tol
        then a UserWarning is raised.
         
    '''
    (desired_lon, desired_lat) = desired_lonlat

    if not isinstance(cubelist_full, iris.cube.CubeList):
        raise UserWarning('cubelist_full needs to be an iris.cube.CubeList object')
    
    cubelist = iris.cube.CubeList()
    
    for cube in cubelist_full:
        all_coord_names = [coord.name() for coord in cube.coords(dim_coords=True)]
        
        (_lat_str, _lon_str) = _get_latlon_str(all_coord_names, latlon_str=latlon_str)
        
        if cube.coord_dims(cube.coords(_lat_str)[0]) != (len(cube.data.shape)-2, ):
            raise UserWarning('Latitude DimCoord of cube should be last-but-one element of data array')
            
        if cube.coord_dims(cube.coords(_lon_str)[0]) != (len(cube.data.shape)-1, ):
            raise UserWarning('Longitude DimCoord of cube should be last element of data array')
            
        ilon = cube.coord(_lon_str).nearest_neighbour_index(desired_lon) 
        ilat = cube.coord(_lat_str).nearest_neighbour_index(desired_lat) 
        # should there be warnings here if ilon,ilat are zero or max?
            
        print('selected point has lon =')
        selected_lon = cube.coord(_lon_str).points[ilon]
        print(selected_lon)
        print('and lat =') 
        selected_lat = cube.coord(_lat_str).points[ilat]
        print(selected_lat)
        if abs(selected_lon - desired_lon) > tol or abs(selected_lat - desired_lat) > tol:
            print(selected_lon, desired_lon)
            print(selected_lat, desired_lat)
            raise UserWarning('Selected point is more than ' + str(tol) + ' away from desired lon, lat')
            
        cubelist.append(cube[..., ilat, ilon])

    return cubelist
    
def write_1point_file(cubelist_full, outfile, desired_vars, desired_lonlat, latlon_str=None):
    ''' 
    Writes a text file suitable for input into a JULES single point run
    
    Args:
    
    * cubelist_full:
        cubelist of data to pick a point from
    * outfile:
        filename to save the result to.
    * desired_vars:
        list of names of variables to extract and save.
    * desired_lonlat:
        A tuple (latitude, longitude) specifying which point is to be picked
    
    Kwargs:
    
    * latlon_str:
        A tuple (lat_str,lon_str) containing the variable names for longitude and latitude in lsmask.
         
    '''

    if not isinstance(cubelist_full, iris.cube.CubeList):
        raise UserWarning('cubelist_full needs to be an iris.cube.CubeList object')
        
    cubelist = pick_1point(cubelist_full, desired_lonlat, latlon_str=latlon_str)
 
    all_var_names = [ cube.var_name for cube in cubelist]
    
    with open(outfile, "w") as outp:
    
        store_number_of_lines = None
        store_indices = []
        store_time_str = None       
    
        outp.write('#')
        
        # check that each desired variable is in cubelist, check that if there's a time dimension,
        # it's consistently named and has same length in all cubes in cubelist
        # then print out first line of file, which gives variable names an number of columns for each variable
        
        for var_name in desired_vars:
            try:
                ind = all_var_names.index(var_name)
            except ValueError:
                print('var_name=', var_name)
                print('desired_vars=', desired_vars)
                raise UserWarning("this var_name is not in the cubelist")
                
            store_indices.append(ind)
            
            all_coord_names = [ coord.name() for coord in cubelist[ind].coords(dim_coords=True) ]
            
            _time_str = _get_time_str(all_coord_names)
            
            if _time_str is None:
                temp_number_of_lines = 1
            else:
                temp_number_of_lines = len(cubelist[ind].coord(_time_str).points)
                          
            if store_number_of_lines is None:
                store_number_of_lines = temp_number_of_lines
                store_time_str = _time_str
            else:
                if store_number_of_lines != temp_number_of_lines:
                    raise UserWarning("variables have time dimensions with different length")
                if store_time_str != _time_str:
                    raise UserWarning("the time coord for each cube in the cubelist need to be have the same name")
                
            if store_time_str is None:
                outp.write(str(var_name)+'('+str(cubelist[ind].data.size)+") ")
            else:    
                indx = [Ellipsis]*cubelist[ind].data.ndim # create a list of :
                
                # which dimension of cubelist[ind].data corresponds to time
                dim_time = cubelist[ind].coord_dims(cubelist[ind].coords(_time_str)[0])[0]
                
                indx[dim_time] = 0 # swap one of the : with 0 so that just the first time is selected
                outp.write(str(var_name)+'('+str(cubelist[ind].data[indx].size)+") ")

        outp.write("\n")

        # now print out data. Each line is a timestep

        for line in range(store_number_of_lines):
            for ind in store_indices:
            
                if store_time_str is None:
                    if cubelist[ind].data.shape == ():
                        outp.write(str(cubelist[ind].data)+" ")
                    else:    
                        outp.write(" ".join(str(x) for x in cubelist[ind].data)+" ")
                else:
                    indx = [Ellipsis]*cubelist[ind].data.ndim # create a list of :
                    
                    # which dimension of cubelist[ind].data corresponds to time
                    dim_time = cubelist[ind].coord_dims(cubelist[ind].coords(_time_str)[0])[0]
                    
                    indx[dim_time] = line # swap one of the : with 0 so that just the first time is selected
                    
                    outp.write(" ".join(str(x) for x in cubelist[ind].data[indx])+" ")
                           
            outp.write("\n")
           

def guess_some_dim_coords(cube, **kwargs):

    raise UserWarning("Deprecated function")
    
            
def add_missing_dims(collapsed_cube, cube):
    ''' Deprecated. Please use iris.util.as_compatible_shape instead.'''

    print("Deprecated function. Please use iris.util.as_compatible_shape instead.")
    
    new_cube = iris.util.Resolve(collapsed_cube, cube)(cube.core_data())

    return new_cube


def pad_latitude(cube_or_cubelist_in, sample_cube_on_input_grid, latlon_str=None):
    '''
    Pads cubes with nans for latitude values that exist in sample_cube_on_input_grid
    but not cube_or_cubelist_in.
    All the cubes in cube_or_cubelist_in and sample_cube_on_input_grid must have lat as the second-to-last dimension.
    Particularly useful if a file has been converted to a points-only and then back to a 2D grid, and you 
    want to make sure it has the same latitude dimension that you started with.
    
    Args:
    
    * cube_or_cubelist_in:
        cube or cubelist to be padded
    * sample_cube_on_input_grid:
        cube with the desired latitude values
    
    Kwargs:
    
    * latlon_str:
        A tuple of strings (lat_str,lon_str) containing additional names for longitude and latitude 
        (if latlon_str=None, the standard names e.g 'longitude', 'lon' are searched for).
    
    '''
      
    if latlon_str == None:
        lat_str='latitude'
        lon_str='longitude'
    else:    
        (lat_str, lon_str) = latlon_str

    if isinstance(cube_or_cubelist_in, iris.cube.Cube):
        cubelist_in = iris.cube.CubeList([cube_or_cubelist_in.copy()])
    elif isinstance(cube_or_cubelist_in, iris.cube.CubeList):
        cubelist_in = iris.cube.CubeList([])
        for cube in cube_or_cubelist_in:
            cubelist_in.append(cube.copy())
    else:
        raise UserWarning('first argument is not an iris Cube or CubeList')       
      
    cubelist_on_input_grid = iris.cube.CubeList()
            
    while cubelist_in:
        small_cube = cubelist_in.pop()
        
        dim_lat = small_cube.ndim-2
        dim_lon = small_cube.ndim-1
        
        if small_cube.coord_dims(lat_str) != (dim_lat,):
            raise UserWarning('lat is not in the expected place (second from last dimension of cube.data)')
        
        new_shape = list(small_cube.shape)
        new_shape[-2] = len(sample_cube_on_input_grid.coord(lat_str).points)
        new_data = np.empty(new_shape)
        
        new_data.fill(np.nan)
        
        list_coords_and_dims = []
        
        for coord in small_cube.coords():
            dims = small_cube.coord_dims(coord)

            if dim_lat in dims: # copy across coords from sample_cube_on_input_grid if they 
                                # depend just on lat or just on lat,lon (if lon is last dim in both data arrays)
                if coord.name() in [c.name() for c in sample_cube_on_input_grid.coords()]:
                                    
                    if dims == (dim_lat,):
                        list_coords_and_dims.append((sample_cube_on_input_grid.coord(coord.name()).copy(), dims))
                        
                    elif dims == (dim_lat, dim_lon):
                        if small_cube.coord_dims(lon_str) == (dim_lon,):
                            if sample_cube_on_input_grid.coord_dims(lon_str) == (dim_lon,): 
                                if np.allclose(small_cube.coord(lon_str).points, sample_cube_on_input_grid.coord(lon_str.points)):
                                    list_coords_and_dims.append((sample_cube_on_input_grid.coord(coord.name()).copy(), dims))
                
                # ignore coords that depend on lat and aren't in sample_cube_on_input_grid, or depend
                # on coords other than lat, lon
                
            else: # append coords that don't depend on lat       
                list_coords_and_dims.append((coord, dims))
            
        new_cube = iris.cube.Cube(new_data)
        new_cube.metadata = copy.deepcopy(small_cube.metadata)
        
        # add the coords to the new cube
        
        for (coord, dims) in list_coords_and_dims:
        
            if small_cube.coords(coord.name(), dim_coords=True):
                new_cube.add_dim_coord( coord, dims )                    
            else:
                new_cube.add_aux_coord( coord, dims )        
        
        # copy the data across to the new cube
        
        ymin = int(np.where(np.isclose(new_cube.coord(lat_str).points, small_cube.coord(lat_str).points[0]))[0])
        ymax = int(np.where(np.isclose(new_cube.coord(lat_str).points, small_cube.coord(lat_str).points[-1]))[0])
        
        if not np.allclose(new_cube.coord(lat_str).points[ymin : ymax+1], small_cube.coord(lat_str).points):
            raise UserWarning('latitudes do not match')
        
        new_cube.data[..., ymin : ymax+1, :] = small_cube.data
          
        cubelist_on_input_grid.append(new_cube)
        
    if isinstance(cube_or_cubelist_in, iris.cube.Cube):
        if len( cubelist_on_input_grid ) == 1:
            result = cubelist_on_input_grid[0]
        else:
            raise UserWarning('CubeList does not have length 1')
    else:
        result = cubelist_on_input_grid
            
    return  result
