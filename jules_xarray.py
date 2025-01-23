"""code to read jules land points to iris cube format"""

import xarray as xr
import pandas as pd
import numpy as np
from ncdata.iris_xarray import (
    cubes_from_xarray,
)  # to convert from xarray to cube only required if need conv_to_cube=True
import dask.array as da
from scipy import stats

DIMS_IN_JULES = [
    "scpool",
    "sclayer",
    "soil",
    "pft",
    "type",
    "tile",
    "snow",
    "scalar",
    "ch4layer",
    "bedrock",
    "seed",
]
# these are the dimension names in the jules file
# need to add more dims if they are not defined here
# (and not land, lat, lon, x, y)


# ####################################################
def shift_time_to_middle_of_bounds(ds, method):
    """
    Shifts the time coordinate to the middle of the time bounds.

    Args:
        ds (xarray.Dataset): The dataset containing the time coordinate and time bounds.
        method (str): The method to use for shifting the time coordinate.

    Returns:
        xarray.Dataset: The dataset with the time coordinate shifted to the middle of the time bounds.
    """
    if "time_bounds" in ds.variables:
        if method is not None:
            if "mean" in method:
                time_bounds = ds["time_bounds"]
                midpoint_time = (
                    time_bounds[0, 1].values - time_bounds[0, 0].values
                ) / 2.0
                ds.coords["time"] = ds.coords["time"] - midpoint_time
                ds = ds.drop_vars(["time_bounds"])
            elif "point" in method:
                return ds
            else:
                raise ValueError("cell_methods not recognised")
    return ds


# ####################################################
def check_cell_methods(ds):
    """
    Check the cell_methods attribute of each variable in the dataset and ensure that there is only one unique cell_methods value.

    Parameters:
    ds (xarray.Dataset): The dataset to check.

    Returns:
    set: A set containing the unique cell_methods values found in the dataset.

    Raises:
    ValueError: If there are multiple unique cell_methods values found in the dataset.
    """
    # Initialize a set to store the unique cell_methods
    cell_methods_set = set()

    for var_name, dataarray in ds.data_vars.items():
        # Skip the 'time_bounds' variable
        if var_name == "time_bounds":
            continue

        cell_methods = dataarray.attrs.get("cell_methods")

        # Add the cell_methods to the set
        if cell_methods is not None:
            method = cell_methods.split(":")[1].strip()
            cell_methods_set.add(method)

    # Check the cell_methods are unique
    if len(cell_methods_set) > 1:
        print(f"WARNING: Different cell_methods in one netcdf file:")
        for cell_methods in cell_methods_set:
            print(f"Cell methods available: {cell_methods}")
        if "mean" in cell_methods_set:
            print("WARNING: Loading variables with mean cell_methods only.")
        else:
            raise ValueError(
                "ERROR: different cell methods in file, and no mean present."
            )
    return cell_methods_set


# ####################################################
def ds_to_sparse_grid(ds):
    """
    Convert a JULES dataset to a grid format.

    This function takes a JULES dataset, sorts the grid, drops the latitude and longitude variables,
    creates a multi-index from the latitude and longitude values, assigns the multi-index to the dataset,
     It also reindexes any dimensions found  in DIMS_IN_JULES to start from 1.

    Parameters:
    ds (xarray.Dataset): The JULES dataset to convert.

    Returns:
    ds (xarray.Dataset): The converted dataset.
    """

    if "y" in ds.dims:
        if len(ds.y.values) == 1 and len(ds.x.values) > 1:
            ds = ds.sel({"y": 0})
        ds = ds.squeeze()

    if ("latitude" in ds.data_vars and "longitude" in ds.data_vars) or (
        "latitude" in ds.coords and "longitude" in ds.coords
    ):
        lats = ds.latitude.values
        lons = ds.longitude.values
        if lats.shape != () and lons.shape != ():
            ds = ds.drop_vars(["latitude", "longitude"])
            multi_index = pd.MultiIndex.from_arrays(
                [lats, lons], names=["latitude", "longitude"]
            )
            if "x" in ds.dims:
                mindex_coords = xr.DataArray(multi_index, dims=["x"])
            elif "land" in ds.dims:
                mindex_coords = xr.DataArray(multi_index, dims=["land"])
            ds = ds.assign_coords(new=mindex_coords)
            ds = ds.unstack()
            ds = ds.drop_vars("new")

    # sort out intial dim values - may well need to add more
    for dim_name in DIMS_IN_JULES:
        if dim_name in ds.dims:
            ds[dim_name] = list(range(1, len(ds[dim_name].values) + 1))

    return ds


# ####################################################
def load(
    filenames, cons=None, shift_time=True, l_check_cell_methods=True, conv_to_cube=True
):
    """
    Load data from one or multiple NetCDF files into an xarray dataset.

    Parameters:
    - filenames (str or list): Path(s) to the NetCDF file(s) to load.
    - cons (str): Variable name to extract from the loaded dataset using iris constraint.
    - shift_time (bool): Whether to shift time values to the middle of the bounds.
    - l_check_cell_methods (bool): Whether to check if cell methods are all the same.

    Returns:
    - xarray.Dataset or iris.cube.CubeList: The loaded dataset or extracted cube(s) based on the provided constraints.
    """
    if isinstance(filenames, list):
        if len(filenames) > 1:
            xrdatasets = xr.open_mfdataset(filenames, chunks={}, use_cftime=True)
        else:
            xrdatasets = xr.open_dataset(filenames[0], chunks={})
    else:
        try:
            xrdatasets = xr.open_dataset(filenames, chunks={}, use_cftime=True)
        except:
            raise ValueError(f"something wrong with one of the files: {filenames}")

    # xrdatasets = xrdatasets.fillna(-999)
    for var_name, dataarray in xrdatasets.data_vars.items():
        if np.issubdtype(dataarray.dtype, np.number):
            xrdatasets[var_name] = dataarray.fillna(-999)

    # check cell methods are all the same
    if l_check_cell_methods:
        cell_methods_set = check_cell_methods(xrdatasets)
        if len(cell_methods_set) > 1:
            for var in list(xrdatasets.variables):
                if "cell_methods" in xrdatasets[var].attrs:
                    if "mean" not in xrdatasets[var].attrs["cell_methods"]:
                        # If it is, add the variable to the dictionary
                        del xrdatasets[var]
    else:
        cell_methods_set = None

    if shift_time:
        xrdatasets = shift_time_to_middle_of_bounds(xrdatasets, cell_methods_set)

    xrdatasets_grid = ds_to_sparse_grid(xrdatasets)
    # xrdatasets_grid.to_netcdf("grid_data.nc")
    xrdatasets.close()

    # Assuming xrdatasets_grid is your xarray Dataset
    coords = ["latitude", "longitude"]
    for coord in coords:
        if coord not in xrdatasets_grid.coords:
            if coord == "latitude":
                if "lat" in xrdatasets_grid.coords:
                    xrdatasets_grid = xrdatasets_grid.rename({"lat": "latitude"})
            if coord == "longitude":
                if "lon" in xrdatasets_grid.coords:
                    xrdatasets_grid = xrdatasets_grid.rename({"lon": "longitude"})

    for coord in coords:
        if coord not in xrdatasets_grid.coords:
            raise ValueError(
                f"{coord} not found in the dataset after attempting to rename 'lat'/'lon'."
            )

    for coord in coords:
        values = xrdatasets_grid[coord].values
        if values.shape == ():
            # If there's only one lat or lon value, skip the rest of the loop
            continue

        # Calculate differences between subsequent coordinates
        diff = np.diff(values)
        # Check if all differences are equal
        is_equally_spaced = np.allclose(diff, diff[0])
        if not is_equally_spaced:
            # Define new equally spaced coordinates
            step_size, _ = stats.mode(diff)
            if step_size == 0:
                raise ValueError(f"ERROR: step_size is 0 for {coord}")
            new_values = np.arange(values.min(), values.max() + step_size, step_size)
            # Reindex the data over the new coordinates
            if coord == "latitude":
                xrdatasets_grid = xrdatasets_grid.reindex(
                    latitude=new_values, longitude=xrdatasets_grid["longitude"].values
                )
            else:  # coord == 'longitude'
                xrdatasets_grid = xrdatasets_grid.reindex(
                    latitude=xrdatasets_grid["latitude"].values, longitude=new_values
                )

    if conv_to_cube:  # convert to iris cube
        cubelist = cubes_from_xarray(xrdatasets_grid)
        if cons is not None:
            try:
                cube = cubelist.extract_cube(cons)
            except:
                raise ValueError(
                    f"variable in iris constraint not found from jules_xarray"
                )
            if "latitude" in [coord.name() for coord in cube.coords()]:
                cube.coord("latitude").units = "degrees"
            if "longitude" in [coord.name() for coord in cube.coords()]:
                cube.coord("longitude").units = "degrees"
            cube_lazy = cube.core_data()
            cube_lazy = da.where(da.isnan(cube_lazy), -9999, cube_lazy)
            cube_lazy = da.ma.masked_less(cube_lazy, -9990)
            cube.data = cube_lazy
            return cube
        else:
            return cubelist
    else:  # return xarray dataset
        return xrdatasets_grid
