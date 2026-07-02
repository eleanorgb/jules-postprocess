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
    "ch4subgrid",
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
def check_cell_methods(ds, cell_method_sel="mean"):
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
        if cell_method_sel in cell_methods_set:
            print(f"WARNING: Loading variables with {cell_method_sel} cell_methods only.")
        else:
            raise ValueError(
                f"ERROR: different cell methods in file, and no {cell_method_sel} present."
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
            elif "points" in ds.dims:
                mindex_coords = xr.DataArray(multi_index, dims=["points"])
            elif "y" in ds.dims:
                mindex_coords = xr.DataArray(multi_index, dims=["y"])
            ds = ds.assign_coords(new=mindex_coords)
            ds = ds.unstack()
            ds = ds.drop_vars("new")

    # sort out intial dim values - may well need to add more
    for dim_name in DIMS_IN_JULES:
        if dim_name in ds.dims:
            ds[dim_name] = list(range(1, len(ds[dim_name].values) + 1))

    return ds


# ####################################################
def grid_to_land_points(ds, lsmask=None, land_dim_name="land", fill_value=-999.0):
    """
    Convert a lat/lon grid xarray Dataset back to JULES land-points format.

    This is the inverse of ds_to_sparse_grid.  The output has a single 'land'
    dimension indexed by an integer, with MultiIndex coordinates 'latitude' and
    'longitude' matching the original land-point layout.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset on a regular lat/lon grid.  Latitude and longitude must be
        present as dimensions or coordinates named 'latitude'/'longitude' (or
        'lat'/'lon', which are renamed automatically).
    lsmask : xarray.DataArray or numpy.ndarray, optional
        Boolean land mask on the same lat/lon grid (True = land).  If None,
        any grid point where *all* variables are finite and not equal to
        fill_value is treated as land.
    land_dim_name : str
        Name to give the output land-point dimension.  Default 'land'.
    fill_value : float
        Value used to mark missing/ocean points (default -999).

    Returns
    -------
    xarray.Dataset
        Dataset with a single land-point dimension instead of lat/lon.
    """
    # ---- normalise coordinate names ----------------------------------------
    rename_map = {}
    if "lat" in ds.dims and "latitude" not in ds.dims:
        rename_map["lat"] = "latitude"
    if "lon" in ds.dims and "longitude" not in ds.dims:
        rename_map["lon"] = "longitude"
    if rename_map:
        ds = ds.rename(rename_map)

    if "latitude" not in ds.dims or "longitude" not in ds.dims:
        raise ValueError(
            "Dataset must have 'latitude' and 'longitude' as dimensions."
        )

    # ---- build land mask ----------------------------------------------------
    if lsmask is not None:
        if isinstance(lsmask, np.ndarray):
            mask_da = xr.DataArray(
                lsmask.astype(bool),
                dims=["latitude", "longitude"],
                coords={"latitude": ds["latitude"], "longitude": ds["longitude"]},
            )
        else:
            mask_da = lsmask.astype(bool)
            # align mask to dataset grid just in case
            mask_da = mask_da.reindex_like(
                ds[["latitude", "longitude"]], method="nearest"
            )
    else:
        # derive mask from first numeric variable that spans lat/lon
        mask_da = None
        for var in ds.data_vars:
            da_var = ds[var]
            if "latitude" in da_var.dims and "longitude" in da_var.dims:
                # reduce any extra dims to get a 2-D (lat, lon) mask
                extra_dims = [d for d in da_var.dims
                              if d not in ("latitude", "longitude")]
                data_2d = da_var.isel({d: 0 for d in extra_dims})
                valid = np.isfinite(data_2d.values) & (
                    data_2d.values != fill_value
                )
                mask_da = xr.DataArray(
                    valid,
                    dims=["latitude", "longitude"],
                    coords={
                        "latitude": ds["latitude"],
                        "longitude": ds["longitude"],
                    },
                )
                break
        if mask_da is None:
            raise ValueError(
                "Could not infer a land mask; please supply one via lsmask."
            )

    # ---- stack to land points -----------------------------------------------
    # Stack lat/lon into a single MultiIndex dimension.
    ds_stacked = ds.stack({land_dim_name: ("latitude", "longitude")})

    # Build the mask aligned to the stacked dimension.
    mask_stacked = mask_da.stack({land_dim_name: ("latitude", "longitude")})
    land_idx = np.where(mask_stacked.values)[0]

    if land_idx.size == 0:
        raise ValueError("Land mask contains no True points.")

    ds_land = ds_stacked.isel({land_dim_name: land_idx})

    # Reset to integer index, preserving lat/lon as plain coordinates.
    lats = ds_land["latitude"].values
    lons = ds_land["longitude"].values

    ds_land = ds_land.drop_vars(["latitude", "longitude", land_dim_name])
    ds_land = ds_land.assign_coords(
        {
            land_dim_name: np.arange(len(land_idx)),
            "latitude": (land_dim_name, lats),
            "longitude": (land_dim_name, lons),
        }
    )

    return ds_land


# ####################################################
def load(
    filenames, cons=None, shift_time=True, l_check_cell_methods=True, cell_method_sel = "mean", conv_to_cube=True
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
    def is_valid_file(filename):
        try:
            ds = xr.open_dataset(filename, chunks={}, use_cftime=True)
            for dim in ds.dims:
                if ds.sizes[dim] == 0:
                    return False
            return True
        except:
            return False

    if isinstance(filenames, list):
        filenames = [f for f in filenames if is_valid_file(f)]
        if not filenames:
            raise ValueError("No valid files to open.")
        if len(filenames) > 1:
            try:
                xrdatasets = xr.open_mfdataset(filenames, chunks={}, use_cftime=True)
            except:
                raise ValueError(f"something wrong with one of the netcdf files")
        else:
            try:
                xrdatasets = xr.open_dataset(filenames[0], chunks={})
            except:
                raise ValueError(f"something wrong with: {filenames[0]}")
    else:
        try:
            xrdatasets = xr.open_dataset(filenames, chunks={}, use_cftime=True)
        except:
            raise ValueError(f"something wrong with: {filenames}")

    # xrdatasets = xrdatasets.fillna(-999)
    for var_name, dataarray in xrdatasets.data_vars.items():
        if np.issubdtype(dataarray.dtype, np.number):
            xrdatasets[var_name] = dataarray.fillna(-999)

    # check cell methods are all the same
    if l_check_cell_methods:
        cell_methods_set = check_cell_methods(xrdatasets, cell_method_sel=cell_method_sel)
        if len(cell_methods_set) > 1:
            for var in list(xrdatasets.variables):
                if "cell_methods" in xrdatasets[var].attrs:
                    if cell_method_sel not in xrdatasets[var].attrs["cell_methods"]:
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
        elif len(values) < 2:
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
