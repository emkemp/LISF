#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
#------------------------------------------------------------------------------
#
# SCRIPT: make_noahmp401_geotiff.py
#
# PURPOSE:  Extracts postprocessed LIS-NoahMP 4.0.1 2d ensemble gridspace data
# and converts to GeoTIFF.
#
# REQUIREMENTS as of 10 Aug 2022:
# * Python 3.9 or higher
# * UNIDATA NetCDF4 Python library (for reading netCDF4 files)
# * GDAL Python library (bundled in osgeo package).
#
# REVISION HISTORY:
# 10 Aug 2022: Eric Kemp/SSAI, first version, liberally borrowing from older
#   make_sm_geotiff.py script.
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import sys

# Third party modules
# NOTE: pylint cannot see the Dataset class in netCDF4 since the latter is
# not written in Python.  We therefore disable a check for this line to
# avoid a known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
from osgeo import gdal, osr

_MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
           "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

_SOIL_LAYERS = ["0-0.1m", "0.1-0.4m",  "0.4-1.0m",  "1.0-2.0m"]

_VARNAMES_2D = ["Albedo_tavg", "AvgSurfT_inst", "AvgSurfT_tavg",
                "CanopInt_inst", "Elevation_inst", "Evap_tavg",
                "Greenness_inst", "LWdown_f_inst", "LWdown_f_tavg",
                "Landcover_inst", "Landmask_inst", "Psurf_f_inst",
                "Psurf_f_tavg", "Qair_f_inst", "Qair_f_tavg",
                "Qg_tavg", "Qh_tavg", "Qle_tavg", "Qs_acc", "Qsb_acc",
                "RHMin_inst", "Tair_f_min", "SWE_inst",
                "SWdown_f_inst", "SWdown_f_tavg", "SnowDepth_inst",
                "Snowcover_inst", "Soiltype_inst", "Tair_f_inst",
                "Tair_f_max", "Tair_f_tavg", "TotalPrecip_acc",
                "Wind_f_inst", "Wind_f_tavg"]

_VARNAMES_3D = ["RelSMC_inst", "SmLiqFrac_inst", "SoilMoist_inst",
                "SoilMoist_tavg", "SoilTemp_inst", "SoilTemp_tavg"]

def _usage():
    """Print command line usage."""
    txt = f"[INFO] Usage: {sys.argv[0]} ldtfile lvtfile yyyymmddhh"
    print(txt)
    print("[INFO]  where:")
    print("[INFO]    ldtfile: LDT parameter file with full lat/lon data")
    print("[INFO]    lvtfile: LVT 2d ensemble gridspace file")
    print("[INFO]    yyyymmddhh: Valid date and time (UTC)")

def _read_cmd_args():
    """Read command line arguments."""

    # Check of argument count is correct.
    if len(sys.argv) != 4:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if LDT parameter file can be opened.
    ldtfile = sys.argv[1]
    ncid_ldt = nc4_dataset(ldtfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_ldt.close()

    # Check of LVT file can be opened
    lvtfile = sys.argv[2]
    ncid_lvt = nc4_dataset(lvtfile, mode='r', format='NETCDF4_CLASSIC')
    ncid_lvt.close()

    yyyymmddhh = sys.argv[3]

    return ldtfile, lvtfile, yyyymmddhh

def _make_geotransform(lon, lat, nxx, nyy):
    """Set affine transformation from image coordinate space to georeferenced
    space.  See https://gdal.org/tutorials/geotransforms_tut.html"""
    xmin = lon.min()
    xmax = lon.max()
    ymin = lat.min()
    ymax = lat.max()
    xres = (xmax - xmin) / float(nxx)
    yres = (ymax - ymin) / float(nyy)
    # Sujay's original code
    #geotransform = (xmax, xres, 0, ymin, 0, -yres)
    # Eric's code...Based on gdal.org/tutorials/geotransforms_tut.html
    # xmin is x-coordinate of upper-left corner of upper-left pixel
    # ymax is y-coordinate of upper-left corner of upper-left pixel
    # Third variable is row rotation, set to zero for north-up image
    # Fourth variable is column rotation, set to zero
    # Sixth variable is n-s pixel resolution (negative for north-up image)
    geotransform = (xmin, xres, 0, ymax, 0, -1*yres)
    return geotransform

def _create_output_raster(outfile, nxx, nyy, geotransform, var1):
    """Create the output raster file (the GeoTIFF), including map projection"""
    output_raster = gdal.GetDriverByName('GTiff').Create(outfile,
                                                         nxx, nyy, 1,
                                                         gdal.GDT_Float32)

    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326) # Corresponds to WGS 84
    output_raster.GetRasterBand(1).SetNoDataValue(-9999)
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(var1)

    return output_raster

def _make_timestring(yyyymmddhh):
    """Create time metadata"""
    validdt = datetime.datetime(year=int(yyyymmddhh[0:4]),
                                month=int(yyyymmddhh[4:6]),
                                day=int(yyyymmddhh[6:8]),
                                hour=int(yyyymmddhh[8:10]))
    time_string = f"Valid {validdt.hour:02d}Z"
    time_string += f" {validdt.day}"
    time_string += f" {_MONTHS[validdt.month-1]}"
    time_string += f" {validdt.year:04d}"
    return time_string

def _set_calc_type(varname):
    """Return string summarizing how variable was calculated"""
    if varname == "RHMin_inst":
        calctype = "Value at time of minimum temperature"
    elif varname[-5:] == "_tavg":
        calctype = "3-hr time average"
    elif varname[-5:] == "_inst":
        calctype = "Instantaneous value"
    elif varname[-4:] == "_acc":
        calctype = "3-hr accumulation"
    elif varname[-4:] == "_min":
        calctype = "3-hr minimum"
    elif varname[-4:] == "_max":
        calctype = "3-hr maximum"
    else:
        print(f"[ERR] Unrecognized variable {varname}")
        sys.exit(1)
    return calctype

def _set_landuse_key():
    """Return string summarizing land use categories"""
    summary = "IGBP-Modified MODIS Land Use Categories:"
    summary += " 1=Evergreen Needleleaf Forest"
    summary += " 2=Evergreen Broadleaf Forest"
    summary += " 3=Deciduous Needleleaf Forest"
    summary += " 4=Deciduous Broadleaf Forest"
    summary += " 5=Mixed Forests"
    summary += " 6=Closed Shrublands"
    summary += " 7=Open Shrublands"
    summary += " 8=Woody Savannas"
    summary += " 9=Savannas"
    summary += " 10=Grasslands"
    summary += " 11=Permanent Wetlands"
    summary += " 12=Croplands"
    summary += " 13=Urban and Built-Up"
    summary += " 14=Cropland/Natural Vegetation Mosaic"
    summary += " 15=Snow and Ice"
    summary += " 16=Barren or Sparsely Vegetated"
    summary += " 17=Water"
    summary += " 18=Wooded Tundra"
    summary += " 19=Mixed Tundra"
    summary += " 20=Barren Tundra"
    return summary

def _set_landmask_key():
    """Return string summarizing land mask"""
    summary = "Land-water mask:"
    summary += " 0=Water"
    summary += " 1=Land"
    return summary

def _set_soiltype_key():
    """Return string summarizing soil categories"""
    summary = "NCAR STATSGO+FAO blended soil texture categories:"
    summary += " 1=Sand"
    summary += " 2=Loamy Sand"
    summary += " 3=Sandy Loam"
    summary += " 4=Silt Loam"
    summary += " 5=Silt"
    summary += " 6=Loam"
    summary += " 7=Sandy Clay Loam"
    summary += " 8=Silty Clay Loam"
    summary += " 9=Clay Loam"
    summary += " 10=Sandy Clay"
    summary += " 11=Silty Clay"
    summary += " 12=Clay"
    summary += " 13=Organic Material"
    summary += " 14=Water"
    summary += " 15=Bedrock"
    summary += " 16=Other (land-ice)"
    return summary

def _proc_soilvar_3d(ncid, varname, longitudes, latitudes, yyyymmddhh):
    """Process all layers of specified variable."""
    for i in range(0, 4): # Loop across four LSM layers
        nrows, ncols = ncid.variables[varname][i,:,:].shape
        geotransform = _make_geotransform(longitudes, latitudes, ncols, nrows)
        filename = f"{varname}.{_SOIL_LAYERS[i]}.{yyyymmddhh}.tif"
        output_raster = \
            _create_output_raster(filename, ncols, nrows,
                                  geotransform,
                                  ncid.variables[varname][i,:,:][::-1,:])
        metadata = {}
        metadata['valid_date_time'] = _make_timestring(yyyymmddhh)
        metadata['varname'] = ncid.variables[varname].getncattr('long_name')
        metadata['soil_layer'] = _SOIL_LAYERS[i]
        metadata['units'] = ncid.variables[varname].getncattr('units')
        metadata['model'] = "NoahMP 4.0.1"
        metadata['processing_note'] = _set_calc_type(varname)
        output_raster.GetRasterBand(1).SetMetadata(metadata)
        output_raster.FlushCache() # Write to disk
        del output_raster

def _proc_soilvar_2d(ncid, varname, longitudes, latitudes, yyyymmddhh):
    """Process specified 2d variable."""
    nrows, ncols = ncid.variables[varname][:,:].shape
    geotransform = _make_geotransform(longitudes, latitudes, ncols, nrows)
    filename = f"{varname}.{yyyymmddhh}.tif"
    output_raster = _create_output_raster(filename, ncols, nrows,
                                          geotransform,
                                          ncid.variables[varname][:,:][::-1,:])
    metadata = {}
    metadata['valid_date_time'] = _make_timestring(yyyymmddhh)
    metadata['varname'] = ncid.variables[varname].getncattr('long_name')
    metadata['units'] = ncid.variables[varname].getncattr('units')
    metadata['model'] = "NoahMP 4.0.1"
    metadata['processing_note'] = _set_calc_type(varname)
    if varname == "Landcover_inst":
        metadata["key"] = _set_landuse_key()
    elif varname == "Landmask_inst":
        metadata["key"] = _set_landmask_key()
    elif varname == "Soiltype_inst":
        metadata["key"] = _set_soiltype_key()
    output_raster.GetRasterBand(1).SetMetadata(metadata)
    output_raster.FlushCache() # Write to disk
    del output_raster

def _main():
    """Main driver"""

    # Process the command line arguments.
    ldtfile, lvtfile, yyyymmddhh = _read_cmd_args()

    # First, fetch latitudes/longitudes.  This is pulled from the LDT parameter
    # file, since LVT output has data voids over water.
    ncid = nc4_dataset(ldtfile, 'r', format='NETCDF4')
    longitudes = ncid.variables["lon"][:,:]
    latitudes = ncid.variables["lat"][:,:]
    ncid.close()

    # Now process all the variables from the LVT file.
    ncid = nc4_dataset(lvtfile, 'r', format='NETCDF4')
    for varname in _VARNAMES_3D:
        _proc_soilvar_3d(ncid, varname, longitudes, latitudes, yyyymmddhh)
    for varname in _VARNAMES_2D:
        _proc_soilvar_2d(ncid, varname, longitudes, latitudes, yyyymmddhh)
    ncid.close()

if __name__ == "__main__":
    _main()