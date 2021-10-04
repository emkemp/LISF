#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: make_s2s_metric_geotiff.py
#
# PURPOSE: Generates GeoTIFF files from S2S metrics netCDF file.
#
# REQUIREMENTS as of 4 Oct 2021:
# * Python 3.8 or higher
# * UNIDATA NetCDF4 Python library
# * GDAL Python library (bundled in osgeo package)
#
# REVISION HISTORY:
# * 4 Oct 2021: Eric Kemp (SSAI), first version.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import sys

# Third party modules
# NOTE: pylint cannot see the Dataset class in netCDF4 since that latter is not
# written in Python. We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
from osgeo import gdal, osr

# Private constants
_VARNAMES = ["RootZone_SM_ANOM", "RootZone_SM_SANOM",
             "Streamflow_ANOM", "Streamflow_SANOM",
             "Surface_SM_ANOM", "Surface_SM_SANOM"]

# Private methods
def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s metricfile" %(sys.argv[0])
    print(txt)
    print("[INFO] where:")
    print("[INFO] metricfile: netcdf file with metric data")

def _read_cmd_args():
    """Read command line arguments."""
    # Check if argument count is correct
    if len(sys.argv) != 2:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if metric file can be opened.
    metricfile = sys.argv[1]
    ncid = nc4_dataset(metricfile, mode="r", format="NETCDF4_CLASSIC")
    ncid.close()

    return metricfile

def _parse_metric_filename(metricfile):
    """Parse the name of the metric file, and same components"""
    components = {}
    component_list = metricfile.split("_")
    for component in component_list:
        key = component.split(".")[0]
        value = component.split(".")[1]
        components[key] = value
    return components

def _get_startdate(metricfile, metric_filename_components):
    """Get start date from name of metric file."""
    try:
        yyyymmdd = metric_filename_components["DP"].split("_")[0]
    except KeyError:
        print("[ERR] Cannot resolve data period from file name %s" \
              %(metricfile))
        sys.exit(1)
    startdate = datetime.datetime(year=int(yyyymmdd[0:4]),
                                  month=int(yyyymmdd[4:6]),
                                  day=int(yyyymmdd[6:8]))
    return startdate

def _get_num_months(metricfile):
    """Get number of months in the metrics file."""
    ncid = nc4_dataset(metricfile, 'r', format="NETCDF4_CLASSIC")
    num_months = ncid.dimensions["lead"].size
    return num_months

def _get_latlon(metricfile):
    """Get latitude and longtitude from metricfile."""
    ncid = nc4_dataset(metricfile, 'r', format="NETCDF4_CLASSIC")
    latitudes = ncid.variables["latitude"][:]
    longitudes = ncid.variables["longitude"][:]
    ncid.close()
    return latitudes, longitudes

def _get_variable(varname, metricfile):
    """Retrieve variable and select metadata from metric file."""
    ncid = nc4_dataset(metricfile, 'r', format="NETCDF4_CLASSIC")
    var = ncid.variables[varname][:,:,:,:]
    units = ncid.variables[varname].units
    long_name = ncid.variables[varname].long_name
    ncid.close()
    return var, units, long_name

def _make_geotransform(lons, lats):
    """Set affine transformation from image coordinate space to georeferenced
    space.  See https://gdal.org/tutorials/geotransforms_tut.html"""
    xmin = lons.min()
    xmax = lons.max()
    ymin = lats.min()
    ymax = lats.max()
    nxx = lons.size
    nyy = lats.size
    xres = (xmax - xmin) / float(nxx)
    yres = (ymax - ymin) / float(nyy)
    # Below is based on gdal.org/tutorials/geotransforms_tut.html
    # xmin is x-coordinate of upper-left corner of upper-left pixel
    # ymax in y-coordinate of upper-left corner of upper-left pixel
    # Third variable is row rotation, set to zero for north-up image
    # Fourth variable is column rotation, set to zero
    # Sixth variable is n-s pixel resolution (negative for north-up image)
    geotransform = (xmin, xres, 0, ymax, 0, -1*yres)
    return geotransform

def _create_output_raster(outfile, nxx, nyy, geotransform, var2d):
    """Create the output raster file (the GeoTIFF), including map projection"""
    output_raster = gdal.GetDriverByName('GTiff').Create(outfile,
                                                         nxx, nyy, 1,
                                                         gdal.GDT_Float32)
    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326) # Corresponds to WGS 84
    output_raster.GetRasterBand(1).SetNoDataValue(-9999)
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(var2d)
    return output_raster

def _set_newdate(date):
    """Set a new date one month in the future."""
    if date.month == 12:
        newdate = datetime.datetime(year=(date.year + 1),
                                    month=1,
                                    day=1)
    else:
        newdate = datetime.datetime(year=date.year,
                                    month=(date.month + 1),
                                    day=1)
    return newdate

def _set_monthly_dates(startdate, num_months):
    """Set start and end of each month of data in the metric file."""
    startdates_month = []
    enddates_month = []
    for imonth in range(0, num_months):
        if imonth == 0:
            startdates_month.append(startdate)
            newdate2 = _set_newdate(startdate)
            enddates_month.append(newdate2)
        else:
            newdate1 = _set_newdate(startdates_month[imonth-1])
            newdate2 = _set_newdate(enddates_month[imonth-1])
            startdates_month.append(newdate1)
            enddates_month.append(newdate2)
    return startdates_month, enddates_month

def _make_geotiff_filename(varname, startdate, enddate,
                           metric_filename_components):
    """Make name of new geotiff file."""
    filename = "PS.%s" %(metric_filename_components["PS"])
    filename += "_SC.%s" %(metric_filename_components["SC"])
    filename += "_DI.%s" %(metric_filename_components["DI"])
    filename += "_GP.%s" %(metric_filename_components["GP"])
    filename += "_GR.%s" %(metric_filename_components["GR"])
    filename += "_AR.%s" %(metric_filename_components["AR"])
    filename += "_PA.LIS-S2S-" %(varname.upper())
    filename += "_DP.%4.4d%2.2d%2.2d-%4.4d%2.2d%2.2d" %(startdate.year,
                                                        startdate.month,
                                                        startdate.day,
                                                        enddate.year,
                                                        enddate.month,
                                                        enddate.day)
    filename += "_TP.%s" %(metric_filename_components["TP"])
    filename += "_DF.TIF"
    return filename

def _sub_driver():
    """First part of driver."""
    metricfile = _read_cmd_args()
    metric_filename_components = _parse_metric_filename(metricfile)
    startdate = _get_startdate(metricfile,
                               metric_filename_components)
    lats, lons = _get_latlon(metricfile)
    geotransform = _make_geotransform(lons, lats) # lons, then lats = x, then y
    num_months = _get_num_months(metricfile)
    startdates_month, enddates_month = \
        _set_monthly_dates(startdate, num_months)
    # We bundle these variables into a dictionary so pylint doesn't complain
    # about too many local variables in the driver method.
    bundle = {
        "metricfile" : metricfile,
        "metric_filename_components" : metric_filename_components,
        "lats" : lats,
        "lons" : lons,
        "geotransform" : geotransform,
        "num_months" : num_months,
        "startdates_month" : startdates_month,
        "enddates_month" : enddates_month,
    }
    return bundle

def _driver():
    """Main driver."""
    # First part of driver is in own method, with return variables saved
    # in a dictionary.  This is to appease pylint.
    var_bundle = _sub_driver()
    metadata = {}
    metadata["generating_process"] = \
        var_bundle["metric_filename_components"]["GP"]
    for varname in _VARNAMES:
        var, units, long_name = _get_variable(varname,
                                              var_bundle["metricfile"])
        metadata["varname"] = varname
        metadata["units"] = units
        metadata["long_name"] = long_name
        for imonth in range(0, var_bundle["num_months"]):
            metadata["start_date"] = "%4.4d%2.2d%2.2d" \
                %(var_bundle["startdates_month"][imonth].year,
                  var_bundle["startdates_month"][imonth].month,
                  var_bundle["startdates_month"][imonth].day)
            metadata["end_date"] = "%4.4d%2.2d%2.2d" \
                %(var_bundle["enddates_month"][imonth].year,
                  var_bundle["enddates_month"][imonth].month,
                  var_bundle["enddates_month"][imonth].day)
            var2d = var[0, imonth, :, :] # For now, just use first ens member
            geotiff_filename = \
                _make_geotiff_filename(varname,
                                       var_bundle["startdates_month"][imonth],
                                       var_bundle["enddates_month"][imonth],
                                      var_bundle["metric_filename_components"])
            output_raster = \
                _create_output_raster(geotiff_filename,
                                      len(var_bundle["lons"]),
                                      len(var_bundle["lats"]),
                                      var_bundle["geotransform"], var2d)
            output_raster.GetRasterBand(1).SetMetadata(metadata)
            output_raster.FlushCache() # Write to disk

# Invoke driver
if __name__ == "__main__":
    _driver()
