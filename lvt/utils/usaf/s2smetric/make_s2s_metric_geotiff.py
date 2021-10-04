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

# Private class
class _MetricGeoTiff:
    """Class for building GeoTIFFs from S2S Metric netCDF file."""

    def __init__(self, metricfile):
        """Constructor"""
        self.metricfile = metricfile

        # Parse the name of the metric file, and save elements.  Assumes
        # filename obeys Air Force Weather file naming convention.
        self.metric_filename_elements = {}
        element_list = self.metricfile.split("_")
        for element in element_list:
            key = element.split(".")[0]
            value = element.split(".")[1]
            self.metric_filename_elements[key] = value

        # Find number of months in metrics file, and save latitudes and
        # longitudes.
        ncid = nc4_dataset(self.metricfile, 'r', format="NETCDF4_CLASSIC")
        self.num_months = ncid.dimensions["lead"].size
        self.latitudes = ncid.variables["latitude"][:]
        self.longitudes = ncid.variables["longitude"][:]
        ncid.close()

        # Set start and end of each month of data in the metric file.
        self.startdates_month = []
        self.enddates_month = []
        for imonth in range(0, self.num_months):
            if imonth == 0:
                self.startdates_month.append(self.get_first_startdate())
                newdate2 = _set_newdate(self.get_first_startdate())
                self.enddates_month.append(newdate2)
            else:
                newdate1 = _set_newdate(self.startdates_month[imonth-1])
                newdate2 = _set_newdate(self.enddates_month[imonth-1])
                self.startdates_month.append(newdate1)
                self.enddates_month.append(newdate2)

    def get_first_startdate(self):
        """Get start date from name of metric file."""
        try:
            yyyymmdd = self.metric_filename_elements["DP"].split("_")[0]
        except KeyError:
            print("[ERR] Cannot resolve data period from file name %s" \
                  %(self.metricfile))
            sys.exit(1)
        return datetime.datetime(year=int(yyyymmdd[0:4]),
                                 month=int(yyyymmdd[4:6]),
                                 day=int(yyyymmdd[6:8]))

    def get_geotransform(self):
        """Set affine transformation from image coordinate space to
        georeferenced space. See
        https://gdal.org/tutorials/geotransforms_tut.html"""
        xmin = self.longitudes.min()
        xmax = self.longitudes.max()
        ymin = self.latitudes.min()
        ymax = self.latitudes.max()
        nxx = self.longitudes.size
        nyy = self.latitudes.size
        xres = (xmax - xmin) / float(nxx)
        yres = (ymax - ymin) / float(nyy)
        # Below is based on gdal.org/tutorials/geotransforms_tut.html
        # xmin is x-coordinate of upper-left corner of upper-left pixel
        # ymax in y-coordinate of upper-left corner of upper-left pixel
        # Third variable is row rotation, set to zero for north-up image
        # Fourth variable is column rotation, set to zero
        # Sixth variable is n-s pixel resolution (negative for north-up image)
        return (xmin, xres, 0, ymax, 0, -1*yres)

    def make_geotiff_filename(self, varname, imonth):
        """Make name of new geotiff file."""
        startdate = self.startdates_month[imonth]
        enddate = self.enddates_month[imonth]
        filename = "PS.%s" %(self.metric_filename_elements["PS"])
        filename += "_SC.%s" %(self.metric_filename_elements["SC"])
        filename += "_DI.%s" %(self.metric_filename_elements["DI"])
        filename += "_GP.%s" %(self.metric_filename_elements["GP"])
        filename += "_GR.%s" %(self.metric_filename_elements["GR"])
        filename += "_AR.%s" %(self.metric_filename_elements["AR"])
        filename += "_PA.LIS-S2S-" %(varname.upper())
        filename += "_DP.%4.4d%2.2d%2.2d-%4.4d%2.2d%2.2d" %(startdate.year,
                                                            startdate.month,
                                                            startdate.day,
                                                            enddate.year,
                                                            enddate.month,
                                                            enddate.day)
        filename += "_TP.%s" %(self.metric_filename_elements["TP"])
        filename += "_DF.TIF"
        return filename

    def create_output_raster(self, outfile, var2d):
        """Create the output raster file (the GeoTIFF), including map
        projection"""
        nxx = self.longitudes.size
        nyy = self.latitudes.size
        geotransform = self.get_geotransform()
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

# Private module methods
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

def _get_variable(varname, metricfile):
    """Retrieve variable and select metadata from metric file."""
    ncid = nc4_dataset(metricfile, 'r', format="NETCDF4_CLASSIC")
    var = ncid.variables[varname][:,:,:,:]
    units = ncid.variables[varname].units
    long_name = ncid.variables[varname].long_name
    ncid.close()
    return var, units, long_name

def _driver():
    """Main driver."""
    metricfile = _read_cmd_args()
    mgt = _MetricGeoTiff(metricfile)
    metadata = {}
    metadata["generating_process"] = mgt.metric_filename_elements["GP"]
    for varname in _VARNAMES:
        var, units, long_name = _get_variable(varname, mgt.metricfile)
        metadata["varname"] = varname
        metadata["units"] = units
        metadata["long_name"] = long_name
        for imonth in range(0, mgt.num_months):
            metadata["start_date"] = "%4.4d%2.2d%2.2d" \
                %(mgt.startdates_month[imonth].year,
                  mgt.startdates_month[imonth].month,
                  mgt.startdates_month[imonth].day)
            metadata["end_date"] = "%4.4d%2.2d%2.2d" \
                %(mgt.enddates_month[imonth].year,
                  mgt.enddates_month[imonth].month,
                  mgt.enddates_month[imonth].day)
            var2d = var[0, imonth, :, :]
            geotiff_filename = \
                mgt.make_geotiff_filename(varname, imonth)
            output_raster = \
                mgt.create_output_raster(geotiff_filename, var2d)
            output_raster.GetRasterBand(1).SetMetadata(metadata)
            output_raster.FlushCache() # Write to disk

# Invoke driver
if __name__ == "__main__":
    _driver()
