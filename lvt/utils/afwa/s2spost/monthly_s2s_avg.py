#!/usr/bin/env python3

"""
#------------------------------------------------------------------------------
#
# SCRIPT: monthly_s2s_avg.py
#
# PURPOSE: Read daily S2S CF-convention netCDF files, calculate monthly
# averages, and write to new CF-convention netCDF file.
#
# REQUIREMENTS as of 16 Sep 2021:
# * Python 3.8 or higher
# * UNIDATA NetCDF4 Python library
# * NumPy Python library.
#
# REVISION HISTORY:
# 16 Sep 2021: Eric Kemp (SSAI), first version.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import sys

# Third-party libraries
import numpy as np
# NOTE: pylint cannot see the Dataset class in netCDF4 since the latter is not
# written in Python.  We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module

# List of variables to average.  This is intended as an internal constant,
# hence the name is prefixed with "_".

_VAR_ACC_LIST = ["Qs_acc", "Qsb_acc", "TotalPrecip_acc"]

_VAR_AVG_LIST = ["Qle_tavg", "Qh_tavg", "Qg_tavg", "Evap_tavg",
                 "AvgSurfT_tavg", "AvgSurfT_inst", "Albedo_tavg",
                 "SWE_inst", "SnowDepth_inst",
                 "SoilMoist_tavg", "SoilMoist_inst",
                 "SoilTemp_tavg", "SoilTemp_inst",
#                 "SmLiqFrac_inst",
                 "CanopInt_inst", "TWS_inst", "GWS_inst", "Snowcover_inst",
                 "Wind_f_tavg", "Wind_f_inst",
                 "Tair_f_tavg", "Tair_f_inst", "Tair_f_min", "Tair_f_max",
                 "Qair_f_tavg", "Qair_f_inst",
                 "Psurf_f_tavg", "Psurf_f_inst",
                 "SWdown_f_tavg", "SWdown_f_inst",
                 "LWdown_f_tavg", "LWdown_f_inst",
                 "Greenness_inst",
#                 "RelSMC_tavg",
                 "RelSMC_inst",
                 "RHMin_inst",
                 "Streamflow_tavg", "FloodedFrac_tavg", "SurfElev_tavg",
                 "SWS_tavg", "RiverStor_tavg", "RiverDepth_tavg",
                 "RiverFlowVelocity_tavg", "FloodStor_tavg",
                 "FloodedArea_tavg"]

_CONST_LIST = ["lat", "lon", "ensemble", "soil_layer",
               "soil_layer_thickness",
               "Landmask_inst", "LANDMASK", "Landcover_inst",
               "Soiltype_inst", "Elevation_inst"]

def _usage():
    """Print command line usage."""
    txt = "[INFO] Usage: %s input_dir output_dir start_yyyymmdd end_yyyymmdd"
    print(txt)
    print("[INFO] where:")
    print("[INFO]  input_dir: directory with daily S2S files in CF convention")
    print("[INFO]  output_dir: directory for output file")
    print("[INFO]  start_yyyymmdd: Starting date/time of daily files")
    print("[INFO]  end_yyyymmdd: Ending date/time of daily files")

def _proc_date(yyyymmdd):
    """Convert YYYYMMDD string to Python date object."""
    if len(yyyymmdd) != 8:
        print("[ERR] Invalid length for YYYYMMDD, must be 8 characters!")
        sys.exit(1)
    year = int(yyyymmdd[0:4])
    month = int(yyyymmdd[4:6])
    day = int(yyyymmdd[6:8])
    try:
        dateobj = datetime.date(year, month, day)
    except ValueError:
        print("[ERR] Invalid YYYYMMDD passed to script!")
        sys.exit(1)
    return dateobj

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 5:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if input directory exists.
    input_dir = sys.argv[1]
    if not os.path.exists(input_dir):
        print("[ERR] Directory %s does not exist!")
        sys.exit(1)

    # Create output directory if it doesn't exist.
    output_dir = sys.argv[2]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get valid starting and ending dates of data.
    start_yyyymmdd = sys.argv[3]
    startdate = _proc_date(start_yyyymmdd)
    end_yyyymmdd = sys.argv[4]
    enddate = _proc_date(end_yyyymmdd)
    if startdate > enddate:
        print("[ERR] Start date is after end date!")
        sys.exit(1)

    return input_dir, output_dir, startdate, enddate

def _create_daily_s2s_filename(input_dir, curdate):
    """Create path to daily S2S netCDF file."""
    name = "%s" %(input_dir)
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += "_GP.LIS-S2S"
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S"
    name += "_DD.%4.4d%2.2d%2.2d" %(curdate.year, curdate.month, curdate.day)
    name += "_DT.0000"
    name += "_DF.NC"
    return name

def _create_monthly_s2s_filename(output_dir, startdate, enddate):
    """Create path to monthly S2S netCDF file."""
    name = "%s" %(output_dir)
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += "_GP.LIS-S2S"
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S"
    name += "_D_.%4.4d%2.2d%2.2d-%4.4d%2.2d%2.2d" \
        %(startdate.year, startdate.month, startdate.day,
          enddate.year, enddate.month, enddate.day)
    name += "_TP.0000-0000"
    name += "_DF.NC"
    return name

def _create_firstguess_monthly_file(infile, outfile):
    """Read daily S2S file, and copy to monthly S2S file with a few
    extra fields.  This allows us to cleanly copy dimensions and all
    attributes.  The numerical values of the arrays in the monthly S2S
    file will be replaced later in the script."""

    if not os.path.exists(infile):
        print("[ERR] %s does not exist!" %(infile))
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    # Create monthly file, copying dimensions and global attributes.
    # NOTE:  We can clean up the global attributes later in the script.
    ncid_out = nc4_dataset(outfile, "w", format='NETCDF4_CLASSIC')
    for dimname in ncid_in.dimensions:
        ncid_out.createDimension(dimname, ncid_in.dimensions[dimname].size)
    for gattrname in ncid_in.__dict__:
        ncid_out.setncattr(gattrname, ncid_in.__dict__[gattrname])

    # Copy the constant and avg fields
    for varname in _CONST_LIST + _VAR_AVG_LIST:
        var_in = ncid_in.variables[varname]
        if "_FillValue" in var_in.__dict__:
            var_out = ncid_out.createVariable(varname, var_in.datatype,
                                              dimensions=var_in.dimensions,
                                              zlib=True,
                                              complevel=1,
                                              shuffle=True,
                                              fill_value=var_in._FillValue)
        elif "missing_value" in var_in.__dict__:
            var_out = ncid_out.createVariable(varname, var_in.datatype,
                                              dimensions=var_in.dimensions,
                                              zlib=True,
                                              complevel=1,
                                              shuffle=True,
                                              fill_value=var_in.missing_value)
        else:
            var_out = ncid_out.createVariable(varname, var_in.datatype,
                                              dimensions=var_in.dimensions,
                                              zlib=True,
                                              complevel=1,
                                              shuffle=True)
        for attrname in var_in.__dict__:
            if attrname == "_FillValue":
                continue
            var_out.setncattr(attrname, var_in.__dict__[attrname])
        if len(var_out.shape) == 4:
            var_out[:,:,:,:] = var_in[:,:,:,:]
        elif len(var_out.shape) == 3:
            var_out[:,:,:] = var_in[:,:,:]
        elif len(var_out.shape) == 2:
            var_out[:,:] = var_in[:,:]
        elif len(var_out.shape) == 1:
            var_out[:] = var_in[:]

    # Copy the var_acc fields, and also create tavg versions.  The values
    # will be overwritten later.
    for orig_varname in _VAR_ACC_LIST:
        var_in = ncid_in.variables[orig_varname]

        varname_tavg = "%s_tavg" %(orig_varname)
        for varname in [orig_varname, varname_tavg]:
            if "_FillValue" in var_in.__dict__:
                var_out = ncid_out.createVariable(varname, var_in.datatype,
                                                  dimensions=var_in.dimensions,
                                                  zlib=True,
                                                  complevel=1,
                                                  shuffle=True,
                                                  fill_value=var_in._FillValue)
            elif "missing_value" in var_in.__dict__:
                var_out = ncid_out.createVariable(varname, var_in.datatype,
                                                  dimensions=var_in.dimensions,
                                                  zlib=True,
                                                  complevel=1,
                                                  shuffle=True,
                                            fill_value=var_in.missing_value)
            else:
                var_out = ncid_out.createVariable(varname, var_in.datatype,
                                                  dimensions=var_in.dimensions,
                                                  zlib=True,
                                                  complevel=1,
                                                  shuffle=True)
            for attrname in var_in.__dict__:
                if attrname == "_FillValue":
                    continue
                var_out.setncattr(attrname, var_in.__dict__[attrname])
            if len(var_out.shape) == 4:
                var_out[:,:,:,:] = var_in[:,:,:,:]
            elif len(var_out.shape) == 3:
                var_out[:,:,:] = var_in[:,:,:]
            elif len(var_out.shape) == 2:
                var_out[:,:] = var_in[:,:]

    ncid_out.close()
    ncid_in.close()

def _read_first_daily_file(infile):
    """Read the first daily S2S file and copy the required variable values to
    appropriate dictionaries."""

    if not os.path.exists(infile):
        print("[ERR] %s does not exist!" %(infile))
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    accs = {}
    avgs = {}
    avgs["counter"] = 1

    # Copy the values of the fields we will average or accumulate
    for varname in _VAR_AVG_LIST:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            avgs[varname] = var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            avgs[varname] = var_in[:,:,:]
        elif len(var_in.shape) == 2:
            avgs[varname] = var_in[:,:]
    for varname in _VAR_ACC_LIST:
        var_in = ncid_in.variables[varname]
        varname_tavg = "%s_tavg" %(varname)
        if len(var_in.shape) == 4:
            avgs[varname_tavg] = var_in[:,:,:,:]
            accs[varname] = var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            avgs[varname_tavg] = var_in[:,:,:]
            accs[varname] = var_in[:,:,:]
        elif len(var_in.shape) == 2:
            avgs[varname_tavg] = var_in[:,:]
            accs[varname] = var_in[:,:]

    ncid_in.close()
    return accs, avgs

def _read_next_daily_file(infile, accs, avgs):
    """Read next daily S2S file and copy the required variable values to
    appropriate dictionaries."""

    if not os.path.exists(infile):
        print("[ERR] %s does not exist!" %(infile))
        sys.exit(1)
    ncid_in = nc4_dataset(infile, 'r', format='NETCDF4_CLASSIC')

    avgs["counter"] += 1

    # Add the values of the fields we will average or accumulate
    for varname in _VAR_AVG_LIST:
        var_in = ncid_in.variables[varname]
        if len(var_in.shape) == 4:
            avgs[varname][:,:,:,:] += var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            avgs[varname][:,:,:] += var_in[:,:,:]
        elif len(var_in.shape) == 2:
            avgs[varname][:,:] += var_in[:,:]
    for varname in _VAR_ACC_LIST:
        var_in = ncid_in.variables[varname]
        varname_tavg = "%s_tavg" %(varname)
        if len(var_in.shape) == 4:
            avgs[varname_tavg][:,:,:,:] += var_in[:,:,:,:]
            accs[varname][:,:,:,:] = var_in[:,:,:,:]
        elif len(var_in.shape) == 3:
            avgs[varname_tavg][:,:,:] = var_in[:,:,:]
            accs[varname][:,:,:] = var_in[:,:,:]
        elif len(var_in.shape) == 2:
            avgs[varname_tavg][:,:] = var_in[:,:]
            accs[varname][:,:] = var_in[:,:]

    ncid_in.close()
    return accs, avgs

def _finalize_avgs(avgs):
    """Finalize averages by dividing by sample size."""
    count = avgs["counter"]
    for varname in avgs:
        if varname == "counter":
            continue
        if len(avgs[varname].shape) == 4:
            avgs[varname][:,:,:,:] /= count
        elif len(avgs[varname].shape) == 3:
            avgs[varname][:,:,:] /= count
        elif len(avgs[varname].shape) == 2:
            avgs[varname][:,:] /= count
    return avgs

def _update_monthly_s2s_values(outfile, accs, avgs):
    """Update the values in the monthly S2S file."""
    ncid = nc4_dataset(outfile, 'a', format='NETCDF4_CLASSIC')
    for dictionary in [accs, avgs]:
        for varname in dictionary:
            if varname == "counter":
                continue
            var = ncid.variables[varname]
            if len(var.shape) == 4:
                var[:,:,:,:] = dictionary[varname][:,:,:,:]
            elif len(var.shape) == 3:
                var[:,:,:] = dictionary[varname][:,:,:]
            elif len(var.shape) == 2:
                var[:,:] = dictionary[varname][:,:]
    ncid.close()

def _driver():
    """Main driver."""

    # Get the directories and dates
    input_dir, output_dir, startdate, enddate = _read_cmd_args()

    # Loop through dates
    curdate = startdate
    delta = datetime.timedelta(days=1)
    while curdate <= enddate:
        print(curdate)
        infile = _create_daily_s2s_filename(input_dir, curdate)
        if curdate == startdate:
            tmp_outfile = "%s/tmp_monthly.nc" %(output_dir)
            _create_firstguess_monthly_file(infile, tmp_outfile)
            accs, avgs = _read_first_daily_file(infile)
        else:
            accs, avgs = _read_next_daily_file(infile, accs, avgs)
        curdate += delta

    # Finalize averages (dividing by number of days).  Then write
    # avgs and accs to file
    avgs = _finalize_avgs(avgs)
    _update_monthly_s2s_values(tmp_outfile, accs, avgs)
    del accs
    del avgs

    # Rename the output file
    outfile = _create_monthly_s2s_filename(output_dir, startdate,
                                           enddate)
    os.rename(tmp_outfile, outfile)

# Invoke driver
if __name__ == "__main__":
    _driver()