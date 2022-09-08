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
# SCRIPT: daily_s2spost_nc.py
#
# PURPOSE:  Merges daily netCDF output from LIS-NoahMP and LIS-HYMAP2 into
# single, CF-compliant netCDF4 file for distribution.  Rewritten to use
# NetCDF4 Python library instead of NCO software, to reduce runtime.
#
# REQUIREMENTS as of 08 Sep 2022:
# * Python 3.9 or higher.
# * UNIDATA NetCDF4 Python library
#
# REFERENCES:
# https://cfconventions.org for specifications of NetCDF Climate and Forecast
#   (CF) Metadata Conventions.
# https://unidata.ucar.edu/software/udunits for documentation on UDUNITS2
#   library, which CF is generally consistent with for unit specifications.
#
# REVISION HISTORY:
# 08 Sep 2022: Eric Kemp/SSAI, first version.
#------------------------------------------------------------------------------
"""

# Standard modules
import configparser
import copy
import datetime
import os
import sys
import time

# Third-party libraries
# NOTE: pylint cannot see the Dataset class in netCDF4 since the latter is not
# written in Python.  We therefore disable a check for this line to avoid a
# known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module

_cell_methods = {
    "Albedo_tavg" : "time: mean area: point where land",
    "AvgSurfT_inst" : "area: point where land",
    "AvgSurfT_tavg" : "time: mean area: point where land",
    "BaseflowStor_tavg" : "time: mean area: point where land",
    "CanopInt_inst" : "area: point where land",
    "Elevation_inst" : "area: point where land",
    "Evap_tavg" : "time: mean area: point where land",
    "FloodedArea_tavg" : "time: mean area: point where land",
    "FloodedFrac_tavg" : "time: mean area: point where land",
    "FloodStor_tavg" : "time: mean area: point where land",
    "GWS_inst" : "area: point where land",
    "Greenness_inst" : "area: point where land",
    "LWdown_f_tavg" : "time: mean",
    "Psurf_f_tavg" : "time: mean",
    "Qair_f_tavg" : "time: mean",
    "Qg_tavg" : "time: mean area: point where land",
    "Qh_tavg" : "time: mean area: point where land",
    "Qle_tavg" : "time: mean area: point where land",
    "Qs_acc" : "time: sum area: point where land",
    "Qsb_acc" : "time: sum area: point where land",
    "RelSMC_inst" : "area: point where land",
    "RHMin_inst" : "area: point where land",
    "RiverDepth_tavg" : "time: mean area: point where land",
    "RiverFlowVelocity_tavg" : "time: mean area: point where land",
    "RiverStor_tavg" : "time: mean area: point where land",
    "SmLiqFrac_inst" : "area: point where land",
    "Snowcover_inst" : "area: point where land",
    "SnowDepth_inst" : "area: point where land",
    "SoilMoist_inst" : "area: point where land",
    "SoilMoist_tavg" : "time: mean area: point where land",
    "SoilTemp_inst" : "area: point where land",
    "SoilTemp_tavg" : "time: mean area: point where land",
    "Streamflow_tavg" : "time: mean area: point where land",
    "SurfElev_tavg" : "time: mean area: point where land",
    "SWdown_f_tavg" : "time: mean",
    "SWE_inst" : "area: point where land",
    "SWS_tavg" : "time: mean area: point where land",
    "Tair_f_max" : "time: maximum",
    "Tair_f_min" : "time: minimum",
    "Tair_f_tavg" : "time: mean",
    "TotalPrecip_acc" : "time: sum",
    "TWS_inst" : "area: point where land",
    "Wind_f_tavg" : "time: mean",
}

# Private methods.
def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} configfile noahmp_file hymap2_file"
    txt += " output_dir YYYYMMDDHH model_forcing"
    print(txt)
    print("[INFO] where:")
    print("[INFO] configfile: Path to s2spost config file")
    print("[INFO] noahmp_file: LIS-NoahMP netCDF file (2d ensemble gridspace)")
    txt = "[INFO] hymap2_file: LIS-HYMAP2 netCDF file (2d ensemble gridspace)"
    print(txt)
    print("[INFO] output_dir: Directory to write merged output")
    print("[INFO] YYYYMMDDHH is valid year,month,day,hour of data (in UTC)")
    print("[INFO] model_forcing: ID for atmospheric forcing for LIS")

def _read_cmd_args():
    """Read command line arguments."""

    # Check if argument count is correct.
    if len(sys.argv) != 7:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    # Check if input files exist.
    configfile = sys.argv[1]
    if not os.path.exists(configfile):
        print(f"[ERR] {configfile} does not exist!")
        sys.exit(1)

    noahmp_file = sys.argv[2]
    if not os.path.exists(noahmp_file):
        print(f"[ERR] {noahmp_file} does not exist!")
        sys.exit(1)

    hymap2_file = sys.argv[3]
    if not os.path.exists(hymap2_file):
        print(f"[ERR] {hymap2_file} does not exist!")
        sys.exit(1)

    # Create output directory if it doesn't exist.
    output_dir = sys.argv[4]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get valid date and time of data.
    yyyymmddhh = sys.argv[5]

    if len(yyyymmddhh) != 10:
        print("[ERR] Invalid length for YYYYMMDDHH, must be 10 characters!")
        sys.exit(1)
    year = int(yyyymmddhh[0:4])
    month = int(yyyymmddhh[4:6])
    day = int(yyyymmddhh[6:8])
    hour = int(yyyymmddhh[8:10])

    try:
        curdt = datetime.datetime(year, month, day, hour)
    except ValueError:
        print("[ERR] Invalid YYYYMMDDHH passed to script!")
        sys.exit(1)

    # Get ID of model forcing
    model_forcing = sys.argv[6]

    return configfile, noahmp_file, hymap2_file, output_dir, curdt, \
        model_forcing

def _read_config(configfile):
    """Read from s2spost config file."""
    config = configparser.ConfigParser()
    config.read(configfile)
    return config

def _create_final_filename(output_dir, curdt, model_forcing):
    """Create final filename, following 557 convention."""
    name = f"{output_dir}"
    name += "/PS.557WW"
    name += "_SC.U"
    name += "_DI.C"
    name += f"_GP.LIS-S2S-{model_forcing}"
    name += "_GR.C0P25DEG"
    name += "_AR.AFRICA"
    name += "_PA.LIS-S2S"
    name += f"_DD.{curdt.year:04d}{curdt.month:02d}{curdt.day:02d}"
    name += f"_DT.{curdt.hour:02d}00"
    name += "_DF.NC"
    if len(os.path.basename(name)) > 128:
        print("[ERR] Output file name is too long!")
        print(f"[ERR] {os.path.basename(name)} exceeds 128 characters!")
        sys.exit(1)
    return name

def _merge_files(ldtfile, noahmp_file, hymap2_file, merge_file):
    """Copy LDT, NoahMP and HYMAP2 fields into same file."""

    #rootgrp_hymap2 = nc4_dataset(hymap2_file, "r")
    #rootgrb_ldt = nc4_dataset(ldtfile, "r")

    # Copy data from NoahMP to Merged file
    src1 = nc4_dataset(noahmp_file, "r")
    src2 = nc4_dataset(hymap2_file, "r")
    src3 = nc4_dataset(ldtfile, "r")
    dst = nc4_dataset(merge_file, "w", format="NETCDF4")

    # Define all dimensions, variables, and attributes from src1
    dimension_dict = {
        "east_west" : "lon",
        "north_south" : "lat",
    }
    for dimname in src1.dimensions:
        dimname1 = dimname
        if dimname in dimension_dict:
            dimname1 = dimension_dict[dimname]
        dst.createDimension(dimname1, src1.dimensions[dimname].size)
    for gattrname in src1.__dict__:
        dst.setncattr(gattrname, src1.__dict__[gattrname])
    for name, variable in src1.variables.items():

        # Special handling for lat
        if name == "lat":
            dst.createVariable(name, variable.datatype, ("lat"))
        elif name == "lon":
            dst.createVariable(name, variable.datatype, ("lon"))
        else:
            # Need to account for new CF dimension names
            dimensions = []
            for dimension in variable.dimensions:
                if dimension in dimension_dict:
                    dimensions.append(dimension_dict[dimension])
                else:
                    dimensions.append(dimension)
            dst.createVariable(name, variable.datatype,
                               dimensions)
        # Extra CF attributes
        attrs = copy.deepcopy(src1[name].__dict__)
        if name == "time":
           attrs["calendar"] = "standard"
           attrs["axis"] = "T"
           attrs["bounds"] = "time_bnds"
        elif name == "lat":
            attrs["axis"] = "Y"
        elif name == "lon":
            attrs["axis"] = "X"
        if name in _cell_methods:
            attrs["cell_methods"] = _cell_methods[name]
        dst[name].setncatts(attrs)

    # Add select variables and attributes from src2
    src2_excludes = ["lat", "lon", "time", "ensemble", "RunoffStor_tavg",
                     "BaseflowStor_tavg"]
    for name, variable in src2.variables.items():
        if name in src2_excludes:
            continue
        dimensions = []
        for dimension in variable.dimensions:
            if dimension in dimension_dict:
                dimensions.append(dimension_dict[dimension])
            else:
                dimensions.append(dimension)
        dst.createVariable(name, variable.datatype,
                           dimensions)
        # Extra CF attributes
        attrs = copy.deepcopy(src2[name].__dict__)
        if name in _cell_methods:
            attrs["cell_methods"] = _cell_methods[name]
        dst[name].setncatts(attrs)

    # Add LANDMASK variable and attributes from src3
    dimensions = []
    for dimension in src3["LANDMASK"].dimensions:
        if dimension in dimension_dict:
            dimensions.append(dimension_dict[dimension])
        else:
            dimensions.append(dimension)
    dst.createVariable("LANDMASK", src3["LANDMASK"].datatype,
                       dimensions)
    dst["LANDMASK"].setncatts(src3["LANDMASK"].__dict__)

    # Add time_bnds variable
    dst.createDimension("nv", 2)
    dst.createVariable("time_bnds", "f4", ("time", "nv"))

    # Add soil_layer and soil_layer_thickness variables
    dst.createDimension("soil_layer",
                        src1.dimensions["SoilMoist_profiles"].size)
    dst.createVariable("soil_layer", "i4", ("soil_layer"))
    attrs = {
        "long_name" : "soil layer level",
        "axis" : "Z",
        "positive" : "down",
    }
    dst["soil_layer"].setncatts(attrs)
    dst.createVariable("soil_layer_thickness", "f4", ("soil_layer"))
    attrs = {
        "long_name" : "soil layer thicknesses",
        "units" : "m",
    }
    dst["soil_layer_thickness"].setncatts(attrs)

    # Write data from src1
    for name, variable in src1.variables.items():
        # Special handling for lat and lon, which should be 1d arrays
        # for lat/lon projection in CF convention
        if name == "lat":
            dst[name][:] = src1[name][:,0]
        elif name == "lon":
            dst[name][:] = src1[name][0,:]
        elif len(variable.dimensions) == 4:
            dst[name][:,:,:,:] = src1[name][:,:,:,:]
        elif len(variable.dimensions) == 3:
            dst[name][:,:,:] = src1[name][:,:,:]
        elif len(variable.dimensions) == 2:
            dst[name][:,:] = src1[name][:,:]
        elif len(variable.dimensions) == 1:
            dst[name][:] = src1[name][:]

    # Write data from src2
    for name, variable in src2.variables.items():
        if name in src2_excludes:
            continue
        if len(variable.dimensions) == 4:
            dst[name][:,:,:,:] = src2[name][:,:,:,:]
        elif len(variable.dimensions) == 3:
            dst[name][:,:,:] = src2[name][:,:,:]
        elif len(variable.dimensions) == 2:
            dst[name][:,:] = src2[name][:,:]
        elif len(variable.dimensions) == 1:
            dst[name][:] = src2[name][:]

    # Write data from src3
    dst["LANDMASK"][:,:] = src3["LANDMASK"][:,:]

    # Write time_bnds
    dst["time_bnds"][0,0] = -1440.
    dst["time_bnds"][0,1] = 0.

    # Write soil layer data
    dst["soil_layer"][0] = 1
    dst["soil_layer"][1] = 2
    dst["soil_layer"][2] = 3
    dst["soil_layer"][3] = 4
    dst["soil_layer_thickness"][0] = 0.1
    dst["soil_layer_thickness"][1] = 0.3
    dst["soil_layer_thickness"][2] = 0.6
    dst["soil_layer_thickness"][3] = 1.0

    src1.close()
    src2.close()
    src3.close()
    dst.close()

# Test driver
if __name__ == "__main__":
    ldtfile = "/discover/nobackup/projects/usaf_lis/emkemp/AFWA/s2spost/work/LDT/lis_input.s2s_global.noahmp401_hymap2.merit.25km.nc"
    noahmp_file = "/discover/nobackup/projects/usaf_lis/emkemp/AFWA/s2spost/work/CCM4/SURFACEMODEL/201912/LIS_HIST_201912310000.d01.nc"
    hymap2_file = "/discover/nobackup/projects/usaf_lis/emkemp/AFWA/s2spost/work/CCM4/ROUTING/201912/LIS_HIST_201912310000.d01.nc"

    merge_file = "/discover/nobackup/projects/usaf_lis/emkemp/AFWA/s2spost/work/output/test.nc"

    _merge_files(ldtfile, noahmp_file, hymap2_file, merge_file)
