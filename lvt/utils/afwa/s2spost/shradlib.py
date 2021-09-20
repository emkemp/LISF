#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: shradlib.py
#
# PURPOSE: Contains function used by S2S postprocessing scripts.  Based on
# Shrad_modules.py written by Shrad Shukula.
#
# REVISION HISTORY:
# 20 Sep 2021: Eric Kemp (SSAI), first version.  Excludes plotmap and
# write_netcdf_files, which do not appear to be used, are difficult to
# refactor to pylint's expectations, and (in plotmap's case) rely on orphaned
# plotting package (Basemap).
#
#------------------------------------------------------------------------------
"""

# Standard libraries
import errno
import os
from subprocess import check_call

# Third-party libraries
# NOTE: pylint cannot see the classes in netCDF4 since the latter is not
# written in Python.  We therefore disable a check for these line to avoid
# a known false alarm.
# pylint: disable=no-name-in-module
from netCDF4 import Dataset as nc4_dataset
# pylint: enable=no-name-in-module
import numpy as np
from scipy.stats import pearsonr

def calc_r(data1, data2):
    """Calculates Pearson's correlation coefficient between observed and
    forecast anomalies.  Note that forecast anomalies are calculated
    using mean forecast values."""
    corr = np.ones((data1.shape[1], data1.shape[2])) * -99
    for i in range(data1.shape[1]):
        for j in range(data1.shape[2]):
            if np.isfinite(data1[0, i, j]) and np.mean(data1[:, i, j] != 0):
                obs_anom = data1[:, i, j]-np.mean(data1[:, i, j])
                fcst_anom = data2[:, i, j]-np.mean(data2[:, i, j])
                ## calculating correlaiton between anomaly
                corr[i, j] = pearsonr(obs_anom, fcst_anom)[0, ]
    return corr

def makedir(path):
    """Create directory, and ensure mode is correct if preexisting."""
    try:
        os.makedirs(path, 493)
    except OSError as err:
        if err.errno == errno.EEXIST:  # file exists error?
            print('warning:', path, ' exists')
        else:
            raise  # re-raise the exception
        # make sure its mode is right
        os.chmod(path, 493)

def read_nc_files(infile, varname):
    """Extracts data from named variable from netCDF file."""
    ncid = nc4_dataset(infile, 'r')
    data = ncid.variables[varname][:]
    ncid.close()
    return data

def run_cdo(command, infile, outfile):
    """Run the cbo binary"""
    check_call(['cdo', command, infile, outfile])

# def plotmap(var, lats, lons, lat_min, lat_max, lat_int,
#             lon_min, lon_max, lon_int,
#             cmap, cbmin, cbmax, plottype, land_color,
#             ocean_color, levels, labels_1, labels_2,
#             col_under, col_higher, extend, area_thresh=10000):
#     """Creates spatial plots"""

#     if plottype == 'mesh':
#         lonres = (lons.max() - lons.min()) / (len(lons) - 1)
#         latres = (lats.max() - lats.min()) / (len(lats) - 1)
#         lats = lats - 0.5*latres
#         lons = lons - 0.5*lonres
#     elif plottype == 'contour':
#         pass
#     elif plottype == 'color':
#         pass
#     else:
#         raise ValueError('{}: Not a valid option for plottype'.
#                          format(plottype))

#     m = Basemap(projection='cyl', llcrnrlat=lat_min, llcrnrlon=lon_min,
#                 urcrnrlat=lat_max, urcrnrlon=lon_max, rsphere=6371200.,
#                 resolution='i',
#                 area_thresh=area_thresh)
#     xi, yi = m(lons, lats)
#     xi, yi = np.meshgrid(xi, yi)

#     if plottype == 'mesh':
#         cs = m.pcolormesh(xi, yi, var, vmax=cbmax, vmin=cbmin, cmap=cmap)
#     elif plottype == 'contour':
#         cmap = cmap
#         cmap.set_under(col_under)
#         cmap.set_over(col_higher)
#         norm = mpl.pyplot.cm.colors.Normalize(vmax=cbmax, vmin=cbmin,
#                                               clip=False)
#         cs = m.contourf(xi, yi, var, levels, cmap=cmap, xnorm=norm,
#                         extend=extend)
#     elif plottype == 'color':
#         cs = m.pcolor(xi, yi, var, vmax=cbmax, vmin=cbmin, cmap=cmap)

#     #m.drawlsmask(land_color=land_color, ocean_color=ocean_color, lakes=True)
#     m.drawlsmask(ocean_color=ocean_color, lakes=True)
#     m.drawparallels(np.arange(-80, 81, lat_int), labels=labels_1, fontsize=8,
#                     linewidth=0.3, rotation=45)
#     m.drawmeridians(np.arange(0, 360, lon_int),
#                     labels=labels_2, fontsize=8, linewidth=0.3)
#     m.drawcoastlines(linewidth=0.75)
#     m.drawcountries(linewidth=0.75)
#     mpl.rc('xtick', labelsize=11)
#     mpl.rc('ytick', labelsize=11)
#     #if (colorbar):
#         #cbar=m.colorbar(cs, location=colorbar_location, shrink=shrink, extend=extend, pad=pad)
#     #if (label_flag):
#      #   cbar.set_label(cbar_label, fontsize=14)
#     return m, cs


# def write_netcdf_files(filename, var, varname,
#                        description, source, var_units, sig_digit,
#                        lons, lats, sdate, interval):
#     """Writes a netCDF file."""
#     rootgrp = nc4_dataset(file, 'w', format='NETCDF4_CLASSIC')
#     rootgrp.createDimension('longitude', len(lons))
#     rootgrp.createDimension('latitude', len(lats))
#     rootgrp.createDimension('time', None)
#     longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
#     latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
#     times = rootgrp.createVariable('time', 'f8', ('time',))
#     # two dimensions unlimited.
#     varname = rootgrp.createVariable(varname, 'f4',
#                                      ('time', 'latitude', 'longitude',),
#                                      fill_value=nc4_default_fillvals['f4'],
#                                      zlib=True,
#                                      least_significant_digit=sig_digit)
#     # varname = \
#     #     rootgrp.createVariable(varname,'f4',('time','latitude','longitude',), \
#     #                            fill_value=nc.default_fillvals['f4'])

#     rootgrp.description = description
#     rootgrp.history = 'Created ' + time.ctime(time.time())
#     rootgrp.source = source
#     latitudes.units = 'degrees_north'
#     longitudes.units = 'degrees_east'
#     varname.units = var_units
#     string_date = datetime.strftime(sdate, "%Y-%m-%d")
#     times.units = 'days since ' + string_date
#     times.calendar = 'gregorian'
#     latitudes[:] = lats
#     longitudes[:] = lons
#     varname[:, :, :] = var
#     # Interval below is number of days between two dates
#     dates = [sdate + n * timedelta(days=interval) \
#              for n in range(varname.shape[0])]
#     times[:] = nc4_date2num(dates, units=times.units, calendar=times.calendar)
#     rootgrp.close()
