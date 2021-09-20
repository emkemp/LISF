#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: postprocess_NMME.py
#
# PURPOSE: Drives postprocessing of LIS forecasts forced by NMME members.
# Based on Postprocess_NMME_job.sh and
# job_run_convert_Dyn_FCST_to_postproc.scr.
#
# REVISION HISTORY:
# 20 Sep 2021: Eric Kemp (SSAI), first version.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import datetime
import os
import subprocess
import sys

# Hardwired system settings.
# Put into config file. FIXME
_RUNDIR = \
  "/discover/nobackup/projects/fame/FORECASTS/GEOS5/BCSD_Test/" + \
  "NMME_FCST_CODES_AF"

_CSYR=1991 # May need to change
_CEYR=2020

_MODELS = ["CFSv2", "CCSM4", "CCM4", "GEOSv2", "GNEMO", "GFDL"]

_BATCH_SCRIPT = "/discover/nobackup/projects/fame/FORECASTS/GEOS5/" + \
    "BCSD_Test/NMME_FCST_CODES_AF/run_Convert_Dyn_FCST_postproc.scr"

_PY_SCRIPTS = ["Convert_Dyn_FCST_to_anom_SURFACEMODEL.py",
               "Convert_Dyn_FCST_to_sanom_SURFACEMODEL.py",
               "Convert_Dyn_FCST_to_anom_ROUTING.py",
               "Convert_Dyn_FCST_to_sanom_ROUTING.py"]

def _driver():
    """Main driver for postprocessing."""

    # Set month and year for prior calendar month.
    curdate = datetime.date.today()
    curmonth = curdate.month
    curyear = curdate.day

    delta = datetime.timedelta(days=1)
    prevdate = datetime.date(year=curyear, month=curmonth, day=1) - delta
    prevmonth = prevdate.month
    prevyear = prevdate.year

    # Process dynamical forecasts
    os.chdir(_RUNDIR)
    for model in _MODELS:
        for py_script in _PY_SCRIPTS:
            cmd = "sbatch %s %s %s %s %s %s %s %s %s %s %s %s" \
                %(_BATCH_SCRIPT, py_script, curmonth, "NOAHMP", 5,
                  "FAME_DOMAIN", curyear, "FirstLook", "NMME_BC",
                  model, _CSYR, _CEYR)
            returncode = subprocess.call(cmd, shell=True)
            if returncode != 0:
                print("[ERR] Problem with %s" %(py_script))
                sys.exit(1)
    print("[INFO] Model processing submitted.")

    # Process preliminary initial conditions.
    cmd = "run_Convert_Prelim_IC_postproc.scr %s %s %s %s %s %s %s" \
        %(prevmonth, "NOAHMP", "FAME_Domain", prevyear, "FirstLook",
          _CSYR, _CEYR)
    returncode = subprocess.call(cmd, shell=True)
    if returncode != 0:
        print("[ERR] Problem running run_Convert_Prelim_IC_postproc.scr")
        sys.exit(1)

# Invoke driver
if __name__ == "__main__":
    _driver()
