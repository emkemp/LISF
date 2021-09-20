#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT: s2spostlib.py
#
# PURPOSE: Contains function used by S2S postprocessing scripts.  Based on
# All_functions.py.
#
# REVISION HISTORY:
# 20 Sep 2021: Eric Kemp (SSAI), first version.
#
#------------------------------------------------------------------------------
"""

# Third-party libraries
from scipy.stats import percentileofscore as pscore

def calc_empirical_pctl(target_val, obs_clim, mean_type):
    """Calculate percentile score."""
    pctl = pscore(obs_clim, target_val, kind=mean_type)
    return pctl

def _sel_var_rootzone_sm(sel_cim_data, model):
    """Selects climatology for RootZone-SM"""
    if model == 'CLSM':
        # CLSM Rootzone Soil Moisture is in 2nd layer
        var_sel_clim_data = \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=1)
    elif model in ['NOAHMP', 'NoahMP']:
        var_sel_clim_data = \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=0) * \
            0.1 + \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=1) * \
            0.3 + \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=2) * 0.6
    return var_sel_clim_data

def _sel_var_total_sm(sel_cim_data, model):
    """Selects climatology for Total-SM"""
    if model == 'CLSM':
        # Total CLSM soil moisture in 3rd layer
        var_sel_clim_data = \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=2)
    else:
        var_sel_clim_data = \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=0) * \
            0.05 + \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=1) * \
            0.15 + \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=2) * \
            0.3 + \
            sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=3) * 0.5

    return var_sel_clim_data

def sel_var (sel_cim_data, var_name, model):
    """Selects climatology from named variable and model."""
    # Now selecting the climatology further for the given variable
    if var_name == 'RootZone-SM':
        var_sel_clim_data = _sel_var_rootzone_sm(sel_cim_data, model)

    elif var_name == 'Total-SM':
        var_sel_clim_data = _sel_var_total_sm(sel_cim_data, model)

    elif var_name == 'Surface-SM':
        if model == 'CLSM':
            var_sel_clim_data = \
                sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=0)
        else:
            # Sum soil moisture from all layers
            var_sel_clim_data = \
                sel_cim_data.SoilMoist_tavg.isel(SoilMoist_profiles=0)

    elif var_name == 'Total-Runoff':
        # Merge total surface and subsurface runoff
        var_sel_clim_data = \
            sel_cim_data.Qs_tavg + \
            sel_cim_data.Qsb_tavg
    elif var_name == 'TWS':
        var_sel_clim_data = sel_cim_data.TWS_tavg
    elif var_name == 'Precip':
        var_sel_clim_data = sel_cim_data.TotalPrecip_tavg
    elif var_name == 'Air-T':
        var_sel_clim_data = sel_cim_data.Tair_f_tavg
    elif var_name == 'ET':
        var_sel_clim_data = sel_cim_data.Evap_tavg
    elif var_name == 'Streamflow':
        var_sel_clim_data = sel_cim_data.Streamflow_tavg

    return var_sel_clim_data


def get_boundary(region_name):
    """Get boundaries of given region name"""
    boundary_ea = (22, 55, -12, 23)
    boundary_wa = (-19, 26, -5, 25)
    ##boundary_sa = (8, 52, -37, 0)
    ##boundary_fame = (-20, 60, -40, 40)
    boundary_sa = (8, 52, -37, 6)
    boundary_fame = (-20, 55, -40, 40)
    boundary_sa1 = (24, 33, -31, -24)
    if region_name == 'EA':
        boundary = boundary_ea
    elif region_name == 'WA':
        boundary = boundary_wa
    elif region_name == 'SA':
        boundary = boundary_sa
    elif region_name == 'SA1':
        boundary = boundary_sa1
    else:
        boundary = boundary_fame
    return boundary
