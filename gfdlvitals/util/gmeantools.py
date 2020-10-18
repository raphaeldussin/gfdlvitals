""" Generic Suite of Utilities """

import os
import pickle
import sqlite3
import sys

import netCDF4 as nc
import numpy as np

__all__ = [
    "get_web_vars_dict",
    "ncopen",
    "mask_latitude_bands",
    "area_mean",
    "legacy_area_mean",
    "cube_sphere_aggregate",
    "write_sqlite_data",
    "parse_cell_measures",
    "extract_metadata",
    "write_metadata",
]


def get_web_vars_dict():
    """ Reads LM3 variable dictionary """
    return pickle.load(
        open(
            "/home/fms/local/opt/fre-analysis/test/eem/code/cm4_web_analysis/"
            + "etc/LM3_variable_dictionary.pkl",
            "rb",
        )
    )


def ncopen(file, action="exit"):
    """ opens NetCDF file """
    if os.path.exists(file):
        return nc.Dataset(file)
    print("WARNING: Unable to open file " + file)
    if action == "exit":
        sys.exit(0)
    else:
        return None


def mask_latitude_bands(var, cell_area, geolat, region="global"):
    """ Masks field by latitude range """
    if region == "tropics":
        var = np.ma.masked_where(np.logical_or(geolat < -30.0, geolat > 30.0), var)
        cell_area = np.ma.masked_where(
            np.logical_or(geolat < -30.0, geolat > 30.0), cell_area
        )
    elif region == "nh":
        var = np.ma.masked_where(np.less_equal(geolat, 30.0), var)
        cell_area = np.ma.masked_where(np.less_equal(geolat, 30.0), cell_area)
    elif region == "sh":
        var = np.ma.masked_where(np.greater_equal(geolat, -30.0), var)
        cell_area = np.ma.masked_where(np.greater_equal(geolat, -30.0), cell_area)
    elif region != "global":
        raise ValueError("Unknown region specified.")
    return var, cell_area


def area_mean(
    var, cell_area, geolat, geolon, region="global", cell_depth=None,
):
    """ computes area mean scalars """
    if cell_depth is not None:
        if var.shape[0] == cell_depth.shape[0]:
            cell_area = np.tile(cell_area[None, :], (cell_depth.shape[0], 1, 1))
            geolat = np.tile(geolat[None, :], (cell_depth.shape[0], 1, 1))
            geolon = np.tile(geolon[None, :], (cell_depth.shape[0], 1, 1))
        else:
            print(
                "Warning: inconsisent dimensions between varname and the cell depth axis.",
                var.shape[0],
                cell_depth.shape[0],
            )
            null_result = np.ma.masked_where(True, 0.0)
            return null_result, null_result
    cell_area = np.ma.array(cell_area)
    cell_area.mask = var.mask
    var, cell_area = mask_latitude_bands(var, cell_area, geolat, region)
    if cell_depth is not None:
        summed = np.ma.sum(
            var
            * cell_area
            * np.tile(cell_depth[:, None, None], (1, var.shape[1], var.shape[2]))
        )
        var = np.ma.average(var, axis=0, weights=cell_depth)
        res = np.ma.sum(var * cell_area) / cell_area.sum()
        return res, summed.sum()
    res = np.ma.sum(var * cell_area) / cell_area.sum()
    return res, cell_area.sum()


def legacy_area_mean(
    var,
    cell_area,
    geolat,
    geolon,
    cell_frac=None,
    soil_frac=None,
    region="global",
    varname=None,
    cell_depth=None,
    component=None,
):
    """ Legacy version of computing area mean scalars"""
    # Land-specific modifications
    if component == "land":
        module_dict = get_web_vars_dict()
        # Read dictionary of keys
        if varname in module_dict.keys():
            module = module_dict[varname]
        elif varname.lower() in module_dict.keys():
            module = module_dict[varname.lower()]
        else:
            module = ""
        # Create a weighting factor
        if module == "vegn":
            cell_area = cell_area * cell_frac * soil_frac
        else:
            cell_area = cell_area * cell_frac
        # Create a 3-D mask if needed
        if cell_depth is not None:
            if var.shape[0] == cell_depth.shape[0]:
                cell_area = np.tile(cell_area[None, :], (cell_depth.shape[0], 1, 1))
                geolat = np.tile(geolat[None, :], (cell_depth.shape[0], 1, 1))
                geolon = np.tile(geolon[None, :], (cell_depth.shape[0], 1, 1))
            else:
                print(
                    "Warning: inconsisent dimensions between varname and the cell depth axis.",
                    var.shape[0],
                    cell_depth.shape[0],
                )
                null_result = np.ma.masked_where(True, 0.0)
                return null_result, null_result
        # Apply data mask to weighting mask
        cell_area.mask = var.mask
    var, cell_area = mask_latitude_bands(var, cell_area, geolat, region)
    # -- Land depth averaging and summation
    if cell_depth is not None:
        summed = np.ma.sum(
            var
            * cell_area
            * np.tile(cell_depth[:, None, None], (1, var.shape[1], var.shape[2]))
        )
        var = np.ma.average(var, axis=0, weights=cell_depth)
        res = np.ma.sum(var * cell_area) / cell_area.sum()
        return res, summed
    res = np.ma.sum(var * cell_area) / cell_area.sum()
    summed = np.ma.sum(var * cell_area)
    return res, summed


def cube_sphere_aggregate(var, tiles):
    """ Aggregrates distributed cubed-sphere tiles """
    return np.ma.concatenate(
        (
            tiles[0].variables[var][:],
            tiles[1].variables[var][:],
            tiles[2].variables[var][:],
            tiles[3].variables[var][:],
            tiles[4].variables[var][:],
            tiles[5].variables[var][:],
        ),
        axis=-1,
    )


def write_sqlite_data(
    sqlfile, varname, fyear, varmean=None, varsum=None, component=None
):
    """ Writes values to sqlite file """
    conn = sqlite3.connect(sqlfile)
    c = conn.cursor()
    if component == "land":
        sql = (
            "create table if not exists "
            + varname
            + " (year integer primary key, sum float, avg float)"
        )
    else:
        sql = (
            "create table if not exists "
            + varname
            + " (year integer primary key, value float)"
        )
    c.execute(sql)
    if component == "land":
        sql = (
            "insert or replace into "
            + varname
            + " values("
            + fyear[:4]
            + ","
            + str(varsum)
            + ","
            + str(varmean)
            + ")"
        )
    else:
        sql = (
            "insert or replace into "
            + varname
            + " values("
            + fyear[:4]
            + ","
            + str(varmean)
            + ")"
        )
    c.execute(sql)
    conn.commit()
    c.close()
    conn.close()


def parse_cell_measures(attr, key):
    """ Determines cell measures """
    if attr is not None:
        ind = attr.split().index(key + ":") + 1
        return attr.split()[ind]
    return None


def extract_metadata(f, varname, attr):
    """ Extracts variable metadata"""
    if attr in f.variables[varname].__dict__.keys():
        return f.variables[varname].__dict__[attr]
    return None


def write_metadata(sqlfile, varname, attr, value):
    """ Writes metadata to a sqlite file """
    if value is None:
        value = str("")
    conn = sqlite3.connect(sqlfile)
    c = conn.cursor()
    sql = (
        "create table if not exists "
        + str(attr)
        + " (var text primary key, value text)"
    )
    c.execute(sql)
    sql = (
        "insert or replace into "
        + str(attr)
        + ' values("'
        + str(varname)
        + '","'
        + str(value)
        + '")'
    )
    c.execute(sql)
    conn.commit()
    c.close()
    conn.close()


def standard_grid_cell_area(lat, lon, rE=6371.0e3):
    """ Computes area file for a standard grid """
    dlat = lat[1] - lat[0]
    dlon = lon[1] - lon[0]
    area = np.empty((len(lat), len(lon)))
    for j, _lat in enumerate(lat):
        for i, _lon in enumerate(lon):
            lon1 = _lon + dlon / 2.0
            lon0 = _lon - dlon / 2.0
            lat1 = _lat + dlat / 2.0
            lat0 = _lat - dlat / 2.0
            area[j, i] = (
                (np.pi / 180.0)
                * rE
                * rE
                * np.abs(np.sin(np.radians(lat0)) - np.sin(np.radians(lat1)))
                * np.abs(lon0 - lon1)
            )
    return area
