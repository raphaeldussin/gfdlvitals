""" cubed-sphere processing utilities """

import multiprocessing
import numpy as np

import gfdlvitals.util.gmeantools as gmeantools

__all__ = ["process_var", "average"]

GS_TILES = None
DATA_TILES = None
FYEAR = None
OUTDIR = None
LABEL = None
GEOLON = None
GEOLAT = None
CELL_AREA = None


def process_var(varname):
    """ routine to process one variable """
    units = gmeantools.extract_metadata(DATA_TILES[0], varname, "units")
    long_name = gmeantools.extract_metadata(DATA_TILES[0], varname, "long_name")
    if len(DATA_TILES[0].variables[varname].shape) == 3:
        var = gmeantools.cube_sphere_aggregate(varname, DATA_TILES)
        var = np.ma.average(
            var, axis=0, weights=DATA_TILES[0].variables["average_DT"][:]
        )
        for reg in ["global", "tropics", "nh", "sh"]:
            result, area_sum = gmeantools.area_mean(
                var, CELL_AREA, GEOLAT, GEOLON, region=reg
            )
            sqlfile = OUTDIR + "/" + FYEAR + "." + reg + "Ave" + LABEL + ".db"
            gmeantools.write_metadata(sqlfile, varname, "units", units)
            gmeantools.write_metadata(sqlfile, varname, "long_name", long_name)
            gmeantools.write_sqlite_data(sqlfile, varname, FYEAR[:4], result)
            gmeantools.write_sqlite_data(sqlfile, "area", FYEAR[:4], area_sum)


def average(gs_tl, da_tl, year, out, lab):
    """ function to do averaging """
    global GS_TILES
    global DATA_TILES
    global FYEAR
    global OUTDIR
    global LABEL

    GS_TILES = gs_tl
    DATA_TILES = da_tl
    FYEAR = year
    OUTDIR = out
    LABEL = lab

    global GEOLON
    global GEOLAT
    global CELL_AREA

    GEOLAT = gmeantools.cube_sphere_aggregate("grid_latt", GS_TILES)
    GEOLON = gmeantools.cube_sphere_aggregate("grid_lont", GS_TILES)
    CELL_AREA = gmeantools.cube_sphere_aggregate("area", GS_TILES)

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(process_var, DATA_TILES[0].variables.keys())
