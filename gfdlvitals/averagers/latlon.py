""" lat-lon grid processing tools """

import multiprocessing
import numpy as np

import gfdlvitals.util.gmeantools as gmeantools

__all__ = ["process_var", "average"]

FS = None
F = None
FYEAR = None
OUTDIR = None
LABEL = None
GEOLON = None
GEOLAT = None
CELL_AREA = None


def process_var(varname):
    """ routine to process one variable """
    units = gmeantools.extract_metadata(F, varname, "units")
    long_name = gmeantools.extract_metadata(F, varname, "long_name")
    if len(F.variables[varname].shape) == 3:
        var = F.variables[varname][:]
        var = np.ma.average(var, axis=0, weights=F.variables["average_DT"][:])
        for reg in ["global", "tropics", "nh", "sh"]:
            result, area_sum = gmeantools.area_mean(
                var, CELL_AREA, GEOLAT, GEOLON, region=reg
            )
            sqlfile = OUTDIR + "/" + FYEAR + "." + reg + "Ave" + LABEL + ".db"
            gmeantools.write_metadata(sqlfile, varname, "units", units)
            gmeantools.write_metadata(sqlfile, varname, "long_name", long_name)
            gmeantools.write_sqlite_data(sqlfile, varname, FYEAR[:4], result)
            gmeantools.write_sqlite_data(sqlfile, "area", FYEAR[:4], area_sum)


def average(f1, f2, year, out, lab):
    """ averaging function """
    global FS
    global F
    global FYEAR
    global OUTDIR
    global LABEL

    FS = f1
    F = f2
    FYEAR = year
    OUTDIR = out
    LABEL = lab

    # geometry
    global GEOLON
    global GEOLAT
    global CELL_AREA

    lat = FS.variables["lat"][:]
    lon = FS.variables["lon"][:]
    GEOLAT = np.tile(lat[:, None], (1, lon.shape[0]))
    GEOLON = np.tile(lon[None, :], (lat.shape[0], 1))
    CELL_AREA = gmeantools.standard_grid_cell_area(lat, lon)

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(process_var, F.variables.keys())
