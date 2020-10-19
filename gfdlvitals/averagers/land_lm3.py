""" legacy land model averaging routines """

import multiprocessing
import numpy as np

import gfdlvitals.util.gmeantools as gmeantools

__all__ = ["process_var", "average"]

FS = None
FDATA = None
FYEAR = None
OUTDIR = None
LABEL = None
GEOLON = None
GEOLAT = None
CELL_AREA = None
CELL_FRAC = None
SOIL_AREA = None
SOIL_FRAC = None


def process_var(varname):
    """ process a single variable """
    varshape = FDATA.variables[varname].shape
    if len(varshape) >= 3:
        var = FDATA[varname][:]
        var = np.ma.average(var, axis=0, weights=FDATA["average_DT"][:])

    if len(varshape) == 3:
        for reg in ["global", "tropics", "nh", "sh"]:
            sqlfile = OUTDIR + "/" + FYEAR + "." + reg + "Ave" + LABEL + ".db"
            avg, summed = gmeantools.legacy_area_mean(
                var,
                CELL_AREA,
                GEOLAT,
                GEOLON,
                cell_frac=CELL_FRAC,
                soil_frac=SOIL_FRAC,
                region=reg,
                varname=varname,
                component="land",
            )
            if not hasattr(avg, "mask"):
                gmeantools.write_sqlite_data(
                    sqlfile,
                    varname,
                    FYEAR[:4],
                    varmean=avg,
                    varsum=summed,
                    component="land",
                )


def average(f1, f2, year, out, lab):
    """ averaging function """
    global FS
    global FDATA
    global FYEAR
    global OUTDIR
    global LABEL

    FS = f1
    FDATA = f2
    FYEAR = year
    OUTDIR = out
    LABEL = lab

    # geometry
    global GEOLON
    global GEOLAT

    lat = FDATA["lat"][:]
    lon = FDATA["lon"][:]
    GEOLON, GEOLAT = np.meshgrid(lon, lat)

    # land areas and fractions
    global CELL_AREA
    global CELL_FRAC
    global SOIL_AREA
    global SOIL_FRAC

    CELL_AREA = FS["land_area"][:]
    CELL_FRAC = FS["land_frac"][:]

    if "soil_area" in FS.variables.keys():
        SOIL_AREA = FS["soil_area"][0]
    elif "soil_area" in FDATA.variables.keys():
        SOIL_AREA = FDATA["soil_area"][0]
    else:
        raise ValueError("Unable to locate soil area field.")
    SOIL_FRAC = np.ma.array(SOIL_AREA / (CELL_AREA * CELL_FRAC))

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(process_var, FDATA.variables.keys())
