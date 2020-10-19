""" tripolar grid averaging utilities """

import multiprocessing
import numpy as np

import gfdlvitals.util.gmeantools as gmeantools

__all__ = ["process_var", "average"]

FGS = None
FDATA = None
FYEAR = None
OUTDIR = None
LABEL = None
GEOLON = None
GEOLAT = None
CELL_AREA = None


def process_var(varname):
    """ function to process a variable """

    units = gmeantools.extract_metadata(FDATA, varname, "units")
    long_name = gmeantools.extract_metadata(FDATA, varname, "long_name")
    ndims = len(FDATA.variables[varname].shape)
    if ndims >= 3:
        if ndims == 3:
            dims = FDATA.variables[varname].dimensions
            if dims[-2] == "yh" and dims[-1] == "xh":
                var = FDATA.variables[varname][:]
            else:
                return
        elif (ndims == 4) and (varname[0:9] == "tot_layer"):
            var = FDATA.variables[varname][:]
            var = np.ma.sum(var, axis=1)
        else:
            return
        var = np.ma.average(var, axis=0, weights=FDATA.variables["average_DT"][:])
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
    """ function to do the averaging """

    global FGS
    global FDATA
    global FYEAR
    global OUTDIR
    global LABEL

    FGS = f1
    FDATA = f2
    FYEAR = year
    OUTDIR = out
    LABEL = lab

    # geometry
    global GEOLON
    global GEOLAT
    global CELL_AREA

    if "geolat" in FGS.variables.keys():
        GEOLAT = FGS.variables["geolat"][:]
    elif "geolat_t" in FGS.variables.keys():
        GEOLAT = FGS.variables["geolat_t"][:]
    else:
        raise ValueError("Unable to determine geolat.")

    if "geolon" in FGS.variables.keys():
        GEOLON = FGS.variables["geolon"][:]
    elif "geolon_t" in FGS.variables.keys():
        GEOLON = FGS.variables["geolon_t"][:]
    else:
        raise ValueError("Unable to determine GEOLON.")

    if "areacello" in FGS.variables.keys():
        CELL_AREA = FGS.variables["areacello"][:]
    elif "area_t" in FGS.variables.keys():
        CELL_AREA = FGS.variables["area_t"][:]
    else:
        raise ValueError("Unable to determine ocean cell area.")

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(process_var, FDATA.variables.keys())
