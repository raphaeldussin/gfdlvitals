""" LM4-class land model averaging utilities """

import multiprocessing
import re

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
AREA_TYPES = None
CELL_DEPTH = None


def process_var(varname):
    """ process a variable """
    varshape = DATA_TILES[0].variables[varname].shape
    units = gmeantools.extract_metadata(DATA_TILES[0], varname, "units")
    long_name = gmeantools.extract_metadata(DATA_TILES[0], varname, "long_name")
    cell_measures = gmeantools.extract_metadata(DATA_TILES[0], varname, "cell_measures")
    area_measure = gmeantools.parse_cell_measures(cell_measures, "area")
    if (area_measure is not None) and (area_measure != "area_ntrl"):
        if len(varshape) >= 3:
            var = gmeantools.cube_sphere_aggregate(varname, DATA_TILES)
            var = np.ma.average(
                var, axis=0, weights=DATA_TILES[0].variables["average_DT"][:]
            )

            if len(varshape) == 3:
                for reg in ["global", "tropics", "nh", "sh"]:
                    result, area_sum = gmeantools.area_mean(
                        var, AREA_TYPES[area_measure], GEOLAT, GEOLON, region=reg
                    )
                    if not hasattr(result, "mask"):
                        sqlfile = (
                            OUTDIR + "/" + FYEAR + "." + reg + "Ave" + LABEL + ".db"
                        )
                        gmeantools.write_metadata(sqlfile, varname, "units", units)
                        gmeantools.write_metadata(
                            sqlfile, varname, "long_name", long_name
                        )
                        gmeantools.write_metadata(
                            sqlfile, varname, "cell_measure", area_measure
                        )
                        gmeantools.write_sqlite_data(
                            sqlfile, varname, FYEAR[:4], result
                        )
                        gmeantools.write_sqlite_data(
                            sqlfile, area_measure, FYEAR[:4], area_sum
                        )

            elif len(varshape) == 4:
                if varshape[1] == CELL_DEPTH.shape[0]:
                    for reg in ["global", "tropics", "nh", "sh"]:
                        result, volume_sum = gmeantools.area_mean(
                            var,
                            AREA_TYPES[area_measure],
                            GEOLAT,
                            GEOLON,
                            region=reg,
                            cell_depth=CELL_DEPTH,
                        )
                        sqlfile = (
                            OUTDIR + "/" + FYEAR + "." + reg + "Ave" + LABEL + ".db"
                        )
                        gmeantools.write_metadata(sqlfile, varname, "units", units)
                        gmeantools.write_metadata(
                            sqlfile, varname, "long_name", long_name
                        )
                        gmeantools.write_metadata(
                            sqlfile,
                            varname,
                            "cell_measure",
                            area_measure.replace("area", "volume"),
                        )
                        gmeantools.write_sqlite_data(
                            sqlfile, varname, FYEAR[:4], result
                        )
                        gmeantools.write_sqlite_data(
                            sqlfile,
                            area_measure.replace("area", "volume"),
                            FYEAR[:4],
                            volume_sum,
                        )


def average(gs_tl, da_tl, year, out, lab):
    """ averaging function """

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

    for f in [DATA_TILES, GS_TILES]:
        if "geolat_t" in f[0].variables:
            GEOLAT = gmeantools.cube_sphere_aggregate("geolat_t", DATA_TILES)
            GEOLON = gmeantools.cube_sphere_aggregate("geolon_t", DATA_TILES)
            break

    global AREA_TYPES

    AREA_TYPES = {}
    for f in [DATA_TILES, GS_TILES]:
        for v in sorted(f[0].variables):
            if re.match(r".*_area", v) or re.match(r"area.*", v):
                # for now, skip the area variables that depend on time
                timedependent = False
                for d in f[0].variables[v].dimensions:
                    timedependent = timedependent or f[0].dimensions[d].isunlimited()
                if not timedependent:
                    if v not in AREA_TYPES.keys():
                        AREA_TYPES[v] = gmeantools.cube_sphere_aggregate(v, f)

    global CELL_DEPTH

    depth = DATA_TILES[0].variables["zhalf_soil"][:]
    CELL_DEPTH = []
    for i in range(1, len(depth)):
        thickness = round((depth[i] - depth[i - 1]), 2)
        CELL_DEPTH.append(thickness)
    CELL_DEPTH = np.array(CELL_DEPTH)

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(process_var, DATA_TILES[0].variables.keys())
