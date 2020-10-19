""" tripolar ice averaging utilities """

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
average_DT = None


def process_var(v):
    """ routine to process a variable """
    if FDATA.variables[v].shape == CELL_AREA.shape:
        units = gmeantools.extract_metadata(FDATA, v, "units")
        long_name = gmeantools.extract_metadata(FDATA, v, "long_name")
        data = FDATA.variables[v][:]
        for reg in ["global", "nh", "sh"]:
            sqlite_out = OUTDIR + "/" + FYEAR + "." + reg + "Ave" + LABEL + ".db"
            _v, _area = gmeantools.mask_latitude_bands(
                data, CELL_AREA, GEOLAT, region=reg
            )
            _v = np.ma.sum((_v * _area), axis=(-1, -2)) / np.ma.sum(
                _area, axis=(-1, -2)
            )
            gmeantools.write_metadata(sqlite_out, v, "units", units)
            gmeantools.write_metadata(sqlite_out, v, "long_name", long_name)
            gmeantools.write_sqlite_data(
                sqlite_out,
                v + "_mean",
                FYEAR[:4],
                np.ma.average(_v, axis=0, weights=average_DT),
            )
            gmeantools.write_sqlite_data(
                sqlite_out, v + "_max", FYEAR[:4], np.ma.max(_v)
            )
            gmeantools.write_sqlite_data(
                sqlite_out, v + "_min", FYEAR[:4], np.ma.min(_v)
            )


def average(f1, f2, year, out, lab):
    """ averaging function """

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

    GEOLON = FGS.variables["GEOLON"][:]
    GEOLAT = FGS.variables["GEOLAT"][:]

    global average_DT
    average_DT = FDATA.variables["average_DT"][:]

    if "CELL_AREA" in FGS.variables.keys():
        rE = 6371.0e3  # Radius of the Earth in 'm'
        CELL_AREA = FGS.variables["CELL_AREA"][:] * (4.0 * np.pi * (rE ** 2))
    elif "area" in FGS.variables.keys():
        CELL_AREA = FGS.variables["area"][:]
    else:
        print("FATAL: unable to determine cell area used in ice model")

    if "siconc" in FDATA.variables.keys():
        concentration = FDATA.variables["siconc"][:]
    elif "CN" in FDATA.variables.keys():
        concentration = np.ma.sum(FDATA.variables["CN"][:], axis=-3)
    else:
        print("FATAL: unable to determine ice concentration")

    GEOLAT = np.tile(GEOLAT[None, :], (concentration.shape[0], 1, 1))
    GEOLON = np.tile(GEOLON[None, :], (concentration.shape[0], 1, 1))
    CELL_AREA = np.tile(CELL_AREA[None, :], (concentration.shape[0], 1, 1))

    for reg in ["global", "nh", "sh"]:
        sqlite_out = OUTDIR + "/" + FYEAR + "." + reg + "Ave" + LABEL + ".db"
        variables = []
        # area and extent in million square km
        _conc, _area = gmeantools.mask_latitude_bands(
            concentration, CELL_AREA, GEOLAT, region=reg
        )
        variables.append(
            ("area", (np.ma.sum((_conc * _area), axis=(-1, -2)) * 1.0e-12))
        )
        variables.append(
            (
                "extent",
                (
                    np.ma.sum(
                        (np.ma.where(np.greater(_conc, 0.15), _area, 0.0)),
                        axis=(-1, -2),
                    )
                    * 1.0e-12
                ),
            )
        )
        for v in variables:
            gmeantools.write_sqlite_data(
                sqlite_out,
                v[0] + "_mean",
                FYEAR[:4],
                np.ma.average(v[1], weights=average_DT),
            )
            gmeantools.write_sqlite_data(
                sqlite_out, v[0] + "_max", FYEAR[:4], np.ma.max(v[1])
            )
            gmeantools.write_sqlite_data(
                sqlite_out, v[0] + "_min", FYEAR[:4], np.ma.min(v[1])
            )

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(process_var, FDATA.variables.keys())
