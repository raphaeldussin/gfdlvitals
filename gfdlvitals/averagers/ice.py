import numpy as np
import multiprocessing

import gfdlvitals.util.gmeantools as gmeantools

__all__ = ["process_var", "average"]


def process_var(v):
    if fdata.variables[v].shape == cellArea.shape:
        units = gmeantools.extract_metadata(fdata, v, "units")
        long_name = gmeantools.extract_metadata(fdata, v, "long_name")
        data = fdata.variables[v][:]
        for reg in ["global", "nh", "sh"]:
            sqlite_out = outdir + "/" + fYear + "." + reg + "Ave" + label + ".db"
            _v, _area = gmeantools.mask_latitude_bands(
                data, cellArea, geoLat, region=reg
            )
            _v = np.ma.sum((_v * _area), axis=(-1, -2)) / np.ma.sum(
                _area, axis=(-1, -2)
            )
            gmeantools.write_metadata(sqlite_out, v, "units", units)
            gmeantools.write_metadata(sqlite_out, v, "long_name", long_name)
            gmeantools.write_sqlite_data(
                sqlite_out,
                v + "_mean",
                fYear[:4],
                np.ma.average(_v, axis=0, weights=average_DT),
            )
            gmeantools.write_sqlite_data(
                sqlite_out, v + "_max", fYear[:4], np.ma.max(_v)
            )
            gmeantools.write_sqlite_data(
                sqlite_out, v + "_min", fYear[:4], np.ma.min(_v)
            )


def average(f1, f2, year, out, lab):
    global fgs
    global fdata
    global fYear
    global outdir
    global label

    fgs = f1
    fdata = f2
    fYear = year
    outdir = out
    label = lab

    # geometry
    global geoLon
    global geoLat
    global cellArea

    geoLon = fgs.variables["GEOLON"][:]
    geoLat = fgs.variables["GEOLAT"][:]

    global average_DT
    average_DT = fdata.variables["average_DT"][:]

    if "CELL_AREA" in fgs.variables.keys():
        rE = 6371.0e3  # Radius of the Earth in 'm'
        cellArea = fgs.variables["CELL_AREA"][:] * (4.0 * np.pi * (rE ** 2))
    elif "area" in fgs.variables.keys():
        cellArea = fgs.variables["area"][:]
    else:
        print("FATAL: unable to determine cell area used in ice model")

    if "siconc" in fdata.variables.keys():
        concentration = fdata.variables["siconc"][:]
    elif "CN" in fdata.variables.keys():
        concentration = np.ma.sum(fdata.variables["CN"][:], axis=-3)
    else:
        print("FATAL: unable to determine ice concentration")

    geoLat = np.tile(geoLat[None, :], (concentration.shape[0], 1, 1))
    geoLon = np.tile(geoLon[None, :], (concentration.shape[0], 1, 1))
    cellArea = np.tile(cellArea[None, :], (concentration.shape[0], 1, 1))

    for reg in ["global", "nh", "sh"]:
        sqlite_out = outdir + "/" + fYear + "." + reg + "Ave" + label + ".db"
        vars = []
        # area and extent in million square km
        _conc, _area = gmeantools.mask_latitude_bands(
            concentration, cellArea, geoLat, region=reg
        )
        vars.append(("area", (np.ma.sum((_conc * _area), axis=(-1, -2)) * 1.0e-12)))
        vars.append(
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
        for v in vars:
            gmeantools.write_sqlite_data(
                sqlite_out,
                v[0] + "_mean",
                fYear[:4],
                np.ma.average(v[1], weights=average_DT),
            )
            gmeantools.write_sqlite_data(
                sqlite_out, v[0] + "_max", fYear[:4], np.ma.max(v[1])
            )
            gmeantools.write_sqlite_data(
                sqlite_out, v[0] + "_min", fYear[:4], np.ma.min(v[1])
            )

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(process_var, fdata.variables.keys())
