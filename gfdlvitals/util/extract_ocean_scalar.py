""" utilities for extracting ocean scalar fields """

from . import gmeantools

__all__ = ["mom6"]


def mom6(fdata, fyear, outdir):
    """ Reads MOM6 ocean scalar file """

    ignore_list = ["time_bounds", "time_bnds", "average_T2", "average_T1", "average_DT"]

    var_dict = fdata.variables.keys()
    var_dict = list(set(var_dict) - set(ignore_list))

    for varname in var_dict:
        if len(fdata.variables[varname].shape) == 2:
            units = gmeantools.extract_metadata(fdata, varname, "units")
            long_name = gmeantools.extract_metadata(fdata, varname, "long_name")
            result = fdata.variables[varname][0, 0]
            sqlfile = outdir + "/" + fyear + ".globalAveOcean.db"
            gmeantools.write_metadata(sqlfile, varname, "units", units)
            gmeantools.write_metadata(sqlfile, varname, "long_name", long_name)
            gmeantools.write_sqlite_data(sqlfile, varname, fyear[:4], result)
