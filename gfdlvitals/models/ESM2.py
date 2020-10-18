""" top-level module for ESM2-class models """

import tarfile
from gfdlvitals import averagers
from gfdlvitals.util.netcdf import extract_from_tar

__all__ = ["routines"]


def routines(infile):
    """ routines to process ESM2-class models """

    # -- Open the tarfile
    tar = tarfile.open(infile)
    members = tar.getnames()
    # -- Set the model year string
    fyear = str(infile.split("/")[-1].split(".")[0])
    print("Processing " + fyear)

    if members[-1][0:2] == "./":
        modifier = "./"
    else:
        modifier = ""

    # -- Land
    label = "Land"
    modules = "land_month"
    if modifier + fyear + ".land_static.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".land_static.nc")
    else:
        fgs = extract_from_tar(tar, modifier + fyear + ".land_month.nc")
    modules = ["land_month"]
    for module in modules:
        fname = modifier + fyear + "." + module + ".nc"
        if fname in members:
            fdata = extract_from_tar(tar, fname)
            print(fname)
            averagers.land_lm3.average(fgs, fdata, fyear, "./", label)
            fdata.close()
    fgs.close()

    # -- Atmosphere
    label = "Atmos"
    modules = ["atmos_month", "atmos_level"]
    for module in modules:
        # Need to figure out how to handle leading slash
        # fname = './'+fyear+'.'+module+'.nc'
        fname = modifier + fyear + "." + module + ".nc"
        if fname in members:
            ds = extract_from_tar(tar, fname)
            print(fname)
            averagers.latlon.average(ds, ds, fyear, "./", label)
            ds.close()

    # -- Ocean
    label = "Ocean"
    modules = ["ocean_month"]
    if modifier + fyear + ".ocean_static.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ocean_static.nc")
    else:
        fgs = extract_from_tar(tar, modifier + fyear + ".ocean_month.nc")
    for module in modules:
        fname = modifier + fyear + "." + module + ".nc"
        if fname in members:
            fdata = extract_from_tar(tar, fname)
            print(fname)
            averagers.tripolar.average(fgs, fdata, fyear, "./", label)
            fdata.close()
    fgs.close()

    # -- TOPAZ
    label = "TOPAZ"
    modules = [
        "ocean_topaz_fluxes",
        "ocean_topaz_misc",
        "ocean_topaz_sfc_100",
        "ocean_topaz_tracers_month_z",
        "ocean_topaz_wc_btm",
    ]
    if modifier + fyear + ".ocean_static.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ocean_static.nc")
    else:
        fgs = extract_from_tar(tar, modifier + fyear + ".ocean_month.nc")
    for module in modules:
        fname = modifier + fyear + "." + module + ".nc"
        if fname in members:
            fdata = extract_from_tar(tar, fname)
            print(fname)
            averagers.tripolar.average(fgs, fdata, fyear, "./", label)
            fdata.close()
    fgs.close()

    # -- Ice
    label = "Ice"
    modules = ["ice_month"]
    if modifier + fyear + ".ice_static.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ice_static.nc")
    else:
        fgs = extract_from_tar(tar, modifier + fyear + ".ice_month.nc")
    for module in modules:
        fname = modifier + fyear + "." + module + ".nc"
        if fname in members:
            fdata = extract_from_tar(tar, fname)
            print(fname)
            averagers.ice.average(fgs, fdata, fyear, "./", label)
            fdata.close()
    fgs.close()

    # -- Close out the tarfile handle
    tar.close()
