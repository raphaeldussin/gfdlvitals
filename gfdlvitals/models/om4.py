""" OM4 processing routines """

import tarfile

from gfdlvitals import averagers
from gfdlvitals import diags
from gfdlvitals.util import extract_ocean_scalar
from gfdlvitals.util.netcdf import extract_from_tar

__all__ = ["routines"]


def routines(args, infile):
    """ processing routines """
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

    # -- COBALT
    label = "COBALT"
    modules = [
        "ocean_cobalt_sfc",
        "ocean_cobalt_misc",
        "ocean_cobalt_tracers_year",
        "ocean_cobalt_tracers_int",
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

    # -- BLING
    label = "BLING"
    modules = [
        "ocean_bling",
        "ocean_bling_cmip6_omip_2d",
        "ocean_bling_cmip6_omip_rates_year_z",
        "ocean_bling_cmip6_omip_sfc",
        "ocean_bling_cmip6_omip_tracers_month_z",
        "ocean_bling_cmip6_omip_tracers_year_z",
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

    # -- Ocean
    label = "Ocean"
    fname = modifier + fyear + ".ocean_scalar_annual.nc"
    if fname in members:
        print(fname)
        fdata = extract_from_tar(tar, fname)
        extract_ocean_scalar.mom6(fdata, fyear, "./")
        fdata.close()

    # -- AMOC
    label = "Ocean"
    if args.gridspec is not None:
        gs_tar = tarfile.open(args.gridspec)
        ocean_hgrid = extract_from_tar(gs_tar, "./ocean_hgrid.nc")
        topog = extract_from_tar(gs_tar, "./ocean_topog.nc")
        fname = modifier + fyear + ".ocean_annual_z.nc"
        if fname in members:
            vh_file = extract_from_tar(tar, fname)
            diags.amoc.mom6(vh_file, ocean_hgrid, topog, fyear, "./", label)
            ocean_hgrid.close()
            topog.close()
            vh_file.close()
        gs_tar.close()

    # -- Close out the tarfile handle
    tar.close()
