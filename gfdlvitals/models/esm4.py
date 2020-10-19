""" GFDL-ESM4 processor """

import os
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

    # -- Atmosphere
    label = "Atmos"
    modules = ["atmos_month", "atmos_co2_month"]
    # -- open gridspec tiles
    gs_tiles = []
    for tile in range(1, 7):
        gs_tiles.append(
            extract_from_tar(
                tar, modifier + fyear + ".grid_spec.tile" + str(tile) + ".nc"
            )
        )
    # -- data tiles
    for module in modules:
        fname = modifier + fyear + "." + module + ".tile1.nc"
        if fname in members:
            data_tiles = []
            for tile in range(1, 7):
                data_tiles.append(
                    extract_from_tar(
                        tar,
                        modifier + fyear + "." + module + ".tile" + str(tile) + ".nc",
                    )
                )
            print(fname)
            averagers.cubesphere.average(gs_tiles, data_tiles, fyear, "./", label)
            for ds in data_tiles:
                ds.close()

    # -- Aerosols
    label = "AtmosAer"
    modules = ["atmos_month_aer"]
    # -- open gridspec tiles
    gs_tiles = []
    for tile in range(1, 7):
        gs_tiles.append(
            extract_from_tar(
                tar, modifier + fyear + ".grid_spec.tile" + str(tile) + ".nc"
            )
        )
    # -- data tiles
    for module in modules:
        fname = modifier + fyear + "." + module + ".tile1.nc"
        if fname in members:
            data_tiles = []
            for tile in range(1, 7):
                data_tiles.append(
                    extract_from_tar(
                        tar,
                        modifier + fyear + "." + module + ".tile" + str(tile) + ".nc",
                    )
                )
            print(fname)
            averagers.cubesphere.average(gs_tiles, data_tiles, fyear, "./", label)
            for ds in data_tiles:
                ds.close()

    # -- Aerosols (CMIP)
    label = "AeroCMIP"
    modules = ["aerosol_month_cmip"]
    # -- open gridspec tiles
    gs_tiles = []
    for tile in range(1, 7):
        gs_tiles.append(
            extract_from_tar(
                tar, modifier + fyear + ".grid_spec.tile" + str(tile) + ".nc"
            )
        )
    # -- data tiles
    for module in modules:
        fname = modifier + fyear + "." + module + ".tile1.nc"
        if fname in members:
            data_tiles = []
            for tile in range(1, 7):
                data_tiles.append(
                    extract_from_tar(
                        tar,
                        modifier + fyear + "." + module + ".tile" + str(tile) + ".nc",
                    )
                )
            print(fname)
            averagers.cubesphere.average(gs_tiles, data_tiles, fyear, "./", label)
            for ds in data_tiles:
                ds.close()

    # -- Land
    label = "Land"
    modules = ["land_month"]
    # -- open gridspec tiles
    gs_tiles = []
    for tile in range(1, 7):
        gs_tiles.append(
            extract_from_tar(
                tar, modifier + fyear + ".land_static.tile" + str(tile) + ".nc"
            )
        )
    # -- data tiles
    for module in modules:
        fname = modifier + fyear + "." + module + ".tile1.nc"
        if fname in members:
            data_tiles = []
            for tile in range(1, 7):
                data_tiles.append(
                    extract_from_tar(
                        tar,
                        modifier + fyear + "." + module + ".tile" + str(tile) + ".nc",
                    )
                )
            print(fname)
            averagers.land_lm4.average(gs_tiles, data_tiles, fyear, "./", label)
            for ds in data_tiles:
                ds.close()

    # -- Ice
    label = "Ice"
    modules = ["ice_month"]
    fgs = None
    if modifier + fyear + ".ice_static.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ice_static.nc")
    elif modifier + fyear + ".ice_month.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ice_month.nc")
    if fgs is not None:
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
    fgs = None
    if modifier + fyear + ".ocean_static.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ocean_static.nc")
    elif modifier + fyear + ".ocean_month.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ocean_month.nc")
    if fgs is not None:
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
    fgs = None
    if modifier + fyear + ".ocean_static.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ocean_static.nc")
    elif modifier + fyear + ".ocean_month.nc" in members:
        fgs = extract_from_tar(tar, modifier + fyear + ".ocean_month.nc")
    if fgs is not None:
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
        ocean_hgrid = extract_from_tar(gs_tar, "ocean_hgrid.nc")
        topog = extract_from_tar(gs_tar, "ocean_topog.nc")
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

    # -- Do performance timing
    infile = infile.replace("/history/", "/ascii/")
    infile = infile.replace(".nc.tar", ".ascii_out.tar")
    label = "Timing"
    if os.path.exists(infile):
        diags.fms.timing(infile, fyear, "./", label)
