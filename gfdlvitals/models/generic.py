""" Generic GFDL model processing routine """

import tarfile

from gfdlvitals import averagers
from gfdlvitals import diags
from gfdlvitals.util import extract_ocean_scalar
from gfdlvitals.util.netcdf import extract_from_tar

__all__ = ["routines"]


def routines(args, infile):
    """ processing routine """
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
    print("Calculating AMOC")
    label = "Ocean"
    if args.gridspec is not None:
        gs_tar = tarfile.open(args.gridspec)
        _contents = gs_tar.getnames()

        hgrid_file = [x for x in _contents if "ocean_hgrid" in x]
        if len(hgrid_file) == 1:
            hgrid_file = hgrid_file[0]
        elif len(hgrid_file) > 1:
            print("Multiple ocean_hgrid files found ... skipping AMOC")
        else:
            print("No ocean_hgrid file found ... skipping AMOC")

        topog_file = [x for x in _contents if "ocean_topog" in x]
        if len(topog_file) == 1:
            topog_file = topog_file[0]
        elif len(topog_file) > 1:
            print("Multiple ocean_topog files found ... skipping AMOC")
        else:
            print("No ocean_topog file found ... skipping AMOC")

        try:
            ocean_hgrid = extract_from_tar(gs_tar, hgrid_file)
            topog = extract_from_tar(gs_tar, topog_file)
            fname = modifier + fyear + ".ocean_annual_z.nc"
            if fname in members:
                vh_file = extract_from_tar(tar, fname)
                diags.amoc.mom6(vh_file, ocean_hgrid, topog, fyear, "./", label)
                ocean_hgrid.close()
                topog.close()
                vh_file.close()
            gs_tar.close()
        except NameError:
            print("  * Unable to process")

    # -- Close out the tarfile handle
    tar.close()

    # -- Do performance timing
    print("Extracting FMS timings")
    infile = infile.replace("/history/", "/ascii/")
    infile = infile.replace(".nc.tar", ".ascii_out.tar")
    label = "Timing"
    try:
        diags.fms.timing(infile, fyear, "./", label)
    except NameError:
        print("  * Unable to process")
