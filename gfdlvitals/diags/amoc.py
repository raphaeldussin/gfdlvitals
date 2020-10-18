""" Computes AMOC """

import sys
import numpy as np
from gfdlvitals.util import gmeantools

# sys.path.append('/nbhome/ogrp/warsaw_201710_MOM6_2017.10.19/'+
#    'OM4p25_IAF_baseline/mom6/tools/analysis/')

from . import m6toolbox

__all__ = ["mom6"]


def mom6(vh_file, f_ocean_hgrid, f_topog, fyear, outdir, label):
    """ MOM6-specific implemention of AMOC calculation """

    if "vmo" in vh_file.variables.keys():
        vh = (vh_file.variables["vmo"][0].filled(0)) * 1.0e-9
        zt = vh_file.variables["z_i"][:]
        yq = vh_file.variables["yq"][:]
    else:
        print("amoc.py FATAL: vmo variable not present in ocean_annual_z.nc")
        sys.exit(0)

    # -- Get grid info from gridspec file
    x = f_ocean_hgrid.variables["x"][1::2, 1::2]
    y = f_ocean_hgrid.variables["y"][1::2, 1::2]
    depth = f_topog.variables["depth"][:]
    code = m6toolbox.gen_basin_masks(x, y, depth)

    # -- Define atlantic/arctic mask
    atlmask = np.where(np.logical_or(code == 2, code == 4), 1.0, 0.0)

    # -- Compute psi
    psi = m6toolbox.moc_psi(vh, vmsk=atlmask)
    maxsfn = np.max(
        psi[np.logical_and(zt > 500, zt < 2500)][:, np.greater_equal(yq, 20)]
    )
    print("  AMOC = " + str(maxsfn))

    # -- Write to sqlite
    gmeantools.write_sqlite_data(
        outdir + "/" + fyear + ".globalAve" + label + ".db",
        "amoc_vh",
        fyear[:4],
        maxsfn,
    )
