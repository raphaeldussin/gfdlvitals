"""
A collection of useful functions...
"""
import tarfile
import numpy as np
from scipy.io import netcdf


def section2quadmesh(x, z, q, representation="pcm"):
    """
  Creates the appropriate quadmesh coordinates to plot a scalar q(1:nk,1:ni) at
  horizontal positions x(1:ni+1) and between interfaces at z(nk+1,ni), using
  various representations of the topography.

  Returns X(2*ni+1), Z(nk+1,2*ni+1) and Q(nk,2*ni) to be passed to pcolormesh.

  TBD: Optionally, x can be dimensioned as x(ni) in which case it will be extraplated as if it had
  had dimensions x(ni+1).

  Optional argument:

  representation='pcm' (default) yields a step-wise visualization, appropriate for
           z-coordinate models.
  representation='plm' yields a piecewise-linear visualization more representative
           of general-coordinate (and isopycnal) models.
  representation='linear' is the aesthetically most pleasing but does not
           represent the data conservatively.

  """

    if x.ndim != 1:
        raise Exception("The x argument must be a vector")
    if z.ndim != 2:
        raise Exception("The z argument should be a 2D array")
    if q.ndim != 2:
        raise Exception("The z argument should be a 2D array")
    qnk, qni = q.shape
    znk, zni = z.shape
    xni = x.size
    if zni != qni:
        raise Exception("The last dimension of z and q must be equal in length")
    if znk != qnk + 1:
        raise Exception(
            "The first dimension of z must be 1 longer than that of q. q has %i levels"
            % qnk
        )
    if xni != qni + 1:
        raise Exception("The length of x must 1 longer than the last dimension of q")

    if isinstance(z, np.ma.core.MaskedArray):
        z[z.mask] = 0
    if isinstance(q, np.ma.core.MaskedArray):
        qmin = np.amin(q)
        q[q.mask] = qmin

    periodic_domain = (
        abs((x[-1] - x[0]) - 360.0) < 1e-6
    )  # Detect if horizontal axis is a periodic domain

    if representation == "pcm":
        X = np.zeros((2 * qni))
        X[::2] = x[:-1]
        X[1::2] = x[1:]
        Z = np.zeros((qnk + 1, 2 * qni))
        Z[:, ::2] = z
        Z[:, 1::2] = z
        Q = np.zeros((qnk, 2 * qni - 1))
        Q[:, ::2] = q
        Q[:, 1::2] = (q[:, :-1] + q[:, 1:]) / 2.0
    elif representation == "linear":
        X = np.zeros((2 * qni + 1))
        X[::2] = x
        X[1::2] = (x[0:-1] + x[1:]) / 2.0
        Z = np.zeros((qnk + 1, 2 * qni + 1))
        Z[:, 1::2] = z
        Z[:, 2:-1:2] = (z[:, 0:-1] + z[:, 1:]) / 2.0
        Z[:, 0] = z[:, 0]
        Z[:, -1] = z[:, -1]
        Q = np.zeros((qnk, 2 * qni))
        Q[:, ::2] = q
        Q[:, 1::2] = q
    elif representation == "plm":
        X = np.zeros((2 * qni))
        X[::2] = x[:-1]
        X[1::2] = x[1:]
        # PLM reconstruction for Z
        dz = np.roll(z, -1, axis=1) - z  # Right-sided difference
        if not periodic_domain:
            dz[:, -1] = 0  # Non-periodic boundary
        d2 = (
            np.roll(z, -1, axis=1) - np.roll(z, 1, axis=1)
        ) / 2.0  # Centered difference
        d2 = (dz + np.roll(dz, 1, axis=1)) / 2.0  # Centered difference
        s = np.sign(d2)  # Sign of centered slope
        s[dz * np.roll(dz, 1, axis=1) <= 0] = 0  # Flatten extrema
        dz = np.abs(dz)  # Only need magnitude from here on
        S = s * np.minimum(
            np.abs(d2), np.minimum(dz, np.roll(dz, 1, axis=1))
        )  # PLM slope
        Z = np.zeros((qnk + 1, 2 * qni))
        Z[:, ::2] = z - S / 2.0
        Z[:, 1::2] = z + S / 2.0
        Q = np.zeros((qnk, 2 * qni - 1))
        Q[:, ::2] = q
        Q[:, 1::2] = (q[:, :-1] + q[:, 1:]) / 2.0
    else:
        raise Exception("Unknown representation!")

    return X, Z, Q


def get_z(rg, depth, var_name):
    """Returns 3d interface positions from netcdf group rg,
       based on dimension data for variable var_name"""
    if "e" in rg.variables:  # First try native approach
        if len(rg.variables["e"]) == 3:
            result = rg.variables["e"][:]
        elif len(rg.variables["e"]) == 4:
            result = rg.variables["e"][0]
        else:
            raise ValueError("Invalid length for e.")
        return result
    if var_name not in rg.variables:
        raise Exception('Variable "' + var_name + '" not found in netcdf file')
    if len(rg.variables[var_name].shape) < 3:
        raise Exception('Variable "' + var_name + '" must have 3 or more dimensions')
    vdim = rg.variables[var_name].dimensions[-3]
    if vdim not in rg.variables:
        raise Exception(
            'Variable "' + vdim + '" should be a [CF] dimension variable but is missing'
        )
    if "edges" in rg.variables[vdim].ncattrs():
        zvar = getattr(rg.variables[vdim], "edges")
    elif "zw" in rg.variables:
        zvar = "zw"
    else:
        raise Exception(
            'Cannot figure out vertical coordinate from variable "' + var_name + '"'
        )
    if not len(rg.variables[zvar].shape) == 1:
        raise Exception('Variable "' + zvar + '" was expected to be 1d')
    zw = rg.variables[zvar][:]
    zmod = np.zeros((zw.shape[0], depth.shape[0], depth.shape[1]))
    for k in range(zw.shape[0]):
        zmod[k] = -np.minimum(depth, abs(zw[k]))
    return zmod


def rho_wright_97(S, T, P=0):
    """
  Returns the density of seawater for the given salinity, potential temperature
  and pressure.

  Units: salinity in PSU, potential temperature in degrees Celsius and pressure in Pascals.
  """
    a = [7.057924e-4, 3.480336e-7 - 1.112733e-7]
    b = [5.790749e8, 3.516535e6 - 4.002714e4, 2.084372e2, 5.944068e5, -9.643486e3]
    c = [1.704853e5, 7.904722e2, -7.984422, 5.140652e-2, -2.302158e2, -3.079464]

    al0 = a[0] + a[1] * T + a[2] * S
    p_0 = b[0] + b[4] * S + T * (b[1] + T * (b[2] + b[3] * T) + b[5] * S)
    Lambda = c[0] + c[4] * S + T * (c[1] + T * (c[2] + c[3] * T) + c[5] * S)
    return (P + p_0) / (Lambda + al0 * (P + p_0))


def ice9(i, j, source, xcyclic=True, tripolar=True):
    """
  An iterative (stack based) implementation of "Ice 9".

  The flood fill starts at [j,i] and treats any positive value of "source" as
  passable. Zero and negative values block flooding.

  xcyclic = True allows cyclic behavior in the last index. (default)
  tripolar = True allows a fold across the top-most edge. (default)

  Returns an array of 0's and 1's.
  """
    wet_mask = 0 * source
    (nj, ni) = wet_mask.shape
    stack = set()
    stack.add((j, i))
    while stack:
        (j, i) = stack.pop()
        if wet_mask[j, i] or source[j, i] <= 0:
            continue
        wet_mask[j, i] = 1
        if i > 0:
            stack.add((j, i - 1))
        elif xcyclic:
            stack.add((j, ni - 1))
        if i < ni - 1:
            stack.add((j, i + 1))
        elif xcyclic:
            stack.add((j, 0))
        if j > 0:
            stack.add((j - 1, i))
        if j < nj - 1:
            stack.add((j + 1, i))
        elif tripolar:
            stack.add((j, ni - 1 - i))  # Tri-polar fold
    return wet_mask


def ice9_wrapper(x, y, depth, xy0):
    """ Calls Ice-9 function"""
    ji = nearest_ji(x, y, xy0)
    return ice9(ji[1], ji[0], depth)


def mask_from_depth(depth, z_cell_top):
    """
  Generates a "wet mask" for a z-coordinate model based on relative location of
  the ocean bottom to the upper interface of the cell.

  depth (2d) is positiveo
  z_cell_top (scalar) is a negative position of the upper interface of the cell..
  """
    wet = 0 * depth
    wet[depth > -z_cell_top] = 1
    return wet


def moc_psi(vh, vmsk=None):
    """Sums 'vh' zonally and cumulatively in the vertical to yield
       an overturning stream function, psi(y,z)."""
    shape = list(vh.shape)
    shape[-3] += 1
    psi = np.zeros(shape[:-1])
    if len(shape) == 3:
        for k in range(shape[-3] - 1, 0, -1):
            if vmsk is None:
                psi[k - 1, :] = psi[k, :] - vh[k - 1].sum(axis=-1)
            else:
                psi[k - 1, :] = psi[k, :] - (vmsk * vh[k - 1]).sum(axis=-1)
    else:
        for n in range(shape[0]):
            for k in range(shape[-3] - 1, 0, -1):
                if vmsk is None:
                    psi[n, k - 1, :] = psi[n, k, :] - vh[n, k - 1].sum(axis=-1)
                else:
                    psi[n, k - 1, :] = psi[n, k, :] - (vmsk * vh[n, k - 1]).sum(axis=-1)
    return psi


def moc_maskedarray(vh, mask=None):
    """ Coumputes MOC from a masked-array"""
    if mask is not None:
        _mask = np.ma.masked_where(np.not_equal(mask, 1.0), mask)
    else:
        _mask = 1.0
    _vh = vh * _mask
    _vh_btm = np.ma.expand_dims(_vh[:, -1, :, :] * 0.0, axis=1)
    _vh = np.ma.concatenate((_vh, _vh_btm), axis=1)
    _vh = np.ma.sum(_vh, axis=-1) * -1.0
    _vh = _vh[:, ::-1]  # flip z-axis so running sum is from ocean floor to surface
    _vh = np.ma.cumsum(_vh, axis=1)
    _vh = _vh[:, ::-1]  # flip z-axis back to original order
    return _vh


def nearest_ji(x, y, xy0):
    """
  Find (j,i) of cell with center nearest to (x0,y0).
  """
    x0, y0 = xy0
    return np.unravel_index(((x - x0) ** 2 + (y - y0) ** 2).argmin(), x.shape)


def read_nc_from_tar(tar, file, var):
    """Reads NetCDF file that is inside a tar file"""
    tfile = tarfile.open(tar, "r")
    member = [m for m in tfile.getmembers() if file in m.name][0]
    nc = netcdf.netcdf_file(tfile.extractfile(member), "r")
    return nc.variables[var]


def south_of(x, y, xy0, xy1):
    """
  Returns 1 for point south/east of the line that passes through xy0-xy1, 0 otherwise.
  """
    x0 = xy0[0]
    y0 = xy0[1]
    x1 = xy1[0]
    y1 = xy1[1]
    dx = x1 - x0
    dy = y1 - y0
    Y = (x - x0) * dy - (y - y0) * dx
    Y[Y >= 0] = 1
    Y[Y <= 0] = 0
    return Y


def gen_basin_masks(x, y, depth, verbose=False):
    """ Ice-9 based arroproach for defining basin masks """

    def vprint(_str, verbose=verbose):
        if verbose is True:
            print(_str)

    vprint("Generating global wet mask ...")
    wet = ice9_wrapper(
        x, y, depth, (0, -35)
    )  # All ocean points seeded from South Atlantic
    vprint("done.")

    code = 0 * wet

    vprint("Finding Cape of Good Hope ...")
    tmp = 1 - wet
    tmp[x < -30] = 0
    tmp = ice9_wrapper(x, y, tmp, (20, -30.0))
    y_cgh = (tmp * y).min()
    vprint("done.", y_cgh)

    vprint("Finding Melbourne ...")
    tmp = 1 - wet
    tmp[x > -180] = 0
    tmp = ice9_wrapper(x, y, tmp, (-220, -25.0))
    y_mel = (tmp * y).min()
    vprint("done.", y_mel)

    vprint("Processing Persian Gulf ...")
    tmp = wet * (1 - south_of(x, y, (55.0, 23.0), (56.5, 27.0)))
    tmp = ice9_wrapper(x, y, tmp, (53.0, 25.0))
    code[tmp > 0] = 11
    wet = wet - tmp  # Removed named points

    vprint("Processing Red Sea ...")
    tmp = wet * (1 - south_of(x, y, (40.0, 11.0), (45.0, 13.0)))
    tmp = ice9_wrapper(x, y, tmp, (40.0, 18.0))
    code[tmp > 0] = 10
    wet = wet - tmp  # Removed named points

    vprint("Processing Black Sea ...")
    tmp = wet * (1 - south_of(x, y, (26.0, 42.0), (32.0, 40.0)))
    tmp = ice9_wrapper(x, y, tmp, (32.0, 43.0))
    code[tmp > 0] = 7
    wet = wet - tmp  # Removed named points

    vprint("Processing Mediterranean ...")
    tmp = wet * (south_of(x, y, (-5.7, 35.5), (-5.7, 36.5)))
    tmp = ice9_wrapper(x, y, tmp, (4.0, 38.0))
    code[tmp > 0] = 6
    wet = wet - tmp  # Removed named points

    vprint("Processing Baltic ...")
    tmp = wet * (south_of(x, y, (8.6, 56.0), (8.6, 60.0)))
    tmp = ice9_wrapper(x, y, tmp, (10.0, 58.0))
    code[tmp > 0] = 9
    wet = wet - tmp  # Removed named points

    vprint("Processing Hudson Bay ...")
    tmp = wet * (
        (
            1
            - (1 - south_of(x, y, (-95.0, 66.0), (-83.5, 67.5)))
            * (1 - south_of(x, y, (-83.5, 67.5), (-84.0, 71.0)))
        )
        * (1 - south_of(x, y, (-70.0, 58.0), (-70.0, 65.0)))
    )
    tmp = ice9_wrapper(x, y, tmp, (-85.0, 60.0))
    code[tmp > 0] = 8
    wet = wet - tmp  # Removed named points

    vprint("Processing Arctic ...")
    tmp = wet * (
        (1 - south_of(x, y, (-171.0, 66.0), (-166.0, 65.5)))
        * (1 - south_of(x, y, (-64.0, 66.4), (-50.0, 68.5)))  # Lab Sea
        + south_of(x, y, (-50.0, 0.0), (-50.0, 90.0))
        * (1 - south_of(x, y, (0.0, 65.5), (360.0, 65.5)))  # Denmark Strait
        + south_of(x, y, (-18.0, 0.0), (-18.0, 65.0))
        * (1 - south_of(x, y, (0.0, 64.9), (360.0, 64.9)))  # Iceland-Sweden
        + south_of(x, y, (20.0, 0.0), (20.0, 90.0))  # Barents Sea
        + (1 - south_of(x, y, (-280.0, 55.0), (-200.0, 65.0)))
    )
    tmp = ice9_wrapper(x, y, tmp, (0.0, 85.0))
    code[tmp > 0] = 4
    wet = wet - tmp  # Removed named points

    vprint("Processing Pacific ...")
    tmp = wet * (
        (1 - south_of(x, y, (0.0, y_mel), (360.0, y_mel)))
        - south_of(x, y, (-257, 1), (-257, 0)) * south_of(x, y, (0, 3), (1, 3))
        - south_of(x, y, (-254.25, 1), (-254.25, 0)) * south_of(x, y, (0, -5), (1, -5))
        - south_of(x, y, (-243.7, 1), (-243.7, 0))
        * south_of(x, y, (0, -8.4), (1, -8.4))
        - south_of(x, y, (-234.5, 1), (-234.5, 0))
        * south_of(x, y, (0, -8.9), (1, -8.9))
    )
    tmp = ice9_wrapper(x, y, tmp, (-150.0, 0.0))
    code[tmp > 0] = 3
    wet = wet - tmp  # Removed named points

    vprint("Processing Atlantic ...")
    tmp = wet * (1 - south_of(x, y, (0.0, y_cgh), (360.0, y_cgh)))
    tmp = ice9_wrapper(x, y, tmp, (-20.0, 0.0))
    code[tmp > 0] = 2
    wet = wet - tmp  # Removed named points

    vprint("Processing Indian ...")
    tmp = wet * (1 - south_of(x, y, (0.0, y_cgh), (360.0, y_cgh)))
    tmp = ice9_wrapper(x, y, tmp, (55.0, 0.0))
    code[tmp > 0] = 5
    wet = wet - tmp  # Removed named points

    vprint("Processing Southern Ocean ...")
    tmp = ice9_wrapper(x, y, wet, (0.0, -55.0))
    code[tmp > 0] = 1
    wet = wet - tmp  # Removed named points

    vprint("Remapping Persian Gulf points to the Indian Ocean for OMIP/CMIP6 ...")
    code[code == 11] = 5

    code[wet > 0] = -9
    inds = np.unravel_index(wet.argmax(), x.shape)
    j = inds[-2]
    i = inds[-1]
    if j:
        vprint("There are leftover points unassigned to a basin code")
        while j:
            # print(x[j,i],y[j,i],[j,i])
            wet[j, i] = 0
            inds = np.unravel_index(wet.argmax(), x.shape)
            j = inds[-2]
            i = inds[-1]
    else:
        vprint("All points assigned a basin code")

    vprint(
        """
Basin codes:
-----------------------------------------------------------
  (1) Southern Ocean      (6) Mediterranean Sea
  (2) Atlantic Ocean      (7) Black Sea
  (3) Pacific Ocean       (8) Hudson Bay
  (4) Arctic Ocean        (9) Baltic Sea
  (5) Indian Ocean       (10) Red Sea

    """
    )

    return code


# Tests
if __name__ == "__main__":

    import matplotlib.pyplot as plt

    # Test data
    _x = np.arange(5)
    _z = (
        np.array(
            [
                [0, 0.2, 0.3, -0.1],
                [1, 1.5, 0.7, 0.4],
                [2, 2, 1.5, 2],
                [3, 2.3, 1.5, 2.1],
            ]
        )
        * -1
    )
    _q = np.matlib.rand(3, 4)
    print("x=", _x)
    print("z=", _z)
    print("q=", _q)

    _X, _Z, _Q = section2quadmesh(_x, _z, _q)
    print("X=", _X)
    print("Z=", _Z)
    print("Q=", _Q)
    plt.subplot(3, 1, 1)
    plt.pcolormesh(_X, _Z, _Q)

    _X, _Z, _Q = section2quadmesh(_x, _z, _q, representation="linear")
    print("X=", _X)
    print("Z=", _Z)
    print("Q=", _Q)
    plt.subplot(3, 1, 2)
    plt.pcolormesh(_X, _Z, _Q)

    _X, _Z, _Q = section2quadmesh(_x, _z, _q, representation="plm")
    print("X=", _X)
    print("Z=", _Z)
    print("Q=", _Q)
    plt.subplot(3, 1, 3)
    plt.pcolormesh(_X, _Z, _Q)

    plt.show()
