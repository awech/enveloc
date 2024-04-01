import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
from scipy.interpolate import interp1d
from scipy import ndimage

np.set_printoptions(precision=3)


def channel_list(st):
    channels = []
    for tr in st:
        tmp = ".".join(
            [tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel]
        )
        channels.append(
            (tmp, tr.stats.coordinates.longitude, tr.stats.coordinates.latitude)
        )
    return channels


def reduced_displacement(loc, XC, st_temp):
    if np.isnan(loc.latitude):
        return np.nan

    dr = []

    channels = [ll[0] for ll in loc.channels]
    for tr in st_temp:
        if tr.id in channels:
            data = tr.data * 100
            rms = np.sqrt(np.mean(data**2))
            a = gps2dist_azimuth(
                loc.latitude,
                loc.longitude,
                tr.stats.coordinates.latitude,
                tr.stats.coordinates.longitude,
            )
            a = a[0] / 1000.0
            r = np.sqrt(a**2 + loc.depth**2) * 1000 * 100

            if "rd_freq" in XC.__dict__:
                wavelength = (XC._v0 / float(XC.rd_freq)) * 1000 * 100
                dr.append(rms * np.sqrt(wavelength * r))
            else:
                dr.append(rms * r)

    return np.median(dr)


def location_scatter(loc, lats, lons, deps):
    dist = list()
    for y, x in zip(lats, lons):
        a = gps2dist_azimuth(loc.latitude, loc.longitude, y, x)
        dist.append(a[0] / float(1000))

    scatter_x = np.median(dist)
    scatter_z = np.median(np.abs(deps - loc.depth))

    return scatter_x, scatter_z


def rotate_coords(x, y, x0, y0, angle):
    angle = np.radians(
        angle
    )  # "+" for counter-clockwise rotation , "-" for clockwise rotation
    R = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

    A = np.array([x, y])
    newX, newY = np.matmul(R, A)
    newX = newX + x0
    newY = newY + y0

    return newX, newY


def create_rotated_grid(rotate_params):
    import utm

    """ Example dictionary for rotate_params:

    rotate_params = { 'x'    : np.arange(0,120,8),
                      'y'    : np.arange(-15,18,3),
                      'z'    : np.arange(1,67,2),
                      'lat0' : 18.8956,
                      'lon0' : -155.6078,
                      'az'   : 120 }

    lat0,lon0 = origin
    az        = rotation angle, counterclockwise from East
    x,y,z     = grid to be rotated. +y=North, +x=East, +z=Increasing depth
    """

    # convert lat/lon origin into utm x-y origin
    x0, y0, zone_num, zone_let = utm.from_latlon(
        rotate_params["lat0"], rotate_params["lon0"]
    )

    # create blank grids
    LATS = np.ones(
        (len(rotate_params["y"]), len(rotate_params["x"]), len(rotate_params["z"]))
    )
    LONS = np.ones_like(LATS)
    DEPS = np.ones_like(LATS)

    # fill out depth grid
    for k, z in enumerate(rotate_params["z"]):
        DEPS[:, :, k] = z

    # fill out LATS and LONS grids with rotated coordinates
    for i, y in enumerate(rotate_params["y"]):
        for j, x in enumerate(rotate_params["x"]):
            newX, newY = rotate_coords(1000 * x, 1000 * y, x0, y0, rotate_params["az"])
            LATS[i, j, :], LONS[i, j, :] = utm.to_latlon(
                newX, newY, zone_number=zone_num, zone_letter=zone_let
            )

    return LONS, LATS, DEPS


def CCtoHypocenter(CC, XC):
    Nseis0 = len(CC["st"])
    NW2 = len(
        np.where(CC["W2"] > 0)[0]
    )  # number of seismogram pairs contributing to solution

    ## Calculate misfit for every location on the grid ##
    #####################################################
    dT = np.zeros(
        (len(XC._grid["LON"].flatten()), Nseis0)
    )  # delay times to each seismogram
    for ksta in np.arange(Nseis0):
        dT.T[ksta] = (
            CC["st"][ksta].TT.flatten(order="F") * CC["st"][ksta].stats.sampling_rate
        )

    if XC._minimize == "correlation":
        CCshift = np.arange(Nseis0 * (Nseis0 - 1) / 2) * (2 * (XC._mlag + 1) - 1)
        CCshift = np.tile(CCshift, len(XC._grid["LON"].flatten())).reshape(
            len(XC._grid["LON"].flatten()), len(CCshift)
        )

        deltaT = dT.T[CC["j1"]].T - dT.T[CC["i1"]].T
        deltaT[np.where(abs(deltaT) > XC._mlag)] = np.nan
        deltaT = deltaT + XC._mlag
        IND = deltaT + CCshift

        obs_corr = CC["maxC"].flatten()[CC["ii"]]
        all_corr = CC["C"].T[np.where(CC["indx"][0][:] != CC["indx"][1][:])[0]].T
        f = interp1d(
            np.arange(all_corr.size), all_corr.flatten(order="F"), kind=XC.lookup_type
        )
        pred_corr = f(IND.flatten(order="F"))
        pred_corr = np.reshape(pred_corr, IND.shape, order="F")
        residual = obs_corr - pred_corr

    elif XC._minimize == "time":
        obs_lag = CC["LAG"].flatten()[CC["ii"]]
        pred_lag = dT.T[CC["j1"]].T - dT.T[CC["i1"]].T
        residual = obs_lag - pred_lag

    if XC._normType == 1:
        residual = np.abs(residual)
    elif XC._normType == 2:
        residual = np.square(residual)

    misfit = np.einsum("...j,j", residual, CC["W2"], optimize=False) / float(NW2)
    ind_min = misfit.argmin()

    lat = XC._grid["LAT"].flatten(order="F")[ind_min]
    lon = XC._grid["LON"].flatten(order="F")[ind_min]
    dep = XC._grid["DEP"].flatten(order="F")[ind_min]

    return lat, lon, dep, misfit


def create_regrid_array(old_array, lat_inds, lon_inds, dep_inds, zoom=5, order=3):
    new_array = old_array[lat_inds, :, :]
    new_array = new_array[:, lon_inds, :]
    new_array = new_array[:, :, dep_inds]

    new_array = ndimage.zoom(new_array, zoom, order=order, mode="nearest")
    return new_array


def regridHypocenter(XC, misfit):
    misfit = np.reshape(misfit, np.shape(XC._grid["LON"]), order="F")
    inds_min = np.unravel_index(misfit.argmin(), misfit.shape)

    if XC.rotation:
        lat_slice = np.unique(
            np.clip(
                np.arange(inds_min[0] - 2, inds_min[0] + 3),
                0,
                len(XC.grid_size["y"]) - 1,
            )
        )
        lon_slice = np.unique(
            np.clip(
                np.arange(inds_min[1] - 2, inds_min[1] + 3),
                0,
                len(XC.grid_size["x"]) - 1,
            )
        )
    else:
        lat_slice = np.unique(
            np.clip(
                np.arange(inds_min[0] - 2, inds_min[0] + 3),
                0,
                len(XC.grid_size["lats"]) - 1,
            )
        )
        lon_slice = np.unique(
            np.clip(
                np.arange(inds_min[1] - 2, inds_min[1] + 3),
                0,
                len(XC.grid_size["lons"]) - 1,
            )
        )

    dep_slice = np.unique(
        np.clip(
            np.arange(inds_min[2] - 2, inds_min[2] + 3),
            0,
            len(XC.grid_size["deps"]) - 1,
        )
    )
    lons_new = create_regrid_array(
        XC._grid["LON"], lat_slice, lon_slice, dep_slice, zoom=XC.regrid, order=1
    )
    lats_new = create_regrid_array(
        XC._grid["LAT"], lat_slice, lon_slice, dep_slice, zoom=XC.regrid, order=1
    )
    deps_new = create_regrid_array(
        XC._grid["DEP"], lat_slice, lon_slice, dep_slice, zoom=XC.regrid, order=1
    )

    misfit_new = create_regrid_array(
        misfit, lat_slice, lon_slice, dep_slice, zoom=XC.regrid, order=3
    )

    new_inds_min = np.unravel_index(misfit_new.argmin(), misfit_new.shape)

    lat_new = lats_new[new_inds_min]
    lon_new = lons_new[new_inds_min]
    dep_new = deps_new[new_inds_min]

    return lat_new, lon_new, dep_new
