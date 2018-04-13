import numpy as np
np.set_printoptions(precision=3)
from obspy.geodetics.base import gps2dist_azimuth
import itertools
from scipy.interpolate import interp1d, griddata


def channel_list(st):
    channels=[]
    for tr in st:
        tmp='.'.join([tr.stats.network,tr.stats.station,tr.stats.location,tr.stats.channel])
        channels.append((tmp,tr.stats.coordinates.longitude,tr.stats.coordinates.latitude))
    return channels


def reduced_displacement(loc, XC, st_temp):
 
    if np.isnan(loc.latitude):
        return np.nan

    dr = []

    channels=[l[0] for l in loc.channels]
    for tr in st_temp:
        if tr.id in channels:
            data=tr.data*100
            rms=np.sqrt(np.mean(data**2))
            a=gps2dist_azimuth(loc.latitude,loc.longitude,tr.stats.coordinates.latitude,tr.stats.coordinates.longitude)
            a=a[0]/1000.
            r=np.sqrt(a**2+loc.depth**2)*1000*100
            
            if 'rd_freq' in XC.__dict__:
                wavelength=(XC._v0/float(XC.rd_freq))*1000*100
                dr.append(rms*np.sqrt(wavelength*r))
            else:
                dr.append(rms*r)
            
    return np.median(dr)


def location_scatter(loc,lats,lons,deps):

    dist=list()
    for y,x in zip(lats,lons):
        a=gps2dist_azimuth(loc.latitude,loc.longitude,y,x)
        dist.append(a[0]/float(1000))

    scatter_x=np.median(dist)
    scatter_z=np.median(np.abs(deps-loc.depth))

    return scatter_x, scatter_z


def CCtoHypocenter(CC,XC):
    Nseis0=len(CC['st'])
    NW2=len(np.where(CC['W2'] > 0)[0])                      # number of seismogram pairs contributing to solution
    CMAXVEC=CC['maxC'].flatten()[CC['ii']]
    CMAXVEC=np.tile(CMAXVEC,len(CC['tC'])).reshape((len(CC['tC']),len(CMAXVEC)))
    CC1=CMAXVEC-CC['C'].T[np.where(CC['indx'][0][:]!=CC['indx'][1][:])[0]].T
    # CCshift = np.arange(Nseis0*(Nseis0-1)/2)*(2*(XC._mlag+1)-1)-1 # Oct-30-2017 change
    CCshift = np.arange(Nseis0*(Nseis0-1)/2)*(2*(XC._mlag+1)-1)     # Oct-30-2017 change

    ## Calculate misfit for every location on the grid ##
    #####################################################
    dT = np.zeros( (len(XC._grid['LON'].flatten()),Nseis0) )   # delay times to each seismogram
    for ksta in np.arange(Nseis0):
        dT.T[ksta]=CC['st'][ksta].TT.flatten(order='F')*CC['st'][ksta].stats.sampling_rate
    # Oct-30-2017 change:
    # deltaT = dT.T[CC['j1']].T - dT.T[CC['i1']].T + XC._mlag +1
    deltaT = dT.T[CC['j1']].T - dT.T[CC['i1']].T
    deltaT[np.where(abs(deltaT)>XC._mlag)]=np.nan
    deltaT = deltaT+XC._mlag
    # End change
    CCshift = np.tile(CCshift,len(XC._grid['LON'].flatten())).reshape(len(XC._grid['LON'].flatten()),len(CCshift))
    IND=deltaT+CCshift

    f = interp1d(np.arange(CC1.size),CC1.flatten(order='F'),kind=XC.lookup_type)
    CC0 = f(IND.flatten(order='F'))
    CC0 = np.reshape(CC0,IND.shape,order='F')
    
    if XC._normType==1:
        misfit=np.einsum('...j,j',CC0, CC['W2'],optimize=False)/float(NW2)
        ind_min=misfit.argmin()
        # testerror=sum(CC0(ind_min,:))/(NW2-3);
        # chisquare=CC0*((1/testerror)*ones(length(W2new),1));
        # chisquare=reshape(chisquare,size(LON1));
    elif XC._normType==2:
        misfit=np.einsum('...j,j',np.square(CC0), CC['W2'],optimize=False)/float(NW2)
        ind_min=misfit.argmin()
        # testerror=sum(CC0(ind_min,:).^2)/(NW2-3);
        # chisquare=(CC0.^2)*((1/testerror)*ones(length(W2new),1));
        # chisquare=reshape(chisquare,size(LON1));

    lat = XC._grid['LAT'].flatten(order='F')[ind_min]
    lon = XC._grid['LON'].flatten(order='F')[ind_min]
    dep = XC._grid['DEP'].flatten(order='F')[ind_min]

    return lat, lon, dep, misfit


def regridHypocenter(CC,XC,misfit):

    n_newgrid = 20
    ind_min=misfit.argmin()
    lat = XC._grid['LAT'].flatten(order='F')[ind_min]
    lon = XC._grid['LON'].flatten(order='F')[ind_min]
    dep = XC._grid['DEP'].flatten(order='F')[ind_min]


    d_LAT=np.diff(XC.grid_size['lats']).mean()
    d_LON=np.diff(XC.grid_size['lons']).mean()

    new_lats = np.linspace(lat-2*d_LAT,lat+2*d_LAT,n_newgrid)
    new_lons = np.linspace(lon-2*d_LAT,lon+2*d_LAT,n_newgrid)

    new_lats = new_lats[np.where(new_lats>XC.grid_size['lats'][ 0])[0]]
    new_lats = new_lats[np.where(new_lats<XC.grid_size['lats'][-1])[0]]
    new_lons = new_lons[np.where(new_lons>XC.grid_size['lons'][ 0])[0]]
    new_lons = new_lons[np.where(new_lons<XC.grid_size['lons'][-1])[0]]

    misfit = np.reshape(misfit,np.shape(XC._grid['LON']),order='F')
    ii=np.unravel_index(misfit.argmin(),np.shape(XC._grid['LON']),order='C')
    misfit_slice=misfit[:,:,ii[2]]
    LONS,LATS = np.meshgrid(new_lons,new_lats)

    points=np.array([XC._grid['LON'][:,:,ii[2]].flatten(order='F'), XC._grid['LAT'][:,:,ii[2]].flatten(order='F')]).T
    misfit_slice=misfit_slice.flatten(order='F')
    grid_new = griddata(points, misfit_slice, (LONS, LATS), method='cubic')
    ii=np.unravel_index(grid_new.argmin(),np.shape(LONS),order='C')

    lat_new = LATS[ii]
    lon_new = LONS[ii]
    dep_new = dep

    return lat_new, lon_new, dep_new