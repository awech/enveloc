import numpy as np
from scipy.signal import correlate
from itertools import combinations, combinations_with_replacement
from obspy.geodetics.base import gps2dist_azimuth
from obspy import Stream


def setup_indices(st):
    """
        ii, i1, j1 are indices:  ii points from the NxN matrix form to a vector that contains
        only the points below the diagonal while (i1,j1) points from the vector form back to the
        i,j components of the matrix.
    """
    Nseis0 = len(st)
    k = 0
    ii = np.zeros( ( int(Nseis0*(Nseis0-1)/2) , ), dtype=int )
    i1 = ii.T.copy()
    j1 = ii.T.copy()
    for j in np.arange(int(Nseis0-1)):
        for i in np.arange(int(j+1),int(Nseis0)):
            ii[k] = i + j*Nseis0
            i1[k] = i
            j1[k] = j
            k = k + 1
    return ii, i1, j1

    
def drop_stas(st,key):
    keep = np.setdiff1d(np.arange(len(st)),key)
    ST = Stream()
    for ind in keep:
        ST+=st[ind]
    return ST


def station_distances(st):

    Delta = np.zeros( (len(st),len(st)) )
    for comb in combinations(np.arange(len(st)),2):
        xy1 = st[comb[0]].stats.coordinates
        xy2 = st[comb[1]].stats.coordinates
        Delta[comb[0],comb[1]] = gps2dist_azimuth(xy1.latitude,xy1.longitude,xy2.latitude,xy2.longitude)[0]/1000.
    return Delta


def station_xcorr(st,mlag):
    tC = np.linspace(-mlag,mlag,2*mlag+1)/st[0].stats.sampling_rate
    C = np.zeros( (len(tC), int(len(st)+len(st)*(len(st)-1)/2)) )
    LAG = np.zeros((len(st),len(st)))
    maxC = np.zeros((len(st),len(st)))

    tmp_i = list()
    tmp_j = list()
    tmp_I = np.zeros( (len(st),len(st)) )
        
    ref=correlate(st[0].data,st[0].data,mode='full',method='fft')
    mid=ref.argmax()

    for i, comb in enumerate(combinations_with_replacement(np.arange(len(st)),2)):
        tmp_i.append(comb[0])
        tmp_j.append(comb[1])
        tmp_I[comb[0],comb[1]]=1

        scale=np.linalg.norm(st[comb[0]].data)*np.linalg.norm(st[comb[1]].data)
        tmp = correlate(st[comb[0]].data,st[comb[1]].data,mode='full',method='fft')/float(scale)
        # tmp = tmp[mid-mlag-1:mid+mlag]
        tmp = tmp[mid-mlag:mid+mlag+1]     # Oct-30-2017 change when switching from mlag=int(2*dTmax) to mlag=int(dTmax)
        C[:,i] = tmp
        maxC[comb[0],comb[1]]=tmp.max()
        maxC[comb[1],comb[0]]=tmp.max()
        LAG[comb[0],comb[1]]=tC[tmp.argmax()]
        LAG[comb[1],comb[0]]=-tC[tmp.argmax()]
           
    indx=np.array([tmp_i,tmp_j,np.where(tmp_I.flatten())[0]])

    return C,tC,indx,LAG,maxC


def trim_CC(drop_key,Nprev,Nseis0,C,maxC,LAG,ii,indx):

    key = np.setdiff1d(np.arange(Nprev),drop_key)
    keep_indx=np.zeros( ( int(len(key)*(len(key)-1)/2) + len(key), ),dtype=int)

    for i, comb in enumerate(combinations_with_replacement(key,2)):
        keep_indx[i] = comb[0]*Nprev + comb[1]

    INDX = np.zeros((3,len(keep_indx)))
    maxc = np.zeros((len(key),len(key)))
    lag  = np.zeros((len(key),len(key)))
    for i, comb in enumerate(combinations_with_replacement(np.arange(len(key)),2)):
        INDX.T[i] = [comb[0],comb[1],comb[0]*Nseis0 +comb[1]]
        maxc[comb[0]][comb[1]] = maxC[key[comb[0]]][key[comb[1]]]
        maxc[comb[1]][comb[0]] = maxC[key[comb[1]]][key[comb[0]]]
        lag[comb[0]][comb[1]]  =  LAG[key[comb[0]]][key[comb[1]]]
        lag[comb[1]][comb[0]]  =  LAG[key[comb[1]]][key[comb[0]]]

    ind_indx, = np.where(np.in1d(indx[2],keep_indx))
    C=C.T[ind_indx].T
    indx  = INDX
    maxC  = maxc
    Chat  = maxC.flatten()[ii]  # vector form of the maximum correlation coeffient for each pair
    LAG = lag

    return C, indx, LAG, maxC, Chat


def station_weights(Chat,lag,Delta,Nseis0,i1,j1,XC):
    W2 = 1/(np.polyval(XC._p,Chat))
    k,  = np.where(Chat<XC.Cmin)
    W2[k] = 0
    k,  = np.where(Chat>XC._Cmax)
    W2[k] = 0
    k,  = np.where(Delta<XC._dx_min)
    W2[k] = 0
    k,  = np.where(np.abs(lag)>(Delta/XC._v0+XC._dt))
    W2[k] = 0
    k,  = np.where(Delta/XC._v0 > XC.dTmax_s)
    if len(k)>0 and XC.output>2:
        print('{:.0f} station pairs with maximum theoretical travel times > allowed lag times'.format(len(k)))
    W2[k] = 0

    staW=np.zeros(Nseis0, )
    for k in np.arange(Nseis0):
        k1, = np.where(i1==k)
        k2, = np.where(j1==k)
        staW[k] = np.sum(W2[k1])/2.+np.sum(W2[k2])/2.

    return staW, W2


def initial_correlation(XC,st):
    drop_key = np.array([])
    count = 1
    while (drop_key>-1).any() or count==1:
        Nseis0 = len(st)
        err = False
        ii, i1, j1 = setup_indices(st)
        Delta = station_distances(st)
        Delta = Delta.flatten()[ii]
        if not (drop_key>-1).any():
            C,tC,indx,LAG,maxC=station_xcorr(st,XC._mlag)
            Chat  = maxC.flatten()[ii]  # vector form of the maximum correlation coeffient for each pair
            count = 0
        else:
            C, indx, LAG, maxC, Chat = trim_CC(drop_key,Nprev,Nseis0,C,maxC,LAG,ii,indx)

        lag   = LAG.flatten()[ii]       # vector form of the best lag time for each station pair
        staW, W2 = station_weights(Chat,lag,Delta,Nseis0,i1,j1,XC)

        drop_key = np.where(staW==0)[0]
        Nprev = Nseis0
        st=drop_stas(st,drop_key)
        
        if W2.max()==0 or np.unique([tr.stats.station for tr in st]).size < XC.sta_min:
            err = True
            drop_key = np.array([])

    if err:
        CC=dict({'err':err,'st':st})
    else:
        CC=dict({'st':st,'C':C,'LAG':LAG,'maxC':maxC,
                'indx':indx,'Delta':Delta,'W2':W2,'staW':staW,
                'err':err,'ii':ii,'i1':i1,'j1':j1,'tC':tC})
    return CC


def remove_stations(CC,XC,drop_key):
    # The idea of this code is to remove selected stations and reindex behind
    # the scenes

    Nprev=len(CC['st'])
    st=drop_stas(CC['st'],drop_key)

    if not (drop_key>-1).any():
        print('No stations to remove')
        return CC

    C    = CC['C']
    maxC = CC['maxC']
    LAG  = CC['LAG']
    indx = CC['indx']
    while (drop_key>-1).any():
        err = False
        Nseis0=len(st)
        ii, i1, j1 = setup_indices(st)
        Delta = station_distances(st)
        Delta = Delta.flatten()[ii]
        C, indx, LAG, maxC, Chat = trim_CC(drop_key,Nprev,Nseis0,C,maxC,LAG,ii,indx)

        lag   = LAG.flatten()[ii]       # vector form of the best lag time for each station pair
        staW, W2 = station_weights(Chat,lag,Delta,Nseis0,i1,j1,XC)

        drop_key = np.where(staW == 0)[0]
        Nprev = Nseis0
        st=drop_stas(st,drop_key)

        if W2.max()==0 or np.unique([tr.stats.station for tr in st]).size < XC.sta_min:
            err = True
            drop_key = np.array([])


    # if err:
    #     CCrm=dict({'err':err,'st':st})
    # else:
    #     CCrm=dict({'st':st,'C':C,'LAG':LAG,'maxC':maxC,
    #             'indx':indx,'Delta':Delta,'W2':W2,'staW':staW,
    #             'err':err,'ii':ii,'i1':i1,'j1':j1,'tC':CC['tC']})

    CCrm=dict({'st':st,'C':C,'LAG':LAG,'maxC':maxC,
        'indx':indx,'Delta':Delta,'W2':W2,'staW':staW,
        'err':err,'ii':ii,'i1':i1,'j1':j1,'tC':CC['tC']})

    return CCrm


def bootstrap_CC(CC,XC):

    Nprev  = len(CC['st'])
    st   = CC['st']

    drop_key=np.array([])
    bstrapFLAG=True    

    while (drop_key>-1).any() | bstrapFLAG:
        err = False
        Nseis0=len(st)
        if (drop_key>-1).any():
            ii, i1, j1 = setup_indices(st)
            Delta = station_distances(st)
            Delta = Delta.flatten()[ii]
            C, indx, LAG, maxC, Chat = trim_CC(drop_key,Nprev,Nseis0,C,maxC,LAG,ii,indx)
        else:
            ii    = CC['ii']
            i1    = CC['i1']
            j1    = CC['j1']
            Delta = CC['Delta']
            C     = CC['C']
            LAG   = CC['LAG']
            indx  = CC['indx']
            n_data= int(np.round(XC.bootstrap_prct*len(CC['ii'])))
            a = np.floor(len(CC['ii'])*np.random.rand(n_data,)).astype('int')
            bstrap_i = ii[a]
            maxC = CC['maxC'].flatten()
            maxC[bstrap_i] = 0
            maxC = np.reshape(maxC,np.shape(LAG))
            Chat = maxC.flatten()[ii]


        lag   = LAG.flatten()[ii]       # vector form of the best lag time for each station pair
        staW, W2 = station_weights(Chat,lag,Delta,Nseis0,i1,j1,XC)

        drop_key = np.where(staW == 0)[0]
        Nprev = Nseis0
        st=drop_stas(st,drop_key)

        if W2.max()==0 or np.unique([tr.stats.station for tr in st]).size < XC.sta_min:
            err = True
            drop_key = np.array([])

        bstrapFLAG=False

    if err:
        CCnew=dict({'err':err,'st':st})
    else:
        CCnew=dict({'st':st,'C':C,'LAG':LAG,'maxC':maxC,
                'indx':indx,'Delta':Delta,'W2':W2,'staW':staW,
                'err':err,'ii':ii,'i1':i1,'j1':j1,'tC':CC['tC']})

    # CCnew=dict({'st':st,'C':C,'LAG':LAG,'maxC':maxC,
    #     'indx':indx,'Delta':Delta,'W2':W2,'staW':staW,
    #     'err':err,'ii':ii,'i1':i1,'j1':j1,'tC':CC['tC']})

    return CCnew