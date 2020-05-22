import example_utils as utils
import XC_loc.XC_main as XC_main
from obspy.clients.fdsn import Client

### get data & preprocess into envelopes ###
############################################
client=Client('IRIS')
sta_list=['AV.ANON..EHZ','AV.ANPK..EHZ','AV.AZAC..EHZ','AV.ANNE..EHZ','AV.ANPB..EHZ','AV.ANNW..EHZ']
t1 = '2017-03-10 13:39'
t2 = '2017-03-10 13:44'
env = utils.get_IRIS_data(sta_list,t1,t2,f1=1.0,f2=8.0,lowpass=0.2,client=client)

"""
	Everything above can be done however you want.
	What's important is that you get a stream of 
	uniformly sampled and cut waveforms and each trace
	tr in env has: tr.stats.coordinates.latitude
	               tr.stats.coordinates.longitude
	see utils.add_latlon() to see how I do this
	using Obspy's AttribDict

	In this example 'env' is the envelope 
	of bandpass-filtered data. Envelopes are typically
	downsampled to 1-5 sps and smoothed at 5-10 seconds.
"""

# Example locating single window interactively (make sure to disable any matplotlib backend)
# Click on envelope traces you want to remove from data and relocate
#
# NOTE:		The interaction part has only been lightly tested. Possible bugs with the UI here.
#			All the processing and location stuff (sans-interaction mode) are fairly well
#			tested though.
#
XC   = XC_main.XCOR(env)
loc  = XC.locate()
# now do the same but bootstrap to estimate location stability
XC.bootstrap=20
loc  = XC.locate()


# Example autolocating over sliding windows
XC   = XC_main.XCOR(env)
XC.plot = False # don't plot after each window...could also set as argument in XC_main.XCOR(env,plot=False)
loc  = XC.locate(window_length=25,step=10)


# Example bootstrapping to estimate location stability:
XC   = XC_main.XCOR(env,bootstrap=20,plot=False,output=2)
loc = XC.locate(window_length=25,step=10)


# Example using parallel processing:
XC   = XC_main.XCOR(env,plot=False,num_processors=10,output=3,bootstrap=20)
loc  = XC.locate(window_length=25,step=2)


# Example using your own grid:
import numpy as np
small=0.0005
mygrid={ 'deps': np.arange(0.5,30+small,2),
		 'lons': np.arange(-158.4,-157.9+small,0.04),
		 'lats': np.arange(  56.7,  57.1+small,0.04)}
XC   = XC_main.XCOR(env,bootstrap=20,grid_size=mygrid)
loc  = XC.locate()


# Example using rotated grid:

rotation={'x'    : np.arange(-20,20,1),
		  'y'    : np.arange(-20,20,1),
		  'z'    : np.arange(1,39,2),
		  'lat0' : 56.919,
		  'lon0' : -158.1737,
		  'az'   : 30 }


XC   = XC_main.XCOR(env,rotation=rotation)
loc  = XC.locate()



# plot the locations
import matplotlib.pyplot as plt

# if just 1 location:
plt.figure()
plt.plot(loc.longitude,loc.latitude,'o')

# if many locations:
plt.figure()
plt.plot(loc.get_lons(),loc.get_lats(),'o')

