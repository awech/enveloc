import XC_loc.utils as utils
import XC_loc.XC_main as XC_main
from obspy.clients.fdsn import Client
# from obspy.clients.earthworm import Client

### get data & preprocess into envelopes ###
############################################
client=Client('IRIS')
# client=Client('pubavo1.wr.usgs.gov',16022)
sta_list=['AV.ANON..EHZ','AV.ANPK..EHZ','AV.AZAC..EHZ','AV.ANNE..EHZ','AV.ANPB..EHZ','AV.ANNW..EHZ']
t1 = '2017-03-10 13:39'
t2 = '2017-03-10 13:44'
st, env = utils.get_IRIS_data(sta_list,t1,t2,f1=1.0,f2=8.0,lowpass=0.2,client=client)

"""
	Everything above can be done however you want.
	What's important is that you get a stream of 
	uniformly sampled and cut waveforms and each trace
	tr in st has: tr.stats.coordinates.latitude
	              tr.stats.coordinates.longitude
	see utils.add_latlon() to see how I do this
	using Obspy's AttribDict

	In this example 'st' are non-envelope traces
	instrument-corrected into displacement, and
	'env' are the envelope form of those traces.
	You only need 'env' ('st' is optional'), but
	providing 'st' as well allows for a reduced
	displacement to be returned with each location

	Envelopes are typically downsampled to 1-5 sps
	and smoothed at 5-10 seconds.
"""

# Example locating single window interactively (make sure to disable any matplotlib backend)
# Click on envelope traces you want to remove from data and relocate
XC   = XC_main.XCOR(env)
loc  = XC.locate()
# now do the same but bootstrap to estimate location stability
XC.bootstrap=20
loc  = XC.locate()



# Example autolocating over sliding windows
XC   = XC_main.XCOR(env,raw_traces=st)
loc  = XC.locate(window_length=25,step=10)



# Example bootstrapping to estimate location stability:
XC   = XC_main.XCOR(env,raw_traces=st, bootstrap=20)
loc = XC.locate(window_length=25,step=10)




# Example using parallel processing:
XC   = XC_main.XCOR(env,raw_traces=st,visual=False,num_processors=10,output=3,bootstrap=20)
loc  = XC.locate(window_length=25,step=10)




# Example using your own grid:
import numpy as np
small=0.0005
mygrid={ 'deps': np.arange(0.5,30+small,2),
		 'lons': np.arange(-158.4,-157.9+small,0.04),
		 'lats': np.arange(  56.7,  57.1+small,0.04)}
XC   = XC_main.XCOR(env,bootstrap=20,grid_size=mygrid)
loc  = XC.locate()