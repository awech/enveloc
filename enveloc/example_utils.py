from obspy import UTCDateTime, Stream
from obspy.signal.filter import envelope
from obspy.core.util import AttribDict
from obspy.clients.fdsn import Client
import numpy as np
from enveloc.core import XCOR

client = Client('IRIS')

dt = 10

def add_latlon(st,client=None):

	if not client:
		client = Client('IRIS')
	for tr in st:
		print(tr)
		inventory = client.get_stations(network=tr.stats.network,
			                            station=tr.stats.station,
			                            location=tr.stats.location,
			                            channel=tr.stats.channel,
										starttime=tr.stats.starttime,
										endtime=tr.stats.endtime)
		tr.stats.coordinates=AttribDict({
	                        'latitude':  inventory[0][0].latitude,
	                        'longitude': inventory[0][0].longitude,
	                        'elevation': inventory[0][0].elevation})
	return st


def get_data(sta_list,t1,t2,fill_value=0,client=None):
	
	st=Stream()

	if not client:
		client = Client('IRIS')
	for channel in sta_list:
		net, sta, loc, chan = channel.split('.')
		try:
			tr=client.get_waveforms(net, sta, loc, chan, UTCDateTime(t1)-dt,UTCDateTime(t2)+dt,attach_response=True)
			if len(tr)>1:
				if fill_value==0 or fill_value==None:
					tr.detrend('demean')
					tr.taper(max_percentage=None,max_length=1)
				for sub_trace in tr:
					# deal with rare error when sub-traces have different dtypes
					if sub_trace.data.dtype.name != 'int32':
						sub_trace.data=sub_trace.data.astype('int32')
					if sub_trace.data.dtype!=dtype('int32'):
						sub_trace.data=sub_trace.data.astype('int32')
					# deal with rare error when sub-traces have different sample rates
					if sub_trace.stats.sampling_rate!=round(sub_trace.stats.sampling_rate):
						sub_trace.stats.sampling_rate=round(sub_trace.stats.sampling_rate)
				print('Merging gappy data...')
				tr.merge(fill_value=fill_value)
		except:
			print('Cannot get data for {}'.format(channel))
			continue
		st+=tr

	return st


def make_env(st,lowpass=0.2):

	st.detrend('demean')
	for tr in st:
		tr.resample(25.0)
		if tr.stats.npts % 2 ==1:
			tr.trim(starttime=tr.stats.starttime,endtime=tr.stats.endtime+1/tr.stats.sampling_rate,pad=True,fill_value=0)
		tr.data = envelope(tr.data)
		tr.resample(5.0)

	st.filter('lowpass',freq=lowpass)

	return st


def get_IRIS_data(sta_list,t1,t2,f1=1.0,f2=8.0,lowpass=0.2,client=None):

	if not client:
		client = Client('IRIS')
	print('Downloading data...\n')
	st = get_data(sta_list,t1,t2,client=client)
	print('Adding lat/lon info...\n')
	st = add_latlon(st,client=client)
	print('Detrending data...\n')
	st.detrend('demean')
	print('Tapering data...\n')
	st.taper(max_percentage=None,max_length=5)
	print('Filtering data...\n')
	st.filter('bandpass',freqmin=f1,freqmax=f2,corners=3,zerophase=True)
	print('Making envelopes...\n')
	env = make_env(st.copy(),lowpass=lowpass)
	print('Final trim...\n')
	env.trim(env[0].stats.starttime+dt,env[0].stats.endtime-dt)

	return env

def interactive_example():

	t1 = '2018-04-28 13:07'
	t2 = '2018-04-28 13:09'

	FREQMIN = 1.0
	FREQMAX = 8.0
	LOWPASS = 0.2

	sta_list=[
				'HV.BYL..HHZ',
				'HV.DEVL..HHZ',
				'HV.HAT..HHZ',
				'HV.KKO..HHZ',
				'HV.NPT..HHZ',
				'HV.OBL..HHZ',
				'HV.PAUD..HHZ',
				'HV.PUHI..HHZ',
				'HV.RIMD..HHZ',
				'HV.SBL..HHZ',
				'HV.SDH..HHZ',
				'HV.UWB..HHZ',
				'HV.UWE..HHZ',
				'HV.WRM..HHZ',
			]

	env = get_IRIS_data(sta_list,t1,t2,f1=FREQMIN,f2=FREQMAX,lowpass=LOWPASS,client=client)

	print('Creating XCOR object....\n')
	XC = XCOR(env)
	print('Locating from envelopes:\n')
	loc = XC.locate()

	return loc, XC


def test():

	t1 = '2018-04-28 13:07'
	t2 = '2018-04-28 13:09'

	FREQMIN = 1.0
	FREQMAX = 8.0
	LOWPASS = 0.2

	sta_list=[
				'HV.BYL..HHZ',
				'HV.DEVL..HHZ',
				'HV.HAT..HHZ',
				'HV.KKO..HHZ',
				'HV.NPT..HHZ',
				'HV.OBL..HHZ',
				'HV.PAUD..HHZ',
				'HV.PUHI..HHZ',
				'HV.RIMD..HHZ',
				'HV.SBL..HHZ',
				'HV.SDH..HHZ',
				'HV.UWB..HHZ',
				'HV.UWE..HHZ',
				'HV.WRM..HHZ',
			]

	print('Downloading some data for test...')
	env = get_IRIS_data(sta_list,t1,t2,f1=FREQMIN,f2=FREQMAX,lowpass=LOWPASS,client=client)

	print('Creating XCOR object....\n')
	XC = XCOR(env,plot=True,interact=False)
	print('Locating from envelopes:\n')
	loc = XC.locate()


def custom_grid_example():

	t1 = '2018-04-28 13:07'
	t2 = '2018-04-28 13:09'

	FREQMIN = 1.0
	FREQMAX = 8.0
	LOWPASS = 0.2

	sta_list=[
				'HV.BYL..HHZ',
				'HV.DEVL..HHZ',
				'HV.HAT..HHZ',
				'HV.KKO..HHZ',
				'HV.NPT..HHZ',
				'HV.OBL..HHZ',
				'HV.PAUD..HHZ',
				'HV.PUHI..HHZ',
				'HV.RIMD..HHZ',
				'HV.SBL..HHZ',
				'HV.SDH..HHZ',
				'HV.UWB..HHZ',
				'HV.UWE..HHZ',
				'HV.WRM..HHZ',
			]

	env = get_IRIS_data(sta_list,t1,t2,f1=FREQMIN,f2=FREQMAX,lowpass=LOWPASS,client=client)

	print('Creating XCOR object....\n')

	mygrid={ 'deps': np.arange( 0, 14 ,0.5),
			 'lons': np.arange(-155.35,-155.2,0.002),
			 'lats': np.arange(  19.35, 19.45,0.002)}

	XC = XCOR(env,grid_size=mygrid,interact=False)
	print('Locating from envelopes:\n')
	loc = XC.locate()

	return loc, XC


def rotate_grid_example():

	t1 = '2018-04-28 13:07'
	t2 = '2018-04-28 13:09'

	FREQMIN = 1.0
	FREQMAX = 8.0
	LOWPASS = 0.2

	sta_list=[
				'HV.BYL..HHZ',
				'HV.DEVL..HHZ',
				'HV.HAT..HHZ',
				'HV.KKO..HHZ',
				'HV.NPT..HHZ',
				'HV.OBL..HHZ',
				'HV.PAUD..HHZ',
				'HV.PUHI..HHZ',
				'HV.RIMD..HHZ',
				'HV.SBL..HHZ',
				'HV.SDH..HHZ',
				'HV.UWB..HHZ',
				'HV.UWE..HHZ',
				'HV.WRM..HHZ',
			]

	env = get_IRIS_data(sta_list,t1,t2,f1=FREQMIN,f2=FREQMAX,lowpass=LOWPASS,client=client)

	print('Creating XCOR object....\n')

	rotation={'x'    : np.arange(-5,5,0.3),
			  'y'    : np.arange(-3.5,3.5,0.3),
			  'z'    : np.arange(0,25,2),
			  'lat0' : 19.403,
			  'lon0' : -155.281,
			  'az'   : 30 }

	XC = XCOR(env,rotation=rotation,interact=False)
	print('Locating from envelopes:\n')
	loc = XC.locate()

	return loc, XC


def auto_example():
	
	t1 = '2017-05-01 13:00'
	t2 = '2017-05-01 15:00'

	FREQMIN = 0.5
	FREQMAX = 10
	LOWPASS = 0.2

	sta_list=[
				'HV.BYL..HHZ',
				'HV.DEVL..HHZ',
				'HV.HAT..HHZ',
				'HV.KKO..HHZ',
				'HV.NPT..HHZ',
				'HV.OBL..HHZ',
				'HV.PAUD..HHZ',
				'HV.PUHI..HHZ',
				'HV.RIMD..HHZ',
				'HV.SBL..HHZ',
				'HV.SDH..HHZ',
				'HV.UWB..HHZ',
				'HV.UWE..HHZ',
				'HV.WRM..HHZ',
			]

	env = get_IRIS_data(sta_list,t1,t2,f1=FREQMIN,f2=FREQMAX,lowpass=LOWPASS,client=client)

	print('Creating XCOR object....\n')
	XC = XCOR(env,bootstrap=20,plot=False,output=2)
	print('Locating from envelopes:\n')
	locs = XC.locate(window_length=240,step=120)

	return locs, env

def cascadia_example():
	

	t1 = '2020-05-24 04:52:30'
	t2 = '2020-05-24 05:07:30'

	FREQMIN = 1.5
	FREQMAX = 6.0
	LOWPASS = 0.1

	sta_list=[
				'UW.MCW.01.EHZ',
				'PB.B011.--.EHZ',
				'CN.SYMB.--.HHZ',
				'CN.PTRF.--.HHZ',
				'UW.RPW.01.EHZ',
				'CN.VGZ.--.HHZ',
				'PB.B004.--.EHZ',
				'UW.JCW.--.EHZ',
				'UW.SQM.--.EHZ',
				'PB.B003.--.EHZ',
				'PB.B006.--.EHZ',
				'PB.B001.--.EHZ',
				'UW.OBC.--.EHZ',
				'PB.B013.--.EHZ',
				'UW.DOSE.--.HHZ',
				'UW.HDW.--.EHZ',
				'UW.GNW.--.HHZ',
				'UW.GMW.--.EHZ',
				'PB.B014.--.EHZ',
				'UW.SMW.--.EHZ',
				'UW.STOR.--.HHZ',
				'UW.TKEY.--.HHZ',
			]

	env = get_IRIS_data(sta_list,t1,t2,f1=FREQMIN,f2=FREQMAX,lowpass=LOWPASS,client=client)

	mygrid = {'lons': np.arange(-125,-121+0.05,0.075),
			  'lats': np.arange(46.5,49.0+0.05,0.075),
			  'deps': np.arange(20,60+0.1,10)}

	print('Creating XCOR object....\n')
	XC = XCOR(env,grid_size=mygrid,regrid=True,bootstrap=30,output=2)
	print('Locating from envelopes:\n')
	loc = XC.locate()

	return loc, XC