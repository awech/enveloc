from obspy import UTCDateTime, Stream
from obspy.core.util import AttribDict

dt = 10

def add_latlon(st,client=None):
	if not client:
		from obspy.clients.fdsn import Client
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
		from obspy.clients.fdsn import Client
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
					# deal with error when sub-traces have different dtypes
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


def match_sample_rate(st):
	sr = min([tr.stats.sampling_rate for tr in st])

	for tr in st:
		if tr.stats.sampling_rate == 4 * sr:
			tr.decimate(2)
		if tr.stats.sampling_rate == 2 * sr:
			tr.decimate(2)
		if tr.stats.sampling_rate != sr:
			tr.resample(sr)

	if sr == 100:
		st.decimate(2)
		st.decimate(2)
	if sr == 50:
		st.decimate(2)
	elif sr == 40:
		st.decimate(2)

	return st


def make_env(st,lowpass=0.2):
	from obspy.signal.filter import envelope

	for tr in st:
	    tr.data = envelope(tr.data)
	st.filter('lowpass',freq=lowpass)

	sr = st[0].stats.sampling_rate

	if sr == 25:
		st.decimate(5)
	elif sr == 20:
		st.decimate(2)
		st.decimate(2)

	return st


def get_IRIS_data(sta_list,t1,t2,f1=1.0,f2=8.0,lowpass=0.2,client=None):
	if not client:
		from obspy.clients.fdsn import Client
		client = Client('IRIS')
	st = get_data(sta_list,t1,t2,client=client)
	st = add_latlon(st,client=client)
	st = match_sample_rate(st)
	st.filter('bandpass',freqmin=f1,freqmax=f2)
	st = st.remove_response(output='DISP')
	env = make_env(st.copy(),lowpass=lowpass)
	st.trim(st[0].stats.starttime+dt,st[0].stats.endtime-dt)
	env.trim(env[0].stats.starttime+dt,env[0].stats.endtime-dt)

	return st, env

