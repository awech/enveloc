from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
import numpy as np
import sys
from enveloc import xcorr_utils
from enveloc import xloc_utils
from copy import copy, deepcopy
from itertools import combinations
from obspy.geodetics.base import gps2dist_azimuth
from datetime import timedelta


def progress(count, total, status=''):
	"""
	Method to display status bar when calculating traveltimes

	Parameters
	----------
	count : int, required
		number of current iteration
	total : int, required
		number of iterations required
	status : str, optional
		string to be printed to screen with status message

	"""
	
	bar_len = 35
	filled_len = int(round(bar_len * count / float(total)))
	percents = round(100.0 * count / float(total), 1)
	bar = '#' * filled_len + ' ' * (bar_len - filled_len)
	sys.stdout.write('Calculating travel times |{}| {}{} {}{}'.format(bar, percents, '%', status,'\r'))
	sys.stdout.flush()  # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)


def checkvalue(value):
	"""
	Function to switch from 'None' to NaN
	
	Parameters
	----------
	value : number, requred

	"""

	if isinstance(value,float) | isinstance(value,int):
		pass
	else:
		value=np.NaN
	return value


def latlon2xy(loc):
	"""
	Function to change lat/lon to x/y for distance calulation in clustering method

	Parameters
	----------
	loc : obj, required
		location object returned from `XC_locate`

	"""
	lats=np.array(loc.detections.get_lats())
	lons=np.array(loc.detections.get_lons())

	y0=np.nanmean(lats)
	x0=np.nanmean(lons)

	x=[]
	y=[]
	for lat,lon in zip(lats,lons):

	    a=gps2dist_azimuth(lat,lon,y0,lon)[0]/1000.
	    if lat<y0:
	    	a=-a
	    y.append(a)
	    a=gps2dist_azimuth(lat,lon,lat,x0)[0]/1000.
	    if lon<x0:
	    	a=-a
	    x.append(a)

	return x, y


def check_edgeproblems(edge_control,st_tmp):
	"""
	Function to remove traces whose peak amplitude occur with 'edge_control' percent
	of the window edge

	"""
	edge=np.round(edge_control*st_tmp[0].stats.npts)
	for tr in st_tmp:
		if tr.data.argmax()<=edge-1 or tr.data.argmax()>=tr.stats.npts-edge:
			st_tmp.remove(tr)
	return st_tmp


def XC_locate(win,XC):
	"""
	Main location function. Called internally from `XCOR` class's `locate()` method

	Parameters
	----------
	win : list, required
		2 element list of obspy UTCDateTime start and end times for which to cut out a window of seismic data
	XC : obj, required
		object created by `XCOR` class

	"""
	
	st_tmp=XC.traces.slice(win[0],win[1],nearest_sample=True)

	if XC.output > 0:
		print('{}  --  {}'.format(win[0].strftime('%Y.%m.%d %H:%M:%S'),win[1].strftime('%Y.%m.%d %H:%M:%S')))

	""" Remove traces with zero-ed out data """
	for tr in st_tmp:
		if np.any(tr.data<1e-15) and not XC._waveform_loc:
			st_tmp.remove(tr)
	if len(st_tmp) < XC.sta_min:
		if XC.output > 0:
			print('WARNING FROM XC_locate: too many traces with no data')
			return location(starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(st_tmp)) 	# return latitude, longitude, depth # return because these envelopes suck.

	""" Remove traces with flagged gaps """
	for tr in st_tmp:
		if np.any(tr.data==XC._gap_value) and not XC._waveform_loc:
			print('Gap in {} | Removing station...'.format(tr.id))
			st_tmp.remove(tr)
	if len(st_tmp) < XC.sta_min:
		if XC.output > 0:
			print('WARNING FROM XC_locate: too many traces with no data')
			return location(starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(st_tmp)) 	# return latitude, longitude, depth # return because these envelopes suck.


	""" Remove traces with triggers within time window """
	for tr in st_tmp:
		if hasattr(tr,'triggers'):
			if np.any(np.logical_and(tr.triggers>=tr.stats.starttime, tr.triggers<=tr.stats.endtime)):
				st_tmp.remove(tr)
	if len(st_tmp) < XC.sta_min:
		if XC.output > 0:
			print('WARNING FROM XC_locate: too many traces with triggers')
			return location(starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(st_tmp)) 	# return latitude, longitude, depth # return because these envelopes suck.


	####################################################
	"""	Check for traces with peaks at or near the edge.
		Remove problem traces and check that
		enough stations remain  """
	####################################################
	if XC._edge_control>0 and not XC._waveform_loc:
		st_tmp=check_edgeproblems(XC._edge_control,st_tmp)
		if np.unique([tr.stats.station for tr in st_tmp]).size < XC.sta_min:
			if XC.output > 0:
				print('WARNING FROM XC_locate: too many traces with peaks at edge of window')
			return location(starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(st_tmp)) 	# return latitude, longitude, depth # return because these envelopes suck.
	####################################################
	####################################################

	####################################################
	""" perform initial cross-correlation of all traces """
	if XC.detrend:
		st_tmp.detrend('demean')
	CC = xcorr_utils.initial_correlation(XC,st_tmp)
	""" return if not enough stations remain """
	if CC['err']:
		if XC.output > 1:
			print('WARNING FROM XC_locate: need at least {:.0f} stations to correlate above {:.2f}'.format(XC.sta_min,XC.Cmin))
		return location(starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(CC['st'])) 	# return latitude, longitude, depth # return because these envelopes suck.	
	####################################################
	####################################################

	""" make working copy of cross-correlation dictionary """
	CCrm = CC.copy()

	count=0
	drop_key=np.array([])
	""" you're going stay in the loop if it is the first time (count==0) or 
	    while there are stations to remove and the data need relocated from
	    interactive mode """
	while (drop_key>-1).any() or count==0:

		bstrap_lat=list()
		bstrap_lon=list()
		bstrap_dep=list()
		
		""" remove stations and values from cc-matrix if stations
		    get removed during interactive location process   """
		if (drop_key>-1).any():
			if XC.output > 1:
				print('Removing stations.')
			CCrm = xcorr_utils.remove_stations(CCrm,XC,drop_key)

		nan_count = 0
		""" Loop for bootstrap iterations. Last iteration is uses all of the data. """
		for bstrap in np.arange(XC.bootstrap):

			####################################################
			""" Check if you've had too many bad (NaN) bootstrap locations """
			if nan_count >= XC.bootstrap/2.:
				if XC.output > 1:
					print('Too many NaN bootstrap locations.')
				return location(starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(CCrm['st']))
			####################################################
			####################################################


			####################################################
			""" make boostrap working copy of cross-correlation dictionary """
			CCnew=CCrm.copy()
			""" If not the last iteration, throw away % of cc-matrix data  """
			if bstrap != np.arange(XC.bootstrap)[-1]:
				CCnew = xcorr_utils.bootstrap_CC(CCrm,XC)
			""" Not enough stations remain? """	
			if CCnew['err']:
				if XC.output > 1:
					print('WARNING FROM XC_locate: need at least {:.0f} stations to correlate above {:.2f}'.format(XC.sta_min,XC.Cmin))				
				""" Continue if bootstrap loop and iterate # of bad bstrap locations """
				if bstrap != np.arange(XC.bootstrap)[-1]:
					nan_count=nan_count+1
					continue
				### Return empty location object if it is last iteration   """
				else:
					return location(starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(CCrm['st'])) 	# return latitude, longitude, depth # return because these envelopes suck.
			####################################################
			####################################################


			""" Do a grid search to get a hypocenter """
			tmp_lat, tmp_lon, tmp_dep, misfit = xloc_utils.CCtoHypocenter(CCnew,XC)
			if 'windows' not in XC.__dict__ and XC._num_processors==1:
				XC.misfit=misfit


			####################################################
			""" Check if location is on grid edge """
			edge_check=False
			if XC.rotation:
				if np.in1d(tmp_lat,XC._grid['LAT'][0,:,0])[0] or np.in1d(tmp_lat,XC._grid['LAT'][-1,:,0])[0]:
					edge_check=True
				if np.in1d(tmp_lat,XC._grid['LON'][0,:,0])[0] or np.in1d(tmp_lat,XC._grid['LON'][-1,:,0])[0]:
					edge_check=True
			else:
				if np.in1d(tmp_lat,XC.grid_size['lats'].take([0,-1]))[0] or np.in1d(tmp_lon,XC.grid_size['lons'].take([0,-1]))[0]:
					edge_check=True
			if edge_check:	
				if XC.output > 1:
					print('On grid edge!! No location for you.')
				if bstrap != np.arange(XC.bootstrap)[-1]:
					nan_count=nan_count+1
					continue
				else:
					return location(starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(CCrm['st']))
			####################################################
			####################################################



			""" Regrid if XC.regrid==True """
			if XC.regrid:
				if XC.rotation:
					tmp_lat, tmp_lon, tmp_dep = xloc_utils.regridHypocenter_rotated(XC,misfit)
				else:
					tmp_lat, tmp_lon, tmp_dep = xloc_utils.regridHypocenter(XC,misfit)


			####################################################
			""" Append bootstrap location info if not last iteration """
			if bstrap != np.arange(XC.bootstrap)[-1]:
				if XC.output > 1:
					print('...bootstrap...Latitude: {:.3f}, Longitude: {:.3f}, Depth {:.1f}'.format(tmp_lat,tmp_lon,tmp_dep))
				bstrap_lat.append(tmp_lat)
				bstrap_lon.append(tmp_lon)
				bstrap_dep.append(tmp_dep)
			### Last iteration. Set location object info with this location. """
			else:
				loc=location(latitude=tmp_lat,longitude=tmp_lon,depth=tmp_dep,starttime=win[0],endtime=win[1],channels=xloc_utils.channel_list(CCrm['st']))
				""" Optionally calculate reduced displacement"""
				if 'raw_traces' in XC.__dict__:
					loc.reduced_displacement=xloc_utils.reduced_displacement(loc,XC,XC.raw_traces.slice(loc.starttime,loc.endtime))
				""" If bootstraping at all, calculate scatter from bootstrap locations """
				if XC.bootstrap > 1:
					dx,dz = xloc_utils.location_scatter(loc,bstrap_lat,bstrap_lon,bstrap_dep)
					loc.horizontal_scatter=dx
					loc.vertical_scatter=dz
					if XC.output > 0:
						print('--Location--\nLatitude: {:.3f}, Longitude: {:.3f}\t +/- {:.1f}\n--Depth--\n{:.1f} +/- {:.1f} km'.format(tmp_lat,tmp_lon,dx,tmp_dep,dz))
				else:
					if XC.output > 0:
						print('--Location--\nLatitude: {:.3f}, Longitude: {:.3f}, Depth: {:.1f} km'.format(tmp_lat,tmp_lon,tmp_dep))
			####################################################
			####################################################


		####################################################
		""" Done with bootstrap loop. Plot & interact """
		if XC.plot and XC._num_processors==1:
			from enveloc.plotting_utils import XC_plot
			""" set up some variables for plotting """
			CC1=1-CCrm['C'][:,CCrm['indx'][0,:]!=CCrm['indx'][1,:]]
			Nseis0=len(CCrm['st'])
			plot_loc=loc.copy()
			plot_loc.bstrap_lat=np.array(bstrap_lat)
			plot_loc.bstrap_lon=np.array(bstrap_lon)
			CCshift = np.arange(Nseis0*(Nseis0-1)/2)*(2*(XC._mlag+1)-1)-1 # index into cross-corrleation array
			""" Call plotting routine. Return plot_opt variable to see if 
			    we need to relocate """
			plot_opt = XC_plot(CCrm,XC,CC1,misfit,plot_loc)
			drop_key = plot_opt['krm']
			""" relocate with some stations removed """
			if plot_opt['relocate']:
				count = -1
			elif plot_opt['restart']:
				""" relocate with original station list """
				count = -1
				CCrm  = CC.copy()
				drop_key = np.array([])
		####################################################
		####################################################

		count = count + 1
	""" done with while loop looking for revomved station """

	return loc


class location(object):
	""" Define a location object """

	def __init__(self,latitude=None,longitude=None,depth=None,dx=None,dz=None,starttime=None,endtime=None,channels=[],reduced_displacement=None):
		self.latitude  = checkvalue(latitude)
		self.longitude = checkvalue(longitude)
		self.depth     = checkvalue(depth)
		self.horizontal_scatter = checkvalue(dx)
		self.vertical_scatter   = checkvalue(dz)
		self.starttime = starttime
		self.endtime  = endtime
		self.channels = channels
		self.reduced_displacement = reduced_displacement

	def __repr__(self):
		a=self.__dict__.keys()
		A=[]
		for item in a:
			A.append(item+'='+str(self.__getattribute__(item)))
		return 'location(' + '\n'.join(A) + ')'

	def station_latlons(self):
		""" return stas, np.array(lats), np.array(lons) for this location """
		lats=[]
		lons=[]
		stas=[]
		for c in self.channels:
			lats.append(c[2])
			lons.append(c[1])
			stas.append(c[0])

		return stas, np.array(lats), np.array(lons)

	def copy(self):
		return(deepcopy(self))


class event_list(object):
	"""
	Define an events object, which is a list of :class:`core.location` objects

	"""

	def __init__(self,events=[]):
		self.events=events

	def __repr__(self):
		return 'event_list object containing {:.0f} events'.format(len(self.events))

	def get_lats(self):
		""" 
		Return numpy array of all location latitudes in the list

		"""
		lats=np.array([a.latitude for a in self.events])
		return lats

	def get_lons(self):
		""" 
		Return numpy array of all location longitudes in the list

		"""
		lons=np.array([a.longitude for a in self.events])
		return lons

	def get_depths(self):
		""" 
		Return numpy array of all location depths in the list 

		"""
		lons=np.array([a.depth for a in self.events])
		return lons

	def get_times(self):
		""" 
		Return two numpy arrays of all location window start and end times in the list 

		"""
		starttimes=np.array([a.starttime for a in self.events])
		endtimes=np.array([a.endtime for a in self.events])
		return starttimes, endtimes

	def get_reduced_displacement(self):
		""" 
		Return all location reduced displacements in the events list as a numpy array 

		"""
		reduced_displacement=np.array([a.reduced_displacement for a in self.events])
		return reduced_displacement

	def get_numchannels(self):
		""" 
		Returns a numpy array of number of channels used in each location 

		"""
		num=np.array([len(a.station_latlons()[0]) for a in self.events])
		return num

	def get_stations(self):
		""" 
		Return a list of all stations used in all locations in event list 

		"""
		stas=[]
		for a in self.events:
			for sta in a.channels:
				stas.append(sta[0])
		return stas

	def tolist(self):
		"""
		returns the attribute `events` as a list

		"""
		return self.events

	def filter(self,min_t=None,max_t=None,min_lat=None,max_lat=None,min_lon=None,max_lon=None,min_horizontal_scatter=None,
			   max_horizontal_scatter=None,min_depth=None,max_depth=None,min_rd=None,max_rd=None,min_vertical_scatter=None,
			   max_vertical_scatter=None,min_num_channels=None,max_num_channels=None,highpass_loc=None,inplace=False):
		""" 
		Filter this event list by location properties and return a new event list object 

		Parameters
		----------
		Mostly self-explanatory

		"""

		NEW=self.copy()
		NEW.events=np.array(NEW.events)
		if min_t:
			NEW.events=NEW.events[NEW.get_times()[0]>min_t]
		if max_t:
			NEW.events=NEW.events[NEW.get_times()[0]<max_t]
		if min_lat:
			NEW.events=NEW.events[NEW.get_lats()>min_lat]
		if max_lat:
			NEW.events=NEW.events[NEW.get_lats()<max_lat]
		if min_lon:
			NEW.events=NEW.events[NEW.get_lons()>min_lon]
		if max_lon:
			NEW.events=NEW.events[NEW.get_lons()<max_lon]
		if min_horizontal_scatter:
			d_xy=np.array([a.horizontal_scatter for a in self.events])
			NEW.events=NEW.events[d_xy>min_horizontal_scatter]
		if max_horizontal_scatter:
			d_xy=np.array([a.horizontal_scatter for a in self.events])
			NEW.events=NEW.events[d_xy<max_horizontal_scatter]
		if min_depth:
			NEW.events=NEW.events[NEW.get_depths()>min_depth]
		if max_depth:
			NEW.events=NEW.events[NEW.get_depths()<max_depth]
		if min_rd:
			NEW.events=NEW.events[NEW.get_reduced_displacement()>min_rd]
		if max_rd:
			NEW.events=NEW.events[NEW.get_reduced_displacement()<max_rd]
		if min_vertical_scatter:
			d_xy=np.array([a.vertical_scatter for a in self.events])
			NEW.events=NEW.events[d_xy>min_vertical_scatter]
		if max_vertical_scatter:
			d_xy=np.array([a.vertical_scatter for a in self.events])
			NEW.events=NEW.events[d_xy<max_vertical_scatter]
		if min_num_channels:
			num=np.array([len(a.channels) for a in self.events])
			NEW.events=NEW.events[num>min_num_channels]
		if max_num_channels:
			num=np.array([len(a.channels) for a in self.events])
			NEW.events=NEW.events[num<max_num_channels]
		if highpass_loc is not None:
			vals=np.array([l.highpass_loc for l in self.events])
			if highpass_loc:
				NEW.events=NEW.events[vals]
			else:
				NEW.events=NEW.events[~vals]
		if inplace:
			setattr(self,'events',NEW.__getattribute__('events'))
		else:
			return NEW
			# return event_list(NEW.events.tolist())

	def calc_reduced_displacement(self,st,XC,num_processors=1):
		""" 
		Method to calculate a reduced displacement for all locations in the events list.
		The method updates the `reduced_displacement` property in each location object. `XCOR` can also do
		this internally while generating a location, but the idea with this method is to apply the calculation
		after obtaining all initial locations and throwing out bad locations via clustering or maximum horizontal
		scatter. That way you are not spending computation time for events you will never use.

		Parameters
		----------
	
		st : stream, required
			Obspy stream of envelopes. Each trace must have stats.coordinates.latitude/longitude
		XC : `XCOR` object, required
			`XCOR` object used to obtain locations within the event list
		num_processors : int, optional
			Optional integer to parallelize when set >1, but this doesn't appear to be faster. Default = 1.
		
		"""
		if num_processors > 1:
			from multiprocessing import Pool
			pool = Pool(processes=num_processors)
			new_windows=[]
			for l in self.events:
				if not np.isnan(l.latitude):
					new_windows.append(l)
			results=[pool.apply_async(xloc_utils.reduced_displacement,args=(win,XC,st.slice(win.starttime,win.endtime))) for win in new_windows]
			pool.close()
			pool.join()
			RD=[p.get() for p in results]
			times=np.array(self.get_times()[0])
			for i,rd in enumerate(RD):
				ind=np.where(times==new_windows[i].starttime)[0][0]
				self.events[ind].reduced_displacement=rd
		else:
			for l in self.events:
				if not np.isnan(l.latitude):
					l.reduced_displacement=xloc_utils.reduced_displacement(l,XC,st.slice(l.starttime,l.endtime))

	def remove(self,max_scatter=3.0,rm_nan_loc=True,rm_nan_err=True,inplace=False):
		"""
		Method to remove locations from the event list. Returns a new copy.

		Parameters
		----------

		max_scatter : float, optional
			maximum horizontal scatter allowed for a give locationDefault = 3.0
		rm_nan_loc : boolean, optional
			`True` will remove all locations with `nan` locations. Default = `True`
		rm_nan_err : boolean, optional
			  `True` will remove all locations with nan horizontal_scatter (not appropriate 
			  if not bootstrapping to obtain horizontal_scatter, eg. bootstrap=1), Default = `True`
		inplace : boolean, optional
			`True` to operate in place. `False` to return a copy. Default = `False`
		
		"""

		errs=np.array([a.horizontal_scatter for a in self.events])
		IND=np.arange(len(errs)).tolist()
		if rm_nan_loc:
			lats=self.get_lats()
			for ind in np.arange(len(errs)):
				if np.isnan(lats[ind]):
					IND.remove(ind)

		if rm_nan_err:
			for ind in np.arange(len(errs)):
				if np.isnan(errs[ind]):
					if ind in IND:
						IND.remove(ind)

		for ind in np.arange(len(errs)):
			if np.isfinite(errs[ind]) and errs[ind]>max_scatter:
				if ind in IND:
					IND.remove(ind)

		events=[self.events[ind] for ind in IND]
		if inplace:
			self.events=events
		else:
			return event_list(events)


	def remove_highpass(self,inplace=False):
		"""
		Method to remove locations that also have a location with a highpass filter
		
		Parameters
		----------
		
		inplace : boolean, optional
			`True` to operate in place or `False` to return a copy. Default = `False`

		"""

		NEW=self.copy()
		NEW.events=np.array(NEW.events)
		vals=np.array([l.highpass_loc for l in self.events])
		NEW.events=NEW.events[~vals]
		if inplace:
			setattr(self,'events',NEW.__getattribute__('events'))
		else:
			return NEW


	def remove_duplicates(self,distance=25.0,inplace=False):
		"""
		Method to remove locations within a list that have identical starttimes
		
		Parameters
		----------
		
		distance : float, optional
			maximum horizontal distance (km) below which the locations of events 
			with identical starttimes are averaged. Default = 25.0.
		inplace : boolean, optional
			`True` to operate in place or `False` to return a copy. Default = `False`
		
		"""
		T = self.get_times()[0]
		idx_sort = np.argsort(T)
		sorted_T_array = T[idx_sort]
		vals, idx_start, count = np.unique(sorted_T_array, return_counts=True,return_index=True)
		res = np.split(idx_sort, idx_start[1:])
		vals = vals[count > 1]
		res = filter(lambda x: x.size > 1, res)

		IND=np.arange(len(T)).tolist()
		for r in res:
			for comb in combinations(r,2):
				xy1  = self.events[comb[0]]
				xy2  = self.events[comb[1]]
				dist = gps2dist_azimuth(xy1.latitude,xy1.longitude,xy2.latitude,xy2.longitude)[0]/1000.
				if dist<distance:
					if len(self.events[comb[0]].channels)>len(self.events[comb[1]].channels):
						ind = comb[0]
					elif len(self.events[comb[0]].channels)<len(self.events[comb[1]].channels):
						ind = comb[1]
					elif len(self.events[comb[0]].channels) == len(self.events[comb[1]].channels):
						if np.isnan(self.events[comb[0]].horizontal_scatter) and not np.isnan(self.events[comb[1]].horizontal_scatter):
							ind = comb[1]
						elif not np.isnan(self.events[comb[0]].horizontal_scatter) and  np.isnan(self.events[comb[1]].horizontal_scatter):
							ind = comb[0]
						elif np.isnan(self.events[comb[0]].horizontal_scatter) and np.isnan(self.events[comb[1]].horizontal_scatter):
							ind = comb[0]
						else:
							scatter=np.array([self.events[comb[0]].horizontal_scatter, self.events[comb[1]].horizontal_scatter])
							ind = comb[scatter.argmin()]
					if ind in IND:
						IND.remove(ind)
		events=[self.events[ind] for ind in IND]
		if inplace:
			self.events=events
		else:
			return event_list(events)


	def cluster(self,dx=10,dt=60,num_events=4):
		""" 
		Cluster locations in the event_list using :meth:`sklearn.cluster.DBSCAN` 

		Parameters
		----------

		dx : float, optional
			Maximum horizontal distance in km. Default = 10.0
		dt : float, optional
		 	Maximum time difference between detections, in minutes. Default = 60
		num_events : int, optional
			Number of events required within `dx` and `dt`. Default = 4.


		The method returns a :class:`core.detections` object, which contains various different lists as attributes: 

		* detections - events from original event list
		* core_clustered - events who all meet the criteria
		* edge_clustered - events within `dx` & `dt` distance of core_clustered' event, but don't themselves have `num_events` within `dx` & `dt` of them
		* noise          - events that don't meet either criteria above
		* all_clustered  - core_clustered + edge_clustered combined for convenience 

		"""
		det=detections(self.events)
		det.cluster(dx=dx,dt=dt,num_events=num_events)
		return det


	def plot_locations(self,XC):
		"""
		Method to plot the locations for this `event_list`

		Parameters
		----------

		XC : `XCOR` object, required
			`XCOR` object used to create `event_list` object

		"""

		from enveloc.plotting_utils import plot_locations
		
		plot_locations(self,XC)

		return

	def __add__(self,other):
		return event_list(self.events+other.events)

	def copy(self):
		""" Return a copy of the events object """
		return(deepcopy(self))


class detections(object):
	""" 
	`detections` object. Contains different 'event_list' objects as properties. This
	object allows for lists all events, clustered events, and unclustered events to 
	exist all in one place and be modified using the same class methods.

	"""

	def __init__(self,detects=[]):
		if detects.__class__.__name__=='list':
			self.detections=event_list(detects)
		elif detects.__class__.__name__=='event_list':
			self.detections=detects

	def __repr__(self):
		a=self.__dict__.keys()
		A=[]
		for item in a:
			A.append(item+': '+str(len(self.__getattribute__(item).events)) +' events')
		return 'DETECTION object with attributes:\n(' + '\n'.join(A) + ')'

	def copy(self):
		return(deepcopy(self))

	def filter(self,min_t=None,max_t=None,min_lat=None,max_lat=None,min_lon=None,max_lon=None,min_horizontal_scatter=None,
		   max_horizontal_scatter=None,min_depth=None,max_depth=None,min_rd=None,max_rd=None,min_vertical_scatter=None,
		   max_vertical_scatter=None,min_num_channels=None,max_num_channels=None):
		"""
		Filter each event_list in the detections object by location properties and return a new detections object
		with modified event_lists

		Parameters
		----------
		Mostly self-explanatory

		"""
		NEW=self.copy()
		a=NEW.__dict__.keys()
		for item in a:
			setattr(NEW,item,NEW.__getattribute__(item).filter(min_t=min_t,max_t=max_t,min_lat=min_lat,max_lat=max_lat,min_lon=min_lon,
											  max_lon=max_lon,min_horizontal_scatter=min_horizontal_scatter,max_horizontal_scatter=max_horizontal_scatter,
											  min_depth=min_depth,max_depth=max_depth,min_rd=min_rd,max_rd=max_rd,min_vertical_scatter=min_vertical_scatter,
		   									  max_vertical_scatter=max_vertical_scatter,min_num_channels=min_num_channels,max_num_channels=max_num_channels))
		return NEW


	def cluster(self,dx=25,dt=None,num_events=4):
		""" 
		Cluster locations in the 'detections' event_list using :meth:`sklearn.cluster.DBSCAN` 

		Parameters
		----------

		dx : float, optional
			Maximum horizontal distance in km. Default = 25.0
		dt : float, optional
		 	Maximum time difference between detections, in minutes. Default = `None`
		num_events : int, optional
			Number of events required within `dx` and `dt`. Default = 4.


		The method leaves the 'detections' `event_list` untouched, but adds additional 
		event_list properties in place: 

		* core_clustered - events who all meet the criteria
		* edge_clustered - events within `dx` & `dt` distance of core_clustered' event, but don't themselves have `num_events` within `dx` & `dt` of them
		* noise          - events that don't meet either criteria above
		* all_clustered  - core_clustered + edge_clustered combined for convenience 

		"""

		from sklearn.cluster import DBSCAN
		x,y = latlon2xy(self)
		
		if not dt:
			X   = np.array([x,y]).T
		else:
			# scale time to match distance
			t   = self.detections.get_times()[0]
			dtime  = np.array([(t0.datetime-t.min().datetime).total_seconds()/60. for t0 in t])
			dtime  = dtime*(dx/np.float(dt))
			# put distance and time together
			X   = np.array([x,y,dtime]).T

		db  = DBSCAN(eps=dx, min_samples=num_events).fit(X)

		all_detects=np.where(db.labels_>-1)[0]
		tmp=[]
		for ind in all_detects:
			tmp.append(self.detections.events[ind])
		self.all_clustered=event_list(tmp)

		tmp=[]
		for ind in db.core_sample_indices_:
			tmp.append(self.detections.events[ind])
		self.core_clustered=event_list(tmp)

		tmp=[]
		for ind in all_detects:
			if ind not in db.core_sample_indices_:
				tmp.append(self.detections.events[ind])
		self.edge_clustered=event_list(tmp)

		tmp=[]
		noise_inds=np.where(db.labels_==-1)[0]
		for ind in noise_inds:
			tmp.append(self.detections.events[ind])
		self.noise=event_list(tmp)

	def uncluster(self):
		""" 
		Method to remove all event lists not named 'detections'. Operates in place.

		"""
		a=self.__dict__.keys()
		for item in a:
			if item != 'detections':
				self.__delattr__(item)

	def remove(self,max_scatter=3.0,rm_nan_loc=True,rm_nan_err=True,inplace=False):
		"""
		Method to remove locations from the 'detections' event list. Returns a new copy.

		Parameters
		----------

		max_scatter : float, optional
			maximum horizontal scatter allowed for a give locationDefault = 3.0
		rm_nan_loc : boolean, optional
			`True` will remove all locations with `nan` locations. Default = `True`
		rm_nan_err : boolean, optional
			  `True` will remove all locations with nan horizontal_scatter (not appropriate 
			  if not bootstrapping to obtain horizontal_scatter, eg. bootstrap=1), Default = `True`
		inplace : boolean, optional
			`True` to operate in place. `False` to return a copy. Default = `False`
		
		"""

		NEW = self.copy()
		for item in NEW.__dict__.keys():
			setattr(NEW,item,NEW.__getattribute__(item).remove(max_scatter=max_scatter,rm_nan_loc=rm_nan_loc,rm_nan_err=rm_nan_err))
		if inplace:
			for item in self.__dict__.keys():
				setattr(self,item,NEW.__getattribute__(item))
		else:
			return NEW

	def remove_duplicates(self,distance=25.0,inplace=False):
		"""
		Method to remove locations within a list that have identical starttimes

		Parameters
		----------

		distance : float, optional
			maximum horizontal distance (km) below which the locations of events with 
			identical starttimes are averaged. Default = 25.0.
		inplace : boolean, optional
			`True` to operate in place or `False` to return a copy. Default = `False`
		
		"""

		NEW = self.copy()
		for item in NEW.__dict__.keys():
			setattr(NEW,item,NEW.__getattribute__(item).remove_duplicates(distance=distance))
		if inplace:
			for item in self.__dict__.keys():
				setattr(self,item,NEW.__getattribute__(item))
		else:
			return NEW

	def plot_locations(self,XC):
		"""
		Method to plot the locations from each `event_list`

		Parameters
		----------

		XC : `XCOR` object, required
			`XCOR` object used to create `detections` object

		"""

		from enveloc.plotting_utils import plot_locations
		
		plot_locations(self,XC)

		return

	def __add__(self,other):
		NEW = self.copy()
		for item in NEW.__dict__.keys():
			try:
				setattr(NEW,item,NEW.__getattribute__(item)+other.__getattribute__(item))
			except:
				print('Cannot add unlike objects')
				return
		return NEW


class XCOR(object):

	"""
	XCOR is the main class object that needs to be set up in order to get a location.
	Most of the paramaters will self-assign, and in many cases that is fine. The only
	essential variable is the obspy stream 'st', but you will likely want to give thought
	and care to the grid and velocity model.
	
	Parameters
	----------
	
	st : stream, required
		Obspy stream of envelopes. Each trace must have stats.coordinates.latitude/longitude
	model : str, optional
		Obspy taup model or model file from which to calculate travel times.
		If this is not provided, a default file will be called.
	model_dir : str, optional
		Directory where model.npz file is place by obspy's taup program
		If not provided, it will default to within the enveloc module directory
	grid_size : dict, optional
		A dict with keys 'lats', 'lons', and 'deps' each of which are 1D monotonic numpy arrays. 
		If unprovided, this will be created internally based on the extent of the input stations.
	detrend : boolean, optional
		True/False to force a demean on obspy stream `st` or not.
		Default = `True`
	regrid : boolean, optional
		Whether to reinterpolate location on finer grid.
		Default = `True`
	phase_types : list, optional
		list of phases for which to calculate travel times. Ultimate travel time
		used for each grid node for each station will be the minimum calculated
		from the phases in the list provided. Also accepts ['Nkmps'] where N is a 
		constant velocity in km/s. This would be useful for inputing a surface wave 
		velocity, for example. See obspy taup documentation for more details on phase 
		nomenclature.
		Default = [ 's', 'S' ]
	normType : int, optional
		Integer 1 or 2 to use an L1 or L2 norm.
		Default = 1
	plot : boolean, optional
		Boolean flag to plot or not.
		Default = `True`
	interact : boolean, optional
		Boolean flag to interact with plot or not.
		Default = `True`
	output : int, optional
		Interger from 0-3 controlling increasing level of messages the code produces
		Default = 1
	Cmin : float, optional
		Minimum normalized cross-correlation coefficient to be considered.
		Default = 0.5
	Cmax : float, optional
		Maximum normalized cross-correlation coefficient to be considered. Not always 
		necessary, but could be for removing correlated noise spikes, for example.
		Default = 0.995
	sta_min : int, optional
		Minimum number of stations needed to obtain a location.
		Default = 3
	dx_min : float, optional
		Minimum horizontal distance (km) between channels required to consider th
		cross-correlation. Useful to exclude correlations between multiple channels
		at the same station.
		Default = 0.1
	bootstrap : int, optional
		Number of iterations for each location. For each iteration, 'bootstrap_prct'
		of the CC values will be zeroed out, and the resulting station contributions
		are adjusted accordingly to obtain a new location. These iterated locations
		provide a measure of location scatter, which is recorded in the location object.
		The final location will use all of the data (no CC values zeroed out). Thus,
		if boootstrap=N, it will bootstrap N times and get a location using all the data
		on the N+1 iteration. If bootstrap==1, the data are only located once with all the
		data, and no horizontal scatter is estimated.
		Default = 1
	bootstrap_prct : float, optional
		Value between 0 and 1 determining the fraction of cc data to throw out in each
		bootstrap iteration.
		Default = 0.04 (4%)
	lookup_type : str, optional
		Type of interpolation method to use when getting a predicted cross-correlation
		value for each channel pair for each grid ('linear', 'nearest', 'zero', 'slinear', 
		'quadratic', 'cubic'). See scipy.interpolate.interp1d for details.
		Default='cubic'
	rotation : dict, optinal
		Dictionary defining desired rotated grid. Needs keys 'x' (x grid nodes), 
		'y' (y grid nodes), 'z' (depth grid nodes), 'lat0' & 'lon0' (origin lat/lon),
		'az', rotation azimuth (counterclockwise from East)
		Default = `None`
	dTmax_s : float, optional
		Maximum cross-correlation shift in seconds. Defaults to the smaller of:
			a) 1/2 of the obspy stream window length
			b) the maximum predicted inter-station differential time + `dt`
	dt : float, optional
		Seconds of additional allowed cross-correlation shifts beyond the maximumpredicted 
		interchannel differential time. 
		Default = 3.0
	rd_freq	: float, optional
		Frequency at which to calculate surface wave reduced displacement using the
		velocity from the top layer of the velocity model. If not set, body wave reduced 
		displacement will be calculated.
		Default = `None`
	raw_traces : stream, optional
		Obpsy stream of non-envelope traces matching the input `st` variable. Used to
		calculate reduced displacement.
		Default = `[]`
	env_hp : stream, optional
		Obspy stream of second envelopes with traces matchin the input `st` variable.
		This would typically be envelopes generated in a different passband than `st`.
		If a location is obtained for st, it will also check env_hp and set a flag.
		This can be useful for weeding out earthquake detections or correlated noise
		spikes when you are only interested in detecting a band-limited signal. It
		essentially gives you a way to flag windows where you have coherent energy in
		two different passbands. Kind of a niche-use case, but necessary for my needs.
		Default = `[]`
	edge_control : float, optional
		Value between 0 and 1 defining a percentage from the window edge to not allow
		the peak trace energy. That is, if any trace's peak amplitude occurs between a
		with a percentage of the window edge, it is removed. Useful when autolocating.
		Default = 0.03
	num_processors : int, optional
		Number of processors to use if you are iterating through sliding windows to 
		get locations. This can also be set in the `locate()` method.
		Default = 1.
	tt_file : str, optional
		Path to a .npz file containing pre-calculated traveltimes for this station set
		on this grid. Calulated using `save_traveltimes()` method.
		Default = `None`
	waveform_loc : boolean, optional
		Boolean flag to optionally locate waveforms rather than envelopes.
		This flag exists to bypass some of the built in envelope quailty control	
		Default = `False`
	gap_value : int, optional
		Value used to flag data gaps identified in preprocessing.
		Default = -123454321

	return : object

	"""

	def __init__(self,st,model=None,grid_size=None,detrend=True,regrid=False,phase_types=['s','S'],
				 	     normType=1,plot=True,interact=True,output=1,
				 	     Cmin=0.5,Cmax=0.995,sta_min=3,dx_min=0.1,dt=3.0,bootstrap=1,
				 	     bootstrap_prct=0.04,lookup_type='cubic',rotation=None,dTmax_s=None,
				 	     rd_freq=None,raw_traces=[],env_hp=[],edge_control=0.03,num_processors=1,
						 tt_file=None,waveform_loc=False,model_dir=None,gap_value=-123454321):

		if output==True:
			output = 1
		elif output==False:
			output = 0
		self.output    = output
		self.traces    = st
		self.detrend   = detrend
		self.rotation  = rotation
		if not tt_file:
			import os
			if not model_dir:
				self._model_dir=os.path.dirname(os.path.realpath(__file__))
			else:
				self._model_dir=model_dir
			if not model:
				print('Warning!  No model provided.')
				print('Default model used. Maybe not ideal...')
				self._model_file=os.path.dirname(os.path.realpath(__file__))+'/default_vel_model.tvel'
				print('Using '+self._model_file)
			else:
				self._model_file=model

			if type(phase_types)==str:
				phase_types=[phase_types]
			self.phase_types = phase_types

			if grid_size:
				self.grid_size=grid_size
			else:
				if self.rotation:
					self.grid_size=dict({'x':self.rotation['x'],'y':self.rotation['y'],'deps':self.rotation['z']})
				else:
					self.grid_size=dict({'lats':[],'lons':[],'deps':[]})

					print('Warning!  No grid provided.')
					print('Making grid based on stations provided...')

					lats = np.array([tr.stats.coordinates.latitude for tr in st])
					lons = np.array([tr.stats.coordinates.longitude for tr in st])
					dlat = 0.33*(lats.max()-lats.min())
					dlon = 0.33*(lons.max()-lons.min())
					self.grid_size['lats'] = np.linspace(lats.min()-dlat,lats.max()+dlat,20)
					self.grid_size['lons'] = np.linspace(lons.min()-dlon,lons.max()+dlon,25)
					tmp = gps2dist_azimuth(self.grid_size['lats'].max(),self.grid_size['lons'].max(),
										   self.grid_size['lats'].min(),self.grid_size['lons'].min())
					max_depth = np.min([tmp[0]/1000, 60.])
					self.grid_size['deps'] = np.linspace(0,max_depth,10)

			self.calculate_traveltimes()

			if 'kmps' in self.phase_types[0]:
				self._v0     = float(self.phase_types[0].split('kmps')[0])
			else:
				self._v0     = self.model.model.s_mod.v_mod.layers[0][4]
		else:
			self.load_traveltimes(tt_file)
			if 'kmps' in self.phase_types[0]:
				self._v0     = float(self.phase_types[0].split('kmps')[0])
			else:
				self._v0     = self._model_layers[0][4]

		self.regrid      = regrid
		self._normType   = normType
		self.plot        = plot
		
		self.interact    = interact
		if not self.plot:
			self.interact = False
		
		self._p          = np.array([4.9650, -18.9378, 26.7329, -16.5927, 3.8438]) # weight polynomial coefficients
		self.Cmin        = Cmin
		self._Cmax       = Cmax
		self._dx_min     = dx_min
		self.sta_min     = sta_min
		self.bootstrap    = bootstrap
		if bootstrap > 1:
			self.bootstrap = bootstrap + 1
		if self.bootstrap < 1 or self.bootstrap-np.round(self.bootstrap) != 0:
			print('Warning: \'bootstrap\' value must be an integer > 0')
			print('Setting bootstrap = 1')
			self.bootstrap = 1

		self.bootstrap_prct = bootstrap_prct	
		if lookup_type not in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']:
			print('\'lookup_type\' needs to be: \'linear\', \'nearest\', \'zero\', \'slinear\', \'quadratic\' or \'cubic\'.')
			print('Setting lookup_type to \'cubic\'.')
			lookup_type = 'cubic'
		self.lookup_type  = lookup_type

		self._dt         = dt
		
		if dTmax_s:
			self.dTmax_s     = dTmax_s
		else:
			seconds=[(tr.stats.npts-1)/tr.stats.sampling_rate for tr in st]
			self.dTmax_s = np.min([.5*np.min(seconds),xcorr_utils.station_distances(st).max()/self._v0+self._dt])

		dTmax 	= np.round(float(self.dTmax_s)*self.traces[0].stats.sampling_rate)	# can't make a total shift of more than dTmax samples	
		# self._mlag         = int(2*dTmax) # change on 0ct-30-2017
		self._mlag         = int(dTmax)
		
		if raw_traces:
			self.raw_traces = raw_traces
			if rd_freq:
				self.rd_freq = rd_freq

		if env_hp:
			self.env_hp = env_hp
			for tr in self.env_hp:
				tr.TT=self.traces.select(id=tr.id)[0].TT

		self._edge_control   = edge_control
		self._num_processors = num_processors
		self._waveform_loc   = waveform_loc
		self._gap_value = gap_value


	def __repr__(self):
		a=self.__dict__.keys()
		A=[]
		for item in a:
			A.append(item+'='+str(self.__getattribute__(item)))
		return 'parameters(' + '\n'.join(A) + ')'


	def calculate_traveltimes(self):
		"""
		Calculate travel times based on the model and grid provided in the `XCOR` object

		"""

		if 'default_vel_model' in self._model_file:
			try:
				self.model = TauPyModel(model=self._model_dir+'/default_vel_model')
			except:
				from obspy.taup.taup_create import build_taup_model
				build_taup_model(self._model_file,output_folder=self._model_dir)
				self._model_file=self._model_dir+'/'+self._model_file.split('.tvel')[0].split('/')[-1]
				self.model = TauPyModel(model=self._model_file)
		else:
			if self._model_file.split('.')[-1]=='tvel':
				from obspy.taup.taup_create import build_taup_model
				build_taup_model(self._model_file,output_folder=self._model_dir)
				self._model_file=self._model_dir+'/'+self._model_file.split('.tvel')[0].split('/')[-1]
		
			self.model = TauPyModel(model=self._model_file)

		self._grid=dict({'LON':[],'LAT':[],'DEP':[]})
		if self.rotation:
			self._grid['LON'],self._grid['LAT'],self._grid['DEP'] = xloc_utils.create_rotated_grid(self.rotation)
		else:
			self._grid['LON'],self._grid['LAT'],self._grid['DEP'] = np.meshgrid(self.grid_size['lons'],self.grid_size['lats'],self.grid_size['deps'])

		if self.output > 1:
			progress(0, len(self.grid_size['deps']), status='complete')
		"""
		if self.tt_calc=='full':
			for i,tr in enumerate(self.traces):
				TIMES=[]
				for j,lon in enumerate(self._grid['LON'].flatten()):
					deg=locations2degrees(tr.stats.coordinates.latitude,tr.stats.coordinates.longitude,self._grid['LAT'].flatten()[j],lon)
					PATHS=self.model.get_ray_paths(source_depth_in_km=self._grid['DEP'].flatten()[j],distance_in_degree=deg,receiver_depth_in_km=0,phase_list=self.phase_types)
					times=np.array([p.time for p in PATHS])
					TIMES.append(times.min())

				tr.TT=np.reshape(TIMES,self._grid['LON'].shape)
				if self.output:
					progress(i+1, len(self.traces), status='complete')
		else: 
		"""
		if self.rotation:
			LON=self._grid['LON'][:,:,0]
			LAT=self._grid['LAT'][:,:,0]
		else:
			LON,LAT = np.meshgrid(self.grid_size['lons'],self.grid_size['lats'])

		lat_corners = LAT[::LAT.shape[0]-1, ::LAT.shape[1]-1]
		lon_corners = LON[::LON.shape[0]-1, ::LON.shape[1]-1]
		dist=[]
		for tr in self.traces:
			tr.TT=np.zeros(LAT.shape+(len(self.grid_size['deps']),))
			for y,x in zip(lat_corners.flatten(),lon_corners.flatten()):
				dist.append(locations2degrees(tr.stats.coordinates.latitude,tr.stats.coordinates.longitude,y,x))

		max_deg=np.array(dist).max()
		DEGS=np.linspace(0,max_deg,50)
		for i,z in enumerate(self.grid_size['deps']):

			TIMES=[]
			for x in DEGS:
				try:
					PATHS=self.model.get_ray_paths(source_depth_in_km=z,distance_in_degree=x,receiver_depth_in_km=0,phase_list=self.phase_types)
				except:
					# an error can occur when a grid node is on a velocity model boundary. When this occurs, add a slight deviation to make sure
					# get_ray_paths runs correctly
					small_perturbation=0.00001
					PATHS=self.model.get_ray_paths(source_depth_in_km=z+small_perturbation,distance_in_degree=x,receiver_depth_in_km=0,phase_list=self.phase_types)
				if self.output > 2:
					print(self.phase_types)
					print(PATHS)
				times=np.array([p.time for p in PATHS])
				TIMES.append(times.min())

			for tr in self.traces:
				deg=np.array([locations2degrees(tr.stats.coordinates.latitude,tr.stats.coordinates.longitude,y,x) 
													for y,x in zip(LAT.flatten(),LON.flatten())])
				inter_times=np.interp(deg,DEGS,np.array(TIMES))
				TIM = np.reshape(inter_times,LAT.shape)
				tr.TT[:,:,i]=TIM
			
			if self.output > 1:
				progress(i+1, len(self.grid_size['deps']), status='complete')
		print('\n')


	def save_traveltimes(self,outfile):
		"""
		Save an .npz file containing the model and grid. The variables should be created 
		using the `calculate_traveltimes()` method, which is called internally  when initiating
		an XCOR object with a model and grid.

		Parameters
		----------
		file : str, required
			The .npz file you wish to save, ideally with full path provided.
		
		"""

		model = self.model.model.s_mod.v_mod.layers
		phase_types = self.phase_types
		deps  = self.grid_size['deps']

		if self.rotation:
			Xs  = self.grid_size['x']
			Ys  = self.grid_size['y']
		else:
			lats  = self.grid_size['lats']
			lons  = self.grid_size['lons']

		stas = {}
		for tr in self.traces:
			stas[tr.id.replace('.','_')]=tr.TT
		
		if self.rotation:
			np.savez(outfile,y=Ys,     x=Xs,     deps=deps,model=model,phase_types=phase_types,stations=stas,**stas)
		else:
			np.savez(outfile,lats=lats,lons=lons,deps=deps,model=model,phase_types=phase_types,stations=stas,**stas)


	def load_traveltimes(self,file):
		""" 
		Load an .npz file containing a model and grid, and assign the appropriate times
		to each trace in the object. The .npz file can be created initially using the
		`save_traveltimes()` method after `XCOR` has run `calculate_traveltimes()` with a model and grid.

		Parameters
		----------
		file : str, required
			The path/filename to the .npz file you wish to load

		"""
		npzfile=np.load(file)
		self._grid=dict({'LON':[],'LAT':[],'DEP':[]})
		if self.rotation:
			self.grid_size={'x'   :npzfile['x'],
							'y'   :npzfile['y'],
							'deps':npzfile['deps']}
			self._grid['LON'],self._grid['LAT'],self._grid['DEP'] = xloc_utils.create_rotated_grid(self.rotation)
		else:
			self.grid_size={'lons':npzfile['lons'],
							'lats':npzfile['lats'],
							'deps':npzfile['deps']}
			self._grid['LON'],self._grid['LAT'],self._grid['DEP'] = np.meshgrid(self.grid_size['lons'],self.grid_size['lats'],self.grid_size['deps'])
		self._model_layers=npzfile['model']
		self.phase_types=npzfile['phase_types'].tolist()
		for tr in self.traces:
			tr.TT=npzfile[tr.id.replace('.','_')]
			if np.shape(tr.TT)!=np.shape(self._grid['LON']):
				print('Error loading file. Dimension mismatch for station '+tr.id)


	def locate(self,window_length=None,step=None,offset=0, include_partial_windows=False, nearest_sample=True, dTmax_s=None,num_processors=None):
		""" 
		Main method to locate data provided. This method is actually a wrapper to call main envelope 
		cross correlation location algorithm `XC_locate`.
		
		If using this to locate a single window of data, simply call locate() with no input parameters.
		All of the relevant details (plot, interact, map_resolution) are properties set in the initial
		object creation. If using this to iterate over windows and locate, the following parameters must be set:

		Parameters
		----------
		window_length : float, optional
			Length of each sub-window in seconds.
			Default = `None`
		step : float, optional
			The step between the start times of two successive windows in seconds. Can be negative if an offset is given.
			Default = `None`


		These parameters are optional:

		Parameters
		----------
		offset : float, optional
			The offset of the first window in seconds relative to the start time of the whole interval.
			Default = `False`
		include_partial_windows : boolean, optional
			Determines if sub-windows that are shorter then 99.9 % of the desired length are returned.
			Default = `False`
		nearest_sample : boolean, optional
			If set to True, the closest sample is selected, if set to False, the inner (next sample for a 
			start time border, previous sample for an end time border) sample containing the time is selected.
			Default = `True`
		dTmax_s : float, optional
			Maximum cross-correlation shift in seconds. Defaults to the smaller of:
				a) 1/2 of the obspy stream sub-window length
				b) the maximum predicted inter-station differential time + 'dt', set during initial object creation
		num_processors : int, optional
			Number of processors to use in parallel
			Default = `None`

		"""
		if num_processors:
			self._num_processors=num_processors
		if not window_length:
			if dTmax_s:
				self.dTmax_s = dTmax_s
				dTmax 	     = np.round(float(self.dTmax_s)*self.traces[0].stats.sampling_rate)	# can't make a total shift of more than dTmax samples	
				# self._mlag         = int(2*dTmax) # change on 0ct-30-2017
				self._mlag         = int(dTmax)
			
			loc=XC_locate((self.traces[0].stats.starttime,self.traces[0].stats.endtime),self)
		else:
			self._windows=True
			if dTmax_s:
				self.dTmax_s = dTmax_s
			else:
				self.dTmax_s = np.min([.5*window_length,xcorr_utils.station_distances(self.traces).max()/self._v0+self._dt])

			dTmax 	   = np.round(float(self.dTmax_s)*self.traces[0].stats.sampling_rate)	# can't make a total shift of more than dTmax samples	
			# self._mlag         = int(2*dTmax) # change on 0ct-30-2017
			self._mlag         = int(dTmax)

			if not step:
				print('Set step size')
				return

			from obspy.core.util.misc import get_window_times
			windows = get_window_times(starttime=self.traces[0].stats.starttime,endtime=self.traces[0].stats.endtime,window_length=window_length,
									   step=step,offset=offset,include_partial_windows=include_partial_windows)

			bp_traces      = self.traces.copy()
			if 'env_hp' in self.__dict__:
				env_hp=True
				hp_traces     = self.env_hp.copy()
				self.env_hp   = []
			else:
				env_hp=False

			if 'raw_traces' in self.__dict__:
				rd_flag    = True
				raw_traces = self.raw_traces.copy()
				self.__delattr__('raw_traces')
			else:
				rd_flag    = False

			if self._num_processors>1:
				from multiprocessing import Pool
				pool = Pool(processes=self._num_processors)
				results=[pool.apply_async(XC_locate,args=(win,self)) for win in windows]

				pool.close()
				pool.join()
				loc=event_list([p.get() for p in results])
				
				self.traces=[]
				if env_hp:
					new_windows=[]
					for i,l in enumerate(loc.events):
						l.highpass_loc=False
						if not np.isnan(l.latitude):
							new_windows.append(windows[i])
					if self.output > 0:
						print('Calculating highpass locations for {:.0f} locations'.format(len(new_windows)))
					tmp_bootstrap  = copy(self.bootstrap)
					tmp_regrid     = copy(self.regrid)
					self.traces    = hp_traces
					self.bootstrap = 10
					self.regrid    = False

					pool = Pool(processes=self._num_processors)
					results=[pool.apply_async(XC_locate,args=(win,self)) for win in new_windows]
					pool.close()
					pool.join()
					# loc_hp=event_list([p.get() for p in results]).remove(max_scatter=10)
					loc_hp=event_list([p.get() for p in results])
					print(len(loc_hp.events))
					loc_hp.remove(max_scatter=5,rm_nan_loc=True,rm_nan_err=True,inplace=True)
					print(len(loc_hp.events))
					times=np.array(loc.get_times()[0])
					for t in loc_hp.get_times()[0]:
						ind=np.where(times==t)[0][0]
						loc.events[ind].highpass_loc=True
					self.traces=[]

				if rd_flag:
					pool = Pool(processes=self._num_processors)
					new_windows=[]
					for i,l in enumerate(loc.events):
						if not np.isnan(l.latitude):
							new_windows.append((windows[i],l))
					if self.output > 0:
						print('Calculating reduced displacement for {:.0f} locations'.format(len(new_windows)))
					results=[pool.apply_async(xloc_utils.reduced_displacement,args=(win[1],self,raw_traces.slice(win[0][0],win[0][1]))) for win in new_windows]
					pool.close()
					pool.join()
					RD=[p.get() for p in results]
					if not env_hp:
						times=np.array(loc.get_times()[0])
					for i,rd in enumerate(RD):
						ind=np.where(times==new_windows[i][0][0])[0][0]
						loc.events[ind].reduced_displacement=rd
					self.raw_traces = raw_traces

			else:
				loc=event_list([XC_locate(win,self) for win in windows])
				if env_hp:
					if self.output > 0:
						print('Calculating highpass locations...')
					tmp_bootstrap  = copy(self.bootstrap)
					tmp_regrid     = copy(self.regrid)
					self.traces    = hp_traces
					self.bootstrap = 10
					self.regrid    = False
					loc_hp=[]
					for i,l in enumerate(loc.events):
						l.highpass_loc=False
						if not np.isnan(l.latitude):
							loc_hp.append(XC_locate(windows[i],self))
					loc_hp=event_list(loc_hp)
					loc_hp.remove(max_scatter=5,rm_nan_loc=True,rm_nan_err=True,inplace=True)
					times=np.array(loc.get_times()[0])
					for t in loc_hp.get_times()[0]:
						ind=np.where(times==t)[0][0]
						loc.events[ind].highpass_loc=True
					self.traces=[]
				if rd_flag:
					new_windows=[]
					for i,l in enumerate(loc.events):
						if not np.isnan(l.latitude):
							new_windows.append((windows[i],l))
					if self.output > 0:
						print('Calculating reduced displacement for {:.0f} locations'.format(len(new_windows)))
					RD=[xloc_utils.reduced_displacement(win[1],self,raw_traces.slice(win[0][0],win[0][1])) for win in new_windows]
					if not env_hp:
						times=np.array(loc.get_times()[0])
					for i,rd in enumerate(RD):
						ind=np.where(times==new_windows[i][0][0])[0][0]
						loc.events[ind].reduced_displacement=rd
					self.raw_traces = raw_traces
			
			self.traces    = bp_traces
			if env_hp:
				self.env_hp    = hp_traces
				self.bootstrap = tmp_bootstrap
				self.regrid    = tmp_regrid

		return loc


	def plot_grid(self):
		from  enveloc import plotting_utils
		plotting_utils.grid_plot(self)
