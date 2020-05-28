from collections import Counter
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.dates as mdates
from matplotlib.widgets import Button
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


trace_ax  = []
h_trace   = []
h_5       = []
h_maptext = []
h_mapdot  = []
fig       = []
krm       = []
sort_key  = []

relocate = False
restart  = False
plot_opt = dict({'krm':      krm,
	             'relocate': relocate,
	             'restart':  restart})

st_dot_color0  = 'y'
st_dot_color1  = 'darkorange'
st_text_color0 = 'darkblue'
st_text_color1 = 'darkred'
trace_color0   = 'k'
trace_color1   = 'r'
xcorr_color0   = 'royalblue'
xcorr_color1   = 'r'

interval = 0.05

def on_click(event):
	global krm
	if event.inaxes == trace_ax:
		ind = int(np.round(event.ydata))
		if ind > len(h_trace)-1:
			return
		if ind<0:
			return
		if ind in krm:
			krm.remove(ind)
			h_trace[ind][0].set_color(trace_color0)
			h_trace[ind][0].set_linewidth(0.2)
			h_mapdot[ind][0].set_color(st_dot_color0)
			h_mapdot[ind][0].set_zorder(100)
			h_maptext[ind].set_color(st_text_color0)
			h_maptext[ind].set_zorder(100)
			# h_maptext[ind].set_weight('normal')
			for t in h_5[ind]:
				t[0].set_color(xcorr_color0)
				t[0].set_linewidth(0.5)
				t[0].set_zorder(100)

		else:
			krm.append(ind)
			h_trace[ind][0].set_color(trace_color1)
			h_trace[ind][0].set_linewidth(0.5)
			h_mapdot[ind][0].set_color(st_dot_color1)
			h_mapdot[ind][0].set_zorder(1000)
			h_maptext[ind].set_color(st_text_color1)
			h_maptext[ind].set_zorder(1000)
			# h_maptext[ind].set_weight('bold')
			for t in h_5[ind]:
				t[0].set_color(xcorr_color1)
				t[0].set_linewidth(1.0)
				t[0].set_zorder(1000)

		fig.canvas.draw()


class Index(object):

	def done(self, event):
		global krm
		global relocate
		global restart
		krm=[]
		restart=False
		relocate=False
		plt.pause(interval)
		plt.close()

	def reloc(self, event):
		global relocate
		global restart
		restart=False
		relocate=True
		plt.pause(interval)
		plt.close()

	def restart(self, event):
		global relocate
		global restart
		restart=True
		relocate=False
		plt.pause(interval)
		plt.close()


def XC_plot(CC,XC,CC1,misfit,loc):
	del krm[:]
	plot_opt = dict({'krm':       [],
	             'relocate':   False,
	             'restart': False})
	import warnings
	warnings.filterwarnings('ignore')

	global fig
	global trace_ax
	global h_trace
	global h_5
	global h_maptext
	global h_mapdot
	global sort_key

	trace_ax  = []
	h_trace   = []
	h_5       = []
	h_maptext = []
	h_mapdot  = []
	fig       = []


	######## calculate travel times ########
	########################################
	st_sort  = CC['st'].copy()
	Nseis0=len(CC['st'])
	dTs = np.zeros((Nseis0,))
	ind_min=misfit.argmin()
	for ksta in range(Nseis0):
		dTs[ksta]=CC['st'][ksta].TT.flatten(order='F')[ind_min]
		st_sort[ksta].dTs=dTs[ksta]
		st_sort[ksta].staW=CC['staW'][ksta]
	dT = dTs*CC['st'][0].stats.sampling_rate
	deltaT = dT[CC['j1']]-dT[CC['i1']] + XC._mlag
	sort_key = np.argsort(dTs)[::-1]
	st_sort.traces.sort(key=lambda x: x.dTs,reverse=True)
	########################################
	########################################


	gs = gridspec.GridSpec(60, 20)
	fig = plt.figure(num='XC loc',figsize=(6,7.2))


	############### plot map ###############
	########################################
	map_ax = fig.add_subplot(gs[:25,:13], projection=ccrs.PlateCarree())
	if XC.rotation:
		latlims=[XC._grid['LAT'].flatten().min(), XC._grid['LAT'].flatten().max()]
		lonlims=[XC._grid['LON'].flatten().min(), XC._grid['LON'].flatten().max()]
	else:
		latlims=[XC.grid_size['lats'][0],XC.grid_size['lats'][-1]]
		lonlims=[XC.grid_size['lons'][0],XC.grid_size['lons'][-1]]
	map_ax.set_extent(lonlims+latlims, crs=ccrs.PlateCarree())
	map_ax.add_feature(cfeature.LAND,color='silver')
	map_ax.add_feature(cfeature.OCEAN,color='lightblue')
	map_ax.add_feature(cfeature.COASTLINE)
	map_ax.add_feature(cfeature.BORDERS, linestyle=':')
	map_ax.add_feature(cfeature.LAKES, color='lightblue', alpha=0.5)
	map_ax.add_feature(cfeature.RIVERS, color='lightblue')
	
	# add lat/lon labels
	gl = map_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='white', alpha=0.8, linestyle=':')
	gl.xformatter = LongitudeFormatter(zero_direction_label=True)
	gl.yformatter = LatitudeFormatter()
	meridians = np.linspace(lonlims[0],lonlims[-1],4)
	parallels = np.linspace(latlims[0],latlims[-1],4)
	gl.xlocator = mticker.FixedLocator(meridians)
	gl.ylocator = mticker.FixedLocator(parallels)
	gl.bottom_labels = False
	gl.right_labels = False
	gl.xlabel_style = {'size': 6}
	gl.ylabel_style = {'size': 6}

	# add stations to map
	stalats=np.array([tr.stats.coordinates.latitude for tr in st_sort])
	stalons=np.array([tr.stats.coordinates.longitude for tr in st_sort])
	stanames=[tr.stats.station for tr in st_sort]
	for i,name in enumerate(stanames):
		h_mapdot.append(
						map_ax.plot(stalons[i],stalats[i],transform=ccrs.PlateCarree(),
								marker='^', color=st_dot_color0, markeredgecolor='k',
								markersize=8,markeredgewidth=0.5)
						)
		h_maptext.append(
						map_ax.text(stalons[i],stalats[i],name,
									fontsize=9,color=st_text_color0)
						)

	# add location(s) to the map
	rand1=(np.random.random_sample((len(loc.bstrap_lon),))-0.5)/50.
	rand2=(np.random.random_sample((len(loc.bstrap_lon),))-0.5)/50.
	map_ax.plot(loc.bstrap_lon+rand1,loc.bstrap_lat+rand2,'.',color='lightgreen',markeredgecolor='k',markersize=2,markeredgewidth=0.5)
	map_ax.plot(loc.longitude,loc.latitude,'*',markerfacecolor='firebrick',markeredgecolor='k',markersize=10,markeredgewidth=0.5)
	
	if np.isnan(loc.horizontal_scatter):
		plt.title('Lat: {:.3f}, Lon: {:.3f}\nDepth {:.1f} km'.format(loc.latitude,
																	 loc.longitude,
																	 loc.depth),fontsize=8)
	else:
		plt.title('Lat: {:.3f}, Lon: {:.3f} (scatter: {:.1f} km)\nDepth {:.1f} km (scatter {:.1f} km)'.format(loc.latitude,
																											  loc.longitude,
																											  loc.horizontal_scatter,
																											  loc.depth,
																											  loc.vertical_scatter),fontsize=8)
	########################################
	########################################



	############ plot envelopes ############
	########################################
	trace_ax = plt.subplot(gs[30:,:9])
	for i,tr in enumerate(st_sort):
		T=np.linspace(mdates.date2num(tr.stats.starttime.datetime),mdates.date2num(tr.stats.endtime.datetime),len(tr.data))
		h_trace.append(plt.plot(mdates.num2date(T),i+tr.data/tr.data.max(),color=trace_color0,linewidth=0.2))
		plt.text(mdates.date2num(tr.stats.endtime.datetime),i,' {:.1f}s, {:.1f}'.format(tr.dTs,tr.staW/len(CC['st'])),fontsize=6)
	plt.text(mdates.date2num(tr.stats.endtime.datetime),i+1,' time, weight',fontsize=5)
	trace_ax.set_yticks(np.arange(23))
	trace_ax.set_yticklabels([tr.stats.station +'.'+tr.stats.channel for tr in st_sort])

	plt.xticks(np.linspace(mdates.date2num(tr.stats.starttime.datetime),mdates.date2num(tr.stats.endtime.datetime),4))
	plt.tick_params(axis='both', which='major', labelsize=6)
	plt.xlim(T[0],T[-1])
	plt.ylim(-2,len(CC['st'])+0.25)
	xfmt = mdates.DateFormatter('%H:%M:%S')
	trace_ax.xaxis.set_major_formatter(xfmt)
	plt.grid(linestyle='--',linewidth=0.2,zorder=0,axis='x')
	plt.title('Envelopes starting at\n{}'.format(CC['st'][0].stats.starttime.strftime('%Y.%m.%d %H:%M:%S')),fontsize=8)
	########################################
	########################################



	########## plot correlations ###########
	########################################
	cc_ax = plt.subplot(gs[30:,11:])
	h_5=list()
	for i in range(len(CC['st'])):
		k0=sort_key[i]
		tmp_list=list()
		ii=np.where(CC['j1']==k0)[0]
		for sub_ind in ii:
			if CC['W2'][sub_ind]>0:
				DT0 = CC['tC'][np.round(deltaT[sub_ind]).astype('int')]
				tmp_list.append(plt.plot(CC['tC']-DT0,np.abs(DT0)+(1-CC1[:,sub_ind])*10,'-',color=xcorr_color0,linewidth=0.5,alpha=0.3))
		ii=np.where(CC['i1']==k0)[0]
		for sub_ind in ii:
			if CC['W2'][sub_ind]>0:
				DT0 = CC['tC'][np.round(deltaT[sub_ind]).astype('int')]
				tmp_list.append(plt.plot(CC['tC']-DT0,np.abs(DT0)+(1-CC1[:,sub_ind])*10,'-',color=xcorr_color0,linewidth=0.5,alpha=0.3))
		h_5.append(tmp_list)

	plt.grid(linestyle='--',linewidth=0.2,zorder=0)
	plt.tick_params(axis='both', which='major', labelsize=6)
	plt.yticks([],[])
	plt.title('Time shifted cross-correlations (s)',fontsize=8)
	########################################
	########################################


	msft = np.reshape(misfit,np.shape(XC._grid['LON']),order='F')
	ii=np.unravel_index(msft.argmin(),np.shape(XC._grid['LON']),order='C')

	########### plot misfit map ############
	########################################
	misfit_ax1 = plt.subplot(gs[:11,14:], projection=ccrs.PlateCarree())
	misfit_ax1.set_extent(lonlims+latlims, crs=ccrs.PlateCarree())
	misfit_ax1.add_feature(cfeature.COASTLINE)
	misfit_ax1.contourf(XC._grid['LON'][:,:,ii[2]], 
						XC._grid['LAT'][:,:,ii[2]], 
						msft[:,:,ii[2]],40,
                		transform=ccrs.PlateCarree(),
                		cmap='viridis_r')
	misfit_ax1.plot(loc.longitude,loc.latitude,marker='*',color='firebrick',markeredgewidth=0.1,markeredgecolor='firebrick',markersize=6)
	misfit_ax1.plot(stalons,stalats,marker='^',color='gray',markeredgewidth=0.1,markeredgecolor='gray',markersize=4,linewidth=0)
	gl2 = misfit_ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=0.2, color='white', alpha=0.8, linestyle=':')
	gl2.xlocator = mticker.FixedLocator(meridians)
	gl2.ylocator = mticker.FixedLocator(parallels)
	gl2.yformatter = LatitudeFormatter()
	gl2.bottom_labels = False
	gl2.xlabels_top = False
	gl2.ylabels_left = False
	gl2.ylabel_style = {'size': 6}
	plt.title('Misfit Function',fontsize=8)
	########################################
	########################################


	########## plot misfit slice ###########
	########################################
	misfit_ax2 = plt.subplot(gs[12:24,14:])
	map_aspect = (meridians[-1]-meridians[0])/((parallels[-1]-parallels[0])*111.1)
	if map_aspect > 1:
		misfit_ax2.set_aspect((meridians[-1]-meridians[0])/((parallels[-1]-parallels[0])*111.1))
	plt.contourf(XC._grid['LON'][ii[0],:,:], -XC._grid['DEP'][ii[0],:,:], msft[ii[0],:,:],40,cmap='viridis_r')
	plt.plot(loc.longitude,-loc.depth,'*',markerfacecolor='firebrick',markeredgecolor='firebrick',markeredgewidth=0.1)
	misfit_ax2.yaxis.tick_right()
	plt.ylabel('Depth',fontsize=8)
	misfit_ax2.yaxis.set_label_position('right')
	plt.gca().set_xticks(meridians[1:3])
	plt.xlim(lonlims)
	misfit_ax2.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))
	plt.xlabel('Longitude',fontsize=8)
	plt.tick_params(axis='y', which='major', labelsize=6)
	plt.tick_params(axis='x', which='major', labelsize=5)
	plt.grid(linestyle=':',linewidth=0.2,color='w')
	########################################
	########################################

	# fig.subplots_adjust(left=0.1,bottom=0.04,right=0.91,top=0.95,wspace=0.0,hspace=0.0)
	# a1=misfit_ax1.get_position()
	# a2=misfit_ax2.get_position()
	# pos = [a1.bounds[0],a2.bounds[1],a1.width,a1.height]
	# misfit_ax2.set_position(pos)


	if XC.interact:
		callback = Index()
		ax_done = plt.axes([0.13, 0.04, 0.06, 0.02])
		b_done = Button(ax_done, label='Done',color='yellow',hovercolor='lightskyblue')
		b_done.label.set_fontsize(6)
		b_done.on_clicked(callback.done)

		ax_reloc = plt.axes([0.24, 0.04, 0.08, 0.02])
		b_reloc = Button(ax_reloc, label='Relocate',color='yellow',hovercolor='lightskyblue')
		b_reloc.label.set_fontsize(6)
		b_reloc.on_clicked(callback.reloc)

		ax_restart = plt.axes([0.365, 0.04, 0.08, 0.02])
		b_restart = Button(ax_restart, label='Restart',color='yellow',hovercolor='lightskyblue')
		b_restart.label.set_fontsize(6)
		b_restart.on_clicked(callback.restart)

		cid=fig.canvas.mpl_connect('button_press_event',on_click)
		plt.show()
		plt.pause(0.5)
		fig.canvas.mpl_disconnect(cid)
		for i,ind in enumerate(krm):
			krm[i]=sort_key[ind]
	else:
		plt.show()

	plot_opt['krm']=np.array(krm)
	plot_opt['relocate']=relocate
	plot_opt['restart']=restart

	return plot_opt


def grid_plot(XC):
	from obspy.geodetics.base import gps2dist_azimuth

	############### plot map ###############
	########################################

	fig = plt.figure(figsize=(6,7.2))

	map_ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
	if XC.rotation:
		latlims=[XC._grid['LAT'].flatten().min(), XC._grid['LAT'].flatten().max()]
		lonlims=[XC._grid['LON'].flatten().min(), XC._grid['LON'].flatten().max()]
	else:
		latlims=[XC.grid_size['lats'][0],XC.grid_size['lats'][-1]]
		lonlims=[XC.grid_size['lons'][0],XC.grid_size['lons'][-1]]
	map_ax.set_extent(lonlims+latlims, crs=ccrs.PlateCarree())
	map_ax.add_feature(cfeature.LAND,color='silver')
	map_ax.add_feature(cfeature.OCEAN,color='lightblue')
	map_ax.add_feature(cfeature.COASTLINE)
	map_ax.add_feature(cfeature.BORDERS, linestyle=':')
	map_ax.add_feature(cfeature.LAKES, color='lightblue', alpha=0.5)
	map_ax.add_feature(cfeature.RIVERS, color='lightblue')
	
	# add lat/lon labels
	gl = map_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='white', alpha=0.8, linestyle=':')
	gl.xformatter = LongitudeFormatter(zero_direction_label=True)
	gl.yformatter = LatitudeFormatter()
	meridians = np.linspace(lonlims[0],lonlims[-1],4)
	parallels = np.linspace(latlims[0],latlims[-1],4)
	gl.xlocator = mticker.FixedLocator(meridians)
	gl.ylocator = mticker.FixedLocator(parallels)
	gl.bottom_labels = False
	gl.right_labels = False
	gl.xlabel_style = {'size': 6}
	gl.ylabel_style = {'size': 6}

	# add grid nodes
	map_ax.plot(XC._grid['LON'][:,:,0].flatten(),XC._grid['LAT'][:,:,0].flatten(),transform=ccrs.PlateCarree(), marker='+', color='g', markersize=3,linewidth=0,)

	# add stations to map
	stalats=np.array([tr.stats.coordinates.latitude for tr in XC.traces])
	stalons=np.array([tr.stats.coordinates.longitude for tr in XC.traces])
	stanames=[tr.stats.station for tr in XC.traces]
	for i,name in enumerate(stanames):
		map_ax.plot(stalons[i],stalats[i],transform=ccrs.PlateCarree(),
								marker='^', color=st_dot_color0, markeredgecolor='k',
								markersize=8,markeredgewidth=0.5)
		map_ax.text(stalons[i],stalats[i],name,
									fontsize=9,color=st_text_color0)

	if XC.rotation:
		dx=XC.rotation['x'][1]-XC.rotation['x'][0]
		dy=XC.rotation['y'][1]-XC.rotation['y'][0]
		dz=XC.rotation['z'][1]-XC.rotation['z'][0]
		plt.title('Stations and Grid Nodes\ndx={:.1f} km, dy={:.1f} km, dz={:.1f} km\nDepth range: {:.1f} - {:.1f} km'.format(dx,dy,dz,XC.rotation['z'][0],XC.rotation['z'][-1]),fontsize=8)
	else:
		dx=gps2dist_azimuth(XC.grid_size['lats'][0],XC.grid_size['lons'][0],
							XC.grid_size['lats'][0],XC.grid_size['lons'][1])[0]/1000
		dy=gps2dist_azimuth(XC.grid_size['lats'][0],XC.grid_size['lons'][0],
							XC.grid_size['lats'][1],XC.grid_size['lons'][0])[0]/1000
		dz=XC.grid_size['deps'][1]-XC.grid_size['deps'][0]
		plt.title('Stations and Grid Nodes\ndx={:.2f} km, dy={:.2f} km, dz={:.2f} km\nDepth range {:.1f} - {:.1f} km'.format(dx,dy,dz,XC.grid_size['deps'][0],XC.grid_size['deps'][-1]),fontsize=8)

	plt.show()
	
	########################################
	########################################
	
	return


def plot_locations(locs,XC):

	gs = gridspec.GridSpec(4,1)
	fig = plt.figure(figsize=(6,6))

	map_ax = fig.add_subplot(gs[:3,0], projection=ccrs.PlateCarree())
	if XC.rotation:
		latlims=[XC._grid['LAT'].flatten().min(), XC._grid['LAT'].flatten().max()]
		lonlims=[XC._grid['LON'].flatten().min(), XC._grid['LON'].flatten().max()]
	else:
		latlims=[XC.grid_size['lats'][0],XC.grid_size['lats'][-1]]
		lonlims=[XC.grid_size['lons'][0],XC.grid_size['lons'][-1]]
	map_ax.set_extent(lonlims+latlims, crs=ccrs.PlateCarree())
	map_ax.add_feature(cfeature.LAND,color='silver')
	map_ax.add_feature(cfeature.OCEAN,color='lightblue')
	map_ax.add_feature(cfeature.COASTLINE)
	map_ax.add_feature(cfeature.BORDERS, linestyle=':')
	map_ax.add_feature(cfeature.LAKES, color='lightblue', alpha=0.5)
	map_ax.add_feature(cfeature.RIVERS, color='lightblue')
	
	# add lat/lon labels
	gl = map_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='white', alpha=0.8, linestyle=':')
	gl.xformatter = LongitudeFormatter(zero_direction_label=True)
	gl.yformatter = LatitudeFormatter()
	meridians = np.linspace(lonlims[0],lonlims[-1],4)
	parallels = np.linspace(latlims[0],latlims[-1],4)
	gl.xlocator = mticker.FixedLocator(meridians)
	gl.ylocator = mticker.FixedLocator(parallels)
	gl.bottom_labels = False
	gl.right_labels = False
	gl.xlabel_style = {'size': 6}
	gl.ylabel_style = {'size': 6}

	stalats=np.array([tr.stats.coordinates.latitude for tr in XC.traces])
	stalons=np.array([tr.stats.coordinates.longitude for tr in XC.traces])
	stanames=[tr.stats.station for tr in XC.traces]

	hist_ax = fig.add_subplot(gs[-1,0])
	hist_ax.set_xlim(-1,len(XC.traces))
	hist_ax.set_title('Station Usage',fontsize=10)

	for i,name in enumerate(stanames):
		map_ax.plot(stalons[i],stalats[i],transform=ccrs.PlateCarree(),
							marker='^', color=st_dot_color0, markeredgecolor='k',
							markersize=8,markeredgewidth=0.5)
		map_ax.text(stalons[i],stalats[i],name,
								fontsize=6,color=st_text_color0)

	if 'event_list' in str(locs.__class__):
		map_ax.scatter(locs.get_lons(),locs.get_lats(),transform=ccrs.PlateCarree(),
			marker='o',s=18,color='red',edgecolors='k',linewidths=0.5)
		stas=locs.get_stations()
		sta_list=dict(Counter(stas).most_common())
		for tr in XC.traces:
			if tr.id not in sta_list.keys():
				sta_list[tr.id] = 0
		plt.bar(sta_list.keys(), sta_list.values(), width=0.5, color='gray')
	
	elif 'detections' in str(locs.__class__):
		# 1st plot noise
		h1=map_ax.scatter(locs.noise.get_lons(),locs.noise.get_lats(),transform=ccrs.PlateCarree(),
			marker='o',s=8,color='lightsteelblue',edgecolors='k',linewidths=0.5)
		# 2nd plot edge
		h2=map_ax.scatter(locs.edge_clustered.get_lons(),locs.edge_clustered.get_lats(),transform=ccrs.PlateCarree(),
			marker='o',s=12,color='darkorange',edgecolors='k',linewidths=0.5)
		# 2nd plot edge
		h3=map_ax.scatter(locs.core_clustered.get_lons(),locs.core_clustered.get_lats(),transform=ccrs.PlateCarree(),
			marker='o',s=18,color='red',edgecolors='k',linewidths=0.5)

		map_ax.legend([h1,h2,h3],tuple(['noise','edge','core']),loc='lower left',ncol=1,fontsize=8)

		stas=locs.detections.get_stations()
		sta_list=dict(Counter(stas).most_common())
		for tr in XC.traces:
			if tr.id not in sta_list.keys():
				sta_list[tr.id] = 0
		
		b1=hist_ax.bar(sta_list.keys(), sta_list.values(), width=0.5, color='lightsteelblue')
		
		stas=locs.core_clustered.get_stations()
		sta_list=dict(Counter(stas).most_common())
		for tr in XC.traces:
			if tr.id not in sta_list.keys():
				sta_list[tr.id] = 0
		b2=hist_ax.bar(sta_list.keys(), sta_list.values(), width=0.5, color='red')

		stas=locs.edge_clustered.get_stations()
		sta_list=dict(Counter(stas).most_common())
		for tr in XC.traces:
			if tr.id not in sta_list.keys():
				sta_list[tr.id] = 0
		b3=hist_ax.bar(sta_list.keys(), sta_list.values(), width=0.5, color='darkorange')

		hist_ax.legend([b1,b3,b2],tuple(['noise','edge','core']),loc='upper right',ncol=1,fontsize=6)
		hist_ax.set_ylabel('counts')


	plt.xticks(rotation=45,fontsize=6,ha='right')

	plt.show()

	return

