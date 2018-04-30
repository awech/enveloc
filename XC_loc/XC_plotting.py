import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.dates as mdates
from matplotlib.widgets import Button
import numpy as np


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

interval = 0.2

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
	map_ax = fig.add_subplot(gs[:25,:13])
	XC.map.drawcoastlines()
	XC.map.drawmapboundary(fill_color='lightblue')
	XC.map.fillcontinents(color='gray',lake_color='lightblue')
	# add lat/lon labels
	parallels=np.linspace(XC.map.boundarylats.min(),XC.map.boundarylats.max(),4)
	meridians=np.linspace(XC.map.boundarylons.min(),XC.map.boundarylons.max(),4)
	XC.map.drawparallels(parallels[1:3],color='w',linewidth=0.5,dashes=[1,4],labels=[True,False,False,True],fontsize=6,labelstyle='+/-',fmt='%.3f')
	XC.map.drawmeridians(meridians[1:3],color='w',linewidth=0.5,dashes=[1,4],labels=[True,False,False,True],fontsize=6,labelstyle='+/-',fmt='%.3f')
	# add stations to map
	stalats=np.array([tr.stats.coordinates.latitude for tr in st_sort])
	stalons=np.array([tr.stats.coordinates.longitude for tr in st_sort])
	stanames=[tr.stats.station for tr in st_sort]
	for i,name in enumerate(stanames):
		x,y=XC.map(stalons[i],stalats[i])
		h_maptext.append(plt.text(x,y,name,fontsize=9,color=st_text_color0))
		h_mapdot.append(XC.map.plot(stalons[i],stalats[i],'^',color=st_dot_color0,latlon=True,markeredgecolor='k',markersize=8,markeredgewidth=0.5))
	# add location(s) to the map
	rand1=(np.random.random_sample((len(loc.bstrap_lon),))-0.5)/50.
	rand2=(np.random.random_sample((len(loc.bstrap_lon),))-0.5)/50.
	XC.map.plot(loc.bstrap_lon+rand1,loc.bstrap_lat+rand2,'.',color='lightgreen',latlon=True,markeredgecolor='k',markersize=2,markeredgewidth=0.5)
	XC.map.plot(loc.longitude,loc.latitude,'g*',latlon=True,markeredgecolor='k',markersize=10,markeredgewidth=0.5)
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
	trace_ax = fig.add_subplot(gs[30:,:9])
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
	cc_ax = fig.add_subplot(gs[30:,11:])
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
	misfit_ax1 = fig.add_subplot(gs[:11,14:])
	XC.map.drawcoastlines()
	XC.map.contourf(XC._grid['LON'][:,:,ii[2]], XC._grid['LAT'][:,:,ii[2]], msft[:,:,ii[2]],40,latlon=True)
	XC.map.plot(XC._grid['LON'][ii[0],ii[1],ii[2]],-XC._grid['DEP'][ii[0],ii[1],ii[2]],'rs',latlon=True)
	XC.map.drawparallels(parallels[1:3],color='w',linewidth=0.2,dashes=[1,4],labels=[False,True,True,False],fontsize=6,labelstyle='+/-',fmt='%.3f')
	plt.title('Misfit Function',fontsize=8)
	########################################
	########################################


	########## plot misfit slice ###########
	########################################
	misfit_ax2 = fig.add_subplot(gs[12:24,14:])
	map_aspect = (meridians[-1]-meridians[0])/((parallels[-1]-parallels[0])*111.1)
	if map_aspect > 1:
		misfit_ax2.set_aspect((meridians[-1]-meridians[0])/((parallels[-1]-parallels[0])*111.1))
	plt.contourf(XC._grid['LON'][ii[0],:,:], -XC._grid['DEP'][ii[0],:,:], msft[ii[0],:,:],40)
	misfit_ax2.yaxis.tick_right()
	plt.ylabel('Depth',fontsize=8)
	misfit_ax2.yaxis.set_label_position('right')
	plt.xlabel('Longitude',fontsize=8)
	plt.tick_params(axis='both', which='major', labelsize=6)
	########################################
	########################################

	fig.subplots_adjust(left=0.1,bottom=0.04,right=0.91,top=0.95,wspace=0.0,hspace=0.0)
	a1=misfit_ax1.get_position()
	a2=misfit_ax2.get_position()
	pos = [a1.bounds[0],a2.bounds[1],a1.width,a1.height]
	misfit_ax2.set_position(pos)


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




