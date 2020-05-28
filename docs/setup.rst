Usage overview
==============

*enveloc* uses |Obspy_link| for handling and manipulating seismic data and metadata. Before using *enveloc*
I strongly recommend familiarizing yourself with Obspy.

At its most basic level *enveloc* needs 3 things, ordered by importance:

 #. Seismic data
 #. A velocity model
 #. A 3D grid

.. |Obspy_link| raw:: html

	<a href="https://docs.obspy.org/" target="_blank">Obspy</a>


.. _seismic data:

Seismic Data
------------
Seismic data is input as an |stream|. 

.. note::
	Data must be pre-processed into envelopes **externally** before being input to *enveloc*. There were too 
	many unknowns to standardize this within *enveloc*. Each application has different needs, and each
	user has preferred methods of filtering, etc.

Assuming you have a Stream of data, *st*, read from *IRIS* or *mseed* or wherever, here is an example
pre-processing from raw data into envelope data:


.. code-block:: python

	from obspy.signal.filter import envelope

	DT = 10       # seconds of padding
	FREQMIN = 1.0 # minimum frequency
	FREQMAX = 8.0 # maximum frequency
	LOWPASS = 0.2 # lowpass for envelope smoothing

	# filter raw waveforms
	st.detrend('demean')
	st.taper(max_percentage=None,max_length=5) # taper 5 seconds off each end
	st.filter('bandpass',freqmin=FREQMIN,freqmax=FREQMAX,corners=3,zerophase=True)

	# convert waveforms envelopes
	st.detrend('demean')
	for tr in st:
		# downsample for efficiency
		tr.resample(25.0)
		# convert to even number of samples for speed
		if tr.stats.npts % 2 == 1:
			tr.trim(starttime=tr.stats.starttime,endtime=tr.stats.endtime+1/tr.stats.sampling_rate,pad=True,fill_value=0)
		tr.data = envelope(tr.data)
		# downsample envelope
		tr.resample(5.0)		# downsample to 5 sps

	# smooth the envelope
	st.filter('lowpass',freq=LOWPASS)

	# trim to remove taper effects
	st.trim(st[0].stats.starttime+DT,st[0].stats.endtime-DT)

Data can either be a small window surrounding a signal of interest or a longer timeseries to be internally
processed as many smaller windows.

.. note::
	#. If you are trying to locate a single window, make sure to trim off any edge effects, such as a pre-filter taper,
	   prior to inputing the Stream into *enveloc*. If all the envelope traces taper to zero on either side, the signal
	   will look like a box function and things will correlate well and with zero lag, resulting in erroneous locations
	   in the center of the network.
	#. Be careful not to make the window length too small. In general the window length should be longer than the largest
	   predicted differential travel time between all station pairs for all grid nodes. That way you can appropriately
	   shift by at least that much during the cross correlation. Internal checks try and force this to be the case, but
	   I've never been fully satisfied with their implementation.

Each |trace|, tr, within the Stream object must also contain its associated latitude and longitude (in decimal degrees)
in the form of an Obspy AttribDict named '*coordinates*', which itself is a key for Trace's AttribDict '*stats*'. This could 
be extracted from an *inventory*:

.. code-block:: python

    from obspy.clients.fdsn import Client
    from obspy.core.util import AttribDict
    
    client = Client('IRIS')
    for tr in st:
        inventory = client.get_stations( network   = tr.stats.network,
                                         station   = tr.stats.station,
                                         location  = tr.stats.location,
                                         channel   = tr.stats.channel,
                                         starttime = tr.stats.starttime,
                                         endtime   = tr.stats.endtime )

        tr.stats.coordinates = AttribDict(
                                           { 'latitude'  : inventory[0][0].latitude,
                                             'longitude' : inventory[0][0].longitude,
                                             'elevation' : inventory[0][0].elevation }
                                         )

or it could be input manually. For example:

.. code-block:: python

    from obspy.core.util import AttribDict

    for tr in st:
        tr.stats.coordinates = AttribDict( 
                                           {'latitude'  : <value>,
                                            'longitude' : <value>,
                                            'elevation' : <value> } 
                                         )

.. |stream| raw:: html

	<a href="https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html" target="_blank">Obspy Stream</a>

.. |trace| raw:: html

	<a href="https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.html" target="_blank">Trace</a>

.. _velocity model:

Velocity Model
--------------
The velocity model is used to calculate traveltimes from each grid node to each station. Strictly speaking a 
velocity model is not required as input in order for *enveloc* to run. If no model is provided, *enveloc* will run 
using a default model, which is taken from *default_vel_model.tvel*:

.. code-block:: bash

	default - P
	default - S
	   0.000	5.1572		2.9775		2.7200
	   0.050	5.1572		2.9775		2.7200
	   0.050	5.1573		2.9773		2.7200
	   4.000	5.1573		2.9773		2.7200
	   4.000	5.4491		3.1461		2.7200
	  10.000	5.4491		3.1461		2.7200
	  10.000	6.0330		3.4831		2.7200
	  15.000	6.0330		3.4831		2.7200
	  15.000	6.7141		3.8764		2.7200
	  20.000	6.7141		3.8764		2.7200
	  20.000	7.2007		4.1573		2.9000
	  25.000	7.2007		4.1573		2.9000
	  25.000	7.4926		4.3258		2.9000
	  33.000	7.4926		4.3258		2.9000
	  33.000	7.6872		4.4382		3.3000
	  47.000	7.6872		4.4382		3.3000
	  47.000	7.8818		4.5506		3.3000
	  65.000	7.8818		4.5506		3.3000
	  65.000	8.0764		4.6629		3.3000
	  65.000	8.0764		4.6629		3.3000
	 210.000	8.3000		4.5180		3.4258
	 210.000	8.3000		4.5220		3.4258
	 260.000	8.4825		4.6090		3.4561
	6371.000	8.2000		4.7000		3.3198

The file has the structure:

.. code-block:: bash

	comment line - generally info about the P velocity model
	comment line - generally info about the S velocity model
	depth pVel sVel Density
	depth pVel sVel Density

where velocities are assumed to be linear between sample points. Because this type of model file doesn’t give 
complete information Obspy makes the following assumptions:

* modelname - from the filename, with ”.tvel” dropped if present
* radius_of_planet - the largest depth in the model
* meanDensity - 5517.0
* G - 6.67e-11

Alternatively, if you want to search using a constant velocity, like a surface wave velocity, you can
provide this constant as a variable to *enveloc*. In this case the default velocity model is used, but 
only for ray-tracing purposes, and the actual velocity and resulting traveltimes are superseded by the 
constant velocity provided. More details provided in :class:`enveloc.XCOR` class (:ref:`xcor class`).

.. _grid section:

Grid
----
Providing a grid as input is optional, but strongly recommended. If no grid is provided, a 
default grid is calculated based on the lat/lon extent of the input seismic stations. By 
default, an evenly spaced 20 x 25 (lat x lon) grid is produced extending 33% beyond both 
the minimum and maximum latitudinal and longitudinal footprint of the stations provided. 
The vertical grid is 10 evenly spaced nodes spanning from 0 km to a maximum depth based on the
horizontal grid extent.

A grid can be provided as an input argument, *grid_size*, a dictionary with keys 
'*lats*', '*lons*', and '*deps*', each of which are 1D monotonic numpy arrays. For example:

.. code-block:: python

	import numpy as np

	grid_size = { 
	              'deps' : np.arange(0.5,30,2),           # in km
	              'lons' : np.arange(-158.4,-157.9, 0.04),
	              'lats' : np.arange(  56.7,  57.1, 0.04)
	            }

Alternatively, *enveloc* can perform a grid search over a rotated grid. For this, you must provide
an input argument, *rotation*, a dictionary with keys '*x*' (x grid nodes, in km), '*y*' (y grid 
nodes, in km), '*z*' (depth grid nodes, in km), '*lat0*' & '*lon0*' (origin lat/lon), '*az*', rotation 
azimuth (degrees counterclockwise from East). For example:

.. code-block:: python

	import numpy as np
	
	rotation = {
	             'x'    : np.arange(-20,20,1), # in km
	             'y'    : np.arange(-10,10,1), # in km
	             'z'    : np.arange(1,39,2),   # in km
	             'lat0' : 56.919,              # origin latitude
	             'lon0' : -158.1737,           # origin longitude
	             'az'   : 30                   # rotation, in degrees, counterclockwise from East
	           }

This example creates a grid rotated 30 degrees north of east about the origin *lat0, lon0* with 
grid nodes every 1 km extending +/- 20 and +/- 10 km in the rotated x and y directions, respectively.
*grid_size* and *rotation* should not both be used.