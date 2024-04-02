enveloc
=======

A Python package to locate emergent seismicity using network of seimsic stations. 
*enveloc* uses envelope cross correlation of seismic traces on a fixed time window and performs
a 3D grid search to maximize signal coherency. The location process can be performed interactively
with on single short time window, or more automatically on longer data sets by internally cutting
data and processing smaller time windows. The latter approach can also optionally be parallelized 
internally to shorten processing time.

Code is stored on |github| and the latest stable release can be found |releases_link| and 
|pypi|.


.. |github| raw:: html

    <a href="https://github.com/awech/enveloc" target="_blank">github</a>

.. |releases_link| raw:: html

	<a href="https://github.com/awech/enveloc/releases" target="_blank">here</a>

.. |pypi| raw:: html

    <a href="https://pypi.org/project/enveloc/" target="_blank">PyPI</a>

*enveloc* uses |Obspy_link| for handling and manipulating seismic data and metadata.

.. |Obspy_link| raw:: html

	<a href="https://docs.obspy.org/" target="_blank">Obspy</a>

This package is written and maintained by Aaron Wech, and is distributed under the 
GNU General Public Licence v3.

Citation
--------
If you use this package in your work, please cite the following |paper|:
Wech, A.G., and K.C. Creager (2008), Automatic detection and location of Cascadia tremor, Geophys. Res. Lett., 35, L20302.

.. |paper| raw:: html

	<a href="https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL035458" target="_blank">paper</a>


Contents:
---------

.. toctree::
    :numbered:
    :maxdepth: 1

    intro
    installation
    setup
    tutorial
    output
    xcor
    event_list
    detections