# enveloc

A Python package to locate emergent seismicity using a network of seismic stations.
*enveloc* uses envelope cross correlation of seismic traces on a fixed time window and performs
a 3D grid search to maximize signal coherency. The location process can be performed interactively
with a single short time window, or more automatically on longer data sets by internally cutting
data and processing smaller time windows. The latter approach can also optionally be parallelized
internally to shorten processing time.

Code is stored on [github](https://github.com/awech/enveloc) and the latest stable release can be
found [here](https://github.com/awech/enveloc/releases) and on
[PyPI](https://pypi.org/project/enveloc/).

*enveloc* uses [Obspy](https://docs.obspy.org/) for handling and manipulating seismic data and metadata.

This package is written and maintained by Aaron Wech, and is distributed under the
GNU General Public Licence v3.

## Citation

If you use this package in your work, please cite the following
[paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL035458):

> Wech, A.G., and K.C. Creager (2008), Automatic detection and location of Cascadia tremor,
> Geophys. Res. Lett., 35, L20302.

## Contents

1. [Introduction](intro.md)
2. [Installation](installation.md)
3. [Usage overview](setup.md)
4. [Tutorial](tutorial.md)
5. [Output](output.md)
6. API Reference
    - [XCOR object](xcor.md)
    - [event_list object](event_list.md)
    - [detections object](detections.md)
