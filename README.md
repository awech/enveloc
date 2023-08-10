# enveloc
Python package to perform envelope cross-correlation location in seismology

### Dependencies
obspy<br>
cartopy<br>
scikit-learn<br>
utm<br>

### Citation
If you use this package in your work, please cite the following paper: Wech, A.G., and K.C. Creager (2008), Automatic detection and location of Cascadia tremor, Geophys. Res. Lett., 35, L20302.

### Installation
Use anaconda to create a new environment, and install Python v3.11:<br>

`conda create -n enveloc python=3.11`<br>

Next activate that environment:<br>

`source activate enveloc`<br>

Then install with pip:<br>

`pip install enveloc`<br>

Test install with:<br>
`from enveloc.example_utils import test`<br>
`test()`<br>

### Documentation

Go [here](https://enveloc.readthedocs.io/) for more information<br>
