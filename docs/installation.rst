enveloc Installation
====================

*enveloc* is availabe on |github| and as a repository on |pypi|.

I strongly recommend installing *enveloc* using conda because this will:

    #. Simplify the install
    #. Separate your *enveloc* install from your system Python so things don't break

If you do not have either a miniconda or anaconda installation you can follow
the |conda-install| instructions. Once anaconda (or miniconda, my preference) is installed,
create a new conda environment with the following:

.. code-block:: bash

    conda create -n enveloc python=3.8 numpy obspy=1.2.2 cartopy scikit-learn utm

Next activate that environment by calling:

.. code-block:: bash

    source activate enveloc

This ensures your enveloc environment is active, so that  when you call pip, it will install packages
into the enveloc environment. Now install *enveloc*:

.. code-block:: bash

    pip install enveloc

To test the install:

.. code-block:: python

    from enveloc.example_utils import test

    test()

.. |github| raw:: html

    <a href="https://github.com/awech/enveloc" target="_blank">github</a>

.. |pypi| raw:: html

    <a href="https://pypi.org/project/enveloc/" target="_blank">PyPI</a>

.. |conda-install| raw:: html

    <a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">conda-install</a>


Python Version
--------------
The code has been tested on Python 3.7 and 3.8. It may work on other (possibly even 2.7),
but no promises.


Dependencies
------------
* numpy
* scipy
* matplotlib
* obspy
* proj
* cartopy
* scikit-learn
* utm


Notes
-----
Most packages are available through PyPI, but I have been unsuccessful in install proj (on which cartopy relies) 
directly through pip, hence the use of anaconda (that and it makes life easier anyway).