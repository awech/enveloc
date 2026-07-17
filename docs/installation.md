# enveloc Installation

*enveloc* is available on [github](https://github.com/awech/enveloc) and as a repository on
[PyPI](https://pypi.org/project/enveloc/).

I strongly recommend installing *enveloc* using conda because this will:

1. Simplify the install
2. Separate your *enveloc* install from your system Python so things don't break

If you do not have either a miniconda or anaconda installation you can follow
the [conda-install](https://docs.conda.io/en/latest/miniconda.html) instructions. Once anaconda
(or miniconda, my preference) is installed, create a new conda environment with Python v3.11 with
the following:

```bash
conda create -n enveloc python=3.11
```

Next activate that environment by calling:

```bash
source activate enveloc
```

This ensures your enveloc environment is active, so that when you call pip, it will install packages
into the enveloc environment. Now install *enveloc*:

```bash
pip install enveloc
```

To test the install:

```python
from enveloc.example_utils import test

test()
```

## Python Version

The code has been tested on Python 3.11.

## Dependencies

* obspy
* cartopy
* scikit-learn
* utm

## Notes

All packages should be available through PyPI. I had troubles with installing *cartopy* using pip with
versions of Python <3.9. If you can get *cartopy* working (e.g. via conda install), *enveloc* should
still work with earlier Python versions 3.7 & 3.8.
