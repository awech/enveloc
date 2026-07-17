# Tutorial

This section provides some basic examples for how to use **enveloc**. There are way too many knobs to
turn to go over all of the parameters and functionality in detail. The goal here isn't to make you a pro
user. This section is meant to help you understand how it works and if/how it could work for your
application. After gaining some comfort and familiarity, I recommend checking out more parameters in the
[XCOR class](xcor.md) reference page.

## First Steps

### Data input

Start by getting some seismic data. As mentioned in the [Seismic Data](setup.md#seismic-data) section,
data needs to be pre-processed before calling enveloc. Let's get a couple minutes of data from Kilauea
volcano during the onset of lava-lake tremor, an emergent signal:

<figure markdown="span">
  ![lava lake tremor](lava_lake_tremor.png){ width="800" }
</figure>

Then filter and convert to smoothed envelopes:

```python
from enveloc import example_utils

t1 = "2018-04-28 13:07"
t2 = "2018-04-28 13:09"

FREQMIN = 1.0
FREQMAX = 8.0
LOWPASS = 0.2

sta_list=[
           "HV.BYL..HHZ",
           "HV.DEVL..HHZ",
           "HV.HAT..HHZ",
           "HV.KKO..HHZ",
           "HV.NPT..HHZ",
           "HV.OBL..HHZ",
           "HV.PAUD..HHZ",
           "HV.PUHI..HHZ",
           "HV.RIMD..HHZ",
           "HV.SBL..HHZ",
           "HV.SDH..HHZ",
           "HV.UWB..HHZ",
           "HV.UWE..HHZ",
           "HV.WRM..HHZ",
         ]
env = example_utils.get_IRIS_data(sta_list, t1, t2, f1=FREQMIN, f2=FREQMAX,lowpass=LOWPASS)
```

### Creating *XCOR* object

With an Obspy stream of envelopes, `env`, in hand, the first step is to create an
[`XCOR`](xcor.md) object. `XCOR` is everything. It calculates the traveltime grid. It cross correlates
all the traces. It handles the grid search. It has a bunch of possible input parameters (see the
[XCOR class](xcor.md) reference), but for now we can create one with just the envelopes:

```python
from enveloc.core import XCOR
XC = XCOR(env)
```

Upon creation, `XCOR` organizes the traces for cross correlation and internally calculates traveltimes
to all stations, which are stored in the object `XC`. In the above case, where no additional input is
provided, a default grid is created (see the [Grid](setup.md#grid) section) and traveltimes are calculated
using the default *S*-wave velocity model (see the [Velocity Model](setup.md#velocity-model) section).

### Locate a signal

Now with the `XC` created, we can try and locate the signal using the built-in location method
`locate()` from `XCOR`:

```python
loc = XC.locate()
```

Alternatively, you could do all of the above steps in one test:

```python
from enveloc import example_utils
loc, XC = example_utils.interactive_example()
```

By default the code will attempt a single location and produce an interactive plot:

<figure markdown="span">
  ![enveloc interactive plot](enveloc_interact.png){ width="500" }
  <figcaption>
  Interactive plot produced by <em>enveloc</em>. Click on a trace to select. Click on trace again to
  de-select. Selected traces are highlighted, as is the station on the map and all associated cross
  correlograms. If '<em>Relocate</em>' is pressed, the code will attempt to relocate with the selected
  traces removed from the algorithm. '<em>Restart</em>' returns all original traces. '<em>Done</em>'
  exits interactive mode and closes the figure.
  </figcaption>
</figure>

!!! note
    The interaction part has only been lightly tested, and there may be possible bugs with the UI here.
    Make sure to disable any matplotlib backend.

    All the processing and location stuff (sans-interaction mode) are fairly well tested though.

See the [Single location](output.md#single-location) output section for details about the output variable
*loc* in this case.

## Grid Inputs

*enveloc* will automatically produce a grid if none is provided, but creating and inputing a grid is
strongly recommended. This can be done either as a lat/lon/depth or x/y/z grid. In both cases, the custom
grid is input as a dictionary with the relevant grid parameters.

### Custom lat/lon grid

Now we can try using the same data as above, but locating on a custom grid:

```python
import numpy as np

mygrid = {
            "deps": np.arange( 0, 14, 0.5),
            "lons": np.arange(-155.35, -155.2, 0.002),
            "lats": np.arange(19.35, 19.45, 0.002)
         }

XC = XCOR(env, grid_size=mygrid, interact=False)
loc = XC.locate()
```

In this case we input the custom grid as an argument to `XCOR`, and we turn off interactive mode with
`interact=False` (this latter step is unrelated, but an example of how one might do so).

### Custom rotated grid

Again, using the same data as above, we can try locating on a custom rotated grid:

```python
import numpy as np

my_rotation = {
                "x"    : np.arange(-5,5, 0.3),
                "y"    : np.arange(-3.5, 3.5, 0.3),
                "z"    : np.arange(0, 25, 2),
                "lat0" : 19.403,
                "lon0" : -155.281,
                "az"   : 30
              }

XC = XCOR(env, rotation=my_rotation, interact=False)
loc = XC.locate()
```

Where the `rotation` argument is set to the dictionary variable `my_rotation`.
You can view the grid by calling `plot_grid()`:

```python
XC.plot_grid()
```

which produces the following:

<figure markdown="span">
  ![rotated grid](rotated_grid.png){ width="500" }
  <figcaption>Grid plot produced by <em>enveloc</em>.</figcaption>
</figure>

It's OK if stations fall outside the grid.

### Regional Example

Here is an example locating tectonic tremor in the Pacific Northwest of the USA.

```python
from enveloc import example_utils
loc, XC = example_utils.cascadia_example()
```

In this example `XC` is created within `cascadia_example()` by the command:

```python
XC = XCOR(env, grid_size=mygrid, regrid=True, bootstrap=30, output=2)
```

which:

1. Uses a custom lat/lon grid
2. Regrids. After finding the minimum misfit grid-node, it relocates on a finer scale grid
   surrounding the original grid node. (not well-tested on rotated grid)
3. Bootstraps. It attempts 30 locations, throwing away a small percentage (Default=4%) of the
   correlations each time to create a cloud of scattered locations, which can be used to estimate location robustness.
4. Increases output. The `output` variable increases how much output is printed to the screen when
   *enveloc* runs. Higher integers (up to 4) means much chattier.

## Auto-locations

Rather than locate a single time window, *enveloc* is designed such that you can input a long time
series to try and locate many time windows.

### Making Windows

All of the above steps are the same: you pre-process the data beforehand, and hand *enveloc* envelopes,
a velocity model and a grid. However, now instead of supplying an Obspy Stream of envelopes spanning a
few minutes, you input hours or days of data.

```python
from enveloc import example_utils

t1 = "2020-05-24 00:00"
t2 = "2020-05-24 08:00"

FREQMIN = 1.5
FREQMAX = 6.0
LOWPASS = 0.1

sta_list = [
                "PB.B011.--.EHZ",
                "CN.SYMB.--.HHZ",
                "CN.PTRF.--.HHZ",
                "CN.VGZ.--.HHZ",
                "UW.JCW.--.EHZ",
                "PB.B003.--.EHZ",
                "PB.B006.--.EHZ",
                "PB.B001.--.EHZ",
                "PB.B013.--.EHZ",
                "UW.DOSE.--.HHZ",
                "UW.HDW.--.EHZ",
                "UW.GNW.--.HHZ",
                "UW.GMW.--.EHZ",
                "PB.B014.--.EHZ",
                "UW.SMW.--.EHZ",
                "UW.STOR.--.HHZ",
                "UW.TKEY.--.HHZ",
           ]

env = example_utils.get_IRIS_data(sta_list, t1, t2, f1=FREQMIN, f2=FREQMAX, lowpass=LOWPASS)
```

Now we can locate windows of length 300 seconds overlapping by 150 seconds:

```python
import numpy as np
from enveloc.core import XCOR

mygrid = {
            "lons": np.arange(-125, -121+0.05, 0.075),
            "lats": np.arange(46.5, 49.0+0.05, 0.075),
            "deps": np.arange(20, 60+0.1, 8)
         }


XC  = XCOR(env, bootstrap=20, plot=False, grid_size=mygrid, output=2)
locs = XC.locate(window_length=300, step=150)
```

!!! note
    Invoking auto-locations will override the `dTmax_s` set when the `XCOR` object was initialized. This
    is because when created, `XC` doesn't yet know that the data will be sliced up or what the window
    lengths will be in those slices. If `window_length` is passed to `locate()` (thus telling `locate()`
    that you are auto-locating on sub-windows), it will re-calculate a new default `dTmax_s` based on the
    `window_length` provided and station spacing. If you want to override/specify `dTmax_s` for
    auto-locations, you need to pass it as an argument in the `locate()` method
    (e.g., `locs = XC.locate(window_length=300, step=150, dTmax_s=25)`).

### Parallel Processing

Processing all these windows takes some time. That took about ~315 seconds on my 2018 MacBook Pro (the
location step, not the `XC = XCOR(...)` step that calculates traveltimes...more on that later). We can
speed that up by using multiple processors. You can do that by changing the number of processors used:

```python
XC = XCOR(
    env,
    bootstrap=20,
    plot=False,
    grid_size=mygrid,
    output=2,
    num_processors=4
)
locs = XC.locate(window_length=300, step=150)
```

Only looping over location windows is parallelized. Increasing to 4 processors reduces the location step
of 191 windows to ~55 seconds on my machine.

See the [Multiple locations](output.md#multiple-locations) output section for more discussion on the
output `locs`.

### Saving traveltimes

You can see that performing the traveltime calculation can take a long time, especially for high-density
grids. The above example took ~32 seconds on my machine, with a relatively small grid. More vertical grid
nodes make for longer computation times. If you are performing this step routinely, this can be a huge
time sink. For this reason, the `XCOR` object has the method `save_traveltimes()`, that allows you to save
the traveltimes as a compressed numpy *.npz* file to be loaded later.

```python
XC.save_traveltimes("example_tt_file.npz")
```

This file can then be loaded in much faster later when creating the `XC` object with the same station/grid
combo, which speeds things up considerably. The filename can be a full path and is passed via the argument
`tt_file`.

```python
XC = XCOR(
    env,
    bootstrap=20,
    plot=False,
    grid_size=mygrid,
    output=2,
    num_processors=4,
    tt_file="example_tt_file.npz"
)
```

Python clocked this at 0.03 seconds, which is a bit faster.

!!! note
    `XCOR` does try to check to make sure the input grid and stations match the data in the input *.npz*
    file of pre-calculated traveltimes. This has been lightly tested but you should be careful nonetheless.
