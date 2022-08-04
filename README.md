# PX915 - Control of Thin Liquid Films

Code for modelling thin liquid films controlled by baseplate actuators. Written for the PX915 individual project by [Oscar Holroyd](https://warwick.ac.uk/fac/sci/hetsys/people/studentscohort3/holroyd/), supervised by [Radu Cimpeanu](https://warwick.ac.uk/fac/sci/maths/people/staff/cimpeanu/) and [Susana Gomes](https://warwick.ac.uk/fac/sci/maths/people/staff/gomes).


## Installation
* The code relies on [Basilisk](<http://basilisk.fr/>) to model the Navier-Stokes equations. See the [installation page](<http://basilisk.fr/src/INSTALL>) for instructions. Note that bview is *not* required.
* The LQR controls require a functioning [LAPACKE](https://netlib.org/lapack/lapacke.html) implementation.

Once these are installed, an executable can be generated:
```bash
git clone https://github.com/OaHolroyd/px915-basilisk
cd px915-basilisk
make source
make
```


## Running the code

For a given set of parameters, the code will print the total cost

$$ \kappa = \int_0^T \int_0^L \mu (h(x)-1)^2 + (1-\mu) F^2 \text{d}x \text{d}t $$

to `stdout`. All other outputs (dimensionless numbers, simulation progress etc.) are printed to `stderr`.

### Requirements
The code requires an input JSON file containing physical parameters, solver settings and output details. An example of such a file (containing parameters corresponding to a thin liquid water film) is included in [params.json](params.json). If you wish to use an alternative file, pass its' path as the first (and only) argument to `film`.

The code outputs 1D and 2D data under the directory 'out', and full variable outputs under 'dump'. These directories must exist for the code to run (and are created by `make`).

### Visualisation
The output data is stored as plain text, timestamped on the first line and subsequently in a format easily read by [gnuplot](http://www.gnuplot.info/). The basic visualisation script [gen_plots.py](gen_plots.py) is included for quick analysis.

### Restarting
Basilisk includes the option to 'dump' the entire simulation to a single file, which can then be restored to continue the simulation from the output time. By default this occurs every 100 time-units. Since stable travelling waves take a long time to develop for most parameter regimes it is *strongly suggested* that a single run without controls is performed to generate a dump-file with a travelling wave before loading it and beginning controls after this point. To do this set the value `"t0"` to correspond to the file at dump/dump-\<time\>.


## Input parameters
All of the parameters are either in SI units or dimensionless. The keys in [params.json](params.json) are hopefully fairly self-explanatory. However, below is a full description of what they all do.

#### Domain parameters
These 6 parameters describe the spacial and temporal domains:
* **`h0`** the thickness of a uniform film - m
* **`Lx`** the ratio of film thickness to domain length - dimensionless
* **`Ly`** the ratio of film thickness to domain height (including air layer) - dimensionless
* **`theta`** the angle of the plate from horizontal - rad
* **`tmax`** the simulation end time - dimensionless
* **`t0`** the start time (either 0 or matching a dump file) - dimensionless
