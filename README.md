# PX915 - Control of Thin Liquid Films

Code for modelling thin liquid films controlled by baseplate actuators. Written for the PX915 individual project by [Oscar Holroyd](https://warwick.ac.uk/fac/sci/hetsys/people/studentscohort3/holroyd/), supervised by [Radu Cimpeanu](https://warwick.ac.uk/fac/sci/maths/people/staff/cimpeanu/) and [Susana Gomes](https://warwick.ac.uk/fac/sci/maths/people/staff/gomes).

- [Installation](#installation)
- [Running the Code](#running-the-code)
- [Input Parameters](#input-parameters)
- [Mathematical Background](#mathematical-background)


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


## Running the Code

For a given set of parameters, the code will print the total cost
```math
\kappa = \int_0^T \int_0^L \mu (h(x)-1)^2 + (1-\mu) f^2 \text{d}x \text{d}t
```
to `stdout`. All other outputs (dimensionless numbers, simulation progress etc.) are printed to `stderr`.

### Requirements
The code requires an input JSON file containing physical parameters, solver settings and output details. An example of such a file (containing parameters corresponding to a thin liquid water film) is included in [params.json](params.json). If you wish to use an alternative file, pass its' path as the first (and only) argument to `film`.

The code outputs 1D and 2D data under the directory 'out', and full variable outputs under 'dump'. These directories must exist for the code to run (and are created by `make`).

### Visualisation
The output data is stored as plain text, timestamped on the first line and subsequently in a format easily read by [gnuplot](http://www.gnuplot.info/). The basic visualisation script [gen_plots.py](gen_plots.py) is included for quick analysis.

### Restarting
Basilisk includes the option to 'dump' the entire simulation to a single file, which can then be restored to continue the simulation from the output time. By default this occurs every 100 time-units. Since stable travelling waves take a long time to develop for most parameter regimes it is *strongly suggested* that a single run without controls is performed to generate a dump-file with a travelling wave before loading it and beginning controls after this point. To do this set the value `"t0"` to correspond to the file at dump/dump-\<time\>.


## Input Parameters
All of the parameters are either in SI units or dimensionless. The keys in [params.json](params.json) are hopefully fairly self-explanatory. However, below is a full description.

#### Domain parameters
* **`h0`** the thickness of a uniform film - m
* **`Lx`** the ratio of film thickness to domain length
* **`Ly`** the ratio of film thickness to domain height (including air layer)
* **`theta`** the angle of the plate from horizontal - rad
* **`tmax`** the simulation end time - dimensionless
* **`t0`** the start time (either 0 or matching a [dump file](#restarting)) - dimensionless

#### Physical parameters
* **`rho_l`** fluid-phase density - kg m^-3
* **`rho_g`** gas-phase density - kg m^-3
* **`mu_l`** fluid-phase dynamic viscosity - kg m^-1 s^-1
* **`mu_g`** gas-phase dynamic viscosity - kg m^-1 s^-1
* **`gamma`** surface tension - N m^-1
* **`grav`** gravitational acceleration - m s^-2

#### Solver parameters
* **`level`** the grid refinement level (resulting in 2^level gridcells)
* **`dtout`** output timestep - dimensionless

#### Control parameters
* **`M`** number of actuators - integer
* **`P`** number of observers - integer
* **`start`** control start time - dimensionless
* **`width`** actuator and observer width parameter. As this approaches 0 the actuators and observers tend towards Dirac delta distributions. If this is too small for the grid to resolve this may cause unexpected results. - dimensionless
* **`alpha`** control strength parameter - dimensionless
* **`del`** actuator/actuator upstream offset - dimensionless
* **`mu`** interface/control cost weighting - dimensionless
* **`rom`** the reduced order model to use for the control strategy - "benney" or "wr"
* **`strategy`** the type of control to use - "pair", "static" or "dynamic"


## Mathematical Background
The purpose of this code is to control a thin liquid film to the flat, Nusselt solution. The Navier-Stokes equations are too complex to apply any established control theoretical results to, and so we instead turn to a hierarchical control method, using reduced order models.

### Reduced order models
We currently consider two ROMs. Instead of describing the evolution of velocity and pressure, they describe the evolution of the film height $h$ and the heigh-averaged flux $q$, with a forcing term $f$ that comes from fluid injection through the base. This results in a mass-conservation equation
```math
h_t + q_x = f.
```

To close the system, the **Benney** equation slaves the flux to the height
```math
q = \frac{h^3}{3}\left(2-2h_x \cot\theta + \frac{h_{xxx}}{Ca}\right) + Re\left(\frac{8h^6h_x}{15} - \frac{2h^4f}{3}\right),
```
and the **weighted-residual** system requires an additional evolution equation for $q$:
```math
\frac{2Reh^2q_t}{5} + q = \frac{h^3}{3}\left(2-2h_x \cot\theta + \frac{h_{xxx}}{Ca}\right) + Re\left(\frac{18q^2h_x}{35} - \frac{34hqq_x}{35} + \frac{hqf}{5}\right).
```

We can then implement control strategies to control the Navier-Stokes system by making the assumption that it is well-approximated by one of these simpler models and applying controls as if we were aiming to control the alternative system rather than the original.

### Control strategies
There are three strategies that we consider.

#### Paired actuators/observers
Here we simply couple an equal number of observers and actuators, with a fixed shift between them. The control is then simply
```math
f_i = -\alpha \left(h(x_i - \delta) - 1\right).
```

#### Static LQR controls
If we relax the restriction of discrete observers, giving the control access to the full interface, we can further simplify the problem by linearising it:
```math
h_t = J(h-1) + \Psi K(h-1).
```
The control operator $K$ that minimises the cost $\kappa$ is provided by the [Linear-Quadratic Regulator](https://en.wikipedia.org/wiki/Linear%E2%80%93quadratic_regulator).

#### Dynamic LQR controls
Although the static LQR controls work well, they are not disadvantages by restricted observations. If we limit ourselves to discrete observers we must use an estimator to derive our controls rather than the interfacial height directly. The estimator approximates the most unstable Fourier modes of the linearised system:
```math
z_t = \left(\tilde{J} + \tilde{\Psi}\tilde{K}\right)z + L\left(\Phi (h-1) - \tilde{\Phi}z\right),
```
and is forced towards the observed system by the operator $L$, which is also set using the LQR.
