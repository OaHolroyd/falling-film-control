/* Standard headers */
#include <stdlib.h>

/* Basilisk headers */
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "embed.h" // embedded boundaries (supersedes mask)
#include "tension.h" // surface tension
#include "heights.h" // interfacial height


/* ========================================================================== */
/*   PHYSICAL PARAMETERS                                                      */
/* ========================================================================== */
/*  All of the parameters are read in with SI units and non-dimensionalised   */
/*  using the following scales:                                               */
/*    velocity - Us                                                           */
/*    length   - h0                                                           */
/*    pressure - mu_l Us / h0                                                 */
/*  This is the scaling used in Thompson 2016 and Cimpeanu 2021.              */

/* Domain parameters */
#define h0 (175e-6) // film thickness
#define nl 64.0 // dimensionless domain length (Lx/h0)
#define nh 8.0 // dimensionless domain height (Ly/h0)
#define theta 1.047197551 // inclination angle

/* Time parameters */
#define tmax 100.0 // final time
#define T0 0 // initial time (0 or an integer corresponding to a dump file)

/* Physical parameters */
#define rho_l 998.0
#define rho_g 1.17 // densities
#define mu_l (8.967e-4)
#define mu_g (1.836e-5) // (dynamic) viscosities
#define gamma 0.072 // surface tension
#define grav 9.807 // acceleration due to gravity

/* Solver parameters */
#define LEVEL_MAX 8 // maximum refinement level
#define tol_f 0.0001 // fluid fraction tolerance
#define tol_u 0.01 // velocity tolerance
#define dtout 1.0 // output step
#define nout (1<<(LEVEL_MAX-1)) // output resolution

/* Dimensionless numbers */
#define Us ((rho_l*grav*sin(theta)*h0*h0)/(2*mu_l)) // Nusselt surface velocity
#define Re ((rho_l*Us*h0)/mu_l) // Reynolds number
#define Ca ((mu_l*Us)/gamma) // capillary number


/* ========================================================================== */
/*   CONTROL PARAMETERS                                                       */
/* ========================================================================== */
#define PI 3.14159265358979323846
#define C_M 5 // number of controls
#define C_START 80000.0 // control start time
double C_loc[C_M]; // control locations
double C_mag[C_M]; // current control magnitudes
double C_norm; // control normaliser
#define C_W 0.01 // control width parameter
#define C_ALPHA 1.0 // control strength
#define C_PHI 1.0 // observer/control displacement


/* ========================================================================== */
/*   OUTPUT FLAGS                                                             */
/* ========================================================================== */
#define LOG_STEP 10 // log every LOG_STEP steps
#define OUTPUT_DAT 1 // whether to output the data

/* harmonic viscosity averaging TODO: check this vs other types */
#ifdef mu
#undef mu
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#endif

face vector av[]; // acceleration vector field
double G[2]; // gravity

int nx; // x-resolution
double dx; // gridspacing


/* ========================================================================== */
/*   Boundary Conditions                                                      */
/* ========================================================================== */
/* top boundary is embedded (free outflow) */
u.n[embed] = neumann(0.0);
p[embed] = dirichlet(0.0);
pf[embed] = dirichlet(0.0);

/* Bottom boundary (no slip with controls) */
u.n[bottom] = dirichlet(control(x)); // add control at the base
u.t[bottom] = dirichlet(0.0);
// p[bottom] = neumann(0.0);
// pf[bottom] = neumann(0.0);


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* discrete actuator function */
double actuator(double x) {

  return C_norm*exp((cos(2*PI*x/nl)-1.0)/(C_W*C_W));
}

/* Set up the domain size */
void init_domain() {
  init_grid(1 << (LEVEL_MAX));

  L0 = nl; // always longer than it is wide
  X0 = 0.0;
  Y0 = 0.0;
}

/* Sets the internal Basilisk parameters */
void set_params() {
  /* model params */
  TOLERANCE = 1e-4; // default 1e-3
  DT = 5e-2;

  /* physical params */
  rho1 = Re;
  rho2 = rho1 * rho_g/rho_l;
  mu1 = 1.0;
  mu2 = mu1 * mu_g/mu_l;
  f.sigma = 1.0/Ca;

  /* acceleration and gravity */
  a = av;
  G[0] = 2.0/Re;
  G[1] = -2.0/tan(theta)/Re;

  /* output domain params */
  nx = 1<<(LEVEL_MAX-1);
  dx = nl/((double)(nx));
}

/* Sets the control variables (call after set_params) */
void set_Cparams() {
  /* control locations */
  double dc = nl/C_M;
  for (int i = 0; i < C_M; i++) {
    C_loc[i] = (i+0.5)*dc;
    C_mag[i] = 0.0;
  } // i end

  /* set control normaliser */
  C_norm = 1.0;
  double I = 0.0;
  for (int i = 0; i < nx; i++) {
    I += actuator(dx*i - nl/2.0);
  } // i end
  C_norm = 1.0/(dx*I);
}

/* Initialises the fluid */
void init_fluid() {
  if (T0 == 0) {
    /* Nusselt film */
    // fraction(f, y < 1);

    /* cosine perturbation */
    fraction(f, 1.0-y+0.05*sin(1.0*(2.0/(nl))*M_PI*(x+10)));

    /* initialise with Nusselt velocity */
    foreach () {
      u.x[] = f[]*y*(2.0-y) + (1.0-f[]);
      u.y[] = 0.0;
    }

    /* initialise with Nusselt pressure */
    foreach () {
      p[] = f[]*2*cos(theta)/sin(theta)*(1-y);
    }

    boundary({u,f,p});
  } else {
    /* restore from a dump file */
    char dump_file[32];
    sprintf(dump_file, "dump/dump-%04d", T0);
    restore(file = dump_file);
  }
}

/* Get interfacial height at a given x-coord */
#define NP 10
double interfacial_height(double xp, vector hei) {
  if (xp < 0) {
    xp += nl;
  } else if (xp > nl) {
    xp -= nl;
  }

  double dh[NP];
  double yh[NP];
  double yp;
  double y0 = 0.0, y1 = 2.0;

  /* try and find range of possible heights */
  for (int i = 0; i < NP; i++) {
    yp = y0 + i*(y1-y0)/(NP-1);
    Point point = locate(xp,yp);

    if (hei.y[] != nodata) {
      yh[i] = y + height(hei.y[])*Delta;
      dh[i] = abs(y-yh[i]);
    } else {
      yh[i] = -1000;
      dh[i] = 1000;
    }
  } // i end

  /* find the closest one */
  int j = 0;
  for (int i = 1; i < NP; i++) {
    if (dh[i] < dh[j]) { j = i; }
  } // i end

  return yh[j];
}

/* returns the baseplate control velocity as a function of x */
double control(double x) {
  double vc = 0.0;

  for (int i = 0; i < C_M; i++) {
    vc += C_mag[i] * actuator(x-C_loc[i]);
  } // i end

  return -C_ALPHA*vc;
}


/* ========================================================================== */
/*   MAIN                                                                     */
/* ========================================================================== */
int main(int argc, char const *argv[]) {
  /* periodic left/right */
  periodic(right);

  /* set up the model */
  init_domain();
  set_params();
  set_Cparams();

  /* sanity check the dimensionless numbers and Nusselt velocity */
  fprintf(stderr, "Us: %.8lf\n", Us);
  fprintf(stderr, "Re: %.8lf\n", Re);
  fprintf(stderr, "Ca: %.8lf\n", Ca);

  run();

  return EXIT_SUCCESS;
}


/* ========================================================================== */
/*   EVENTS                                                                   */
/* ========================================================================== */
/* set up the fluid */
event init(i=0) {
  /* use solid rather than mask */
  solid(cs, fs, nh + y);
  init_fluid();
}

/* impose acceleration due to gravity */
event acceleration (i++) {
  foreach_face (x) {
    av.x[] += f[]*G[0];
  }
  foreach_face (y) {
    av.y[] += f[]*G[1];
  }
}

/* static grid refinement */
event refinement(i=0) {
  /* refine the entire grid to maximum */
  refine(level < LEVEL_MAX);

  /* unrefine at higher layers */
  unrefine(y > 3 && level > LEVEL_MAX-2);
  unrefine(y > 5 && level > LEVEL_MAX-3);
  unrefine(y > 7 && level > LEVEL_MAX-4);
}

/* set the control magnitudes */
event controls(i++) {
  if (t > C_START) {
    vector dh[];
    heights(f,dh);

    for (int i = 0; i < C_M; i++) {
      C_mag[i] = interfacial_height(C_loc[i] - C_PHI, dh) - 1;
    } // i end

    boundary({u.x, u.y}); // update boundary velocities
  }
}

/* print progress to stderr */
#if LOG_STEP
int first_log = 1;
event output_log(i=0; t<=tmax; i+=LOG_STEP) {
  /* print column headers */
  if (first_log) {
    fprintf(stderr, "  model t        dt      iter         N    elap t  \n");
    first_log = 0;
  }

  /* print data */
  fprintf(stderr, "\r %8.2lf  %8.5lf  %8d  %8ld  %8.0lf ",t,dt,i,grid->n,perf.t);

  /* force output */
  fflush (stderr);
}
#endif

/* output 1D and 2D data */
#if OUTPUT_DAT
int datcount = 0;
char fname[64];
event output_dat(t=0.0; t<=tmax; t += dtout) {
  /* refinement level */
  scalar l[];
  foreach () {
    l[] = level;
  }

  /* vorticity */
  scalar omega[];
  vorticity(u, omega);

  /* speed */
  scalar u_mag[];
  foreach () {
    u_mag[] = sqrt(u.x[]*u.x[] + u.y[]*u.y[]);
  }

  /* velocity components */
  scalar u_x[], u_y[];
  foreach () {
    u_x[] = u.x[];
    u_y[] = u.y[];
  }

  /* distance to interface */
  vector dh[];
  heights(f,dh);

  /* interfacial height */
  scalar yh[];
  foreach () {
    yh[] = interfacial_height(x, dh);
  }

  /* 2D fields */
  sprintf(fname, "out/data-2-%010d.dat", datcount);
  FILE *fp = fopen(fname, "w");
  fprintf(fp, "# t: %lf\n", t);
  output_field({f, l, omega, u_mag, u_x, u_y, p, yh}, fp, box = {{0.0,0.0},{nl,nh}}, n = nout);
  fclose(fp);

  /* 1D interface */
  sprintf(fname, "out/data-1-%010d.dat", datcount);
  fp = fopen(fname, "w");
  int i;
  fprintf(fp, "# t: %lf\n", t);
  for (i = 0; i < nx+1; i++) {
    fprintf(fp, "%lf %lf %lf\n", i*dx, interfacial_height(i*dx, dh), control(i*dx));
  } // i end
  fclose(fp);

  datcount++;
}
#endif

/* TODO: doesn't work, try embed */
event dump_xxx(t=0.0; t+=100) {
  char dump_file[32];
  sprintf(dump_file, "dump/dump-%04.0lf", t);
  dump(file = dump_file);
}

/* finish */
event stop(t=tmax) {
  char dump_file[32];
  sprintf(dump_file, "dump/dump-%04.0lf", t);
  dump(file = dump_file);

  #if LOG_STEP
  fprintf(stderr, "\n");
  #endif
  fprintf(stderr, "final time: %lf\n", t);
}
