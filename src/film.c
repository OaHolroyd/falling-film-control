/* Standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Local headers */
#include "c-utils.h"
#include "parallel.h"
#include "params.h"
#include "b-utils.h"
#include "control.h"

/* Basilisk headers */
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "embed.h" // embedded boundaries (supersedes mask)
#include "tension.h" // surface tension
#include "heights.h" // interfacial height


double *H; // film height
face vector av[]; // acceleration vector field
double G[2]; // gravity
double Ccost; // control cost


/* ========================================================================== */
/*   OUTPUT PARAMETERS                                                        */
/* ========================================================================== */
#define NOUT (1<<(LEVEL-1)) // output resolution
#define LOG_STEP 10 // log every LOG_STEP steps
#define OUTPUT_DAT 1 // dimension of data to output (1 or 2, 0 for no output)
#define DUMP 0 // how often to dump (for restarting)


/* ========================================================================== */
/*   BOUNDARY CONDITIONS                                                      */
/* ========================================================================== */
/* top boundary is embedded (free outflow) */
u.n[embed] = neumann(0.0);
p[embed] = neumann(0.0);
pf[embed] = neumann(0.0);

/* Bottom boundary (no slip with controls) */
u.n[bottom] = dirichlet(control(x)); // add control at the base
u.t[bottom] = dirichlet(0.0);
// p[bottom] = neumann(0.0);
// pf[bottom] = neumann(0.0);


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* Set up the domain size */
void init_domain() {
  init_grid(1 << (LEVEL));

  L0 = LX; // always longer than it is wide
  X0 = 0.0;
  Y0 = 0.0;
}

/* Sets the internal Basilisk parameters */
void set_params() {
  /* model params */
  TOLERANCE = 1e-4; // default 1e-3
  DT = 5e-2;

  /* physical params */
  rho1 = RE;
  rho2 = rho1 * RHO_G/RHO_L;
  mu1 = 1.0;
  mu2 = mu1 * MU_G/MU_L;
  f.sigma = 1.0/CA;

  /* acceleration and gravity */
  // TODO: remove these globals if possible
  a = av;
  G[0] = 2.0/RE;
  G[1] = -2.0/tan(THETA)/RE;
}

/* Initialises the fluid */
void init_fluid() {
  if ((int)T0 == 0) {
    /* Nusselt film */
    // fraction(f, y < 1);

    /* cosine perturbation */
    fraction(f, 1.0-y+0.05*sin(1.0*(2.0/(LX))*M_PI*(x+10)));

    /* initialise with Nusselt velocity */
    foreach () {
      u.x[] = f[]*y*(2.0-y) + (1.0-f[]);
      u.y[] = 0.0;
    }

    /* initialise with Nusselt pressure */
    foreach () {
      p[] = f[]*2*cos(THETA)/sin(THETA)*(1-y);
    }

    boundary({u,f,p});
  } else {
    /* restore from a dump file */
    char dump_file[32];
    sprintf(dump_file, "dump/dump-%04d", (int)T0);
    restore(file = dump_file);
  }

  /* compute heights */
  heights(f,hei);

  /* allocate film */
  H = malloc(N*sizeof(double));

  Ccost = 0.0;
}


/* ========================================================================== */
/*   MAIN                                                                     */
/* ========================================================================== */
int main(int argc, char const *argv[]) {
  /* periodic left/right */
  periodic(right);

  /* if an argument has been submitted, check this input file */
  int err = 0;
  if (argc == 2) {
    err = read_params(argv[1]);
  } else if (argc == 1) {
    err = read_params("params.json");
  } else {
    ABORT("only accepts one argument: the path to a parameter file");
  }

  if (err) {
    ABORT("failed to read parameter file (returned %d)", err);
  }

  /* set up the model */ // TODO allow setting of ct and rt via JSON file
  init_domain();
  set_params();
  control_set(C_STRAT, C_ROM, C_M, C_P, C_W, C_ALPHA, C_MU, C_DEL,
              LX, N, RE, CA, THETA);

  /* sanity check the dimensionless numbers and Nusselt velocity */
  fprintf(stderr, "Us: %.8lf\n", US);
  fprintf(stderr, "Re: %.8lf\n", RE);
  fprintf(stderr, "Ca: %.8lf\n", CA);

  run();

  control_free();
  free(H);

  return EXIT_SUCCESS;
}


/* ========================================================================== */
/*   EVENTS                                                                   */
/* ========================================================================== */
/* set up the fluid */
event init(i=0) {
  /* use solid rather than mask */
  solid(cs, fs, LY + y);
  init_fluid();
}

/* static grid refinement */
event refinement(i=0) {
  /* refine the entire grid to maximum */
  refine(level < LEVEL);

  /* unrefine at higher layers */
  unrefine(y > 3 && level > LEVEL-2);
  unrefine(y > 5 && level > LEVEL-3);
  unrefine(y > 7 && level > LEVEL-4);
}

/* impose acceleration due to gravity TODO: check the air */
event acceleration(i++) {
  foreach_face (x) {
    av.x[] += f[]*G[0];
  }
  foreach_face (y) {
    av.y[] += f[]*G[1];
  }
}

/* set the control magnitudes */
event controls(i++) {
  heights(f,hei);

  /* compute film height */
  for (int i = 0; i < N; i++) {
    H[i] = interfacial_height(ITOX(i));
  } // i end

  if (t >= C_START) {
    control_step(dt, H);
    Ccost += dt * control_cost(H);

    boundary({u.x, u.y}); // update boundary velocities
  }
}

/* print progress to stderr */
#if LOG_STEP
event output_log(i=0; t<=TMAX; i+=LOG_STEP) {
  /* print column headers */
  static int first_log = 1;
  if (first_log) {
    fprintf(stderr, "  model t        dt      iter         N    elap t      cost\n");
    first_log = 0;
  }

  /* print data */
  fprintf(stderr, "\r %8.2lf  %8.5lf  %8d  %8ld  %8.0lf  %8.2lf",
                      t,      dt,    i, grid->n, perf.t, Ccost);

  /* force output */
  fflush (stderr);
}
#endif

/* output 1D and 2D data */
#if OUTPUT_DAT
event output_dat(t=0.0; t<=TMAX; t += DTOUT) {
  /* shift datcount to prevent file overwrites */
  static int datcount = 0;
  char fname[64];

  if (datcount == 0) {
    /* check if there are any files in ./out */
    FILE *fp;
    int i = 0;
    double t0;
    sprintf(fname, "out/data-1-%010d.dat", i);
    while ((fp = fopen(fname, "r"))) {
      /* check if the time is before the current time */
      if (fscanf(fp, "# t: %lf\n", &t0) != 1) { ABORT("missing timestamp"); }
      fclose(fp);
      if (t0 >= t) {
        break;
      }

      /* try the next file */
      i++;
      sprintf(fname, "out/data-1-%010d.dat", i);
    }
    datcount = i;
  }

  FILE *fp;

  /* only compute 2D outputs if required */
  if (OUTPUT_DAT > 1) {
    scalar l[];
    scalar u_mag[];
    scalar u_x[], u_y[];
    scalar omega[];

    /* iterated scalars */
    foreach () {
      l[] = level; // refinement level
      u_mag[] = sqrt(u.x[]*u.x[] + u.y[]*u.y[]); // speed
      u_x[] = u.x[]; // horizontal velocity
      u_y[] = u.y[]; // vertical velocity
    }

    /* vorticity */
    vorticity(u, omega);

    /* 2D fields */
    sprintf(fname, "out/data-2-%010d.dat", datcount);
    fp = fopen(fname, "w");
    if (!fp) { ABORT("'%s' could not be opened", fname); }
    fprintf(fp, "# t: %lf\n", t);
    output_field({f, l, omega, u_mag, u_x, u_y, p}, fp, box = {{0.0,0.0},{LX,LY}}, n = NOUT);
    fclose(fp);
  }

  /* 1D interface */
  sprintf(fname, "out/data-1-%010d.dat", datcount);
  fp = fopen(fname, "w");
  if (!fp) { ABORT("'%s' could not be opened", fname); }
  fprintf(fp, "# t: %lf\n", t);
  double dx = LX/((double)(NOUT));
  for (int i = 0; i < NOUT+1; i++) {
    fprintf(fp, "%lf %lf %lf %lf\n", i*dx, interfacial_height(i*dx), control(i*dx), estimator(i*dx));
  } // i end
  fclose(fp);

  datcount++;
}
#endif

/* dump output every 100 time units */
#if DUMP
event dump_xxx(t=0.0; t+=DUMP) {
  char dump_file[32];
  sprintf(dump_file, "dump/dump-%04.0lf", t);
  dump(file = dump_file);
}
#endif

/* finish */
event stop(t=TMAX) {
  #if DUMP
  char dump_file[32];
  sprintf(dump_file, "dump/dump-%04.0lf-end", t);
  dump(file = dump_file);
  #endif

  #if LOG_STEP
  fprintf(stderr, "\n");
  #endif
  fprintf(stderr, "final time: %lf\n", t);
  fprintf(stderr, "total cost: %lf\n", Ccost);
}
