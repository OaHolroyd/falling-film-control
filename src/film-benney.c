/* Standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* Local headers */
#include "c-utils.h"
#include "parallel.h"
#include "params.h"
#include "linalg.h"
#include "control.h"


/* convert from index to location */
#define ITOX(i) (DX*(i+0.5))

/* convert from location to index */
#define XTOI(x) ((int)(x/DX - 0.5))


/* ========================================================================== */
/*   VARIABLES                                                                */
/* ========================================================================== */
double *x; // x coordinates
double *h, *h_1, *h_2; // interface
double *q; // flux
double *f; // control magnitudes
double *H, *Q, *S, U; // target state
double *hx, *hxxx, *qf; // interface/flux derivatives
double *res; // residual
double **J; // Jacobian
double **qh; // qh
double **CM; // control matrix

#define DT 0.05 // maximum timestep
#define ITERMAX 20 // max iterations
double Cost; // control cost


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* initialises variables */
void init_variables(void) {
  /* allocate memory */
  x = malloc(N*sizeof(double));
  h = malloc(N*sizeof(double));
  h_1 = malloc(N*sizeof(double));
  h_2 = malloc(N*sizeof(double));
  q = malloc(N*sizeof(double));
  f = malloc(N*sizeof(double));
  H = malloc(N*sizeof(double));
  Q = malloc(N*sizeof(double));
  S = malloc(N*sizeof(double));
  hx = malloc(N*sizeof(double));
  hxxx = malloc(N*sizeof(double));
  qf = malloc(N*sizeof(double));
  res = malloc(2*N*sizeof(double));
  J = malloc_f2d(2*N, 2*N);
  CM = malloc_f2d(N, N);
  qh = malloc_f2d(N, N);

  Cost = 0.0;
}

/* sets the initial condition */
void initial_condition(void) {
  for (int i = 0; i < N; i++) {
    x[i] = ITOX(i);
    h[i] = 1.0 + 0.05*sin((2.0/(LX))*M_PI*(x[i]+10.0));
    h_1[i] = h[i]; h_2[i] = h[i];
    q[i] = 0.0;
    f[i] = 0.0;
  } // i end
}

/* set the target state */
void init_target(void) {
  /* target state */
  U = 0.0;
  for (int i = 0; i < N; i++) {
    H[i] = 1.0;
    Q[i] = 2.0/3.0;
    S[i] = 0.0;
  } // i end
}

/* frees variables */
void free_variables(void) {
  free(x);
  free(h);
  free(h_1);
  free(h_2);
  free(q);
  free(f);
  free(H);
  free(Q);
  free(S);
  free(hx);
  free(hxxx);
  free(qf);
  free(res);
  free_2d(J);
  free_2d(CM);
  free_2d(qh);
}

/* outputs the data */
void output(double t) {
  static int datcount = 0;

  if (OUTPUT >= 0) {
    /* L2 interfacial deviation */
    double dh = 0.0;
    for (int i = 0; i < N; i++) {
      dh += (H[i]-h[i])*(H[i]-h[i]);
    } // i end
    dh = sqrt(DX*dh);

    /* L2 control magnitude */
    double dc = 0.0;
    for (int i = 0; i < N; i++) {
      dc += f[i]*f[i];
    } // i end
    dc = sqrt(DX*dc);

    /* L2 estimator deviation */
    double de = 0.0;
    for (int i = 0; i < N; i++) {
      de += (h[i]-H[i]-estimator(ITOX(i)))*(h[i]-H[i]-estimator(ITOX(i)));
    } // i end
    de = sqrt(DX*de);

    /* append data to file */
    FILE *fp = fopen("out/benney-0.dat", "a");
    fprintf(fp, "%lf %lf %lf %lf %lf\n", t-C_START, dh, de, dc, Cost);
    fclose(fp);
  }

  /* output 1D/2D data */
  if (OUTPUT >= 1) {
    /* shift datcount to prevent file overwrites */
    char fname[64];

    if (datcount == 0) {
      /* check if there are any files in ./out */
      FILE *fp;
      int i = 0;
      double t0;
      sprintf(fname, "out/benney-1-%010d.dat", i);
      while ((fp = fopen(fname, "r"))) {
        /* check if the time is before the current time */
        if (fscanf(fp, "# t: %lf\n", &t0) != 1) { ABORT("missing timestamp"); }
        fclose(fp);
        if (t0 >= t) {
          break;
        }

        /* try the next file */
        i++;
        sprintf(fname, "out/benney-1-%010d.dat", i);
      }
      datcount = i;
    }

    /* 1D interface */
    sprintf(fname, "out/benney-1-%010d.dat", datcount);
    FILE *fp = fopen(fname, "w");
    if (!fp) { ABORT("'%s' could not be opened", fname); }
    fprintf(fp, "# t: %lf\n", t);
    for (int i = 0; i < N; i++) {
      fprintf(fp, "%lf %lf %lf %lf %lf\n", x[i], h[i], f[i], estimator(x[i]), q[i]);
    } // i end
    fclose(fp);
  }

  datcount++;
}

/* outputs dimensionless numbers and other details */
void output_numbers(void) {
  FILE *fp = fopen("out/numbers.dat", "w");
  fprintf(fp, "%.8lf\n", RE);
  fprintf(fp, "%.8lf\n", CA);
  fprintf(fp, "%.8lf\n", C_START);
  fclose(fp);
}

/* prints the log */
void print_log(int i, double t) {
  /* print column headers */
  static int first_log = 1;
  static time_t t0 = 0;
  if (first_log) {
    fprintf(stderr, "  model t        dt      iter         N    elap t      cost\n");
    first_log = 0;
    t0 = time(NULL);
  }

  /* print data */
  fprintf(stderr, "\r %8.2lf  %8.5lf  %8d  %8d  %8ld  %8.2lf",
                      t,      DT,    i, N, time(NULL)-t0, Cost);

  /* force output */
  fflush (stderr);
}

/* computes q, hx, hxxx (required for residual and Jacobian) */
void compute_precursors(void) {
  double h3;
  for (int i = 2; i < N-2; i++) {
    h3 = h[i]*h[i]*h[i];
    hx[i] = 0.5/DX*(h[i+1]-h[i-1]);
    hxxx[i] = 0.5/(DX*DX*DX) * (h[i+2]-2.0*h[i+1]+2.0*h[i-1]-h[i-2]);
    q[i] = (h3/3.0)*( 2.0 - 2.0*hx[i]/tan(THETA) + hxxx[i]/CA + RE*(h3*hx[i]*(8.0/5.0) - h[i]*f[i]*2.0) );
    qf[i] = (-2.0*RE/3.0)*h[i]*h3;
  } // i end

  h3 = h[0]*h[0]*h[0];
  hx[0] = 0.5/DX*(h[1]-h[N-1]);
  hxxx[0] = 0.5/(DX*DX*DX) * (h[2]-2.0*h[1]+2.0*h[N-1]-h[N-2]);
  q[0] = (h3/3.0)*( (2.0 - 2.0*hx[0]/tan(THETA) + hxxx[0]/CA) + RE*(h3*hx[0]*(8.0/5.0) - h[0]*f[0]*2.0) );

  h3 = h[1]*h[1]*h[1];
  hx[1] = 0.5/DX*(h[2]-h[0]);
  hxxx[1] = 0.5/(DX*DX*DX) * (h[3]-2.0*h[2]+2.0*h[0]-h[N-1]);
  q[1] = (h3/3.0)*( (2.0 - 2.0*hx[1]/tan(THETA) + hxxx[1]/CA) + RE*(h3*hx[1]*(8.0/5.0) - h[1]*f[1]*2.0) );

  h3 = h[N-2]*h[N-2]*h[N-2];
  hx[N-2] = 0.5/DX*(h[N-1]-h[N-3]);
  hxxx[N-2] = 0.5/(DX*DX*DX) * (h[0]-2.0*h[N-1]+2.0*h[N-3]-h[N-4]);
  q[N-2] = (h3/3.0)*( (2.0 - 2.0*hx[N-2]/tan(THETA) + hxxx[N-2]/CA) + RE*(h3*hx[N-2]*(8.0/5.0) - h[N-2]*f[N-2]*2.0) );

  h3 = h[N-1]*h[N-1]*h[N-1];
  hx[N-1] = 0.5/DX*(h[0]-h[N-2]);
  hxxx[N-1] = 0.5/(DX*DX*DX) * (h[1]-2.0*h[0]+2.0*h[N-2]-h[N-3]);
  q[N-1] = (h3/3.0)*( (2.0 - 2.0*hx[N-1]/tan(THETA) + hxxx[N-1]/CA) + RE*(h3*hx[N-1]*(8.0/5.0) - h[N-1]*f[N-1]*2.0) );
}

/* compute the residual */
void compute_residual(int use_CM) {
  /* residual */
  for (int i = 1; i < N-1; i++) {
    double ht = (h[i] - 4.0*h_1[i]/3.0 + h_2[i]/3.0)/(DT*2.0/3.0);
    double qx = 0.5/DX * (q[i+1] - q[i-1]);
    res[i] = ht - f[i] + qx;
  } // i end

  double ht = (h[0] - 4.0*h_1[0]/3.0 + h_2[0]/3.0)/(DT*2.0/3.0);
  double qx = 0.5/DX * (q[1] - q[N-1]);
  res[0] = ht - f[0] + qx;

  ht = (h[N-1] - 4.0*h_1[N-1]/3.0 + h_2[N-1]/3.0)/(DT*2.0/3.0);
  qx = 0.5/DX * (q[0] - q[N-2]);
  res[N-1] = ht - f[N-1] + qx;

  if (use_CM) {
    // f+CM*(h-H)-S
    for (int i = 0; i < N; i++) {
      res[N+i] = f[i]-S[i];
      for (int k = 0; k < N; k++) {
        res[N+i] += CM[i][k] * (h[k]-H[k]);
      } // k end
    } // i end
  } else {
    // f+0*(h-H)-S
    for (int i = 0; i < N; i++) {
      res[N+i] = f[i]-S[i];
    } // i end
  }
}

/* computes the Jacobian (via qh) */
void compute_jacobian(int use_CM) {
  double h3, c;

  /* construct qh */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      qh[i][j] = 0.0;
    } // j end
  } // i end

  /* central diagonal */
  for (int i = 0; i < N; i++) {
    double h2 = h[i]*h[i];
    qh[i][i] = h2*(2.0-2.0*hx[i]/tan(THETA) + hxxx[i]/CA + (16.0*RE/5.0)*h2*h[i]*hx[i] - (8.0*RE/3.0)*h[i]*f[i]);
  } // i end

  /* 1st diagonals */
  for (int i = 1; i < N-1; i++) {
    h3 = h[i]*h[i]*h[i];
    c = h3*(-2.0/(3.0*tan(THETA)) + (8.0*RE/15.0)*h3);
    qh[i][i-1] = -0.5/DX * c;
    qh[i][i+1] = 0.5/DX * c;
  } // i end

  h3 = h[0]*h[0]*h[0];
  c = h3*(-2.0/(3.0*tan(THETA)) + (8.0*RE/15.0)*h3);
  qh[0][N-1] = -0.5/DX * c;
  qh[0][1] = 0.5/DX * c;

  h3 = h[N-1]*h[N-1]*h[N-1];
  c = h3*(-2.0/(3.0*tan(THETA)) + (8.0*RE/15.0)*h3);
  qh[N-1][N-2] = -0.5/DX * c;
  qh[N-1][0] = 0.5/DX * c;

  /* second diagonals */
  for (int i = 2; i < N-2; i++) {
    c = h[i]*h[i]*h[i]/(3.0*CA*DX*DX*DX);
    qh[i][i-2] = -0.5*c;
    qh[i][i-1] += c;
    qh[i][i+1] += -c;
    qh[i][i+2] = 0.5*c;
  } // i end

  c = h[0]*h[0]*h[0]/(3.0*CA*DX*DX*DX);
  qh[0][N-2] = -0.5*c;
  qh[0][N-1] += c;
  qh[0][1] += -c;
  qh[0][2] = 0.5*c;

  c = h[1]*h[1]*h[1]/(3.0*CA*DX*DX*DX);
  qh[1][N-1] = -0.5*c;
  qh[1][0] += c;
  qh[1][2] += -c;
  qh[1][3] = 0.5*c;

  c = h[N-2]*h[N-2]*h[N-2]/(3.0*CA*DX*DX*DX);
  qh[N-2][N-4] = -0.5*c;
  qh[N-2][N-3] += c;
  qh[N-2][N-1] += -c;
  qh[N-2][0] = 0.5*c;

  c = h[N-1]*h[N-1]*h[N-1]/(3.0*CA*DX*DX*DX);
  qh[N-1][N-3] = -0.5*c;
  qh[N-1][N-2] += c;
  qh[N-1][0] += -c;
  qh[N-1][1] = 0.5*c;


  /* set J */
  /* J00 (I/(2*DT/3)+Dx*qh) */
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < N; j++) {
      J[i][j] = (0.5/DX)*(qh[i+1][j] - qh[i-1][j]);
    } // j end
    J[i][i] += 3.0/(2.0*DT);
  } // i end

  for (int j = 0; j < N; j++) {
    J[0][j] = (0.5/DX)*(qh[1][j] - qh[N-1][j]);
  } // j end
  J[0][0] += 3.0/(2.0*DT);

  for (int j = 0; j < N; j++) {
    J[N-1][j] = (0.5/DX)*(qh[0][j] - qh[N-2][j]);
  } // j end
  J[N-1][N-1] += 3.0/(2.0*DT);

  /* J01 (-I+Dx*qf) */
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < N; j++) {
      J[i][N+j] = 0.0;
    } // j end
    J[i][N+i-1] = -0.5/DX*qf[i-1];
    J[i][N+i] = -1.0;
    J[i][N+i+1] = 0.5/DX*qf[i+1];
  } // i end

  for (int j = 0; j < N; j++) {
    J[0][N+j] = 0.0;
  } // j end
  J[0][2*N-1] = -0.5/DX*qf[N-1];
  J[0][N] = -1.0;
  J[0][N+1] = 0.5/DX*qf[1];

  for (int j = 0; j < N; j++) {
    J[N-1][N+j] = 0.0;
  } // j end
  J[N-1][2*N-2] = -0.5/DX*qf[N-2];
  J[N-1][2*N-1] = -1.0;
  J[N-1][N] = 0.5/DX*qf[0];

  /* J10 (CM) */
  if (use_CM) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        J[N+i][j] = CM[i][j];
      } // j end
    } // i end
  } else {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        J[N+i][j] = 0.0;
      } // j end
    } // i end
  }

  /* J11 (I) */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      J[N+i][N+j] = 0.0;
    } // j end
    J[N+i][N+i] = 1.0;
  } // i end
}


/* ========================================================================== */
/*   MAIN                                                                     */
/* ========================================================================== */
int main(int argc, char const *argv[]) {
  /* ========== */
  /*   SET UP   */
  /* ========== */
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

  /* ensure that T0 is 0.0 */
  if (T0 != 0.0) {
    T0 = 0.0;
    fprintf(stderr, "WARNING: Benney solver does not support dump/restore\n");
  }

  /* set up the model */
  init_variables();
  initial_condition();
  init_target();

  /* compute the controls */
  control_set(C_STRAT, C_ROM, C_M, C_P, C_W, C_ALPHA, C_MU, C_DEL,
              LX, N, RE, CA, THETA);
  control_output();
  control_matrix(CM);

  /* sanity check the dimensionless numbers and Nusselt velocity */
  fprintf(stderr, "Us: %.8lf\n", US);
  fprintf(stderr, "Re: %.8lf\n", RE);
  fprintf(stderr, "Ca: %.8lf\n", CA);
  output_numbers();


  /* ============ */
  /*   TIMESTEP   */
  /* ============ */
  double t = 0.0;
  int nsteps = 0;
  output(t);
  while (t < TMAX) {
    /* log */
    if (nsteps%LOG_STEP == 0) {
      print_log(nsteps, t);
    }

    /* turn off controls before C_START */
    int use_CM = t>=C_START;

    /* iterate to a solution */
    for (int k = 0; k < ITERMAX; k++) {
      /* compute q, hx, hxxx, res */
      compute_precursors();
      compute_residual(use_CM);

      /* stop if residual is small enough */
      double nres = 0.0;
      for (int i = 0; i < 2*N; i++) {
        nres += res[i]*res[i];
      } // i end
      if (nres < 1e-14) { break; }

      /* iterate to improve h and f */
      compute_jacobian(use_CM);
      dsv(J, res, 2*N);

      for (int i = 0; i < N; i++) {
        h[i] -= res[i];
      } // i end
      for (int i = 0; i < N; i++) {
        f[i] -= res[N+i];
      } // i end
    } // k end

    /* check for blowup */
    if (isnan(h[0])) {
      fprintf(stderr, "\nBlowup at t = %lf\n", t);
      break;
    }

    /* copy backwards */
    for (int i = 0; i < N; i++) {
      h_2[i] = h_1[i];
      h_1[i] = h[i];
    } // i end


    /* step time forward */
    nsteps++;
    t += DT;

    /* output data */
    if (fabs(fmod(t, DTOUT)) < 0.5*DT || fabs(fmod(t, DTOUT)-DTOUT) < 0.5*DT) {
      output(t);
    }

    /* contribute to cost */
    if (t > C_START) {
      for (int i = 0; i < N; i++) {
        Cost += DT*DX * (C_MU*(h[i]-H[i])*(h[i]-H[i]) + (1.0-C_MU)*f[i]*f[i]);
      } // i end
    }
  }

  /* clean up */
  control_free();
  free_variables();
  fprintf(stderr, "\n");

  return EXIT_SUCCESS;
}
