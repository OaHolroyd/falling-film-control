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

#include "debug.h"


/* convert from index to location */
#define ITOX(i) (DX*(i+0.5))

/* convert from location to index */
#define XTOI(x) ((int)(x/DX - 0.5))


/* ========================================================================== */
/*   VARIABLES                                                                */
/* ========================================================================== */
double *x; // x coordinates
double *h, *h_1, *h_2; // interface
double *q, *q_1, *q_2; // flux
double *f; // control magnitudes
double *H, *Q, *S, U; // target state
double *hx, *hxxx, *qx; // interface/flux derivatives
double *res; // residual
double **J; // Jacobian
double **Rq_th; // Rq_th
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
  q_1 = malloc(N*sizeof(double));
  q_2 = malloc(N*sizeof(double));
  f = malloc(N*sizeof(double));
  H = malloc(N*sizeof(double));
  Q = malloc(N*sizeof(double));
  S = malloc(N*sizeof(double));
  hx = malloc(N*sizeof(double));
  hxxx = malloc(N*sizeof(double));
  qx = malloc(N*sizeof(double));
  res = malloc(3*N*sizeof(double));
  J = malloc_f2d(3*N, 3*N);
  CM = malloc_f2d(N, N);
  Rq_th = malloc_f2d(N, N);

  Cost = 0.0;
}

/* sets the initial condition */
void initial_condition(void) {
  for (int i = 0; i < N; i++) {
    x[i] = ITOX(i);
    h[i] = 1.0 + 0.05*sin((2.0/(LX))*M_PI*(x[i]+10.0));
    h_1[i] = h[i]; h_2[i] = h[i];
    q[i] = 0.0;
    q_1[i] = q[i]; q_2[i] = q[i];
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
  free(q_1);
  free(q_2);
  free(f);
  free(H);
  free(Q);
  free(S);
  free(hx);
  free(hxxx);
  free(qx);
  free(res);
  free_2d(J);
  free_2d(CM);
  free_2d(Rq_th);
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
    FILE *fp = fopen("out/wr-0.dat", "a");
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
      sprintf(fname, "out/wr-1-%010d.dat", i);
      while ((fp = fopen(fname, "r"))) {
        /* check if the time is before the current time */
        if (fscanf(fp, "# t: %lf\n", &t0) != 1) { ABORT("missing timestamp"); }
        fclose(fp);
        if (t0 >= t) {
          break;
        }

        /* try the next file */
        i++;
        sprintf(fname, "out/wr-1-%010d.dat", i);
      }
      datcount = i;
    }

    /* 1D interface */
    sprintf(fname, "out/wr-1-%010d.dat", datcount);
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
  for (int i = 2; i < N-2; i++) {
    hx[i] = 0.5/DX*(h[i+1]-h[i-1]);
    hxxx[i] = 0.5/(DX*DX*DX) * (h[i+2]-2.0*h[i+1]+2.0*h[i-1]-h[i-2]);
    qx[i] = 0.5/DX*(q[i+1]-q[i-1]);
  } // i end

  hx[0] = 0.5/DX*(h[1]-h[N-1]);
  hxxx[0] = 0.5/(DX*DX*DX) * (h[2]-2.0*h[1]+2.0*h[N-1]-h[N-2]);
  qx[0] = 0.5/DX*(q[1]-q[N-1]);

  hx[1] = 0.5/DX*(h[2]-h[0]);
  hxxx[1] = 0.5/(DX*DX*DX) * (h[3]-2.0*h[2]+2.0*h[0]-h[N-1]);
  qx[1] = 0.5/DX*(q[2]-q[0]);

  hx[N-2] = 0.5/DX*(h[N-1]-h[N-3]);
  hxxx[N-2] = 0.5/(DX*DX*DX) * (h[0]-2.0*h[N-1]+2.0*h[N-3]-h[N-4]);
  qx[N-2] = 0.5/DX*(q[N-1]-q[N-3]);

  hx[N-1] = 0.5/DX*(h[0]-h[N-2]);
  hxxx[N-1] = 0.5/(DX*DX*DX) * (h[1]-2.0*h[0]+2.0*h[N-2]-h[N-3]);
  qx[N-1] = 0.5/DX*(q[0]-q[N-2]);
}

/* compute the residual */
void compute_residual(int use_CM) {
  /* Rh (ht - U*hx - f + qx) */
  for (int i = 0; i < N; i++) {
    double ht = (h[i] - 4.0*h_1[i]/3.0 + h_2[i]/3.0)/(DT*2.0/3.0);
    res[i] = ht - U*hx[i] - f[i] + qx[i];
  } // i end

  /* Rq (qt - U*Dx*q - Rq_t/RE) */
  for (int i = 0; i < N; i++) {
    double h2 = h[i]*h[i];
    double q2 = q[i]*q[i];
    double qt = (q[i] - 4.0*q_1[i]/3.0 + q_2[i]/3.0)/(DT*2.0/3.0);
    double Rq_t = -5.0*q[i]/(2.0*h2)+(5.0*h[i]/6.0)*(2.0-2.0*hx[i]/tan(THETA)+hxxx[i]/CA)
                  + RE*(9.0*q2*hx[i]/(7.0*h2)-17.0*q[i]*qx[i]/(7.0*h[i])+f[i]*q[i]/(2.0*h[i]));
    res[N+i] = qt - U*qx[i] - Rq_t/RE;
  } // i end

  /* Rf (-CM*(h-H)+S-f) */
  // TODO: check if this could be made negative
  if (use_CM) {
    // -CM*(h-H)+S-f
    for (int i = 0; i < N; i++) {
      res[2*N+i] = S[i]-f[i];
      for (int k = 0; k < N; k++) {
        res[2*N+i] -= CM[i][k] * (h[k]-H[k]);
      } // k end
    } // i end
  } else {
    // S-f
    for (int i = 0; i < N; i++) {
      res[2*N+i] = S[i]-f[i];
    } // i end
  }
}

/* computes the Jacobian (via Rq_th) */
void compute_jacobian(int use_CM) {
  double c;

  /* construct Rq_th */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Rq_th[i][j] = 0.0;
    } // j end
  } // i end

  /* central diagonal */
  for (int i = 0; i < N; i++) {
    double h2 = h[i]*h[i];
    double q2 = q[i]*q[i];
    Rq_th[i][i] = 5.0*q[i]/(h2*h[i]) + (5.0/6.0)*(2.0-2.0*hx[i]/tan(THETA) + hxxx[i]/CA)
                  - 18.0*RE*q2*hx[i]/(7.0*h2*h[i])+17.0*RE*q[i]*qx[i]/(7.0*h2)-0.5*RE*f[i]*q[i]/h2;
  } // i end

  /* 1st diagonals */
  for (int i = 1; i < N-1; i++) {
    c = -5.0*h[i]/(3.0*tan(THETA)) + RE*9.0*q[i]*q[i]/(7.0*h[i]*h[i]);
    Rq_th[i][i-1] = -0.5/DX * c;
    Rq_th[i][i+1] = 0.5/DX * c;
  } // i end

  c = -5.0*h[0]/(3.0*tan(THETA)) + RE*9.0*q[0]*q[0]/(7.0*h[0]*h[0]);
  Rq_th[0][N-1] = -0.5/DX * c;
  Rq_th[0][1] = 0.5/DX * c;

  c = -5.0*h[N-1]/(3.0*tan(THETA)) + RE*9.0*q[N-1]*q[N-1]/(7.0*h[N-1]*h[N-1]);
  Rq_th[N-1][N-2] = -0.5/DX * c;
  Rq_th[N-1][0] = 0.5/DX * c;

  /* second diagonals */
  for (int i = 2; i < N-2; i++) {
    c = 5.0*h[i]/(6.0*CA*DX*DX*DX);
    Rq_th[i][i-2] = -0.5*c;
    Rq_th[i][i-1] += c;
    Rq_th[i][i+1] += -c;
    Rq_th[i][i+2] = 0.5*c;
  } // i end

  c = 5.0*h[0]/(6.0*CA*DX*DX*DX);
  Rq_th[0][N-2] = -0.5*c;
  Rq_th[0][N-1] += c;
  Rq_th[0][1] += -c;
  Rq_th[0][2] = 0.5*c;

  c = 5.0*h[1]/(6.0*CA*DX*DX*DX);
  Rq_th[1][N-1] = -0.5*c;
  Rq_th[1][0] += c;
  Rq_th[1][2] += -c;
  Rq_th[1][3] = 0.5*c;

  c = 5.0*h[N-2]/(6.0*CA*DX*DX*DX);
  Rq_th[N-2][N-4] = -0.5*c;
  Rq_th[N-2][N-3] += c;
  Rq_th[N-2][N-1] += -c;
  Rq_th[N-2][0] = 0.5*c;

  c = 5.0*h[N-1]/(6.0*CA*DX*DX*DX);
  Rq_th[N-1][N-3] = -0.5*c;
  Rq_th[N-1][N-2] += c;
  Rq_th[N-1][0] += -c;
  Rq_th[N-1][1] = 0.5*c;


  /* set J */
  /* J00 (I/(2*DT/3)-U*Dx) */
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < N; j++) {
      J[i][j] = 0.0;
    } // j end
    J[i][i-1] = U*0.5/DX;
    J[i][i] = 3.0/(2.0*DT);
    J[i][i+1] = -U*0.5/DX;
  } // i end

  for (int j = 0; j < N; j++) {
    J[0][j] = 0.0;
  } // j end
  J[0][N-1] = U*0.5/DX;
  J[0][0] = 3.0/(2.0*DT);
  J[0][1] = -U*0.5/DX;

  for (int j = 0; j < N; j++) {
    J[N-1][j] = 0.0;
  } // j end
  J[N-1][N-2] = U*0.5/DX;
  J[N-1][N-1] = 3.0/(2.0*DT);
  J[N-1][0] = -U*0.5/DX;


  /* J01 (Dx) */
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < N; j++) {
      J[i][N+j] = 0.0;
    } // j end
    J[i][N+i-1] = -0.5/DX;
    J[i][N+i+1] = 0.5/DX;
  } // i end

  for (int j = 0; j < N; j++) {
    J[0][N+j] = 0.0;
  } // j end
  J[0][2*N-1] = -0.5/DX;
  J[0][N+1] = 0.5/DX;

  for (int j = 0; j < N; j++) {
    J[N-1][N+j] = 0.0;
  } // j end
  J[N-1][2*N-2] = -0.5/DX;
  J[N-1][N] = 0.5/DX;


  /* J02 (-I) */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      J[i][2*N+j] = 0.0;
    } // j end
    J[i][2*N+i] = -1.0;
  } // i end


  /* J10 (-Rq_th/RE) */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      J[N+i][j] = -Rq_th[i][j]/RE;
    } // j end
  } // i end


  /* J11 (I/(2*DT/3)-RE*U*Dx-Rq_tq/RE) */
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < N; j++) {
      J[N+i][N+j] = 0.0;
    } // j end
    double Rq_tq = -2.5/(h[i]*h[i]) + RE*18.0*hx[i]*q[i]/(7.0*h[i]*h[i])
                   - 17.0*RE*qx[i]/(7.0*h[i]) +RE*f[i]/(2*h[i]);
    J[N+i][N+i-1] = RE*U*0.5/DX;
    J[N+i][N+i] = 3.0/(2.0*DT)-Rq_tq/RE;
    J[N+i][N+i+1] = -RE*U*0.5/DX;
  } // i end

  for (int j = 0; j < N; j++) {
    J[N][N+j] = 0.0;
  } // j end
  double Rq_tq = -2.5/(h[0]*h[0]) + RE*18.0*hx[0]*q[0]/(7.0*h[0]*h[0])
                 - 17.0*RE*qx[0]/(7.0*h[0]) +RE*f[0]/(2*h[0]);
  J[N][2*N-1] = RE*U*0.5/DX;
  J[N][N] = 3.0/(2.0*DT)-Rq_tq/RE;
  J[N][N+1] = -RE*U*0.5/DX;

  for (int j = 0; j < N; j++) {
    J[2*N-1][N+j] = 0.0;
  } // j end
  Rq_tq = -2.5/(h[N-1]*h[N-1]) + RE*18.0*hx[N-1]*q[N-1]/(7.0*h[N-1]*h[N-1])
          - 17.0*RE*qx[N-1]/(7.0*h[N-1]) +RE*f[N-1]/(2*h[N-1]);
  J[2*N-1][2*N-2] = RE*U*0.5/DX;
  J[2*N-1][2*N-1] = 3.0/(2.0*DT)-Rq_tq/RE;
  J[2*N-1][N] = -RE*U*0.5/DX;


  /* J12 (-Rq_tf/RE) */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      J[N+i][2*N+j] = 0.0;
    } // j end
    J[N+i][2*N+i] = -0.5*q[i]/h[i];
  } // i end


  /* J20 (-CM) */
  if (use_CM) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        J[2*N+i][j] = -CM[i][j];
      } // j end
    } // i end
  } else {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        J[2*N+i][j] = 0.0;
      } // j end
    } // i end
  }


  /* J21 (0) */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      J[2*N+i][N+j] = 0.0;
    } // j end
  } // i end


  /* J22 (-I) */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      J[2*N+i][2*N+j] = 0.0;
    } // j end
    J[2*N+i][2*N+i] = -1.0;
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
      for (int i = 0; i < 3*N; i++) {
        nres += res[i]*res[i];
      } // i end
      if (nres < 1e-14) { break; }

      /* iterate to improve h and f */
      compute_jacobian(use_CM);
      dsv(J, res, 3*N);

      for (int i = 0; i < N; i++) {
        h[i] -= res[i];
      } // i end
      for (int i = 0; i < N; i++) {
        q[i] -= res[N+i];
      } // i end
      for (int i = 0; i < N; i++) {
        f[i] -= res[2*N+i];
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
      q_2[i] = q_1[i];
      q_1[i] = q[i];
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
