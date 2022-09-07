#ifndef CONTROL_CORE_H
#define CONTROL_CORE_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "c-utils.h"
#include "control.h"


/* convert from index to location */
#define ITOX(i) (DX*(i+0.5))

/* convert from location to index */
#define XTOI(x) ((int)(x/DX - 0.5))


/* ========================================================================== */
/*   GLOBAL DEFINITIONS                                                       */
/* ========================================================================== */
/* film constants */
static int N; // number of gridcells
static double LX; // domain length
static double DX; // gridspacing
static double RE; // Reynolds' number
static double CA; // capillary number
static double THETA; // plate angle

/* control constants */
static int M; // number of actuators
static int P; // number of observers
static double W; // width parameter
static double ALPHA; // control strength parameter
static double MU; // control cost parameter
static double NORM; // normalising constant
static double DEL; // observer offset (upstream)
static rom_t RT; // type of reduced order model

/* location arrays */
static double *Aloc; // actuator locations
static double *Oloc; // observer locations

/* variable arrays */
static double *Amag; // actuator magnitudes


/* ========================================================================== */
/*   FUNCTION DECLARATIONS                                                    */
/* ========================================================================== */
/* forcing matrix (N-by-M) */
void forcing_matrix(double **F);

/* Jacobian (N-by-N) */
void benney_jacobian(double **J);

/* Actuator (N-by-M) */
void benney_actuator(double **Psi);

/* (the transpose of the) Observer (N-by-P) */
void benney_observer(double **Phi);

/* Jacobian (2N-by-2N) */
void wr_jacobian(double **J);

/* Actuator (2N-by-M) */
void wr_actuator(double **Psi);

/* (the transpose of the) Observer (2N-by-2P) */
void wr_observer(double **Phi);


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* interpolates the discretised film height h at a position x */
double interp(double x, double *h) {
  /* account for periodicity */
  if (x < 0) {
    x += LX;
  } else if (x > LX) {
    x -= LX;
  }

  /* get neighbouring indices */
  double i = x/DX - 0.5;
  double _i0 = floor(i);
  double _i1 = ceil(i);
  int i0 = (int) _i0;
  int i1 = (int) _i1;


  /* cover the case that x is on a gridpoint */
  if (i0 == i1) {
    return h[i0];
  }

  /* linear interpolation */
  return h[i0]*(i1-i) + h[i1]*(i-i0);
}

/* (periodic) actuator function (only valid on [-3Lx/2, 3Lx/2]) */
double actuator(double x) {
  /* account for periodicity */
  if (x < -LX/2) {
    x += LX;
  } else if (x > LX/2) {
    x -= LX;
  }

  return NORM*exp((cos(2*M_PI*x/LX)-1.0)/(W*W));
}

/* sets the common control parameters and allocates common memory */
void internal_control_set(rom_t rt, int m, int p, double w, double alpha, double mu, double del, double lx, int n, double re, double ca, double theta) {
  /* set constants */
  N = n;
  LX = lx;
  DX = LX/N;
  RE = re;
  CA = ca;
  THETA = theta;

  M = m;
  P = p;
  W = w;
  ALPHA = alpha;
  MU = mu;
  DEL = del;
  RT = rt;

  /* actuator locations/magnitudes */
  Aloc = malloc(M*sizeof(double));
  Amag = malloc(M*sizeof(double));
  for (int i = 0; i < M; i++) {
    Aloc[i] = (i+0.5) * LX / M;
    Amag[i] = 0.0;
  } // i end

  /* observer locations */
  Oloc = malloc(P*sizeof(double));
  for (int i = 0; i < P; i++) {
    Oloc[i] = (i+0.5) * LX / P - DEL;

    /* wrap via periodicity */
    // TODO: this will break if |DEL| > LX
    if (Oloc[i] < 0) {
      Oloc[i] += LX;
    } else if (Oloc[i] > LX) {
      Oloc[i] -= LX;
    }
  } // i end

  /* control normaliser */
  NORM = 1.0;
  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    integral += actuator(DX*i - LX/2);
  } // i end
  NORM = 1.0/(DX*integral);
}

/* frees the common memory */
void internal_control_free(void) {
  free(Aloc);
  free(Oloc);
  free(Amag);
}

/* outputs common arrays */
void internal_control_output(void) {
  output_d1d("out/Aloc.dat", Aloc, M);
  output_d1d("out/Oloc.dat", Oloc, P);

  double **F = malloc_f2d(N, M);
  forcing_matrix(F);
  output_d2d("out/F.dat", F, N, M);

  double **Psi, **Phi;
  switch (RT) {
    case BENNEY:
      /* actuator matrix */
      Psi = malloc_f2d(N, M);
      benney_actuator(Psi);

      /* observer matrix (actually the transpose) */
      Phi = malloc_f2d(N, P);
      benney_observer(Phi);

      output_d2d("out/Psi.dat", Psi, N, M);
      output_d2d("out/Phi.dat", Phi, N, P);
      break;
    case WR:
      /* actuator matrix */
      Psi = malloc_f2d(2*N, M);
      wr_actuator(Psi);

      /* observer matrix (actually the transpose) */
      Phi = malloc_f2d(2*N, 2*P);
      wr_observer(Phi);

      output_d2d("out/Psi.dat", Psi, 2*N, M);
      output_d2d("out/Phi.dat", Phi, 2*N, 2*P);
      break;
  }

  free_2d(F);
  free_2d(Psi);
  free_2d(Phi);
}


/* ======================= */
/*  ROM Matrix Generators  */
/* ======================= */
/* forcing matrix (N-by-M) */
void forcing_matrix(double **F) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      F[i][j] = actuator(ITOX(i)-Aloc[j]);
    } // j end
  } // i end
}

/* Jacobian (N-by-N) */
void benney_jacobian(double **J) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      J[i][j] = 0.0;
    } // j end
  } // i end

  double c0 = -1/(3*CA) * (1/(DX*DX*DX*DX));
  double c1 = 1/DX + (2/tan(THETA)/3 - 8*RE/15) * (1/(DX*DX)) + 1/(3*CA) * (4/(DX*DX*DX*DX));
  double c2 = (2/tan(THETA)/3 - 8*RE/15) * (-2/(DX*DX)) - 1/(3*CA) * (6/(DX*DX*DX*DX));
  double c3 = -1/DX + (2/tan(THETA)/3 - 8*RE/15) * (1/(DX*DX)) + 1/(3*CA) * (4/(DX*DX*DX*DX));
  double c4 = -1/(3*CA) * (1/(DX*DX*DX*DX));
  for (int i = 2; i < N-2; i++) {
    J[i][i-2] = c0;
    J[i][i-1] = c1;
    J[i][i] = c2;
    J[i][i+1] = c3;
    J[i][i+2] = c4;
  } // i end

  J[0][N-2] = c0;
  J[0][N-1] = c1;
  J[0][0] = c2;
  J[0][1] = c3;
  J[0][2] = c4;

  J[1][N-1] = c0;
  J[1][0] = c1;
  J[1][1] = c2;
  J[1][1+1] = c3;
  J[1][1+2] = c4;

  J[N-2][N-4] = c0;
  J[N-2][N-3] = c1;
  J[N-2][N-2] = c2;
  J[N-2][N-1] = c3;
  J[N-2][0] = c4;

  J[N-1][N-3] = c0;
  J[N-1][N-2] = c1;
  J[N-1][N-1] = c2;
  J[N-1][0] = c3;
  J[N-1][1] = c4;
}

/* Actuator (N-by-M) */
void benney_actuator(double **Psi) {
  /* forcing matrix */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  /* actuator matrix */
  for (int i = 1; i < N-1; i++) {
    for (int j = 0; j < M; j++) {
      Psi[i][j] = F[i][j] + (RE/(3*DX)) * (F[i+1][j]-F[i-1][j]);
    } // j end
  } // i end
  for (int j = 0; j < M; j++) {
    Psi[0][j] = F[0][j] + (RE/(3*DX)) * (F[1][j]-F[N-1][j]);
    Psi[N-1][j] = F[N-1][j] + (RE/(3*DX)) * (F[0][j]-F[N-2][j]);
  } // j end

  free_2d(F);
}

/* (the transpose of the) Observer (N-by-P) */
void benney_observer(double **Phi) {
  /* P < N use approximate delta-functions */
  if (P != N) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < P; j++) {
        Phi[i][j] = DX*actuator(ITOX(i)-Oloc[j]);
      } // j end
    } // i end
  }

  /* if P == N then just pass all the information */
  else {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        Phi[i][j] = 0.0;
      } // j end
      Phi[i][i] = 1.0;
    } // i end
  }
}

/* Jacobian (2N-by-2N) */
void wr_jacobian(double **J) {
  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < 2*N; j++) {
      J[i][j] = 0.0;
    } // j end
  } // i end

  double c0, c1, c2, c3, c4;

  /* top right */
  c1 = 1/(2*DX);
  c3 = -1/(2*DX);
  for (int i = 1; i < N-1; i++) {
    J[i][N+i-1] = c1;
    J[i][N+i+1] = c3;
  } // i end

  J[0][N+N-1] = c1;
  J[0][N+1] = c3;

  J[N-1][N+N-2] = c1;
  J[N-1][N+0] = c3;

  /* bottom left */
  c0 = 1/(3*CA) * (-1/(2*DX*DX*DX));
  c1 = (8*RE/35 - 2/tan(THETA)/3) * (-1/(2*DX)) + 1/(3*CA) * (1/(DX*DX*DX));
  c2 = 2;
  c3 = (8*RE/35 - 2/tan(THETA)/3) * (1/(2*DX)) + 1/(3*CA) * (-1/(DX*DX*DX));
  c4 = 1/(3*CA) * (1/(2*DX*DX*DX));
  for (int i = 2; i < N-2; i++) {
    J[N+i][i-2] = c0;
    J[N+i][i-1] = c1;
    J[N+i][i] = c2;
    J[N+i][i+1] = c3;
    J[N+i][i+2] = c4;
  } // i end

  J[N+0][N-2] = c0;
  J[N+0][N-1] = c1;
  J[N+0][0] = c2;
  J[N+0][1] = c3;
  J[N+0][2] = c4;

  J[N+1][N-1] = c0;
  J[N+1][0] = c1;
  J[N+1][1] = c2;
  J[N+1][1+1] = c3;
  J[N+1][1+2] = c4;

  J[N+N-2][N-4] = c0;
  J[N+N-2][N-3] = c1;
  J[N+N-2][N-2] = c2;
  J[N+N-2][N-1] = c3;
  J[N+N-2][0] = c4;

  J[N+N-1][N-3] = c0;
  J[N+N-1][N-2] = c1;
  J[N+N-1][N-1] = c2;
  J[N+N-1][0] = c3;
  J[N+N-1][1] = c4;

  /* bottom right */
  c1 = 6*RE/105 * 1/(2*DX);
  c2 = -1;
  c3 = 6*RE/105 * -1/(2*DX);
  for (int i = 1; i < N-1; i++) {
    J[N+i][N+i-1] = c1;
    J[N+i][N+i] = c2;
    J[N+i][N+i+1] = c3;
  } // i end

  J[N+0][N+N-1] = c1;
  J[N+0][N+0] = c2;
  J[N+0][N+1] = c3;

  J[N+N-1][N+N-2] = c1;
  J[N+N-1][N+N-1] = c2;
  J[N+N-1][N+0] = c3;
}

/* Actuator (2N-by-M) */
void wr_actuator(double **Psi) {
  /* forcing matrix (TODO: see if rotating this makes it faster) */
  double **F = malloc_f2d(N, M);
  forcing_matrix(F);

  /* actuator matrix */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      Psi[i][j] = F[i][j];
      Psi[i+N][j] = 2*RE/15 * F[i][j];
    } // j end
  } // i end

  free_2d(F);
}

/* (the transpose of the) Observer (2N-by-2P) */
void wr_observer(double **Phi) {
  /* P < N use approximate delta-functions */
  if (P != N) {
    for (int i = 0; i < 2*N; i++) {
      for (int j = 0; j < 2*P; j++) {
        Phi[i][j] = 0.0;
      } // j end
    } // i end

    /* observe interfacial height */
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < P; j++) {
        Phi[i][j] = DX*actuator(ITOX(i)-Oloc[j]);
      } // j end
    } // i end

    /* observe flux */ // TODO: should the q = 2/3 h assumption be here?
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < P; j++) {
        Phi[N+i][P+j] = DX*actuator(ITOX(i)-Oloc[j]);
      } // j end
    } // i end
  }

  /* if P == N then just pass all the information */
  else {
    for (int i = 0; i < 2*N; i++) {
      for (int j = 0; j < 2*N; j++) {
        Phi[i][j] = 0.0;
      } // j end
      Phi[i][i] = 1.0;
    } // i end
  }
}


#endif
