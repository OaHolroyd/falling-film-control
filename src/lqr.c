#include "lqr.h"

#include <math.h>

#include <lapacke.h>

#include <stdio.h> // TODO: remove


#define COMPLEX double complex


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* computes the optimal gain matrix K for the system (A, B) with quadratic cost
   given by (u, v). n is the system size and m is the control dimension. */
   // TODO: add a lqr_work variant for repeated computations
int dlqr(double **A, double **B, double u, double v, int n, int m, double **K) {
  int i, j, k;
  double c = -1/v;

  /* generate the Hamiltonian */
  COMPLEX *H = malloc(2*n*2*n*sizeof(COMPLEX));

  k = 0;
  for (i = 0; i < n; i++) { // top left
    for (j = 0; j < n; j++) {
      H[2*n*i+j] = A[i][j];
    } // j end
  } // i end

  for (i = 0; i < n; i++) { /* top right (-PHI * V^-1 * PSI^T) */
    for (j = i; j < n; j++) {
      H[2*n*i+(j+n)] = 0.0;
      for (k = 0; k < m; k++) {
        H[2*n*i+(j+n)] += B[i][k] * B[j][k];
      } // k end
      H[2*n*i+(j+n)] *= c;
      H[2*n*j+(i+n)] = H[2*n*i+(j+n)]; // use symmetry of PHI * PSI^T
    } // j end
  } // i end

  for (i = 0; i < n; i++) { // bottom left
    for (j = 0; j < n; j++) {
      H[2*n*(i+n)+j] = 0.0;
    } // j end
    H[2*n*(i+n)+i] = -u;
  } // i end

  for (i = 0; i < n; i++) { // bottom right
    for (j = 0; j < n; j++) {
      H[2*n*(i+n)+(j+n)] = -A[j][i];
    } // j end
  } // i end


  // FILE *fp = fopen("out/H.dat", "w");
  // for (int i = 0; i < 2*n; i++) {
  //   for (int j = 0; j < 2*n; j++) {
  //     fprintf(fp, "%20.10lf ", creal(H[2*n*i+j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute eigenvalues/vectors */
  COMPLEX *w = malloc(2*n*sizeof(COMPLEX));
  COMPLEX *vr = malloc(2*n*2*n*sizeof(COMPLEX));
  int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', 2*n, H,
                           2*n, w, NULL, 2*n, vr, 2*n);


  /* extract stable eigenvectors */ // TODO: remove the need for E
  COMPLEX *E = malloc(2*n*2*n*sizeof(COMPLEX));
  COMPLEX *Q = malloc(n*n*sizeof(COMPLEX));
  COMPLEX *P = malloc(n*n*sizeof(COMPLEX));
  k = 0;
  for (i = 0; i < 2*n; i++) {
    if (creal(w[i]) < 0) {
      for (j = 0; j < 2*n; j++) {
        E[2*n*j+k] = vr[2*n*j+i];
      } // j end
      k++;
    }
  } // i end

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      Q[n*i+j] = E[2*n*i+j];
      P[n*i+j] = E[2*n*(i+n)+j];
    } // j end
  } // i end


  /* compute P */
  int *ipiv = malloc(n*sizeof(int));
  info = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, n, Q, n, ipiv, P, n); // use COL_MAJOR to solve XA=B


  /* compute K */
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      K[i][j] = 0.0;
      for (k = 0; k < n; k++) {
        K[i][j] += B[k][i] * creal(P[k*n+j]);
      } // k end
      K[i][j] *= c;
    } // j end
  } // i end

  // fp = fopen("out/K.dat", "w");
  // for (int i = 0; i < m; i++) {
  //   for (int j = 0; j < n; j++) {
  //     fprintf(fp, "%20.10lf ", K[i][j]);
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/A.dat", "w");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < n; j++) {
  //     fprintf(fp, "%20.10lf ", A[i][j]);
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/B.dat", "w");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < m; j++) {
  //     fprintf(fp, "%20.10lf ", B[i][j]);
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* free workspace */
  free(H);
  free(E);
  free(Q);
  free(P);
  free(w);
  free(vr);
  free(ipiv);

  return info;
}

int zlqr(COMPLEX **A, COMPLEX **B, double u, double v, int n, int m, COMPLEX **K) {
  int i, j, k;
  double c = -1/v;

  /* generate the Hamiltonian */
  COMPLEX *H = malloc(2*n*2*n*sizeof(COMPLEX));

  k = 0;
  for (i = 0; i < n; i++) { // top left
    for (j = 0; j < n; j++) {
      H[2*n*i+j] = A[i][j];
    } // j end
  } // i end

  for (i = 0; i < n; i++) { /* top right (-PHI * V^-1 * PSI^T) */
    for (j = i; j < n; j++) {
      H[2*n*i+(j+n)] = 0.0;
      for (k = 0; k < m; k++) {
        H[2*n*i+(j+n)] += B[i][k] * conj(B[j][k]);
      } // k end
      H[2*n*i+(j+n)] *= c;
      H[2*n*j+(i+n)] = H[2*n*i+(j+n)]; // use symmetry of PHI * PSI^T
    } // j end
  } // i end

  for (i = 0; i < n; i++) { // bottom left
    for (j = 0; j < n; j++) {
      H[2*n*(i+n)+j] = 0.0;
    } // j end
    H[2*n*(i+n)+i] = -u;
  } // i end

  for (i = 0; i < n; i++) { // bottom right
    for (j = 0; j < n; j++) {
      H[2*n*(i+n)+(j+n)] = -conj(A[j][i]);
    } // j end
  } // i end


  // FILE *fp = fopen("out/Hr.dat", "w");
  // for (int i = 0; i < 2*n; i++) {
  //   for (int j = 0; j < 2*n; j++) {
  //     fprintf(fp, "%20.10lf ", creal(H[2*n*i+j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/Hi.dat", "w");
  // for (int i = 0; i < 2*n; i++) {
  //   for (int j = 0; j < 2*n; j++) {
  //     fprintf(fp, "%20.10lf ", cimag(H[2*n*i+j]));
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* compute eigenvalues/vectors */
  COMPLEX *w = malloc(2*n*sizeof(COMPLEX));
  COMPLEX *vr = malloc(2*n*2*n*sizeof(COMPLEX));
  int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', 2*n, H,
                           2*n, w, NULL, 2*n, vr, 2*n);


  /* extract stable eigenvectors */ // TODO: remove the need for E
  COMPLEX *E = malloc(2*n*2*n*sizeof(COMPLEX));
  COMPLEX *Q = malloc(n*n*sizeof(COMPLEX));
  COMPLEX *P = malloc(n*n*sizeof(COMPLEX));
  k = 0;
  for (i = 0; i < 2*n; i++) {
    if (creal(w[i]) < 0) {
      for (j = 0; j < 2*n; j++) {
        E[2*n*j+k] = vr[2*n*j+i];
      } // j end
      k++;
    }
  } // i end

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      Q[n*i+j] = E[2*n*i+j];
      P[n*i+j] = E[2*n*(i+n)+j];
    } // j end
  } // i end


  /* compute P */
  int *ipiv = malloc(n*sizeof(int));
  info = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, n, Q, n, ipiv, P, n); // use COL_MAJOR to solve XA=B


  /* compute K */
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      K[i][j] = 0.0;
      for (k = 0; k < n; k++) {
        K[i][j] += conj(B[k][i]) * creal(P[k*n+j]);
      } // k end
      K[i][j] *= c;
    } // j end
  } // i end

  // fp = fopen("out/K.dat", "w");
  // for (int i = 0; i < m; i++) {
  //   for (int j = 0; j < n; j++) {
  //     fprintf(fp, "%20.10lf ", K[i][j]);
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/A.dat", "w");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < n; j++) {
  //     fprintf(fp, "%20.10lf ", A[i][j]);
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);

  // fp = fopen("out/B.dat", "w");
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < m; j++) {
  //     fprintf(fp, "%20.10lf ", B[i][j]);
  //   } // j end
  //   fprintf(fp, "\n");
  // } // i end
  // fclose(fp);


  /* free workspace */
  free(H);
  free(E);
  free(Q);
  free(P);
  free(w);
  free(vr);
  free(ipiv);

  return info;
}


// TODO: make this better
int LQR_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb) {
  return LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
}

int LQR_zgesv(int n, int nrhs, COMPLEX *a, int lda, int *ipiv, COMPLEX *b, int ldb) {
  return LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
}