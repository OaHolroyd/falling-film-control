#ifndef PARAMS_H
#define PARAMS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h> // seems to be required for jsmn

#include "jsmn.h"

/* ========================================================================== */
/*   PHYSICAL PARAMETERS                                                      */
/* ========================================================================== */
/*  All of the parameters are read in with SI units and non-dimensionalised   */
/*  using the following scales:                                               */
/*    velocity - US                                                           */
/*    length   - h0                                                           */
/*    pressure - mu_l US / h0                                                 */
/*  This is the scaling used in Thompson 2016 and Cimpeanu 2021.              */

/* Domain parameters */
double H0 = 100.0e-6; // film thickness
double LX = 64.0; // dimensionless domain length (Lx/h0)
double LY = 8.0; // dimensionless domain height (Ly/h0)
double THETA = 1.0; // inclination angle
double TMAX = 1.0; // final time
double T0 = 0.0; // initial time (0 or an integer corresponding to a dump file)

/* Physical parameters */
double RHO_L = 1000.0;
double RHO_G = 1.0;  // densities
double MU_L = 1.0e-3;
double MU_G = 1.0e-5; // (dynamic) viscosities
double GAMMA = 0.1;  // surface tension
double GRAV = 10;  // acceleration due to gravity

/* Solver parameters */
int LEVEL = 8; // maximum refinement level
double DTOUT = 1.0; // output step


/* ========================================================================== */
/*   CONTROL PARAMETERS                                                       */
/* ========================================================================== */
int C_M = 5; // number of controls
int C_P = 5; // number of observers
double C_START = 100.0; // control start time
double C_W = 0.01; // control width parameter
double C_ALPHA = 1.0; // control strength
double C_PHI = 1.0; // observer/control displacement

double C_norm; // control normaliser
double *C_loc; // control locations
double *C_mag; // current control magnitudes


/* checks if a token and a string match */
static int jsoneq(const char *json, jsmntok_t *tok, const char *s) {
  if (tok->type == JSMN_STRING && (int)strlen(s) == tok->end - tok->start &&
      strncmp(json + tok->start, s, tok->end - tok->start) == 0) {
    return 0;
  }
  return -1;
}


/* read params from file */
int read_params(char *fname) {
  /* fname should not be NULL */
  if (!fname) { return -1; }

  /* attempt to open file */
  FILE *fp = fopen(fname, "r");
  if (!fp) { return -2; }

  /* try and read entire file into the string */
  fseek(fp, 0, SEEK_END);
  long len = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  char *s = malloc(len+1);
  if (!s) { return -3; }
  fread(s, 1, len, fp);
  s[len] = '\0';
  fclose(fp);

  /* initialise parser */
  jsmn_parser p;
  jsmntok_t t[128];
  jsmn_init(&p);
  int r = jsmn_parse(&p, s, strlen(s), t, 128);

  /* check for simple matches */
  int n, m;
  for (int i = 0; i < r-1; i++) {
    /* read the parameters in chunks */
    if (jsoneq(s, &t[i], "DOMAIN") == 0) {
      n = 6;
      m = n;
      for (int j = i+2; j < i+2+2*m; j++) {
        if (jsoneq(s, &t[j], "h0") == 0) {
          j++;
          H0 = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "Lx") == 0) {
          j++;
          LX = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "Ly") == 0) {
          j++;
          LY = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "theta") == 0) {
          j++;
          THETA = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "tmax") == 0) {
          j++;
          TMAX = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "t0") == 0) {
          j++;
          T0 = strtod(s+t[j].start, NULL);
          n--;
        } else {
          fprintf(stderr, "bad token\n");
        }
      } // j end
      if (n) {
        fprintf(stderr, "ERROR: missing domain parameter\n");
        abort();
      }
    }

    else if (jsoneq(s, &t[i], "PHYSCIAL") == 0) {
      n = 6;
      m = n;
      for (int j = i+2; j < i+2+2*m; j++) {
        if (jsoneq(s, &t[j], "rho_l") == 0) {
          j++;
          RHO_L = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "rho_g") == 0) {
          j++;
          RHO_G = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "mu_l") == 0) {
          j++;
          MU_L = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "mu_g") == 0) {
          j++;
          MU_G = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "gamma") == 0) {
          j++;
          GAMMA = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "grav") == 0) {
          j++;
          GRAV = strtod(s+t[j].start, NULL);
          n--;
        } else {
          fprintf(stderr, "bad token\n");
        }
      } // j end
      if (n) {
        fprintf(stderr, "ERROR: missing physical parameter\n");
        abort();
      }
    }

    else if (jsoneq(s, &t[i], "SOLVER") == 0) {
      n = 2;
      m = n;
      for (int j = i+2; j < i+2+2*m; j++) {
        if (jsoneq(s, &t[j], "level") == 0) {
          j++;
          LEVEL = atoi(s+t[j].start);
          n--;
        } else if (jsoneq(s, &t[j], "dtout") == 0) {
          j++;
          DTOUT = strtod(s+t[j].start, NULL);
          n--;
        } else {
          fprintf(stderr, "bad token\n");
        }
      } // j end
      if (n) {
        fprintf(stderr, "ERROR: missing solver parameter\n");
        abort();
      }
    }

    else if (jsoneq(s, &t[i], "CONTROL") == 0) {
      n = 6;
      m = n;
      for (int j = i+2; j < i+2+2*m; j++) {
        if (jsoneq(s, &t[j], "M") == 0) {
          j++;
          C_M = atoi(s+t[j].start);
          n--;
        } else if (jsoneq(s, &t[j], "P") == 0) {
          j++;
          C_P = atoi(s+t[j].start);
          n--;
        } else if (jsoneq(s, &t[j], "start") == 0) {
          j++;
          C_START = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "width") == 0) {
          j++;
          C_W = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "alpha") == 0) {
          j++;
          C_ALPHA = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "phi") == 0) {
          j++;
          C_PHI = strtod(s+t[j].start, NULL);
          n--;
        } else {
          fprintf(stderr, "bad token\n");
        }
      } // j end
      if (n) {
        fprintf(stderr, "ERROR: missing control parameter\n");
        abort();
      }
    }
  } // i end

  return 0;
}

#endif
