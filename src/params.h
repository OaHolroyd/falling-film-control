#ifndef PARAMS_H
#define PARAMS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h> // seems to be required for jsmn

#include "c-utils.h"
#include "parallel.h"
#include "control.h"
#include "jsmn.h"


/* ========================================================================== */
/*   OUTPUT PARAMETERS                                                        */
/* ========================================================================== */
#define NOUT (1<<(LEVEL-1)) // output resolution
#define LOG_STEP 10 // log every LOG_STEP steps
#define OUTPUT_DAT 1 // whether to output data
#define DUMP 100 // how often to dump (for restarting)


/* ========================================================================== */
/*   DIMENSIONLESS NUMBERS                                                    */
/* ========================================================================== */
#define US ((RHO_L*GRAV*sin(THETA)*H0*H0)/(2*MU_L)) // Nusselt surface velocity
#define RE ((RHO_L*US*H0)/MU_L) // Reynolds number
#define CA ((MU_L*US)/GAMMA) // capillary number


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
int OUTPUT = 0; // output type
#define N (1<<LEVEL) // max refinement grid count
#define DX (LX/N) // max refinement grid spacing


/* ========================================================================== */
/*   CONTROL PARAMETERS                                                       */
/* ========================================================================== */
int C_M = 5; // number of controls
int C_P = 5; // number of observers
double C_START = 100.0; // control start time
double C_W = 0.01; // control width parameter
double C_ALPHA = 1.0; // control strength
double C_DEL = 1.0; // observer/control displacement
double C_MU = 0.1; // control cost parameter
rom_t C_ROM = BENNEY; // reduced order model
control_t C_STRAT = STATIC; // control strategy


/* ========================================================================== */
/*   FUNCTION DEFINITIONS                                                     */
/* ========================================================================== */
/* checks if a token and a string match */
static int jsoneq(const char *json, jsmntok_t *tok, const char *s) {
  if (tok->type == JSMN_STRING && (int)strlen(s) == tok->end - tok->start &&
      strncmp(json + tok->start, s, tok->end - tok->start) == 0) {
    return 0;
  }
  return -1;
}

/* read params from file */
int read_params(const char *fname) {
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
  if (fread(s, 1, len, fp) != len) { ABORT("failed to load json"); }
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
        ABORT("missing domain parameter");
      }
    }

    else if (jsoneq(s, &t[i], "PHYSICAL") == 0) {
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
        ABORT("missing physical parameter");
      }
    }

    else if (jsoneq(s, &t[i], "SOLVER") == 0) {
      n = 3;
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
        } else if (jsoneq(s, &t[j], "output") == 0) {
          j++;
          OUTPUT = atoi(s+t[j].start);
          n--;
        } else {
          fprintf(stderr, "bad token\n");
        }
      } // j end
      if (n) {
        ABORT("missing solver parameter");
      }
    }

    else if (jsoneq(s, &t[i], "CONTROL") == 0) {
      n = 9;
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
        } else if (jsoneq(s, &t[j], "del") == 0) {
          j++;
          C_DEL = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "mu") == 0) {
          j++;
          C_MU = strtod(s+t[j].start, NULL);
          n--;
        } else if (jsoneq(s, &t[j], "rom") == 0) {
          j++;
          if (!strncmp(s+t[j].start, "benney", t[j].end-t[j].start)) {
            C_ROM = BENNEY;
          } else if (!strncmp(s+t[j].start, "wr", t[j].end-t[j].start)) {
            C_ROM = WR;
          } else {
            ABORT("invalid ROM type");
          }
          n--;
        } else if (jsoneq(s, &t[j], "strategy") == 0) {
          j++;
          if (!strncmp(s+t[j].start, "pair", t[j].end-t[j].start)) {
            C_STRAT = PAIR;
          } else if (!strncmp(s+t[j].start, "lqr", t[j].end-t[j].start)) {
            C_STRAT = LQR;
          } else if (!strncmp(s+t[j].start, "static", t[j].end-t[j].start)) {
            C_STRAT = STATIC;
          } else if (!strncmp(s+t[j].start, "dynamic", t[j].end-t[j].start)) {
            C_STRAT = DYNAMIC;
          } else {
            ABORT("invalid control strategy");
          }
          n--;
        } else {
          fprintf(stderr, "bad token\n");
        }
      } // j end
      if (n) {
        ABORT("missing control parameter");
      }
    }
  } // i end

  free(s);
  return 0;
}


#endif
