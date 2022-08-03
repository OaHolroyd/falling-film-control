#ifndef B_UTILS_H
#define B_UTILS_H

/* Standard headers */
#include <stdlib.h>
#include <math.h>

/* Local headers */
#include "parallel.h"
#include "params.h"

/* Basilisk headers */
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "embed.h" // embedded boundaries (supersedes mask)
#include "tension.h" // surface tension
#include "heights.h" // interfacial height


/* harmonic viscosity averaging TODO: check this vs other types */
#ifdef mu
#undef mu
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#endif

/* convert from index to location */
#define ITOX(i) (DX*(i+0.5))

/* convert from location to index */
#define XTOI(x) ((int)(x/DX - 0.5))

/* get interfacial height at a given x-coord */
#define NP 10
vector hei[]; // heights
double interfacial_height(double xp) {
  if (xp < 0) {
    xp += LX;
  } else if (xp > LX) {
    xp -= LX;
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
      dh[i] = fabs(y-yh[i]);
      return yh[i];
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


#endif
