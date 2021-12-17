#include "constraint.h"
#include <stdio.h>
#include <string.h>

#define DIMENSION 10

#define PHI (45. / 180. * M_PI)
static double y [DIMENSION];
static double *last_x = NULL; /* the last x seen for which y is valid */

static double rotmatrix [DIMENSION][DIMENSION];

static void compute_rotation (void)
{
    int done = 0;
    int i, j;
    double g1 [DIMENSION];
    double g2 [DIMENSION];
    double v  [DIMENSION][DIMENSION];
    double w  [DIMENSION][DIMENSION];
    double I  [DIMENSION][DIMENSION];

    if (done) {
        return;
    }

    for (i=0; i<DIMENSION; i++) {
        g1 [i] = g2 [i] = 0;
        for (j=0; j<DIMENSION; j++) {
            v [i][j] = w [i][j] = I [i][j] = 0;
            if (i==j) {
                I [i][j] = 1;
            }
        }
    }
    g1 [0] = g1 [1] = g2 [2] = g2 [3] = pow (0.5, 0.5);
    for (i=0; i<DIMENSION; i++) {
        for (j=0; j<DIMENSION; j++) {
            v [i][j] = (g1 [i] * g1 [j] + g2 [i] * g2 [j]) * (cos (PHI) - 1);
            w [i][j] = (g1 [i] * g2 [j] - g2 [i] * g1 [j]) * sin (PHI);
        }
    }
    for (i=0; i<DIMENSION; i++) {
        for (j=0; j<DIMENSION; j++) {
            rotmatrix [i][j] = I [i][j] + v [i][j] + w [i][j];
        }
    }
    done = 1;
}

static void rotate (double *x)
{
    int i, j;
    if (x == last_x) {
        return;
    }
    compute_rotation ();
    for (i=0; i<DIMENSION; i++) {
        y [i] = 0;
        for (j=0; j<DIMENSION; j++) {
            y [i] += rotmatrix [i][j] * x [j];
        }
    }
}

static double f1 (double *x)
{
    rotate (x);
    return y [0];
}

static double g (double *x)
{
    int i;
    double s = 0;
    rotate (x);
    for (i=1; i<rotated.dimension; i++) {
        s += pow (y [i], 2) - 10 * cos (4 * M_PI * y [i]);
    }
    return 1. + 10 * (rotated.dimension - 1) + s;
}

static double f2 (double *x)
{
    double rg = g (x);
    return rg * exp (-y [0] / rg);
}

/* Limit the y [0] range to [-0.3, 0.3] */

static double g1 (double *x)
{
    rotate (x);
    return y [0] - 0.3;
}

static double g2 (double *x)
{
    rotate (x);
    return -0.3 - y [0];
}

struct multi_problem rotated =
{ .dimension      = DIMENSION
, .nfunc          = 4
, .nconstraint    = 2
, .lower          = (double [])
                    { -.3, -.3, -.3, -.3, -.3, -.3, -.3, -.3, -.3, -.3 }
, .upper          = (double [])
                    {  .3,  .3,  .3,  .3,  .3,  .3,  .3,  .3,  .3,  .3 }
, .enforce_bounds = 1
, .popsize        = 50
, .generations    = 2500
, .f              = { &f1, &f2, &g1, &g2 }
, .name           = "Rotated Problem"
};
