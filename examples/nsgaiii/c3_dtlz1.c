#include "optimize.h"
#include <stdio.h>

/* Defaults */
#define NOBJ 3
#define DIM  7

static double g (double *x, int k)
{
    int i;
    double s = 0;
    for (i=0; i<k; i++) {
        s += pow (x [i] - 0.5, 2) - cos (20 * M_PI * (x [i] - 0.5));
    }
    return 100 * (k + s);
}

static void f (double *x, int nx, double *y, int ny)
{
    int i, j;
    int nobj = ny - NOBJ;
    int k = (nx - nobj + 1);
    double gv;
    assert (nobj > 0);
    assert (k > 0);
    assert (nobj >= NOBJ);
    gv = g (x + nx - k, k);
    for (i=0; i<nobj; i++) {
        double p = 0.5;
        for (j=0; j<nobj-i-1; j++) {
            p *= x [j];
        }
        if (j < nobj-1) {
            p *= (1 - x [j]);
        }
        y [i] = p * (1 + gv);
    }
    for (j=0; j<NOBJ; j++) {
        double c = y [j];
        for (i=0; i<nobj; i++) {
            if (i != j) {
                c += y [i] / 0.5;
            }
        }
        y [nobj + j] = 1 - c;
    }
}

struct multi_problem c3_dtlz1 =
{ .dimension      = DIM
, .n_obj          = NOBJ
, .nconstraint    = NOBJ
, .lower          = 0
, .upper          = 1
, .popsize        = 0
, .generations    = 750
, .f              = f
, .name           = "C3-DTLZ1"
};
