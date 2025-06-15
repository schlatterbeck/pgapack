#include "optimize.h"
#include <stdio.h>

#define NOBJ  3
#define DIM  12

static double g (double *x, int k)
{
    int i;
    double s = 0;
    for (i=0; i<k; i++) {
        s += pow (x [i] - 0.5, 2);
    }
    return s;
}

static void f (double *x, int nx, double *y, int ny)
{
    int i, j;
    int nobj = ny - NOBJ;
    int k = nx - nobj + 1;
    double gv;
    assert (nobj > 0);
    assert (k > 0);
    assert (nobj >= NOBJ);
    gv = g (x + nx - k, k);
    for (i=0; i<nobj; i++) {
        double p = 1;
        for (j=0; j<nobj-i-1; j++) {
            p *= cos (pow (x [j], 100) * M_PI / 2.0);
        }
        if (j < nobj-1) {
            p *= sin (pow (x [j], 100) * M_PI / 2.0);
        }
        y [i] = p * (1 + gv);
    }
    for (j=0; j<NOBJ; j++) {
        double c = pow (y [j], 2) / 4;
        for (i=0; i<nobj; i++) {
            if (i != j) {
                c += pow (y [i], 2);
            }
        }
        y [nobj + j] = 1 - c;
    }
}

struct multi_problem c3_dtlz4 =
{ .dimension      = DIM
, .n_obj          = NOBJ
, .nconstraint    = NOBJ
, .lower          = 0
, .upper          = 1
, .popsize        = 0
, .generations    = 750
, .f              = f
, .name           = "C3-DTLZ4"
};
