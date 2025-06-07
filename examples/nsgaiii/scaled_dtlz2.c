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
    int k = nx - ny + 1;
    double gv;
    assert (k > 0);
    gv = g (x + nx - k, k);
    for (i=0; i<ny; i++) {
        double p = 1;
        for (j=0; j<ny-i-1; j++) {
            p *= cos (x [j] * M_PI / 2.0);
        }
        if (j < ny-1) {
            p *= sin (x [j] * M_PI / 2.0);
        }
        y [i] = p * (1 + gv) * pow (10, i);
    }
}

struct multi_problem scaled_dtlz2 =
{ .dimension      = DIM
, .n_obj          = NOBJ
, .nconstraint    = 0
, .lower          = 0
, .upper          = 1
, .popsize        = 0
, .generations    = 250
, .f              = f
, .name           = "Scaled DTLZ2"
};
