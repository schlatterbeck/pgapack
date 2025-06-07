#include "optimize.h"
#include <stdio.h>

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
    int k = (nx - ny + 1);
    double gv;
    assert (k > 0);
    gv = g (x + nx - k, k);
    for (i=0; i<ny; i++) {
        double p = 0.5;
        for (j=0; j<ny-i-1; j++) {
            p *= x [j];
        }
        if (j < ny-1) {
            p *= (1 - x [j]);
        }
        y [i] = p * (1 + gv) * pow (10, i);
    }
}

struct multi_problem scaled_dtlz1 =
{ .dimension      = DIM
, .n_obj          = NOBJ
, .nconstraint    = 0
, .lower          = 0
, .upper          = 1
, .popsize        = 0
, .generations    = 0
, .f              = f
, .name           = "Scaled DTLZ1"
};
