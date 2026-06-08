#include "optimize.h"
#include <stdio.h>

#define NOBJ  3
#define DIM  22

static double g (double *x, unsigned int k)
{
    unsigned int i;
    double s = 0;
    for (i=0; i<k; i++) {
        s += x [i];
    }
    return 1 + s * 9.0 / k;
}

static void f (double *x, unsigned int nx, double *y, unsigned int ny)
{
    unsigned int i;
    unsigned int k = nx - ny + 1;
    double gv;
    double h;
    assert (k > 0);
    gv = g (x + nx - k, k);
    h = ny;
    for (i=0; i<ny - 1; i++) {
        y [i] = x [i];
        h -= y [i] / (1 + gv) * (1 + sin (3 * M_PI * y [i]));
    }
    y [ny - 1] = 1 + gv * h;
}

struct multi_problem dtlz7 =
{ .dimension      = DIM
, .n_obj          = NOBJ
, .nconstraint    = 0
, .lower          = 0
, .upper          = 1
, .popsize        = 500
, .f              = f
, .name           = "DTLZ7"
};
