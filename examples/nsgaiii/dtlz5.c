#include "optimize.h"
#include <stdio.h>

#define NOBJ  3
#define DIM  12

static double g (double *x, unsigned int k)
{
    unsigned int i;
    double s = 0;
    for (i=0; i<k; i++) {
        s += pow (x [i] - 0.5, 2);
    }
    return s;
}

static double theta (double *x, unsigned int i, double gv)
{
    if (i == 0) {
        return x [i];
    }
    /* from paper with factor pi/2 included again
     * return M_PI / (4 * (1 + gv)) * (1 + 2 * gv * x [i]);
     */
    /* Seems the factor pi/2 is included twice, remove it here */
    return 1 / (2 * (1 + gv)) * (1 + 2 * gv * x [i]);
}

static void f (double *x, unsigned int nx, double *y, unsigned int ny)
{
    unsigned int i, j;
    unsigned int k = nx - ny + 1;
    double gv;
    assert (k > 0);
    gv = g (x + nx - k, k);
    for (i=0; i<ny; i++) {
        double p = 1;
        for (j=0; j<ny-i-1; j++) {
            double th = theta (x, j, gv);
            p *= cos (th * M_PI / 2.0);
        }
        if (j < ny-1) {
            double th = theta (x, j, gv);
            p *= sin (th * M_PI / 2.0);
        }
        y [i] = p * (1 + gv);
    }
}

struct multi_problem dtlz5 =
{ .dimension      = DIM
, .n_obj          = NOBJ
, .nconstraint    = 0
, .lower          = 0
, .upper          = 1
, .popsize        = 50
, .f              = f
, .name           = "DTLZ5"
};
