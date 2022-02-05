#include "optimize.h"
#include <stdio.h>

#define NOBJ  3
#define DIM  12
#define K (DIM - NOBJ + 1)

static double g (double *x)
{
    int i;
    double s = 0;
    for (i=0; i<K; i++) {
        s += pow (x [i] - 0.5, 2);
    }
    return s;
}

static void f (double *x, double *y)
{
    int i, j;
    double gv = g (x + DIM - K);
    for (i=0; i<NOBJ; i++) {
        double p = 1;
        for (j=0; j<NOBJ-i-1; j++) {
            p *= cos (pow (x [j], 100) * M_PI / 2.0);
        }
        if (j < NOBJ-1) {
            p *= sin (pow (x [j], 100) * M_PI / 2.0);
        }
        y [i] = p * (1 + gv);
    }
    for (j=0; j<NOBJ; j++) {
        double c = pow (y [j], 2) / 4;
        for (i=0; i<NOBJ; i++) {
            if (i != j) {
                c += pow (y [i], 2);
            }
        }
        y [NOBJ + j] = 1 - c;
    }
}

struct multi_problem c3_dtlz4 =
{ .dimension      = DIM
, .nfunc          = NOBJ * 2
, .nconstraint    = NOBJ
, .lower          = (double []){ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
, .upper          = (double []){ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
, .popsize        = 0
, .generations    = 750
, .f              = f
, .name           = "C3-DTLZ4"
};
