#include "constraint.h"
#include <stdio.h>

#define NOBJ  3
#define DIM  12
#define K (DIM - NOBJ + 1)

/* Multipliers for constraint sum2, starting with NOBJ = 3 */
static double c1_mul [] = {9, 12.5, 12.5, 15, 15};

static double g (double *x)
{
    int i;
    double s = 0;
    for (i=0; i<K; i++) {
        s += pow (x [i] - 0.5, 2) - cos (20 * M_PI * (x [i] - 0.5));
    }
    return 100 * (K + s);
}

static void f (double *x, double *y)
{
    int i, j;
    double gv = g (x + DIM - K);
    double s = 0;
    for (i=0; i<NOBJ; i++) {
        double p = 1;
        for (j=0; j<NOBJ-i-1; j++) {
            p *= cos (x [j] * M_PI / 2.0);
        }
        if (j < NOBJ-1) {
            p *= sin (x [j] * M_PI / 2.0);
        }
        y [i] = p * (1 + gv);
        s  += pow (y [i], 2);
    }
    y [NOBJ] = -(s - 16) * (s - c1_mul [NOBJ - 3]);
}

struct multi_problem c1_dtlz3 =
{ .dimension      = DIM
, .nfunc          = NOBJ + 1
, .nconstraint    = 1
, .lower          = (double []){ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
, .upper          = (double []){ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
, .popsize        = 0
, .generations    = 1000
, .f              = f
, .name           = "C1-DTLZ3"
};
