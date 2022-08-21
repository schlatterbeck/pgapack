#include "optimize.h"
#include <stdio.h>

static double f1 (double *x)
{
    return 1 - exp (-4 * x [0]) * pow (sin (6 * M_PI * x [0]), 6);
}

static double g (double *x)
{
    int i;
    double s = 0;
    for (i=1; i<zdt6.dimension; i++) {
        s += x [i];
    }
    return 1. + 9 * pow (s / 9., 0.25);
}

static double f2 (double *x)
{
    double rg = g (x);
    return rg * (1 - pow (f1 (x) / rg, 2));
}

struct multi_problem zdt6 =
{ .dimension      = 10
, .nfunc          = 2
, .nconstraint    = 0
, .lower          = (double []){ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
, .upper          = (double []){ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
, .enforce_bounds = 1
, .generations    = 500
, .crossover_prob = 0.000001
, .f              = { &f1, &f2 }
, .name           = "Zitzler et. al. (ZDT6)"
};
