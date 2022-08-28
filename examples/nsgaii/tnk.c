#include "optimize.h"
#include <stdio.h>

static double f1 (double *x)
{
    return x [0];
}

static double f2 (double *x)
{
    return x [1];
}

static double g1 (double *x)
{
    return -pow (x [0], 2) - pow (x [1], 2) + 1
           + 0.1 * cos (16 * atan (x [0] / x [1]));
}

static double g2 (double *x)
{
    return pow (x [0] - 0.5, 2) + pow (x [1] - 0.5, 2) - 0.5;
}

struct multi_problem tnk =
{ .dimension      = 2
, .nfunc          = 4
, .nconstraint    = 2
, .lower          = (double []){    0,    0 }
, .upper          = (double []){ M_PI, M_PI }
, .enforce_bounds = 1
, .f              = { &f1, &f2, &g1, &g2 }
, .name           = "Tanaka et. al. (TNK)"
};
