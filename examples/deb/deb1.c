#include "constraint.h"
#include <stdio.h>

static double f (double *x)
{
    return pow (pow (x [0], 2) + x [1] - 11, 2)
         + pow (x [0] + pow (x [1], 2) - 7, 2);
}

static double g1 (double *x)
{
    return pow (x [0] - 0.05, 2) + pow (x [1] - 2.5, 2) - 4.84;
}

static double g2 (double *x)
{
    return 4.84 - pow (x [0], 2) - pow (x [1] - 2.5, 2);
}

struct constrained_problem deb_1 =
{ .dimension      = 2
, .nfunc          = 3
, .lower          = (double []){ 0,  0 }
, .upper          = (double []){ 6,  6 }
, .enforce_bounds = 1
, .popsize        = 4
, .generations    = 150
, .f              = { &f, &g1, &g2 }
};
