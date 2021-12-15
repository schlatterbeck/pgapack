#include "constraint.h"
#include <stdio.h>

static double f1 (double *x)
{
    return x [0];
}

static double f2 (double *x)
{
    return (1 + x [1]) / x [0];
}

static double g1 (double *x)
{
    return 6 - (x [1] + 9 * x [0]);
}

static double g2 (double *x)
{
    return 1 + x [1] - 9 * x [0];
}

struct multi_problem constr =
{ .dimension      = 2
, .nfunc          = 4
, .nconstraint    = 2
, .lower          = (double []){ 0.1, 0 }
, .upper          = (double []){ 1,   5 }
, .enforce_bounds = 1
, .popsize        = 10
, .generations    = 10
, .f              = { &f1, &f2, &g1, &g2 }
, .name           = "CONSTR"
};
