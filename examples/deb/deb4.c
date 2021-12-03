#include "constraint.h"
#include <stdio.h>

static double f (double *x)
{
    return x [0] + x [1] + x [2];
}

static double g1 (double *x)
{
    return 0.0025 * (x [3] + x [5]) - 1;
}

static double g2 (double *x)
{
    return 0.0025 * (x [4] + x [6] - x [3]) - 1;
}

static double g3 (double *x)
{
    return 0.01 * (x [7] - x [4]) - 1;
}

static double g4 (double *x)
{
    return 833.33252 * x [3] + 100 * x [0] - x [0] * x [5] - 83333.333;
}

static double g5 (double *x)
{
    return 1250 * x [4] + x [1] * x [3] - x [1] * x [6] - 1250 * x [3];
}

static double g6 (double *x)
{
    return x [2] * x [4] + 1.25e6 - x [2] * x [7] - 2.5e3 * x [4];
}

struct constrained_problem deb_4 =
{ .dimension      = 8
, .nfunc          = 7
, .lower          = (double []){ 1e2, 1e3, 1e3, 10,  10,  10,  10,  10 }
, .upper          = (double []){ 1e4, 1e4, 1e4, 1e3, 1e3, 1e3, 1e3, 1e3 }
, .enforce_bounds = 1
, .generations    = 5000
, .popsize        = 5
, .f              = { &f, &g1, &g2, &g3, &g4, &g5, &g6 }
};
