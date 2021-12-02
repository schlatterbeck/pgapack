#include "constraint.h"
#include <stdio.h>

static double f (double *x)
{
    double r = 0;
    int i;
    for (i=0; i<4; i++) {
        r += (1 - x [i]) * x [i];
    }
    r *= 5;
    for (i=4; i<13; i++) {
        r -= x [i];
    }
    return r;
}

static double g1 (double *x)
{
    return 2 * (x [0] + x [1]) + x [9] + x [10] - 10;
}

static double g2 (double *x)
{
    return 2 * (x [0] + x [2]) + x [9] + x [11] - 10;
}

static double g3 (double *x)
{
    return 2 * (x [1] + x [2]) + x [10] + x [11] - 10;
}

static double g4 (double *x)
{
    return x [9] - 8 * x [0];
}

static double g5 (double *x)
{
    return x [10] - 8 * x [1];
}

static double g6 (double *x)
{
    return x [11] - 8 * x [2];
}

static double g7 (double *x)
{
    return x [9] - 2 * x [3] - x [4];
}

static double g8 (double *x)
{
    return x [10] - 2 * x [5] - x [6];
}

static double g9 (double *x)
{
    return x [11] - 2 * x [7] - x [8];
}

struct constrained_problem deb_3 =
{ .dimension      = 13
, .nfunc          = 10
, .lower          = (double []){ 0, 0, 0, 0, 0, 0, 0, 0, 0,   0,   0,   0, 0 }
, .upper          = (double []){ 1, 1, 1, 1, 1, 1, 1, 1, 1, 100, 100, 100, 1 }
, .enforce_bounds = 1
, .iterations     = 1000
, .f              = { &f, &g1, &g2, &g3, &g4, &g5, &g6, &g7, &g8, &g9 }
};
