#include "optimize.h"
#include <stdio.h>

static double f1 (double *x)
{
    return -(106780.37 * (x [1] + x [2]) + 61704.67);
}

static double f2 (double *x)
{
    return -(3000. * x [0]);
}

static double f3 (double *x)
{
    return -(305700. * 2289. * x [1] / pow (0.06 * 2289., 0.65));
}

static double f4 (double *x)
{
    return -(250. * 2289. * exp (-39.75 * x [1] + 9.9 * x [2] + 2.74));
}

static double f5 (double *x)
{
    return -(25. * (1.39 / (x [0] * x [1]) + 4940. * x [2] - 80.));
}

static double g1 (double *x)
{
    return 0.00139 / (x [0] * x [1]) + 4.94 * x [2] - 0.08 - 1;
}

static double g2 (double *x)
{
    return 0.000306 / (x [0] * x [1]) + 1.082 * x [2] - 0.0986 - 1;
}

static double g3 (double *x)
{
    return 12.307 / (x [0] * x [1]) + 49408.24 * x [2] + 4051.02 - 5e4;
}

static double g4 (double *x)
{
    return 2.098 / (x [0] * x [1]) + 8046.33 * x [2] - 696.71 - 16e3;
}

static double g5 (double *x)
{
    return 2.138 / (x [0] * x [1]) + 7883.39 * x [2] - 705.04 - 1e4;
}

static double g6 (double *x)
{
    return 0.417 / (x [0] * x [1]) + 1721.26 * x [2] - 136.54 - 2e3;
}

static double g7 (double *x)
{
    return 0.164 / (x [0] * x [1]) + 631.13 * x [2] - 54.48 - 550;
}

struct multi_problem water_m =
{ .dimension      = 3
, .nfunc          = 12
, .nconstraint    = 7
, .maximize       = 1
, .lower          = (double []){ 0.01, 0.01, 0.01 }
, .upper          = (double []){ 0.45, 0.10, 0.10 }
, .enforce_bounds = 1
, .crossover_prob = 0.0000000001
, .f              = { &f1, &f2, &f3, &f4, &f5
                    , &g1, &g2, &g3, &g4, &g5, &g6, &g7
                    }
, .name           = "Ray, Tai, Seow (WATER)"
};
