#include "constraint.h"
#include <stdio.h>

static double f1 (double *x)
{
    return pow (x [0] - 2, 2) + pow (x [1] - 1, 2) + 2;
}

static double f2 (double *x)
{
    return 9 * x [0] - pow (x [1] - 1, 2);
}

static double g1 (double *x)
{
    return pow (x [0], 2) + pow (x [1], 2) - 225;
}

static double g2 (double *x)
{
    return 10 + x [0] - 3 * x [1];
}

struct multi_problem srn =
{ .dimension      = 2
, .nfunc          = 4
, .nconstraint    = 2
, .lower          = (double []){ -20, -20 }
, .upper          = (double []){  20,  20 }
, .enforce_bounds = 1
, .popsize        = 10
, .generations    = 10
, .f              = { &f1, &f2, &g1, &g2 }
, .name           = "SRN"
};
