#include "constraint.h"
#include <stdio.h>

static double f1 (double *x)
{
    return pow (x [0], 2);
}

static double f2 (double *x)
{
    return pow (x [0] - 2, 2);
}

struct multi_problem sch =
{ .dimension      = 1
, .nfunc          = 2
, .nconstraint    = 0
, .lower          = (double []){ -1e3,  -1e3 }
, .upper          = (double []){  1e3,   1e3 }
, .enforce_bounds = 0
, .popsize        = 10
, .generations    = 150
, .f              = { &f1, &f2 }
, .name           = "Schaffer's study (SCH)"
};
