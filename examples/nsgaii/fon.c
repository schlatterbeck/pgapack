#include "optimize.h"
#include <stdio.h>

static double f1 (double *x)
{
    int i;
    double s = 0;
    for (i=0; i<fon.dimension; i++) {
        s -= pow (x [i] - 1 / sqrt (3), 2);
    }
    return 1 - exp (s);
}

static double f2 (double *x)
{
    int i;
    double s = 0;
    for (i=0; i<fon.dimension; i++) {
        s -= pow (x [i] + 1 / sqrt (3), 2);
    }
    return 1 - exp (s);
}

struct multi_problem fon =
{ .dimension      = 3
, .nfunc          = 2
, .nconstraint    = 0
, .lower          = (double []){ -4, -4, -4 }
, .upper          = (double []){  4,  4,  4 }
, .enforce_bounds = 0
, .popsize        = 10
, .generations    = 150
, .f              = { &f1, &f2 }
, .name           = "Fonseca and Fleming's study (FON)"
};
