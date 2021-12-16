#include "constraint.h"
#include <stdio.h>

static double f1 (double *x)
{
    return x [0];
}

static double g (double *x)
{
    int i;
    double s = 0;
    for (i=1; i<zdt4.dimension; i++) {
        s += pow (x [i], 2) - 10 * cos (4 * M_PI * x [i]);
    }
    return 1. + 10. * (zdt4.dimension - 1.) + s;
}

static double f2 (double *x)
{
    double rg = g (x);
    return 1 - sqrt (x [0] / rg);
}

struct multi_problem zdt4 =
{ .dimension      = 10
, .nfunc          = 2
, .nconstraint    = 0
, .lower          = (double []){ 0, -5, -5, -5, -5, -5, -5, -5, -5, -5
  }
, .upper          = (double []){ 1,  5,  5,  5,  5,  5,  5,  5,  5,  5
  }
, .enforce_bounds = 1
, .popsize        = 100
, .generations    = 150
, .f              = { &f1, &f2 }
, .name           = "Zitzler et. al. (ZDT4)"
};
