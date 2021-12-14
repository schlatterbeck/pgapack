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
    for (i=1; i<zdt2.dimension; i++) {
        s += x [i];
    }
    s *= 9. / (zdt2.dimension - 1.0);
    return 1. + s;
}

static double f2 (double *x)
{
    double rg = g (x);
    return rg * (1 - pow (x [0] / rg, 2));
}

struct multi_problem zdt2 =
{ .dimension      = 30
, .nfunc          = 2
, .nconstraint    = 0
, .lower          = (double [])
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  }
, .upper          = (double [])
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  }
, .enforce_bounds = 1
, .popsize        = 10
, .generations    = 1000
, .f              = { &f1, &f2 }
, .name           = "Zitzler et. al. (ZDT2)"
};
