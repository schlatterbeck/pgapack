#include "optimize.h"
#include <stdio.h>

static double f1 (double *x)
{
    return x [0];
}

static double g (double *x)
{
    int i;
    double s = 0;
    for (i=1; i<zdt3.dimension; i++) {
        s += x [i];
    }
    s *= 9. / (zdt3.dimension - 1);
    return 1. + s;
}

static double f2 (double *x)
{
    double rg = g (x);
    return rg * (1 - sqrt (x [0] / rg) - x [0] / rg * sin (10 * M_PI * x [0]));
}

struct multi_problem zdt3 =
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
, .generations    = 500
, .f              = { &f1, &f2 }
, .name           = "Zitzler et. al. (ZDT3)"
};
