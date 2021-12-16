#include "constraint.h"
#include <stdio.h>

static double f1 (double *x)
{
    int i;
    double s = 0;
    for (i=0; i<kur.dimension - 1; i++) {
        s -= 10 * exp (-0.2 * sqrt (pow (x [i], 2) + pow (x [i+1], 2)));
    }
    return s;
}

static double f2 (double *x)
{
    int i;
    double s = 0;
    for (i=0; i<kur.dimension; i++) {
        s += pow (fabs (x [i]), 0.8) + 5 * sin (pow (x [i], 3));
    }
    return s;
}

struct multi_problem kur =
{ .dimension      = 3
, .nfunc          = 2
, .nconstraint    = 0
, .lower          = (double []){ -5, -5, -5 }
, .upper          = (double []){  5,  5,  5 }
, .enforce_bounds = 0
, .popsize        = 10
, .generations    = 150
, .f              = { &f1, &f2 }
, .name           = "Kursawe's study (KUR)"
};
