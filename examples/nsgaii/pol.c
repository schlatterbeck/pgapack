#include "constraint.h"
#include <stdio.h>

static double f1 (double *x)
{
    double a1 = 0.5 * sin (1) - 2 * cos (1) +     sin (2) - 1.5 * cos (2);
    double a2 = 1.5 * sin (1) -     cos (1) + 2 * sin (2) - 0.5 * cos (2);
    double b1 = 0.5 * sin (x [0]) - 2   * cos (x [0])
              +       sin (x [1]) - 1.5 * cos (x [1]);
    double b2 = 1.5 * sin (x [0]) -       cos (x [0])
              + 2   * sin (x [1]) - 0.5 * cos (x [1]);
    return (1 + pow (a1 - b1, 2) + pow (a2 - b2, 2));
}

static double f2 (double *x)
{
    return (pow (x [0] + 3, 2) + pow (x [1] + 1, 2));
}

struct multi_problem pol =
{ .dimension      = 2
, .nfunc          = 2
, .nconstraint    = 0
, .lower          = (double []){ -M_PI, -M_PI }
, .upper          = (double []){  M_PI,  M_PI }
, .enforce_bounds = 0
, .popsize        = 10
, .generations    = 150
, .f              = { &f1, &f2 }
, .name           = "Poloni's study (POL)"
};
