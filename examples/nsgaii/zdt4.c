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
, .lower          = (double []){ 0, -5, -5, -5, -5, -5, -5, -5, -5, -5 }
, .upper          = (double []){ 1,  5,  5,  5,  5,  5,  5,  5,  5,  5 }
, .enforce_bounds = 1
, .dither         = 0.5
, .jitter         = 0.005
, .crossover_prob = 0.0000001
, .f              = { &f1, &f2 }
, .name           = "Zitzler et. al. (ZDT4)"
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][10] =
        { { 1,    0, 0, 0, 0, 0, 0, 0, 0, 0 }
        , { 0.75, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
        , { 0.5,  0, 0, 0, 0, 0, 0, 0, 0, 0 }
        , { 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
        };
    size_t sz = sizeof (xx) / (zdt4.dimension * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double f;
        double s = 0;
        printf ("%d\n", i);
        for (j=0; j<zdt4.nfunc; j++) {
            f = zdt4.f [j] (x);
            printf ("f%d: %e\n", j, f);
        }
    }
}
#endif
