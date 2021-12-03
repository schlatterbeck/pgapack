#include "constraint.h"
#include <stdio.h>

#define eps 1e-2

static double f (double *x)
{
    return exp (x [0] * x [1] * x [2] * x [3] * x [4]);
}

static inline double eq1 (double *x)
{
    return x [0] * x [0] + x [1] * x [1] + x [2] * x [2] + x [3] * x [3]
         + x [4] * x [4] - 10;
}

static double g1 (double *x)
{
    return -(eq1 (x) + eps);
}

static double g2 (double *x)
{
    return  eq1 (x) - eps;
}

static inline double eq2 (double *x)
{
    return x [1] * x [2] - 5 * x [3] * x [4];
}

static double g3 (double *x)
{
    return -(eq2 (x) + eps);
}

static double g4 (double *x)
{
    return eq2 (x) - eps;
}

static inline double eq3 (double *x)
{
    return pow (x [0], 3) + pow (x [1], 3) + 1;
}

static double g5 (double *x)
{
    return -(eq3 (x) + eps);
}

static double g6 (double *x)
{
    return eq3 (x) - eps;
}

struct constrained_problem deb_7 =
{ .dimension      = 5
, .nfunc          = 7
, .lower          = (double []){ -2.3, -2.3, -3.2, -3.2, -3.2 }
, .upper          = (double []){  2.3,  2.3,  3.2,  3.2,  3.2 }
, .enforce_bounds = 1
, .generations    = 2000
, .popsize        = 4
, .f              = { &f, &g1, &g2, &g3, &g4, &g5, &g6 }
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][7] =
        { { -1.717143,   1.595709,  1.827247, -0.7636413, -0.7636450 }
        , { -0.8560713, -0.7196526, 2.817094,  0.6613182, -0.6131167 }
        };
    size_t sz = sizeof (xx) / (deb_7.dimension * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double g;
        double s = 0;
        for (j=0; j<deb_7.nfunc - 1; j++) {
            g = deb_7.f [j+1] (x);
            printf ("c: %e\n", g);
            if (g > 0) {
                s += g;
            }
        }
        printf ("constraints: %e\n", s);
        printf ("evaluation:  %e\n", deb_7.f [0] (x));
        printf ("e1: %e\n", eq1 (x));
        printf ("e2: %e\n", eq2 (x));
        printf ("e3: %e\n", eq3 (x));
    }
}
#endif
