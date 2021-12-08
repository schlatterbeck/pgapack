#include "constraint.h"
#include <stdio.h>

static double f (double *x)
{
    return 5.3578547 * pow (x [2], 2) + 0.8356891 * x [0] * x [4]
         + 37.293239 * x [0] - 40792.141;
}

static double g1 (double *x)
{
    return 2.2053e-3 * x [2] * x [4] - 85.334407 - 5.6858e-3 * x [1] * x [4]
         - 6.262e-4 * x [0] * x [3];
}

static double g2 (double *x)
{
    return 85.334407 - 92 + 5.6858e-3 * x [1] * x [4] + 6.262e-4 * x [0] * x [3]
         - 2.2053e-3 * x [2] * x [4];
}

static double g3 (double *x)
{
    return 90 - ( 80.51249 + 7.1317e-3 * x [1] * x [4]
                + 2.9955e-3 * x [0] * x [1] + 2.1813e-3 * pow (x [2], 2)
                );
}

static double g4 (double *x)
{
    return 80.51249 + 7.1317e-3 * x [1] * x [4] + 2.9955e-3 * x [0] * x [1]
         + 2.1813e-3 * pow (x [2], 2) - 110;
}

static double g5 (double *x)
{
    return 20 - ( 9.300961 + 4.7026e-3 * x [2] * x [4]
                + 1.2547e-3 * x [0] * x [2] + 1.9085e-3 * x [2] * x [3]
                );
}

static double g6 (double *x)
{
    return 9.300961 + 4.7026e-3 * x [2] * x [4] + 1.2547e-3 * x [0] * x [2]
         + 1.9085e-3 * x [2] * x [3] - 25;
}

struct constrained_problem deb_6 =
{ .dimension      = 5
, .nfunc          = 7
, .lower          = (double []){  78, 33, 27, 27, 27 }
, .upper          = (double []){ 102, 45, 45, 45, 45 }
, .enforce_bounds = 1
, .generations    = 500
, .popsize        = 5
, .f              = { &f, &g1, &g2, &g3, &g4, &g5, &g6 }
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][7] =
        { { 78,    33,    29.99526, 45,    36.77581 }
        , { 80.49, 35.07, 32.05,    40.33, 33.34    }
        };
    size_t sz = sizeof (xx) / (deb_6.dimension * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double g;
        double s = 0;
        for (j=0; j<deb_6.nfunc - 1; j++) {
            g = deb_6.f [j+1] (x);
            printf ("c: %e\n", g);
            if (g > 0) {
                s += g;
            }
        }
        printf ("constraints: %e\n", s);
        printf ("evaluation:  %11.5f\n", deb_6.f [0] (x));
    }
}
#endif
