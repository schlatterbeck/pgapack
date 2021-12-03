#include "constraint.h"
#include <stdio.h>

static double f (double *x)
{
    return pow (x [0], 2) + pow (x [1], 2) + x [0] * x [1] - 14 * x [0]
         - 16 * x [1] + pow (x [2] - 10, 2) + 4 * pow (x [3] - 5, 2)
         + pow (x [4] - 3, + 2) + 2 * pow (x [5] - 1, 2) + 5 * pow (x [6], 2)
         + 7 * pow (x [7] - 11, 2) + 2 * pow (x [8] - 10, 2)
         + pow (x [9] - 7, 2) + 45;
}

static double g1 (double *x)
{
    return 4 * x [0] + 5 * x [1] - 3 * x [6] + 9 * x [7] - 105;
}

static double g2 (double *x)
{
    return 10 * x [0] - 8 * x [1] - 17 * x [6] + 2 * x [7];
}

static double g3 (double *x)
{
    return -8 * x [0] + 2 * x [1] + 5 * x [8] - 2 * x [9] - 12;
}

static double g4 (double *x)
{
    return 3 * pow (x [0] - 2, 2) + 4 * pow (x [1] - 3, 2)
         + 2 * pow (x [2], 2) - 7 * x [3] - 120;
}

static double g5 (double *x)
{
    return 5 * pow (x [0], 2) + 8 * x [1] + pow (x [2] - 6, 2)
         - 2 * x [3] - 40;
}

static double g6 (double *x)
{
    return pow (x [0], 2) + 2 * pow (x [1] - 2, 2) - 2 * x [0] * x [1]
         + 14 * x [4] - 6 * x [5];
}

static double g7 (double *x)
{
    return 0.5 * pow (x [0] - 8, 2) + 2 * pow (x [1] - 4, 2)
         + 3 * pow (x [4], 2) - x [5] - 30;
}

static double g8 (double *x)
{
    return -3 * x [0] + 6 * x [1] + 12 * pow (x [8] - 8, 2) - 7 * x [9];
}

struct constrained_problem deb_8 =
{ .dimension      = 10
, .nfunc          = 9
, .lower          = (double []){-10, -10, -10, -10, -10, -10, -10, -10, -10, -10}
, .upper          = (double []){ 10,  10,  10,  10,  10,  10,  10,  10,  10,  10}
, .enforce_bounds = 1
, .iterations     = 5000
, .popsize        = 30
, .f              = { &f, &g1, &g2, &g3, &g4, &g5, &g6, &g7, &g8 }
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][10] =
        { { 2.171996, 2.363683, 8.773926, 5.095984, 0.9906548
          , 1.430574, 1.321644, 9.828726, 8.280092, 8.375927
          }
        };
    size_t sz = sizeof (xx) / (deb_8.dimension * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double g;
        double s = 0;
        for (j=0; j<deb_8.nfunc - 1; j++) {
            g = deb_8.f [j+1] (x);
            printf ("c: %e\n", g);
            if (g > 0) {
                s += g;
            }
        }
        printf ("constraints: %e\n", s);
        printf ("evaluation:  %e\n", deb_8.f [0] (x));
    }
}
#endif
