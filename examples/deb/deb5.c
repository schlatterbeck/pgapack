#include "constraint.h"
#include <stdio.h>

static double f (double *x)
{
    return pow (x [0] - 10, 2) + 5 * pow (x [1] - 12, 2) + pow (x [2], 4)
         + 3 * pow (x [3] - 11, 2) + 10 * pow (x [4], 6) + 7 * pow (x [5], 2)
         + pow (x [6], 4) - 4 * x [5] * x [6] - 10 * x [5] - 8 * x [6];
}

static double g1 (double *x)
{
    return 2 * pow (x [0], 2) + 3 * pow (x [1], 4) + x [2]
         + 4 * pow (x [3], 2) + 5 * x [4] - 127;
}

static double g2 (double *x)
{
    return 7 * x [0] + 3 * x [1] + 10 * pow (x [2], 2) + x [3] - x [4] - 282;
}

static double g3 (double *x)
{
                        
    return 23 * x [0] + pow (x [1], 2) + 6 * pow (x [5], 2) - 8 * x [6] - 196;
}

static double g4 (double *x)
{
    return 4 * pow (x [0], 2) + pow (x [1], 2) - 3 * x [0] * x [1]
         + 2 * pow (x [2], 2) + 5 * x [5] - 11 * x [6];
}

static double g5 (double *x)
{
    return 1250 * x [4] + x [1] * x [3] - x [1] * x [6] - 1250 * x [3];
}

struct constrained_problem deb_5 =
{ .dimension      = 7
, .nfunc          = 5
, .lower          = (double []){ -10, -10, -10, -10, -10, -10, -10 }
, .upper          = (double []){  10,  10,  10,  10,  10,  10,  10 }
, .enforce_bounds = 1
, .iterations     = 1000
, .f              = { &f, &g1, &g2, &g3, &g4, &g5 }
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][7] =
        { { 2.330499, 1.951372, -0.4775414, 4.365726, -0.624487
          , 1.038131, 1.594227
          }
        };
    size_t sz = sizeof (xx) / (deb_5.dimension * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double g;
        double s = 0;
        for (i=0; i< deb_5.nfunc - 1; i++) {
            g = deb_5.f [i+1] (x);
            printf ("c: %e\n", g);
            if (g > 0) {
                s += g;
            }
        }
        printf ("constraints: %11.7f\n", s);
        printf ("evaluation:  %11.7f\n", deb_5.f [0] (x));
    }
}
#endif
