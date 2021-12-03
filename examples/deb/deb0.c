#include "constraint.h"
#include <stdio.h>

static double f (double *x)
{
    return pow (x [0] - 0.8, 2) + pow (x [1] - 0.3, 2);
}

static double g1 (double *x)
{
    return (pow (x [0] - 0.2, 2) + pow (x [1] - 0.5, 2)) / 0.16 - 1;
}

static double g2 (double *x)
{
    return 1 - (pow (x [0] + 0.5, 2) + pow (x [1] - 0.5, 2)) / 0.81;
}

struct constrained_problem deb_0 =
{ .dimension      = 2
, .nfunc          = 3
, .lower          = (double []){ -5, -5 }
, .upper          = (double []){  5,  5 }
, .enforce_bounds = 0
, .popsize        = 10
, .f              = { &f, &g1, &g2 }
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][2] =
        { {1.43975,             0.2857709}
        , {1.5519118010997772,  2.5159401595592499} /* p 54 */
        , {3.6450952291488647,  4.5980691909790039} /* p  0 */
        , {0.83757007122039795, 8.1951549649238586} /* p 0 old */
        , {-2.102253e+16,       1.643393e+17}
        , {0.9321836,           0.5648792}
        };
    size_t sz = sizeof (xx) / (deb_0.dimension * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double g;
        double s = 0;
        for (j=0; j<deb_0.nfunc - 1; j++) {
            g = deb_0.f [j+1] (x);
            printf ("c: %e\n", g);
            if (g > 0) {
                s += g;
            }
        }
        printf ("constraints: %e\n", s);
        printf ("evaluation:  %e\n", deb_0.f [0] (x));
    }
}
#endif
