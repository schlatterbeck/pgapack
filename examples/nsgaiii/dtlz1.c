#include "optimize.h"
#include <stdio.h>

#define NOBJ 3
#define DIM  7
#define K (DIM - NOBJ + 1)

static double g (double *x)
{
    int i;
    double s = 0;
    for (i=0; i<K; i++) {
        s += pow (x [i] - 0.5, 2) - cos (20 * M_PI * (x [i] - 0.5));
    }
    return 100 * (K + s);
}

static void f (double *x, double *y)
{
    int i, j;
    double gv = g (x + DIM - K);
    for (i=0; i<NOBJ; i++) {
        double p = 0.5;
        for (j=0; j<NOBJ-i-1; j++) {
            p *= x [j];
        }
        if (j < NOBJ-1) {
            p *= (1 - x [j]);
        }
        y [i] = p * (1 + gv);
    }
}

struct multi_problem dtlz1 =
{ .dimension      = DIM
, .nfunc          = NOBJ
, .nconstraint    = 0
, .lower          = (double []){ 0, 0, 0, 0, 0, 0, 0 }
, .upper          = (double []){ 1, 1, 1, 1, 1, 1, 1 }
, .popsize        = 0
, .generations    = 0
, .f              = f
, .name           = "DTLZ1"
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][DIM] =
        { {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}
        , {.75, .85, 0.5, 0.5, 0.5, 0.5, 0.5}
        , {  1,   1, 0.5, 0.5, 0.5, 0.5, 0.5}
        , {  0,   0, 0.5, 0.5, 0.5, 0.5, 0.5}
        };
    double yy [NOBJ];
    size_t sz = sizeof (xx) / (sizeof (double) * DIM);
    for (i=0; i<sz; i++) {
        double s = 0;
        f (xx [i], yy);
        printf ("G: %e\n", g (xx [i] + DIM - K));
        for (j=0; j<NOBJ; j++) {
            printf ("F %d %e\n", j, yy [j]);
            s += yy [j];
        }
        printf ("\ns: %e\n", s);
    }
}
#endif /* DEBUG_EVAL */
