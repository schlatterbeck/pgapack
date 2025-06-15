#include "optimize.h"
#include <stdio.h>

/* Defaults */
#define NOBJ 3
#define DIM  7

static double g (double *x, int k)
{
    int i;
    double s = 0;
    for (i=0; i<k; i++) {
        s += pow (x [i] - 0.5, 2) - cos (20 * M_PI * (x [i] - 0.5));
    }
    return 100 * (k + s);
}

static void f (double *x, int nx, double *y, int ny)
{
    int i, j;
    int k = nx - ny + 1;
    double gv;
    assert (k > 0);
    gv = g (x + nx - k, k);
    for (i=0; i<ny; i++) {
        double p = 0.5;
        for (j=0; j<ny-i-1; j++) {
            p *= x [j];
        }
        if (j < ny-1) {
            p *= (1 - x [j]);
        }
        y [i] = p * (1 + gv);
    }
}

/* DTLZ suggest k = 5, i.e. for 3 objectives n = 3 + 5 - 1 = 7 */
struct multi_problem dtlz1 =
{ .dimension      = DIM
, .n_obj          = NOBJ
, .nconstraint    = 0
, .lower          = 0
, .upper          = 1
, .popsize        = 0
, .generations    = 0
, .f              = f
, .name           = "DTLZ1"
};

#ifdef DEBUG_EVAL
#include <stdio.h>
#define NOBJ 3
#define DIM  7
#define K (DIM - NOBJ + 1)
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
        f (xx [i], DIM, yy, NOBJ);
        printf ("G: %e\n", g (xx [i] + DIM - K, K));
        for (j=0; j<NOBJ; j++) {
            printf ("F %d %e\n", j, yy [j]);
            s += yy [j];
        }
        printf ("\ns: %e\n", s);
    }
}
#endif /* DEBUG_EVAL */
