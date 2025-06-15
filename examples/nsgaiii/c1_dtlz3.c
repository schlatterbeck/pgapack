#include "optimize.h"
#include <stdio.h>

#define NOBJ  3
#define DIM  12

/* Multipliers for constraint sum2, starting with NOBJ = 3
 * Higher indeces use last value
 */
static double c1_mul [] = {9, 12.5, 12.5, 15, 15};
#define C1_SIZE (sizeof (c1_mul) / sizeof (*c1_mul))

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
    int nobj = ny - 1; /* 1 constraint */
    int k = nx - nobj + 1;
    size_t c1_idx = nobj - 3;
    double gv;
    double s = 0;
    assert (k > 0);
    gv = g (x + nx - k, k);
    if (c1_idx >= C1_SIZE) {
        c1_idx = C1_SIZE - 1;
    }
    for (i=0; i<nobj; i++) {
        double p = 1;
        for (j=0; j<nobj-i-1; j++) {
            p *= cos (x [j] * M_PI / 2.0);
        }
        if (j < nobj-1) {
            p *= sin (x [j] * M_PI / 2.0);
        }
        y [i] = p * (1 + gv);
        s  += pow (y [i], 2);
    }
    y [nobj] = -(s - 16) * (s - c1_mul [c1_idx]);
}

struct multi_problem c1_dtlz3 =
{ .dimension      = DIM
, .n_obj          = NOBJ
, .nconstraint    = 1
, .lower          = 0
, .upper          = 1
, .popsize        = 0
, .generations    = 1000
, .f              = f
, .name           = "C1-DTLZ3"
};

#ifdef DEBUG_EVAL
#include <stdio.h>
#define K (DIM - NOBJ + 1)
int main ()
{
    int i, j;
    static double xx [][DIM] =
        { { 0.4757814, 0.2537747, 0.5, 0.5, 0.5, 0.5
          , 0.5,       0.5000003, 0.5, 0.5, 0.5, 0.5
          }
        , { 0.5911212, 0.703804,  0.5, 0.5, 0.5, 0.5
          , 0.5,       0.5000001, 0.5, 0.5, 0.5, 0.5
          }
        , { 0.8452794, 0.2859174, 0.5, 0.5, 0.5, 0.5
          , 0.5,       0.5,       0.5, 0.5, 0.5, 0.5
          }
        , { 0.5911212, 0.2872442, 0.5, 0.5, 0.5, 0.5
          , 0.5,       0.5,       0.5, 0.5, 0.5, 0.5
          }
        };
    double yy [NOBJ + 1];
    size_t sz = sizeof (xx) / (sizeof (double) * DIM);
    for (i=0; i<sz; i++) {
        double s = 0;
        f (xx [i], yy);
        printf ("G: %e\n", g (xx [i] + DIM - K, K));
        if (yy [NOBJ] <= 0) {
            for (j=0; j<NOBJ+1; j++) {
                printf ("F %d %e\n", j, yy [j]);
                s += yy [j];
            }
        } else {
            printf ("F %d %e\n", NOBJ, yy [NOBJ]);
        }
        //printf ("\ns: %e\n", s);
    }
}
#endif /* DEBUG_EVAL */
