#include "constraint.h"
#include <stdio.h>
#include <string.h>

#define DIMENSION 10

#define PHI (45. / 180. * M_PI)
static double y [DIMENSION];
static double *last_x = NULL; /* the last x seen for which y is valid */

static double rotmatrix [DIMENSION][DIMENSION];

static void compute_rotation (void)
{
    int done = 0;
    int i, j;
    double g1 [DIMENSION];
    double g2 [DIMENSION];
    double v  [DIMENSION][DIMENSION];
    double w  [DIMENSION][DIMENSION];
    double I  [DIMENSION][DIMENSION];

    if (done) {
        return;
    }

    for (i=0; i<DIMENSION; i++) {
        g1 [i] = g2 [i] = 0;
        for (j=0; j<DIMENSION; j++) {
            v [i][j] = w [i][j] = I [i][j] = 0;
            if (i==j) {
                I [i][j] = 1;
            }
        }
    }
    g1 [0] = g1 [1] = g2 [2] = g2 [3] = pow (0.5, 0.5);
    for (i=0; i<DIMENSION; i++) {
        for (j=0; j<DIMENSION; j++) {
            v [i][j] = (g1 [i] * g1 [j] + g2 [i] * g2 [j]) * (cos (PHI) - 1);
            w [i][j] = (g1 [i] * g2 [j] - g2 [i] * g1 [j]) * sin (PHI);
        }
    }
    for (i=0; i<DIMENSION; i++) {
        for (j=0; j<DIMENSION; j++) {
            rotmatrix [i][j] = I [i][j] + v [i][j] + w [i][j];
        }
    }
    done = 1;
}

static void rotate (double *x)
{
    int i, j;
    if (x == last_x) {
        return;
    }
    compute_rotation ();
    for (i=0; i<DIMENSION; i++) {
        y [i] = 0;
        for (j=0; j<DIMENSION; j++) {
            y [i] += rotmatrix [i][j] * x [j];
        }
    }
}

static double f1 (double *x)
{
    rotate (x);
    return y [0];
}

static double g (double *x)
{
    int i;
    double s = 0;
    rotate (x);
    for (i=1; i<rotated.dimension; i++) {
        s += pow (y [i], 2) - 10 * cos (4 * M_PI * y [i]);
    }
    return 1. + 10 * (rotated.dimension - 1) + s;
}

static double f2 (double *x)
{
    double rg = g (x);
    return rg * exp (-y [0] / rg);
}

/* Limit the y [0] range to [-0.3, 0.3] */

static double g1 (double *x)
{
    rotate (x);
    return y [0] - 0.3;
}

static double g2 (double *x)
{
    rotate (x);
    return -0.3 - y [0];
}

struct multi_problem rotated =
{ .dimension      = DIMENSION
, .nfunc          = 4
, .nconstraint    = 2
, .lower          = (double [])
                    { -.3, -.3, -.3, -.3, -.3, -.3, -.3, -.3, -.3, -.3 }
, .upper          = (double [])
                    {  .3,  .3,  .3,  .3,  .3,  .3,  .3,  .3,  .3,  .3 }
, .enforce_bounds = 1
, .popsize        = 50
, .generations    = 2500
, .f              = { &f1, &f2, &g1, &g2 }
, .name           = "Rotated Problem"
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][DIMENSION] =
        { { -0.25606602, 0.04393398, -0.10606602, -0.10606602
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.23899495, 0.04100505, -0.09899495, -0.09899495
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.22192388,  0.03807612, -0.09192388, -0.09192388
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.20485281,  0.03514719, -0.08485281, -0.08485281
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.18778175,  0.03221825, -0.07778175, -0.07778175
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.17071068,  0.02928932, -0.07071068, -0.07071068
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.15363961,  0.02636039, -0.06363961, -0.06363961
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.13656854,  0.02343146, -0.05656854, -0.05656854
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.11949747,  0.02050253, -0.04949747, -0.04949747
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.10242641,  0.01757359, -0.04242641, -0.04242641
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.08535534,  0.01464466, -0.03535534, -0.03535534
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.06828427,  0.01171573, -0.02828427, -0.02828427
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.0512132,  0.0087868, -0.0212132, -0.0212132
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.03414214,  0.00585786, -0.01414214, -0.01414214
          , 0, 0, 0, 0, 0, 0
          }
        , { -0.01707107,  0.00292893, -0.00707107, -0.00707107
          , 0, 0, 0, 0, 0, 0
          }
        , { 0, 0, 0, 0
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.01707107, -0.00292893,  0.00707107,  0.00707107
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.03414214, -0.00585786,  0.01414214,  0.01414214
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.0512132, -0.0087868,  0.0212132,  0.0212132
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.06828427, -0.01171573,  0.02828427,  0.02828427
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.08535534, -0.01464466,  0.03535534,  0.03535534
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.10242641, -0.01757359,  0.04242641,  0.04242641
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.11949747, -0.02050253,  0.04949747,  0.04949747
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.13656854, -0.02343146,  0.05656854,  0.05656854
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.15363961, -0.02636039,  0.06363961,  0.06363961
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.17071068, -0.02928932,  0.07071068,  0.07071068
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.18778175, -0.03221825,  0.07778175,  0.07778175
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.20485281, -0.03514719,  0.08485281,  0.08485281
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.22192388, -0.03807612,  0.09192388,  0.09192388
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.23899495, -0.04100505,  0.09899495,  0.09899495
          , 0, 0, 0, 0, 0, 0
          }
        , { 0.25606602, -0.04393398,  0.10606602,  0.10606602
          , 0, 0, 0, 0, 0, 0
          }
        };
    size_t sz = sizeof (xx) / (rotated.dimension * sizeof (double));
    printf ("Example: rotated optimum\n");
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double f;
        double s = 0;
        printf ("%d\n", i);
        for (j=0; j<rotated.nfunc; j++) {
            f = rotated.f [j] (x);
            printf ("F %d  %e\n", j, f);
        }
    }
}
#endif
