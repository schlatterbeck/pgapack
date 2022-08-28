#include "constraint.h"
#include <stdio.h>

/*
 * This is the welded beam design from the beginning of the paper.
 * We map the variables h, l, t, b to x [0-3].
 * The cache holds tau_s, tau_ss, tau, p_c, delta, sigma in that order.
 */

static double cache [6];
static double *last_x; /* the last x seen for which the cache is valid */
#define TAU_S  (cache [0])
#define TAU_SS (cache [1])
#define TAU    (cache [2])
#define P_C    (cache [3])
#define DELTA  (cache [4])
#define SIGMA  (cache [5])

#define H (x [0])
#define L (x [1])
#define T (x [2])
#define B (x [3])

static inline double tau_s (double *x)
{
    return 6000. / (sqrt (2) * H * L);
}

static inline double tau_ss (double *x)
{
    return 6000. * (14 + 0.5 * L)
         * sqrt (0.25 * (pow (L, 2) + pow (H + T, 2)))
         / ( 2 * 0.707 * H * L
           * (pow (L, 2) / 12 + 0.25 * pow (H + T, 2))
           );
}

/* This already uses tau_s and tau_ss from the cache */
static inline double tau (double *x)
{
    return sqrt
        ( pow (TAU_S, 2) + pow (TAU_SS, 2)
        + ( L * TAU_S * TAU_SS
          / sqrt (0.25 * (pow (L, 2) + pow (H + T, 2)))
          )
        );
}

static inline double sigma (double *x)
{
    return 504000.0 / (pow (T, 2) * B);
}

static inline double p_c (double *x)
{
    return 64746.022 * (1 - 0.0282346 * T) * T * pow (B, 3);
}

static inline double delta (double *x)
{
    return 2.1952 / (pow (T, 3) * B);
}

static void update_cache (double *x)
{
    if (x == last_x) {
        return;
    }
    cache [0] = tau_s (x);
    cache [1] = tau_ss (x);
    cache [2] = tau (x);
    cache [3] = p_c (x);
    cache [4] = delta (x);
    cache [5] = sigma (x);
}

static double f (double *x)
{
    update_cache (x);
    return 1.10471 * pow (x [0], 2) * x [1]
         + 0.04811 * x [2] * x [3] * (14.0 + x [1]);
}

static double g1 (double *x)
{
    update_cache (x);
    return TAU - 13600;
}

static double g2 (double *x)
{
    update_cache (x);
    return SIGMA - 30000.0;
}

static double g3 (double *x)
{
    return x [0] - x [3];
}

static double g4 (double *x)
{
    update_cache (x);
    return 6000.0 - P_C;
}

static double g5 (double *x)
{
    update_cache (x);
    return DELTA - 0.25;
}

struct constrained_problem deb_9 =
{ .dimension      = 4
, .nfunc          = 6
, .lower          = (double []){0.125, 0.1, 0.1, 0.1}
, .upper          = (double []){   10,  10,  10,  10}
, .enforce_bounds = 1
, .popsize        = 10
, .generations    = 640
, .f              = { &f, &g1, &g2, &g3, &g4, &g5 }
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][4] =
        { { 0.2444,    6.2187,   8.2915,   0.2444 }
        , { 0.2489,    6.1730,   8.1789,   0.2533 }
        , { 0.2443689, 6.218606, 8.291473, 0.2443689 }
        };
    size_t sz = sizeof (xx) / (deb_9.dimension * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double g;
        double s = 0;
        for (j=0; j<deb_9.nfunc - 1; j++) {
            g = deb_9.f [j+1] (x);
            printf ("c: %e\n", g);
            if (g > 0) {
                s += g;
            }
        }
        printf ("constraints: %e\n", s);
        printf ("evaluation:  %e\n", deb_9.f [0] (x));
    }
}
#endif
