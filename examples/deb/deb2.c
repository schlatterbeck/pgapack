#include "constraint.h"
#include <stdio.h>

/* Pre-computation relies on g1 being called first */
static const double a [] =
    { 17.505, 11.275, 214.228, 7.458, 0.961, 1.612, 0.146
    , 107.99, 922.693, 926.832, 18.766, 1072.163, 8961.448, 0.063
    , 71084.33, 2802713.0
    };
static const double b [] =
    { 1053.6667, 35.03, 665.585, 584.463, 265.916, 7.046
    , 0.222, 273.366, 1286.105, 1444.046, 537.141, 3247.039, 26844.086
    , 0.386, 14e4, 12146108.0
    };
static double y [17];
static double c [17];

static void init_y_c (double *x)
{
    /* Dependencies on X[i] after formulas */
    /* Deb has a typo here, he has x1 + x2 which would translate to
     * x [0] + x [1] in this implementation
     */
    y  [0] = x [1] + x [2] + 41.6;                            /* 0 1       */
    c  [0] = 0.024 * x [3] - 4.62;                            /*       3   */
    y  [1] = 12.5 / c [0] + 12.0;                             /*       3   */
    c  [1] = (3.535e-4 * x [0] + 0.5311) * x [0]              /* 0     3   */
           + 0.08705 * y [1] * x [0];                    
    c  [2] = 0.052 * x [0] + 78.0 + 0.002377 * y [1] * x [0]; /* 0     3   */
    y  [2] = c [1] / c [2];                                   /* 0     3   */
    y  [3] = 19.0 * y [2];                                    /* 0     3   */
    c  [3] = 0.04782 * (x [0] - y [2])                        /* 0 1   3   */
           + 0.1956 * pow (x [0] - y [2], 2) / x [1]
           + 0.6376 * y [3] + 1.594 * y [2];
    c  [4] = 100.0 * x [1];                                   /* 1         */
    c  [5] = x [0] - y [2] - y [3];                           /* 0     3   */
    c  [6] = 0.95 - c [3] / c [4];                            /* 0 1   3   */
    y  [4] = c [5] * c [6];                                   /* 0 1   3   */
    y  [5] = x [0] - y [4] - y [3] - y [2];                   /* 0 1   3   */
    c  [7] = 0.995 * (y [3] + y [4]);                         /* 0 1   3   */
    y  [6] = c [7] / y [0];                                   /* 0 1   3   */
    y  [7] = c [7] / 3798.0;                                  /* 0 1   3   */
    c  [8] = y [6] - 0.0663 * y [6] / y [7] - 0.3153;         /* 0 1   3   */
    y  [8] = 96.82 / c [8] + 0.321 * y [0];                   /* 0 1   3   */
    y  [9] = 1.29 * y [4] + 1.258 * y [3] + 2.29 * y [2]      /* 0 1   3   */
           + 1.71 * y [5];
    y [10] = 1.71 * x [0] - 0.452 * y [3] + 0.58 * y [2];     /* 0     3   */
    c  [9] = 12.3 / 752.3;                                    /*           */
    c [10] = 1.75 * y [1] * 0.995 * x [0];                    /* 0     3   */
    c [11] = 0.995 * y [9] + 1998.0;                          /* 0 1   3   */
    y [11] = c [9] * x [0] + c [10] / c [11];                 /* 0 1   3   */
    y [12] = c [11] - 1.75 * y [1];                           /* 0 1   3   */
    y [13] = 3623.0 + 64.4 * x [1] + 58.4 * x [2]             /* 0 1 2 3 4 */
           + 146312.0 / (y [8] + x [4]);
    c [12] = 0.995 * y [9] + 60.8 * x [1] + 48.0 * x [3]      /* 0 1 2 3 4 */
           - 0.1121 * y [13] - 5095.0;
    y [14] = y [12] / c [12];                                 /* 0 1 2 3 4 */
    y [15] = 148e3 - 331e3 * y [14] + 40 * y [12]             /* 0 1 2 3 4 */
           - 61.0 * y [14] * y [12];
    c [13] = 2324.0 * y [9] - 2874e4 * y [1];                 /* 0 1   3   */
    y [16] = 1413e4 - 1328.0 * y [9] - 531.0 * y [10]         /* 0 1   3   */
           + c [13] / c [11];
    c [14] = y [12] / y [14] - y [12] / 0.52;                 /* 0 1 2 3 4 */
    c [15] = 1.104 - 0.72 * y [14];                           /* 0 1 2 3 4 */
    c [16] = y [8] + x [4];                                   /* 0 1   3 4 */
}

static double f (double *x)
{
    /* 0 1 2 3 4 */

    /* Deb has 0.004324 as factor before y [4] while Schittkowski has
     * 0.00423, we're using the original. Looks like Deb found a better
     * evaluation due to that typo.
     */
    return 0.1365 - 5.843e-7 * y [16] + 1.17e-4 * y [13] + 2.358e-5 * y [12]
         + 1.502e-6 * y [15] + 0.0321 * y [11] + 0.004324 * y [4]
         + 1e-4 * (c [14] / c [15]) + 37.48 * (y [1] / c [11]);
}

static double g1 (double *x)
{
    init_y_c (x);
    return x [2] - 1.5 * x [1];
}

static double g2 (double *x)
{
    /* y0 */
    return 213.1 - y [0];
}

static double g3 (double *x)
{
    /* y0 */
    return y [0] - 405.23;
}

static inline double g_j2 (double *x, int j)
{
    return a [j - 2] - y [j - 1];
}

static inline double g_j18 (double *x, int j)
{
    return y [j - 1] - b [j - 2];
}

static inline double g36 (double *x)
{
    /* 0 1   3   */
    return 0.28 / 0.72 * y [4] - y [3];
}

static inline double g37 (double *x)
{
    /* 0 1   3   */
    return 3496.0 * y [1] / c [11] - 21;
}

static inline double g38 (double *x)
{
    /* 0 1   3 4 */
    return 110.6 + y [0] - 62212.0 / c [16];
}

#define g_j2_instatiate(j) \
static double g_j2_ ## j (double *x) \
{ \
    return g_j2 (x, j); \
}
g_j2_instatiate(2)
g_j2_instatiate(3)
g_j2_instatiate(4)
g_j2_instatiate(5)
g_j2_instatiate(6)
g_j2_instatiate(7)
g_j2_instatiate(8)
g_j2_instatiate(9)
g_j2_instatiate(10)
g_j2_instatiate(11)
g_j2_instatiate(12)
g_j2_instatiate(13)
g_j2_instatiate(14)
g_j2_instatiate(15)
g_j2_instatiate(16)
g_j2_instatiate(17)

#define g_j18_instatiate(j) \
static double g_j18_ ## j (double *x) \
{ \
    return g_j18 (x, j); \
}
g_j18_instatiate(2)
g_j18_instatiate(3)
g_j18_instatiate(4)
g_j18_instatiate(5)
g_j18_instatiate(6)
g_j18_instatiate(7)
g_j18_instatiate(8)
g_j18_instatiate(9)
g_j18_instatiate(10)
g_j18_instatiate(11)
g_j18_instatiate(12)
g_j18_instatiate(13)
g_j18_instatiate(14)
g_j18_instatiate(15)
g_j18_instatiate(16)
g_j18_instatiate(17)

struct constrained_problem deb_2 =
{ .dimension      = 5
, .nfunc          = 39
, .lower          = (double []){ 704.4148,   68.6,    0.0,  193.0,    25.0 }
, .upper          = (double []){ 906.3855,  288.88, 134.75, 287.0966, 84.1988}
, .enforce_bounds = 1
, .popsize        = 10
, .generations    = 600
, .f              = { f, g1, g2, g3
                    , g_j2_2, g_j2_3, g_j2_4, g_j2_5, g_j2_6, g_j2_7
                    , g_j2_8, g_j2_9 , g_j2_10, g_j2_11, g_j2_12, g_j2_13
                    , g_j2_14, g_j2_15 , g_j2_16, g_j2_17
                    , g_j18_2, g_j18_3, g_j18_4, g_j18_5, g_j18_6, g_j18_7
                    , g_j18_8, g_j18_9, g_j18_10, g_j18_11, g_j18_12
                    , g_j18_13, g_j18_14, g_j18_15, g_j18_16, g_j18_17
                    , g36, g37, g38
                    }
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    /*
     * According to Schittkowski, Deb
     * Expected values -1.9051338, -1.914595
     * But Deb has constraint violations. Probably because of a typo.
     * Note that the best value obtained with Schittkowski's optimizer is
     * -1.9051553 which is *below* the value of the given solution here.
     */
    static double xx [][5] =
//        { { 705.1803,   68.60005,  102.90001,  282.324999, 37.5850413 }
        { { 0.705180328772e+03, 0.686000529425e+02, 0.102900013236e+03
          , 0.282324998587e+03, 0.375850413432e+02
          }
        , { 707.337769, 68.600273, 102.900146, 282.024841, 84.198792  }
        };
    size_t sz = sizeof (xx) / (5 * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double s = 0;
        for (j=1; j<deb_2.nfunc; j++) {
            double val = deb_2.f [j] (x);
            printf ("C %d: %e\n", j, val);
            if (val > 0) {
                s += val;
            }
        }
        printf ("Constraints: %e\n", s);
        printf ("Evaluation:  %e\n", deb_2.f [0] (x));
    }
    for (i=0; i<17; i++) {
        printf ("Y [%d]: %e\n", i, y [i]);
    }
    for (i=0; i<17; i++) {
        printf ("C [%d]: %e\n", i, c [i]);
    }
    for (i=0; i<16; i++) {
        printf ("A [%d]: %e\n", i, a [i]);
    }
    for (i=0; i<16; i++) {
        printf ("B [%d]: %e\n", i, b [i]);
    }
}
#endif
