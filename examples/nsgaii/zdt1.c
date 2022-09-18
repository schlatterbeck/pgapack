#include "optimize.h"
#include <stdio.h>

static double f1 (double *x)
{
    return x [0];
}

static double g (double *x)
{
    int i;
    double s = 0;
    for (i=1; i<zdt1.dimension; i++) {
        s += x [i];
    }
    s *= 9. / (zdt1.dimension - 1.0);
    return 1. + s;
}

static double f2 (double *x)
{
    double rg = g (x);
    return rg * (1 - sqrt (x [0] / rg));
}

struct multi_problem zdt1 =
{ .dimension      = 30
, .nfunc          = 2
, .nconstraint    = 0
, .lower          = (double [])
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  }
, .upper          = (double [])
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  }
, .enforce_bounds = 1
, .f              = { &f1, &f2 }
, .name           = "Zitzler et. al. (ZDT1)"
};

#ifdef DEBUG_EVAL
#include <stdio.h>
int main ()
{
    int i, j;
    static double xx [][30] =
        { {  0.08322629,  0.49475,    0.9394394,  0.4029122,  0.380806
          ,  0.9450716,   0.5945786,  0.3606442,  0.4050519,  0.4009566
          ,  0.3002905,   0.4911561,  0.837816,   1.214143,   0.8872529
          ,  0.4067645,   0.8997738,  0.7513624,  0.6462557,  0.3775246
          ,  0.4569765,   1.287615,   0.8355889,  0.7239505,  0.7427183
          ,  0.6415458,   0.5772695,  0.2390083,  0.2099733,  0.0463388
          }
        };
    size_t sz = sizeof (xx) / (zdt1.dimension * sizeof (double));
    for (i=0; i<sz; i++) {
        double *x = xx [i];
        double f;
        double s = 0;
        printf ("%d\n", i);
        for (j=0; j<zdt1.nfunc; j++) {
            f = zdt1.f [j] (x);
            printf ("f%d: %e\n", j, f);
        }
    }
}
#endif
