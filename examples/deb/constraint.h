#include <math.h>

/* The constrained_problem contains one specific problem:
 * We store the dimensionality of the problem (the number of variables)
 * and the number of functions. The first function is always the
 * objective function. The rest are constraints. Constraints are --
 * contrary to Deb -- defined as being minimized, so
 * constraints must be <= 0, i.e.,   gn(X) <= 0 for all n
 */
struct constrained_problem
{
    int dimension;              /* The dimension of the problem */
    int nfunc;                  /* Number of constraints + 1 (objective f) */
    double (*lower);            /* Init ranges lower bounds */
    double (*upper);            /* Init ranges upper bounds */
    int enforce_bounds;         /* Enforce bounds on init range */
    int iterations;             /* Number of iterations if 0 default = 100 */
    int popsize;                /* Population size if 0 default = 60 */
    double (*f [])(double *);   /* Functions: 0th is objective function */
};

extern struct constrained_problem deb_0;
extern struct constrained_problem deb_1;
extern struct constrained_problem deb_2;
extern struct constrained_problem deb_3;
extern struct constrained_problem deb_4;
extern struct constrained_problem deb_5;
extern struct constrained_problem deb_6;
extern struct constrained_problem deb_7;
extern struct constrained_problem deb_8;
extern struct constrained_problem deb_9;
