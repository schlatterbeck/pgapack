#include <math.h>
#include "pgapack.h"

/*
 * We store the dimensionality of the problem (the number of variables)
 * and the number of functions. There are multiple evaluations and
 * optionally some constraints. The number of evaluations is
 * nfunc - nconstraint
 * Constraints are defined as being minimized, so
 * constraints must be <= 0, i.e.,   gn(X) <= 0 for all n
 */
struct multi_problem
{
    int dimension;              /* The dimension of the problem */
    int nfunc;                  /* Number of evaluation functions */
    int nconstraint;            /* Number of constraints */
    int maximize;               /* Maximization problem ? */
    double (*lower);            /* Init ranges lower bounds */
    double (*upper);            /* Init ranges upper bounds */
    int generations;            /* Number of generations if 0 default = 100 */
    int popsize;                /* Population size if 0 default = 60 */
    char *name;                 /* Name of the example */
    void (*f)(double *x, double *y); /* input: x, output: y */
};

extern struct multi_problem dtlz1;
extern struct multi_problem dtlz2;
extern struct multi_problem dtlz3;
extern struct multi_problem dtlz4;
extern struct multi_problem scaled_dtlz1;
extern struct multi_problem scaled_dtlz2;
extern struct multi_problem convex_dtlz2;
extern struct multi_problem neg_dtlz2;
extern struct multi_problem c1_dtlz1;
extern struct multi_problem c1_dtlz3;
extern struct multi_problem c2_dtlz2;
extern struct multi_problem c2_convex_dtlz2;
extern struct multi_problem c3_dtlz1;
extern struct multi_problem c3_dtlz4;
