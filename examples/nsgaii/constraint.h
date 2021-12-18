#include <math.h>

/* The constrained_problem contains one specific problem:
 * We store the dimensionality of the problem (the number of variables)
 * and the number of functions. The first function is always the
 * objective function. The rest are constraints. Constraints are --
 * contrary to Deb -- defined as being minimized, so
 * constraints must be <= 0, i.e.,   gn(X) <= 0 for all n
 */
struct multi_problem
{
    int dimension;              /* The dimension of the problem */
    int nfunc;                  /* Number of constraints + 1 (objective f) */
    int nconstraint;            /* Number of constraints */
    double (*lower);            /* Init ranges lower bounds */
    double (*upper);            /* Init ranges upper bounds */
    int enforce_bounds;         /* Enforce bounds on init range */
    int generations;            /* Number of generations if 0 default = 100 */
    int popsize;                /* Population size if 0 default = 60 */
    char *name;                 /* Name of the example */
    double dither;              /* Dither for differential evolution */
    double jitter;              /* Jitter for differential evolution */
    double crossover_prob;      /* Crossover probability for DE */
    double (*f [])(double *);   /* Functions: 0th is objective function */
};

extern struct multi_problem fon;
extern struct multi_problem kur;
extern struct multi_problem pol;
extern struct multi_problem sch;
extern struct multi_problem zdt1;
extern struct multi_problem zdt2;
extern struct multi_problem zdt3;
extern struct multi_problem zdt4;
extern struct multi_problem zdt6;
extern struct multi_problem constr;
extern struct multi_problem srn;
extern struct multi_problem tnk;
extern struct multi_problem water;
extern struct multi_problem rotated;
