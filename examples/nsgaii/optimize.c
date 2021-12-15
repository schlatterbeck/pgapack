
/*  Constrained function optimizer
 *  Functions taken from Deb et. al., 2002, see README.rst
 */
#include <pgapack.h>
#include "constraint.h"

static struct multi_problem *problems [] =
{ &sch
, &fon
, &pol
, &kur
, &zdt1
, &zdt2
, &zdt3
, &zdt4
, &zdt6
, &constr
, &srn
, &tnk
, &water
};
static const int nproblems =
    sizeof (problems) / sizeof (struct multi_problem *);

static struct multi_problem *problem;

double evaluate (PGAContext *ctx, int p, int pop, double *aux)
{
    int i;
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    double *params = (double *)ind->chrom;

    for (i=0; i<problem->nfunc - 1; i++) {
        aux [i] = problem->f [i + 1] (params);
    }
    return problem->f [0] (params);
}

int main (int argc, char **argv)
{
    PGAContext *ctx;
    int popsize = 100;
    int fidx = 0;
    int maxiter = 250;
    int sum_constraints = PGA_FALSE;

    if (argc > 1) {
        fidx = atoi (argv [1]);
        if (fidx < 0 || fidx > nproblems - 1) {
            fprintf
                ( stderr, "Usage: %s [f-index]\nIndex in range 0-%d\n"
                , argv [0], nproblems - 1
                );
            exit (1);
        }
    }
    problem = problems [fidx];
    if (problem->generations > maxiter) {
        maxiter = problem->generations;
    }
    if (0 && problem->popsize > popsize) {
        popsize = problem->popsize;
    }
    if (problem->nconstraint > 0 && argc > 2) {
        sum_constraints = atoi (argv [2]) ? PGA_TRUE : PGA_FALSE;
    }
    printf ("Example: %s", problem->name);
    if (problem->nconstraint > 0) {
        printf (" sum constraints: %s", sum_constraints ? "yes" : "no");
    }
    printf ("\n");
    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, problem->dimension, PGA_MINIMIZE);
    
    PGASetRandomSeed         (ctx, 1);
    PGASetPopSize            (ctx, popsize);
    PGASetNumReplaceValue    (ctx, popsize);
    PGASetSelectType         (ctx, PGA_SELECT_LINEAR);
    PGASetPopReplaceType     (ctx, PGA_POPREPL_NSGA_II);
    PGASetMutationOnlyFlag   (ctx, PGA_TRUE);
    PGASetMutationType       (ctx, PGA_MUTATION_DE);
    PGASetDECrossoverProb    (ctx, 0.8);
    PGASetDECrossoverType    (ctx, PGA_DE_CROSSOVER_BIN);
    PGASetDEVariant          (ctx, PGA_DE_VARIANT_RAND);
    PGASetDEScaleFactor      (ctx, 0.85);
    PGASetRealInitRange      (ctx, problem->lower, problem->upper);
    PGASetMaxGAIterValue     (ctx, maxiter);
    PGASetNumAuxEval         (ctx, problem->nfunc - 1);
    PGASetNumConstraint      (ctx, problem->nconstraint);
    PGASetSumConstraintsFlag (ctx, sum_constraints);
    PGASetNoDuplicatesFlag   (ctx, PGA_TRUE);
    if (problem->enforce_bounds) {
        PGASetMutationBounceBackFlag (ctx, PGA_TRUE);
    };
    
    PGASetUp   (ctx);
    PGARun     (ctx, evaluate);
    PGADestroy (ctx);
    return 0;
}
