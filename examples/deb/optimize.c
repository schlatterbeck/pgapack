
/*  Constrained function optimizer
 *  Functions taken from Deb, 2000, see README.rst
 */
#include <pgapack.h>
#include "constraint.h"

static struct constrained_problem *problems [] =
{ &deb_0
, &deb_1
, &deb_2
, &deb_3
, &deb_4
, &deb_5
, &deb_6
};
static const int nproblems =
    sizeof (problems) / sizeof (struct constrained_problem *);

static struct constrained_problem *problem;

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
    int popsize = 60;
    int fidx = 0;
    int maxiter = 100;

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
    if (problem->iterations > maxiter) {
        maxiter = problem->iterations;
    }
    if (problem->popsize > popsize) {
        popsize = problem->popsize;
    }
    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, problem->dimension, PGA_MINIMIZE);
    
    PGASetRandomSeed(ctx, 1);

    PGASetPopSize          (ctx, popsize);
    PGASetNumReplaceValue  (ctx, popsize);
    PGASetSelectType       (ctx, PGA_SELECT_LINEAR);
    PGASetPopReplaceType   (ctx, PGA_POPREPL_PAIRWISE_BEST);
    PGASetMutationOnlyFlag (ctx, PGA_TRUE);
    PGASetMutationType     (ctx, PGA_MUTATION_DE);
    PGASetDECrossoverProb  (ctx, 0.8);
    PGASetDECrossoverType  (ctx, PGA_DE_CROSSOVER_BIN);
    PGASetDEVariant        (ctx, PGA_DE_VARIANT_RAND);
    PGASetDEScaleFactor    (ctx, 0.85);
    if (problem->enforce_bounds) {
        PGASetMutationBounceBackFlag (ctx, PGA_TRUE);
    };

    PGASetRealInitRange    (ctx, problem->lower, problem->upper);
    PGASetMaxGAIterValue   (ctx, maxiter);
    PGASetNumAuxEval       (ctx, problem->nfunc - 1);
    
    PGASetUp(ctx);

    PGARun (ctx, evaluate);

    PGADestroy(ctx);
    return 0;
}

