
/*  Constrained function optimizer
 *  Functions taken from Deb, 2000, see README.rst
 */
#include <pgapack.h>
#include "constraint.h"

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
    int maxiter = 100;

    problem = &deb_0;
    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, problem->dimension, PGA_MINIMIZE);
    
    PGASetRandomSeed(ctx, 1);

    PGASetPopSize          (ctx, 60);
    PGASetSelectType       (ctx, PGA_SELECT_LINEAR);
    PGASetNumReplaceValue  (ctx, 60);
    PGASetPopReplaceType   (ctx, PGA_POPREPL_PAIRWISE_BEST);
    PGASetMutationOnlyFlag (ctx, PGA_TRUE);
    PGASetMutationType     (ctx, PGA_MUTATION_DE);
    PGASetDECrossoverProb  (ctx, 0.8);
    PGASetDECrossoverType  (ctx, PGA_DE_CROSSOVER_BIN);
    PGASetDEVariant        (ctx, PGA_DE_VARIANT_RAND);
    PGASetDEScaleFactor    (ctx, 0.85);

    PGASetRealInitRange    (ctx, deb_0.lower, deb_0.upper);
    PGASetMaxGAIterValue   (ctx, maxiter);
    PGASetNumAuxEval       (ctx, problem->nfunc - 1);
    
    PGASetUp(ctx);

    PGARun (ctx, evaluate);

    PGADestroy(ctx);
    return 0;
}

