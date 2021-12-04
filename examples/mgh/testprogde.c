#include <stdio.h>
#include <assert.h>
#include "pgapack.h"

/* The fortran subroutine */
extern void objfcn_ (int *, double *, double *, int *);

double x [100];
int nprob = 0;

static const int popsizes [] =
/*    1/10 2/11   3/12  4/13 5/14 6/15 7/16 8/17 9/18 */
    {    6,  50,  2000,   60,  60,   6,   4,  40,   8
    ,    4,   4,     4,   10,  16,   4,   4,   4,  10
    };

/* Initialization ranges: Many badly scaled functions don't work at all
 * if we do not provide a suitable initialization range. The default for
 * is the interval [0, 1].
 */
static const double * init_lower [] =
    { NULL
    , (double []){ 0, 0, 0, 0, 0, 0 }
    , NULL
    , (double []){ 0, 0 }
    , (double []){ -10, -10, -10 }
    , NULL
    , (double []){ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
                 , -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
                 , -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
                 , -1
                 }
    , NULL
    , NULL
    , (double []){ 0, 0 }
    , (double []){ 0, 0, 0, 0 }
    , (double []){ 0, 0, 0 }
    , NULL
    , NULL
    , NULL
    , NULL
    , NULL
    , NULL
    };

static const double * init_upper [] =
    { NULL
    , (double []){ 20, 20, 20, 20, 20, 20 }
    , NULL
    , (double []){ 20, 20 }
    , (double []){ 10, 10, 10 }
    , NULL
    , (double []){ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
                 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
                 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
                 , 1
                 }
    , NULL
    , NULL
    , (double []){ 1e8, 1 }
    , (double []){ 10, 10, 10, 10 }
    , (double []){ 100, 100, 100 }
    , NULL
    , NULL
    , NULL
    , NULL
    , NULL
    , NULL
    };

double evaluate (PGAContext *ctx, int p, int pop)
{
    int i;
    int dimension = PGAGetStringLength (ctx);
    double e;

    for (i=0; i<dimension; i++) {
        x [i] = PGAGetRealAllele (ctx, p, pop, i);
    }
    objfcn_ (&dimension, x, &e, &nprob);
    return e;
}

/*
 * We get 3 arguments:
 * - The problem number
 * - The dimension of the problem
 * - The number of iterations
 */
int main (int argc, char **argv)
{
    int dimension = 0;
    int maxiter   = 0;
    PGAContext *ctx;
    if (argc != 4) {
        fprintf ( stderr
                , "Usage: %s <problem> <dimension> <iterations>\n"
                , argv [0]
                );
        exit (1);
    }
    nprob     = atoi (argv [1]);
    dimension = atoi (argv [2]);
    maxiter   = atoi (argv [3]);
    if (nprob < 1 || nprob > 18) {
        fprintf (stderr, "Problem number 1 <= p <= 18\n");
        exit (1);
    }
    if (dimension < 2 || dimension > 100) {
        fprintf (stderr, "Dimension 2 <= D <= 100\n");
        exit (1);
    }
    if (maxiter < 10) {
        fprintf (stderr, "Iterations >= 10\n");
        exit (1);
    }
    ctx = PGACreate (&argc, argv, PGA_DATATYPE_REAL, dimension, PGA_MINIMIZE);
    PGASetMaxGAIterValue   (ctx, maxiter);
    PGASetPopSize          (ctx, popsizes [nprob - 1]);
    PGASetNumReplaceValue  (ctx, popsizes [nprob - 1]);
    if (init_lower [nprob - 1] != NULL) {
        assert (init_upper [nprob - 1] != NULL);
        PGASetRealInitRange (ctx, init_lower [nprob-1], init_upper [nprob-1]);
    }
    PGASetRandomSeed       (ctx, 1);
    PGASetSelectType       (ctx, PGA_SELECT_LINEAR);
    PGASetPopReplaceType   (ctx, PGA_POPREPL_PAIRWISE_BEST);
    PGASetMutationOnlyFlag (ctx, PGA_TRUE);
    PGASetMutationType     (ctx, PGA_MUTATION_DE);
    PGASetDECrossoverProb  (ctx, 0.2);
    PGASetDECrossoverType  (ctx, PGA_DE_CROSSOVER_BIN);
    PGASetDEVariant        (ctx, PGA_DE_VARIANT_RAND);
    PGASetDEScaleFactor    (ctx, 0.85);
    PGASetUp               (ctx);
    PGARun                 (ctx, evaluate);
    PGADestroy             (ctx);
}
