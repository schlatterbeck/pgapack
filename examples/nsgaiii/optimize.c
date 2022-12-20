
/*  Constrained function optimizer
 *  Functions taken from Deb and Jain, 2014 (both papers) see README.rst
 */
#include "optimize.h"

static struct multi_problem *problems [] =
{ &dtlz1
, &dtlz2
, &dtlz3
, &dtlz4
, &scaled_dtlz1
, &scaled_dtlz2
, &convex_dtlz2
, &neg_dtlz2
, &c1_dtlz1
, &c1_dtlz3
, &c2_dtlz2
, &c2_convex_dtlz2
, &c3_dtlz1
, &c3_dtlz4
};
static const int nproblems =
    sizeof (problems) / sizeof (struct multi_problem *);

static struct multi_problem *problem;

double evaluate (PGAContext *ctx, int p, int pop, double *aux)
{
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    double *params = (double *)ind->chrom;
    double result [problem->nfunc];

    problem->f (params, result);
    memcpy (aux, result + 1, sizeof (double) * (problem->nfunc - 1));

    return result [0];
}

int main (int argc, char **argv)
{
    PGAContext *ctx;
    int popsize = 100;
    int fidx = 0;
    int maxiter = 400;
    void *p = NULL;
    int np = LIN_dasdennis (3, 12, &p, 0, 1, NULL);
    int seed = 1;
    int direction;
    int refdir = 0;
    double directions [][3] = {{0.5, 5, 50}, {0.8, 2, 20}};
    #define NDIR (sizeof (directions) / (3 * sizeof (double)))
    MPI_Comm comm;

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
    if (argc > 2) {
        seed = atoi (argv [2]);
    }
    if (argc > 3) {
        refdir = 1;
    }
    problem = problems [fidx];
    if (problem->generations != 0) {
        maxiter = problem->generations;
    }
    if (problem->popsize > popsize) {
        popsize = problem->popsize;
    }
    direction = problem->maximize ? PGA_MAXIMIZE : PGA_MINIMIZE;
    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, problem->dimension, direction);
    
    PGASetRandomSeed                (ctx, seed);
    PGASetPopSize                   (ctx, popsize);

    PGASetNumReplaceValue         (ctx, popsize);
    PGASetSelectType              (ctx, PGA_SELECT_LINEAR);
    PGASetPopReplaceType          (ctx, PGA_POPREPL_NSGA_II);
    PGASetPopReplaceType          (ctx, PGA_POPREPL_NSGA_III);
    PGASetMutationOnlyFlag        (ctx, PGA_TRUE);
    PGASetMutationType            (ctx, PGA_MUTATION_DE);
    PGASetDECrossoverProb         (ctx, 0.0);
    PGASetDECrossoverType         (ctx, PGA_DE_CROSSOVER_BIN);
    PGASetDEVariant               (ctx, PGA_DE_VARIANT_RAND);
    PGASetDEScaleFactor           (ctx, 0.40);
    PGASetDEJitter                (ctx, 0.30);
    /* DOES NOT WORK! PGASetDEDither (ctx, 0.01); */
    PGASetRealInitRange           (ctx, problem->lower, problem->upper);
    PGASetMaxGAIterValue          (ctx, maxiter);
    PGASetNumAuxEval              (ctx, problem->nfunc - 1);
    PGASetNumConstraint           (ctx, problem->nconstraint);
    PGASetSumConstraintsFlag      (ctx, PGA_FALSE);
    PGASetNoDuplicatesFlag        (ctx, PGA_TRUE);
    PGASetMutationBounceBackFlag  (ctx, PGA_TRUE);
    if (refdir) {
        PGASetReferenceDirections (ctx, NDIR, directions, 7, 0.05);
    } else {
        PGASetReferencePoints     (ctx, np, p);
    }
    if (fidx == 11) {
        /* Avoid regression test failing due to rounding error on different
         * architectures. We reduce the precision printed by one.
         */
        PGASetMultiObjPrecision (ctx, 13);
    }
    
    PGASetUp   (ctx);
    comm = PGAGetCommunicator (ctx);
    if (PGAGetRank (ctx, comm) == 0) {
        printf ("Example: %s\n", problem->name);
    }
    PGARun     (ctx, evaluate);
    PGADestroy (ctx);
    return 0;
}
