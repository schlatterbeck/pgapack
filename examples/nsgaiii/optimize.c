
/*  Constrained function optimizer
 *  Functions taken from Deb and Jain, 2014 (both papers) see README.rst
 */
#include <unistd.h>
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

void usage (char *name, int nproblems)
{
    fprintf
        ( stderr
        , "Usage: %s [-2] [-d] [-g maxiter] [-p popsize] [-r seed] "
          "[-s] [f-index]\n"
          "-2: Use NSGA-II not NSGA-III\n"
          "-d: Use reference directions\n"
          "-g: Maximum number of generations\n"
          "-n: Set old (quadratic) non-dominated sorting algorithm\n"
          "-p: Population size\n"
          "-r: Random seed (uppercase -R is also accepted)\n"
          "-s: Use SBX crossover and Polynomial mutation instead of DE\n"
          "f-index is the function to call in range 0-%d\n"
        , name, nproblems - 1
        );
}

int main (int argc, char **argv)
{
    PGAContext *ctx;
    int popsize = 100;
    int popsize_seen = 0;
    int fidx = 0;
    int maxiter = 400;
    int maxiter_seen = 0;
    void *p = NULL;
    int np = LIN_dasdennis (3, 12, &p, 0, 1, NULL);
    int seed = 1;
    int direction;
    int refdir = 0;
    double directions [][3] = {{0.5, 5, 50}, {0.8, 2, 20}};
    #define NDIR (sizeof (directions) / (3 * sizeof (double)))
    MPI_Comm comm;
    int opt;
    int repl_type = PGA_POPREPL_NSGA_III;
    int use_de = 1;
    int nondom = PGA_NDSORT_JENSEN;

    while ((opt = getopt (argc, argv, "2dg:r:R:s")) != -1) {
        switch (opt) {
        case '2':
            repl_type = PGA_POPREPL_NSGA_II;
            break;
        case 'd':
            refdir = 1;
            break;
        case 'g':
            maxiter = atoi (optarg);
            maxiter_seen = 1;
            break;
        case 'n':
            nondom = PGA_NDSORT_NSQUARE;
            break;
        case 'p':
            popsize = atoi (optarg);
            popsize_seen = 1;
            break;
        case 'r':
        case 'R':
            seed = atoi (optarg);
            break;
        case 's':
            /* Use SBX crossover/Poly mutation instead of DE */
            use_de = 0;
            break;
        default:
            usage (argv [0], nproblems);
            exit (1);
        }
    }

    if (optind < argc) {
        fidx = atoi (argv [optind]);
        if (fidx < 0 || fidx > nproblems - 1) {
            usage (argv [0], nproblems);
            exit (1);
        }
    }
    problem = problems [fidx];
    if (!maxiter_seen && problem->generations != 0) {
        maxiter = problem->generations;
    }
    if (!popsize_seen && problem->popsize > popsize) {
        popsize = problem->popsize;
    }
    direction = problem->maximize ? PGA_MAXIMIZE : PGA_MINIMIZE;
    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, problem->dimension, direction);
    
    PGASetRandomSeed              (ctx, seed);
    PGASetPopSize                 (ctx, popsize);
    if (use_de) {
        PGASetSelectType          (ctx, PGA_SELECT_LINEAR);
        PGASetMixingType          (ctx, PGA_MIX_MUTATE_ONLY);
        PGASetMutationType        (ctx, PGA_MUTATION_DE);
        PGASetDECrossoverProb     (ctx, 0.0);
        PGASetDECrossoverType     (ctx, PGA_DE_CROSSOVER_BIN);
        PGASetDEVariant           (ctx, PGA_DE_VARIANT_RAND);
        PGASetDEScaleFactor       (ctx, 0.40);
        PGASetDEJitter            (ctx, 0.30);
        PGASetNumReplaceValue     (ctx, popsize);
        /* DOES NOT WORK! PGASetDEDither (ctx, 0.01); */
    } else {
        PGASetMutationType        (ctx, PGA_MUTATION_POLY);
        PGASetCrossoverType       (ctx, PGA_CROSSOVER_SBX);
        PGASetDECrossoverProb     (ctx, 0.9);
        PGASetCrossoverSBXEta     (ctx, 20);
        PGASetMutationPolyEta     (ctx, 20);
        // FIXME
        PGASetNumReplaceValue     (ctx, popsize);
    }

    PGASetSortND                  (ctx, nondom);
    PGASetPopReplaceType          (ctx, repl_type);
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
    /* p is freed by PGADestroy if we called PGASetReferencePoints */
    if (refdir) {
        free (p);
    }
    return 0;
}
