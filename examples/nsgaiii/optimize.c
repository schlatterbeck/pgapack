
/*  Constrained function optimizer
 *  Functions taken from Deb and Jain, 2014 (both papers) see README.rst
 */
#include <unistd.h>
#include <time.h>
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
static int nfunc = 0;
static int dimension = 0;
static clock_t t_ranking = 0;
static unsigned int (*SortND)
    (PGAContext *, PGAIndividual **, size_t, int) = NULL;

double evaluate (PGAContext *ctx, int p, int pop, double *aux)
{
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    double *params = (double *)ind->chrom;
    double result [nfunc];

    problem->f (params, dimension, result, nfunc);
    memcpy (aux, result + 1, sizeof (double) * (nfunc - 1));

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
          "-D: Dimension (number of decision variables)\n"
          "-g: Maximum number of generations\n"
          "-n: Set old (quadratic) non-dominated sorting algorithm\n"
          "-o: Number of objectives\n"
          "-p: Population size\n"
          "-r: Random seed (uppercase -R is also accepted)\n"
          "-s: Use SBX crossover and Polynomial mutation instead of DE\n"
          "-t: Measure timing of non-dominated sorting algo and write to file\n"
          "f-index is the function to call in range 0-%d\n"
        , name, nproblems - 1
        );
}

unsigned int ranking_plugin
    (PGAContext *ctx, PGAIndividual **start, size_t n, int goal)
{
    unsigned int retval = 0;
    clock_t t = clock ();
    retval = SortND (ctx, start, n, goal);
    t_ranking += clock () - t;
    return retval;
}

int main (int argc, char **argv)
{
    PGAContext *ctx;
    int i;
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
    int n_obj = 0;
    double *lower = NULL;
    double *upper = NULL;
    char *timing_file = NULL;

    while ((opt = getopt (argc, argv, "2dD:g:no:p:r:R:st:")) != -1) {
        switch (opt) {
        case '2':
            repl_type = PGA_POPREPL_NSGA_II;
            break;
        case 'd':
            refdir = 1;
            break;
        case 'D':
            dimension = atoi (optarg);
            break;
        case 'g':
            maxiter = atoi (optarg);
            maxiter_seen = 1;
            break;
        case 'n':
            nondom = PGA_NDSORT_NSQUARE;
            break;
        case 'o':
            n_obj = atoi (optarg);
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
        case 't':
            timing_file = optarg;
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
    if (n_obj <= 0) {
        n_obj = problem->n_obj;
    } else if (dimension <= 0) {
        /* DTLZ suggest k = 5, i.e. dim = n_obj + 5 - 1
         * But this differs by problem, compute from problem
         */
        dimension = n_obj + problem->dimension - problem->n_obj;
    }
    if (dimension <= 0) {
        dimension = problem->dimension;
    }
    if (dimension <= n_obj) {
        fprintf (stderr, "Need dimension >= n_obj\n");
        exit (1);
    }
    nfunc = n_obj + problem->nconstraint;

    if (!maxiter_seen && problem->generations != 0) {
        maxiter = problem->generations;
    }
    if (!popsize_seen && problem->popsize > popsize) {
        popsize = problem->popsize;
    }
    direction = problem->maximize ? PGA_MAXIMIZE : PGA_MINIMIZE;
    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, dimension, direction);

    lower = malloc (sizeof (*lower) * dimension);
    upper = malloc (sizeof (*upper) * dimension);
    if (lower == NULL || upper == NULL) {
        PGAFatalPrintf (ctx, "Out of memory allocating bounds");
    }
    for (i=0; i<dimension; i++) {
        lower [i] = problem->lower;
        upper [i] = problem->upper;
    }
    
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
    if (repl_type == PGA_POPREPL_NSGA_III) {
        if (refdir) {
            PGASetReferenceDirections (ctx, NDIR, directions, 7, 0.05);
        } else {
            PGASetReferencePoints     (ctx, np, p);
        }
    }

    PGASetSortND                  (ctx, nondom);
    PGASetPopReplaceType          (ctx, repl_type);
    PGASetRealInitRange           (ctx, lower, upper);
    PGASetMaxGAIterValue          (ctx, maxiter);
    PGASetNumAuxEval              (ctx, nfunc - 1);
    PGASetNumConstraint           (ctx, problem->nconstraint);
    PGASetSumConstraintsFlag      (ctx, PGA_FALSE);
    PGASetNoDuplicatesFlag        (ctx, PGA_TRUE);
    PGASetMutationBounceBackFlag  (ctx, PGA_TRUE);
    PGASetCrossoverBounceBackFlag (ctx, PGA_TRUE);
    if (fidx == 11) {
        /* Avoid regression test failing due to rounding error on different
         * architectures. We reduce the precision printed by one.
         */
        PGASetMultiObjPrecision (ctx, 13);
    }
    
    PGASetUp   (ctx);
    if (timing_file != NULL) {
        SortND = ctx->cops.SortND;
        ctx->cops.SortND = ranking_plugin;
    }
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
    free (lower);
    free (upper);
    if (timing_file != NULL) {
        FILE *tf = fopen (timing_file, "w");
        if (tf == NULL) {
            fprintf (stderr, "Cannot open timing file\n");
        }
        fprintf
            ( tf
            , "Timing: %.4Lf\n"
            , (long double)(t_ranking) / (long double)CLOCKS_PER_SEC
            );
        fclose (tf);
    }
    return 0;
}
