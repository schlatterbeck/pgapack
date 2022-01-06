#include <assert.h>
#include "pgapack.h"

/*
 * Ackley's function used in [1]
 * [1] Kalyanmoy Deb and Debayan Deb. Analysing mutation schemes for
 *     real-parameter genetic algorithms. International Journal of
 *     Artificial Intelligence and Soft Computing, 4(1):1–28, February 2014.
 * Note that this is not a constrained problem, it's included here for
 * testing the simulated binary crossover (SBX) and polynomial mutation.
 */

double evaluate (PGAContext *ctx, int p, int pop, double *dummy)
{
    int i;
    int dim = PGAGetStringLength (ctx);
    double s1 = 0, s2 = 0;
    assert (dim == 15);

    for (i=0; i<dim; i++) {
        double xi = PGAGetRealAllele (ctx, p, pop, i);
        s1 += xi * xi;
        s2 += cos (2 * M_PI * xi);
    }
    s1 /= dim;
    s2 /= dim;
    return 20 + M_E - 20 * exp (-0.2 * sqrt (s1)) - exp (s2);
}

int stop_cond (PGAContext *ctx)
{
    double best = PGAGetBestReport (ctx, PGA_OLDPOP, 0);
    if (best <= 0.01) {
        return PGA_TRUE;
    }
    return PGA_FALSE;
}

double lower [] = {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5};
double upper [] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

int main( int argc, char **argv )
{
    PGAContext *ctx; 
    int seed = 0;

    if (argc > 1 && atoi (argv [1]) > 0) {
        seed = atoi (argv [1]);
    }

    ctx = PGACreate (&argc, argv, PGA_DATATYPE_REAL, 15, PGA_MINIMIZE);
    PGASetPopSize                   (ctx, 150);
    PGASetNumReplaceValue           (ctx, 100);
    if (seed) {
        PGASetRandomSeed            (ctx, 1);
    }
    PGASetRealInitRange             (ctx, lower, upper);
    PGASetCrossoverType             (ctx, PGA_CROSSOVER_SBX);
    PGASetCrossoverProb             (ctx, 1.0);
    PGASetUniformCrossoverProb      (ctx, 0.9);
    PGASetCrossoverSBXEta           (ctx, 2);
    PGASetCrossoverBounceBackFlag   (ctx, PGA_TRUE);
    /* This may be a good idea for non-separable functions */
    /* PGASetCrossoverSBXOncePerString (ctx, PGA_TRUE); */
    PGASetMutationType              (ctx, PGA_MUTATION_POLY);
    PGASetMutationProb              (ctx, 1.0 / 15);
    PGASetMutationPolyEta           (ctx, 20);
    PGASetMutationAndCrossoverFlag  (ctx, PGA_TRUE);
    PGASetMutationBounceBackFlag    (ctx, PGA_TRUE);
    PGASetUserFunction              (ctx, PGA_USERFUNCTION_STOPCOND, stop_cond);
    PGASetUp                        (ctx);
    PGARun                          (ctx, evaluate);
    PGADestroy                      (ctx);
    return 0;
}
