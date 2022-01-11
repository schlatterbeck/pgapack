
/*  Constrained function optimizer
 *  Functions taken from Deb and Jain, 2014 (both papers) see README.rst
 */
#include "constraint.h"

static struct multi_problem *problems [] =
{ &dtlz1
, &dtlz2
//, &dtlz3
//, &dtlz4
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
    memcpy (aux, result + 1, sizeof (double) * problem->nfunc - 1);

    return result [0];
}

int main (int argc, char **argv)
{
    PGAContext *ctx;
    int popsize = 100;
    int fidx = 0;
    int maxiter = 300;
    int sum_constraints = PGA_FALSE;
    void *p = NULL;
    int np = LIN_dasdennis (3, 12, &p, 0, 1, NULL);

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
    printf ("Example: %s", problem->name);
    if (problem->nconstraint > 0) {
        printf (" sum constraints: %s", sum_constraints ? "yes" : "no");
    }
    printf ("\n");
    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, problem->dimension, PGA_MINIMIZE);
    
    //PGASetRandomSeed                (ctx, 1);
    PGASetPopSize                   (ctx, popsize);
# if 0
    PGASetNumReplaceValue           (ctx, 60);
    //PGASetSelectType                (ctx, PGA_SELECT_TOURNAMENT);
    PGASetSelectType             (ctx, PGA_SELECT_LINEAR);
    //PGASetTournamentWithReplacement (ctx, PGA_FALSE);
    //PGASetTournamentSize            (ctx, 1.5);
    PGASetPopReplaceType            (ctx, PGA_POPREPL_NSGA_II);
    PGASetMutationOnlyFlag          (ctx, PGA_FALSE);
    //PGASetMutationAndCrossoverFlag  (ctx, PGA_TRUE);
    PGASetMutationType              (ctx, PGA_MUTATION_POLY);
    PGASetMutationProb              (ctx, 1.0 / problem->dimension);
    PGASetMutationPolyEta           (ctx, 20);
    PGASetMutationBounceBackFlag    (ctx, PGA_TRUE);
    PGASetCrossoverType             (ctx, PGA_CROSSOVER_SBX);
    PGASetCrossoverProb             (ctx, 1.0);
    PGASetUniformCrossoverProb      (ctx, 1.0 / problem->dimension);
    PGASetUniformCrossoverProb      (ctx, 0.2);
    PGASetCrossoverBounceBackFlag   (ctx, PGA_TRUE);
    PGASetCrossoverSBXEta           (ctx, 30);
    PGASetCrossoverSBXEta           (ctx, 15);
    //PGASetCrossoverSBXOncePerString (ctx, PGA_TRUE);
    PGASetRealInitRange             (ctx, problem->lower, problem->upper);
    PGASetMaxGAIterValue            (ctx, maxiter);
    PGASetNumAuxEval                (ctx, problem->nfunc - 1);
    PGASetNumConstraint             (ctx, problem->nconstraint);
    //PGASetSumConstraintsFlag        (ctx, sum_constraints);
    PGASetNoDuplicatesFlag          (ctx, PGA_TRUE);
    PGASetReferencePoints           (ctx, np, p);
# endif

# if 1
    PGASetNumReplaceValue        (ctx, popsize);
    PGASetSelectType             (ctx, PGA_SELECT_LINEAR);
    PGASetPopReplaceType         (ctx, PGA_POPREPL_NSGA_II);
    PGASetPopReplaceType         (ctx, PGA_POPREPL_NSGA_III);
    PGASetMutationOnlyFlag       (ctx, PGA_TRUE);
    PGASetMutationType           (ctx, PGA_MUTATION_DE);
    PGASetDECrossoverProb        (ctx, 0.0);
    PGASetDECrossoverType        (ctx, PGA_DE_CROSSOVER_BIN);
    PGASetDEVariant              (ctx, PGA_DE_VARIANT_RAND);
    PGASetDEScaleFactor          (ctx, 0.40);
    PGASetDEJitter               (ctx, 0.30);
    /* DOES NOT WORK! PGASetDEDither (ctx, 0.01); */
    PGASetRealInitRange          (ctx, problem->lower, problem->upper);
    PGASetMaxGAIterValue         (ctx, maxiter);
    PGASetNumAuxEval             (ctx, problem->nfunc - 1);
    PGASetNumConstraint          (ctx, problem->nconstraint);
    PGASetSumConstraintsFlag     (ctx, sum_constraints);
    PGASetNoDuplicatesFlag       (ctx, PGA_TRUE);
    PGASetMutationBounceBackFlag (ctx, PGA_TRUE);
    PGASetReferencePoints        (ctx, np, p);
# endif
    
    PGASetUp   (ctx);
    PGARun     (ctx, evaluate);
    PGADestroy (ctx);
    return 0;
}
