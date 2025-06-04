
/*  Constrained function optimizer
 *  Functions taken from Deb et. al., 2002, see README.rst
 */
#include <pgapack.h>
#include "optimize.h"

void as_c (PGAContext *ctx, FILE *fp, int p, int pop)
{
    int i, j;
    int popsize = PGAGetPopSize (ctx);
    int dim = PGAGetNumAuxEval (ctx) + 1;
    const double *aux = NULL;
    double eval;
    
    fprintf (fp, "double pop [][%d] = {\n", dim);
    for (j=0; j<popsize; j++) {
        eval = PGAGetEvaluation (ctx, j, pop, &aux);
        
        fprintf (fp, "{%e", eval);
        for (i=0; i<dim-1; i++) {
            fprintf (fp, ", %e", aux [i]);
        }
        fprintf (fp, "}");
        if (j < popsize-1) {
            fprintf (fp, ",\n");
        } else {
            fprintf (fp, "\n");
        }
    }
    fprintf (fp, "};\n");
    fflush (fp);
}

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
, &rotated
, &deb7
, &water_m
, &zdt1_m
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
    int direction;
    double crossover_prob = 0.8;
    int epsilon_generation = 0;
    MPI_Comm comm;
    int c_output = 0;

    if (argc > 1 && argv [1][0] == '-' && argv [1][1] == 'C') {
        c_output = 1;
        argc--;
        argv++;
    }

    if (argc > 1) {
        fidx = atoi (argv [1]);
        if (fidx < 0 || fidx > nproblems - 1) {
            fprintf
                ( stderr, "Usage: %s [-C] [f-index]\nIndex in range 0-%d\n"
                , argv [0], nproblems - 1
                );
            exit (1);
        }
    }
    problem = problems [fidx];
    if (problem->generations) {
        maxiter = problem->generations;
    }
    if (problem->popsize) {
        popsize = problem->popsize;
    }
    if (problem->nconstraint > 0) {
        if (argc > 2) {
            sum_constraints = atoi (argv [2]) ? PGA_TRUE : PGA_FALSE;
        }
        if (argc > 3) {
            epsilon_generation = atoi (argv [3]);
        }
    }
    if (problem->crossover_prob > 0) {
        crossover_prob = problem->crossover_prob;
    }
    direction = problem->maximize ? PGA_MAXIMIZE : PGA_MINIMIZE;
    ctx = PGACreate
        (&argc, argv, PGA_DATATYPE_REAL, problem->dimension, direction);
    
    PGASetRandomSeed         (ctx, 1);
    PGASetPopSize            (ctx, popsize);
    PGASetNumReplaceValue    (ctx, popsize);
    PGASetEpsilonGeneration  (ctx, epsilon_generation);
    PGASetSelectType         (ctx, PGA_SELECT_LINEAR);
    PGASetPopReplaceType     (ctx, PGA_POPREPL_NSGA_II);
    PGASetMutationOnlyFlag   (ctx, PGA_TRUE);
    PGASetMutationType       (ctx, PGA_MUTATION_DE);
    PGASetDECrossoverProb    (ctx, crossover_prob);
    PGASetDECrossoverType    (ctx, PGA_DE_CROSSOVER_BIN);
    PGASetDEVariant          (ctx, PGA_DE_VARIANT_RAND);
    PGASetDEScaleFactor      (ctx, 0.85);
    PGASetRealInitRange      (ctx, problem->lower, problem->upper);
    PGASetMaxGAIterValue     (ctx, maxiter);
    PGASetNumAuxEval         (ctx, problem->nfunc - 1);
    PGASetNumConstraint      (ctx, problem->nconstraint);
    PGASetSumConstraintsFlag (ctx, sum_constraints);
    PGASetNoDuplicatesFlag   (ctx, PGA_TRUE);
    PGASetMultiObjPrecision  (ctx, problem->precision ? problem->precision:14);
    if (problem->dither) {
        PGASetDEDither       (ctx, 0.5);
    }
    if (problem->jitter) {
        PGASetDEJitter       (ctx, 0.005);
    }
    if (problem->enforce_bounds) {
        PGASetMutationBounceBackFlag (ctx, PGA_TRUE);
    };
    
    if (c_output) {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_PRINTSTRING, as_c);
        PGASetPrintFrequencyValue (ctx, 1);
        PGASetPrintOptions (ctx, PGA_REPORT_STRING);
    }
    
    PGASetUp   (ctx);
    comm = PGAGetCommunicator (ctx);
    if (PGAGetRank (ctx, comm) == 0) {
        printf ("Example: %s", problem->name);
        if (problem->nconstraint > 0) {
            printf (" sum constraints: %s", sum_constraints ? "yes" : "no");
        }
        printf ("\n");
    }
    PGARun     (ctx, evaluate);
    PGADestroy (ctx);
    return 0;
}
