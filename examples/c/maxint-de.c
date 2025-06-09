/*
 *  This is a test program for PGAPack.  The objective is to maximize each
 *  allele.  The evaluation function sums all allele values.
 */
#include <pgapack.h>

int maximum = 0;

double evaluate (PGAContext *ctx, int p, int pop, double *dummy)
{
    int  stringlen, i;
    double sum = 0;

    stringlen = PGAGetStringLength (ctx);
     
    for (i=stringlen-1; i>=0; i--) {
        sum += (double)PGAGetIntegerAllele (ctx, p, pop, i);
    }

    return (double)sum;
}

void end_of_gene (PGAContext *ctx)
{
    PGAPrintPopulation (ctx, stdout, PGA_NEWPOP);
}

void pre_eval (PGAContext *ctx, int pop)
{
    static int seen = 0;
    if (!seen) {
        PGAPrintPopulation (ctx, stdout, pop);
        seen = 1;
    }
}

int check_stop (PGAContext *ctx)
{
    double v = 0;
    double m = maximum * PGAGetStringLength (ctx);
    int best = PGAGetBestIndex (ctx, PGA_OLDPOP);
    v = PGAGetEvaluation (ctx, best, PGA_OLDPOP);
    if (v >= m) {
        return PGA_TRUE;
    }
    return PGACheckStoppingConditions (ctx);
}

/*  Get an integer parameter from the user.  Since this is
 *  typically a parallel program, we must only do I/O on the
 *  "master" process -- process 0.  Once we read the parameter,
 *  we broadcast it to all the other processes, then every 
 *  process returns the correct value.
 */
int GetIntegerParameter (char *query)
{
    int  rank, tmp = 0;
    char buf [100];

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf ("%s", query);
        if (fgets (buf, sizeof (buf) - 1, stdin) != NULL) {
            tmp = atoi (buf);
        }
    }
    MPI_Bcast (&tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return tmp;
}

int GetYN (char *query)
{
    int rank, tmp = 0;
    char buf [100];
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf ("%s", query);
        if (fgets (buf, sizeof (buf) - 1, stdin) != NULL) {
            if (buf [0] == 'y' || buf [0] == 'Y') {
                tmp = 1;
            }
        }
    }
    MPI_Bcast (&tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return tmp;
}

int main (int argc, char **argv)
{
    PGAContext *ctx;
    int         len, maxiter, popsize, nodup = 0, verbose = 0;
    int        *lower, *upper;
    int         i;

    MPI_Init (&argc, &argv);

    popsize = GetIntegerParameter ("Population size?\n");
    len     = GetIntegerParameter ("String length?\n");
    maximum = GetIntegerParameter ("Maximum?\n");
    maxiter = GetIntegerParameter ("How many iterations?\n");
    nodup   = GetYN ("Avoid duplicates?\n");
    verbose = GetYN ("Verbose reporting?\n");

    ctx = PGACreate (&argc, argv, PGA_DATATYPE_INTEGER, len, PGA_MAXIMIZE);
    lower = malloc (sizeof (int) * len);
    if (lower == NULL) {
        fprintf (stderr, "Cannot allocate lower bound\n");
        exit (1);
    }
    memset (lower, 0, sizeof (int) * len);
    upper = malloc (sizeof (int) * len);
    if (upper == NULL) {
        fprintf (stderr, "Cannot allocate upper bound\n");
        exit (1);
    }
    for (i=0; i<len; i++) {
        upper [i] = maximum;
    }

    PGASetRandomSeed             (ctx, 1);
    PGASetPopSize                (ctx, popsize);
    PGASetIntegerInitPermute     (ctx, 1, len);
    PGASetMaxGAIterValue         (ctx, maxiter);
    PGASetIntegerInitRange       (ctx, lower, upper);
    PGASetNumReplaceValue        (ctx, popsize);
    PGASetMixingType             (ctx, PGA_MIX_MUTATE_ONLY);
    PGASetMutationType           (ctx, PGA_MUTATION_DE);
    PGASetDECrossoverProb        (ctx, 0.2);
    PGASetDECrossoverType        (ctx, PGA_DE_CROSSOVER_BIN);
    PGASetDEVariant              (ctx, PGA_DE_VARIANT_RAND);
    PGASetDEScaleFactor          (ctx, 0.75);
    PGASetDEDither               (ctx, 0.75);
    PGASetDEDitherPerIndividual  (ctx, PGA_TRUE);
    PGASetMutationBounceBackFlag (ctx, PGA_TRUE);
    PGASetUserFunction           (ctx, PGA_USERFUNCTION_STOPCOND, check_stop);

    if (nodup) {
        PGASetNoDuplicatesFlag (ctx, PGA_TRUE);
    }
    if (verbose) {
        PGASetUserFunction (ctx, PGA_USERFUNCTION_PRE_EVAL, pre_eval);
        PGASetUserFunction (ctx, PGA_USERFUNCTION_ENDOFGEN, end_of_gene);
    }

    PGASetUp (ctx);

    PGARun (ctx, evaluate);
    PGADestroy (ctx);
    free (lower);
    free (upper);

    MPI_Finalize ();

    return 0;
}
