/*
 *  This is a test program for PGAPack.  The objective is to maximize each
 *  allele.  The evaluation function sums all allele values.
 */
#include <pgapack.h>

int myMutation (PGAContext *ctx, int p, int pop, double mr)
{
    int         stringlen, i, v, count;

    stringlen = PGAGetStringLength (ctx);
    count     = 0;

    for (i=stringlen-1; i>=0; i--) {
        if (PGARandomFlip (ctx, mr)) {
            v = PGARandomInterval (ctx, 1, stringlen);
            PGASetIntegerAllele (ctx, p, pop, i, v);
            count++;
        }
    }
    return count;
}

double evaluate (PGAContext *ctx, int p, int pop, double *dummy)
{
    int  stringlen, i, sum;

    stringlen = PGAGetStringLength (ctx);
    sum       = 0;
     
    for (i=stringlen-1; i>=0; i--) {
        sum += PGAGetIntegerAllele (ctx, p, pop, i);
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
        printf (query);
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
        printf (query);
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
    int         len, maxiter, rtr = 0, nam = 0, nodup = 0, verbose = 0;

    MPI_Init (&argc, &argv);

    len     = GetIntegerParameter ("String length?\n");
    maxiter = GetIntegerParameter ("How many iterations?\n");
    rtr     = GetYN ("Use restricted tournament replacement?\n");
    nam     = GetYN ("Use negative assortative mating (NAM)?\n");
    nodup   = GetYN ("Avoid duplicates?\n");
    verbose = GetYN ("Verbose reporting?\n");

    ctx = PGACreate (&argc, argv, PGA_DATATYPE_INTEGER, len, PGA_MAXIMIZE);

    PGASetRandomSeed (ctx, 1);
    PGASetUserFunction (ctx, PGA_USERFUNCTION_MUTATION, (void *)myMutation);
    PGASetIntegerInitPermute (ctx, 1, len);

    PGASetMaxGAIterValue (ctx, maxiter);
    PGASetNumReplaceValue (ctx, 90);
    PGASetMutationAndCrossoverFlag (ctx, PGA_TRUE);
    PGASetPrintOptions (ctx, PGA_REPORT_AVERAGE);
    if (rtr) {
       PGASetPopReplaceType (ctx, PGA_POPREPL_RTR);
    }
    if (nam) {
       PGASetNAMWindowSize (ctx, 10);
    }
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

    MPI_Finalize ();

    return 0;
}
