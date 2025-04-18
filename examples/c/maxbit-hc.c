/*
 *  This is a test program for PGAPack.  The objective is to maximize the
 *  number of 1-bits in a chromosome.
 */

#include <pgapack.h>

double NumberOfSetBits (PGAContext *, int, int, double *);

/* User defined hillclimber: We get a random bit index and if the bit is
 * 0 se set it to 1
 */
/*******************************************************************
*               user defined hillclimber                           *
*   ctx - contex variable                                          *
*   p   - chromosome index in population                           *
*   pop - which population to refer to                             *
*   We get a random bit index and set the bit to 1                 *
*   Most hillclimbers make use of problem knowledge: We know a 1   *
*   is better.                                                     *
*******************************************************************/
void hillclimb (PGAContext *ctx, int p, int pop)
{
    int len = PGAGetStringLength (ctx);
    int idx = PGARandomInterval (ctx, 0, len - 1);
    PGASetBinaryAllele (ctx, p, pop, idx, 1);
}

/*******************************************************************
*               user defined stopping check                        *
*   ctx - contex variable                                          *
*******************************************************************/
int stopcond (PGAContext *ctx)
{
    int best = PGAGetBestIndex (ctx, PGA_OLDPOP);
    if (PGAGetEvaluation (ctx, best, PGA_OLDPOP) == PGAGetStringLength (ctx)) {
        return PGA_TRUE;
    }
    return PGACheckStoppingConditions (ctx);
}

/*******************************************************************
*                   user main program                              *
*******************************************************************/
int main (int argc, char **argv)
{
    PGAContext *ctx;
    int full_report = 0;

    if (argc > 1 && atoi (argv [1])) {
        full_report = 1;
    }
    
    ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 256, PGA_MAXIMIZE);
    PGASetRandomSeed (ctx, 1);
    PGASetUserFunction (ctx, PGA_USERFUNCTION_HILLCLIMB, hillclimb);
    PGASetUserFunction (ctx, PGA_USERFUNCTION_STOPCOND, stopcond);
    PGASetRandomDeterministic (ctx, PGA_TRUE);
    /* Extended reporting */
    if (full_report) {
        PGASetPrintOptions (ctx, PGA_REPORT_ONLINE);
        PGASetPrintOptions (ctx, PGA_REPORT_OFFLINE);
        PGASetPrintOptions (ctx, PGA_REPORT_GENE_DISTANCE);
        PGASetPrintOptions (ctx, PGA_REPORT_STRING);
        PGASetPrintOptions (ctx, PGA_REPORT_WORST);
        PGASetPrintOptions (ctx, PGA_REPORT_AVERAGE);
    }
    
    PGASetUp (ctx);
    PGARun (ctx, NumberOfSetBits);
    PGADestroy (ctx);
    
    return 0;
}


/*******************************************************************
*               user defined evaluation function                   *
*   ctx - contex variable                                          *
*   p   - chromosome index in population                           *
*   pop - which population to refer to                             *
*******************************************************************/
double NumberOfSetBits (PGAContext *ctx, int p, int pop, double *dummy)
{
    int i, nbits, stringlen;

    stringlen = PGAGetStringLength (ctx);
    
    nbits = 0;
    for (i=0; i<stringlen; i++) {
	if (PGAGetBinaryAllele (ctx, p, pop, i)) {
	    nbits++;
        }
    }
    
    return (double) nbits;
}
