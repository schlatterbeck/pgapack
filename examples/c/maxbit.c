/*
 *  This is a test program for PGAPack.  The objective is to maximize the
 *  number of 1-bits in a chromosome.
 */

#include <pgapack.h>

double NumberOfSetBits (PGAContext *, int, int, double *);

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
