/*
COPYRIGHT

The following is a notice of limited availability of the code, and disclaimer
which must be included in the prologue of the code and in all source listings
of the code.

(C) COPYRIGHT 2008 University of Chicago

Permission is hereby granted to use, reproduce, prepare derivative works, and
to redistribute to others. This software was authored by:

D. Levine
Mathematics and Computer Science Division 
Argonne National Laboratory Group

with programming assistance of participants in Argonne National 
Laboratory's SERS program.

GOVERNMENT LICENSE

Portions of this material resulted from work developed under a
U.S. Government Contract and are subject to the following license: the
Government is granted for itself and others acting on its behalf a paid-up,
nonexclusive, irrevocable worldwide license in this computer software to
reproduce, prepare derivative works, and perform publicly and display
publicly.

DISCLAIMER

This computer code material was prepared, in part, as an account of work
sponsored by an agency of the United States Government. Neither the United
States, nor the University of Chicago, nor any of their employees, makes any
warranty express or implied, or assumes any legal liability or responsibility
for the accuracy, completeness, or usefulness of any information, apparatus,
product, or process disclosed, or represents that its use would not infringe
privately owned rights.
*/

/*****************************************************************************
*     FILE: utility.c: This file contains routines that perform utility
*                      functions.
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
*****************************************************************************/

#include <pgapack.h>

/*U****************************************************************************
  PGAMean - calculates the mean value of an array of elements

  Category: Utility

  Inputs:
    ctx  - context variable
    a    - array to take the mean of
    n    - number of elements in array a

  Outputs:
    The mean of the n elements in array a

  Example:
    PGAContext *ctx;
    double a[100], mean;
    :
    mean = PGAMean(ctx, a, 100);

****************************************************************************U*/
double PGAMean ( PGAContext *ctx, double *a, int n)
{
    int i;
    double result;

    PGADebugEntered("PGAMean");

    result = 0.;
    for( i=n-1; i>=0; i-- )
        result += a[i];

    PGADebugExited("PGAMean");

    return (result / (double) n );
}


/*U****************************************************************************
  PGAStddev - calculates the standard deviation of an array of elements

  Category: Utility

  Inputs:
    ctx  - context variable
    a    - array to take the standard deviation of
    n    - number of elements in array a
    mean - the mean of the elements in array a

  Outputs:
    The standard deviation of the n elements in array a

  Example:
    PGAContext *ctx;
    double a[100], mean, sigma;
    :
    mean  = PGAMean(ctx, a, 100);
    sigma = PGAStddev(ctx, a, 100, mean);

****************************************************************************U*/
double PGAStddev ( PGAContext *ctx, double *a, int n, double mean)
{
    int i;
    double result;

    PGADebugEntered("PGAStddev");

    result = 0;
    for(i=n-1; i>=0; i--)
        result += (a[i] - mean) * (a[i] - mean);
    result  = sqrt(result/n);

    PGADebugExited("PGAStddev");
    return (result);
}

/*U****************************************************************************
  PGARound - Mathematically round a double to an integer, using 0.5 as the
  cutoff value.

  Category: Utility

  Inputs:
     ctx - context variable
     x - the number to be rounded

  Outputs:
     The rounded number.

  Example:
     PGAContext *ctx;
     int y;
     y = PGARound(ctx, -78.6);

****************************************************************************U*/
int PGARound(PGAContext *ctx, double x)
{
    double ipart, frac;

    PGADebugEntered("PGARound");

    frac = modf(x, &ipart);
    if (frac <= -0.5)
      ipart--;
    else if (frac >= 0.5)
      ipart++;

    PGADebugExited("PGARound");

    return ((int)ipart);
}

/*U****************************************************************************
  PGACopyIndividual - copies string p1 in population pop1 to position p2 in
  population pop2

  Category: Generation

  Inputs:
     ctx  - context variable
     p1   - string to copy
     pop1 - symbolic constant of population containing string p1
     p2   - string to copy p1 to
     pop2 - symbolic constant of population containing string p2

  Outputs:
     String p2 is an exact copy of string p1.

  Example:
    PGAContext *ctx;
    int i,j;
    :
    PGACopyIndividual(ctx, i, PGA_OLDPOP, j, PGA_NEWPOP);

****************************************************************************U*/
void PGACopyIndividual( PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAIndividual *source, *dest;

    PGADebugEntered ("PGACopyIndividual");

    source = PGAGetIndividual (ctx, p1, pop1);
    dest   = PGAGetIndividual (ctx, p2, pop2);

    dest->evalfunc         = source->evalfunc;
    dest->fitness          = source->fitness;
    dest->evaluptodate     = source->evaluptodate;
    dest->auxtotal         = source->auxtotal;
    dest->auxtotaluptodate = source->auxtotaluptodate;
    if (ctx->ga.NumAuxEval) {
        memcpy (dest->auxeval, source->auxeval, ctx->ga.NumAuxEval);
    } else {
        dest->auxeval = NULL;
    }

    (*ctx->cops.CopyString)(ctx, p1, pop1, p2, pop2);

    PGADebugExited("PGACopyIndividual");
}

/*U****************************************************************************
  PGACheckSum - maps a string to a number to be used a verification check
  PGA_DATATYPE_USER is not supported.

  Category: Utility

  Inputs:
     ctx    - context variable
     p      - string index
     pop    - symbolic constant for the population

  Outputs:
     An integer representing the "value" of the string.

  Example:
     PGAContext *ctx;
     int p, sum;
     :
     sum = PGACheckSum(ctx, p, PGA_NEWPOP);

****************************************************************************U*/
int PGACheckSum(PGAContext *ctx, int p, int pop)
{
    long stringlen, totalchars, charbits, i, j, checksum, totalbytes, out_bit;
    unsigned char *message, specimen;

    PGADebugEntered("PGACheckSum");

    stringlen = PGAGetStringLength(ctx);
    switch (ctx->ga.datatype) {
      case PGA_DATATYPE_BINARY:
        totalbytes = ctx->ga.tw * sizeof(PGABinary);
        break;
      case PGA_DATATYPE_INTEGER:
        totalbytes = stringlen * sizeof(PGAInteger);
        break;
      case PGA_DATATYPE_REAL:
        totalbytes = stringlen * sizeof(PGAReal);
        break;
      case PGA_DATATYPE_CHARACTER:
        totalbytes = stringlen * sizeof(PGACharacter);
        break;
      default:
        totalbytes = 0;
        PGAError(ctx, "PGACheckSum: User datatype checksum may be invalid.",
                 PGA_WARNING, PGA_VOID, NULL);
        break;
      }

    message = (unsigned char *)PGAGetIndividual(ctx, p, pop)->chrom;
    totalchars = totalbytes / sizeof(unsigned char);
    charbits = sizeof(unsigned char) * 8;
    checksum = 0;
    for (i = 0; i < totalchars; i++) {
      specimen = *(message + i);
      for (j = 0; j < charbits; j++) {
        out_bit = (checksum & 0x80000000);
        checksum = (checksum << 1) + ((specimen & 0x80) == 0x80);
        if (out_bit)
          checksum ^= 0x04c11db7;
        specimen <<= 1;
      }
    }

    PGADebugExited("PGACheckSum");

    return (checksum);
}

/*U***************************************************************************
  PGAGetWorstIndex - returns the index of the string with the worst
  evaluation function value in population pop

  Category: Utility

  Inputs:
     ctx - context variable
     pop - symbolic constant of the population to find the worst string in

  Outputs:
     Index of the string with the worst evaluation function value

  Example:
     PGAContext *ctx;
     int worst;
     :
     worst = PGAGetWorstIndex(ctx,PGA_OLDPOP);

***************************************************************************U*/
int PGAGetWorstIndex(PGAContext *ctx, int pop)
{
    int     p, worst_indx = 0;
    double  eval, worst_eval;

    PGADebugEntered("PGAGetWorstIndex");

    for (p = 0; p < ctx->ga.PopSize; p++)
	if (!PGAGetEvaluationUpToDateFlag(ctx, p, pop))
	    PGAError(ctx, "PGAGetWorstIndex: Evaluate function not up to "
		     "date:", PGA_FATAL, PGA_INT, (void *) &p);

    worst_eval = PGAGetEvaluation(ctx, 0, pop);
    switch (PGAGetOptDirFlag(ctx)) {
    case PGA_MAXIMIZE:
	for (p = 1; p < ctx->ga.PopSize; p++) {
	    eval = PGAGetEvaluation(ctx, p, pop);
	    if (eval < worst_eval) {
		worst_indx = p;
		worst_eval = eval;
	    }
	}
	break;

    case PGA_MINIMIZE:
	for (p = 1; p < ctx->ga.PopSize; p++) {
	    eval = PGAGetEvaluation(ctx, p, pop);
	    if (eval > worst_eval) {
		worst_indx = p;
		worst_eval = eval;
	    }
	}
	break;    
    }

    PGADebugExited("PGAGetWorstIndex");

    return (worst_indx);
}

/*U***************************************************************************
  PGAGetBestIndex - returns the index of the string with the best evaluation
  function value in population pop

  Category: Utility

  Inputs:
     ctx - context variable
     pop - symbolic constant of the population to find the best string in

  Outputs:
     Index of the string with the best evaluation function value

  Example:
     PGAContext *ctx;
     int best;
     :
     best = PGAGetBestIndex(ctx,PGA_OLDPOP);

***************************************************************************U*/
int PGAGetBestIndex(PGAContext *ctx, int pop)
{
    int     p, Best_indx = 0;

    PGADebugEntered("PGAGetBestIndex");
     
    for (p = 0; p < ctx->ga.PopSize; p++)
        if (!PGAGetEvaluationUpToDateFlag(ctx, p, pop))
	    PGAError(ctx, "PGAGetBestIndex: Evaluate function not up to "
                     "date:", PGA_FATAL, PGA_INT, (void *) &p);

    for (p=1; p<ctx->ga.PopSize; p++) {
        if (PGAEvalCompare (ctx, p, pop, Best_indx, pop) > 0) {
            Best_indx = p;
        }
    }
     
    PGADebugExited("PGAGetBestIndex");

    return (Best_indx);
}


/*I****************************************************************************
  PGAGetIndividual - translate string index and population symbolic constant
  into pointer to an individual

  Inputs:
     ctx - context variable
     p   - string index
     pop - symbolic constant of the population the string is in

  Outputs:
     Address of the PGAIndividual structure for string p in population pop

  Example:
    PGAIndividual *source;
    PGAContext *ctx;
    int p;
    :
    source = PGAGetIndividual ( ctx, p, PGA_NEWPOP );

****************************************************************************I*/
PGAIndividual *PGAGetIndividual ( PGAContext *ctx, int p, int pop)
{
    PGAIndividual *ind;

    PGADebugEntered("PGAGetIndividual");

#ifdef OPTIMIZE
    ind = (pop == PGA_OLDPOP) ? ctx->ga.oldpop : ctx->ga.newpop;

    if (p>=0)
      ind += p;
    else
      ind += (p == PGA_TEMP1) ? ctx->ga.PopSize : ctx->ga.PopSize+1;
#else
    if (pop == PGA_OLDPOP)
      ind = ctx->ga.oldpop;
    else
      if (pop == PGA_NEWPOP)
	ind = ctx->ga.newpop;
      else
        PGAError(ctx, "PGAGetIndividual: Invalid value of pop:",
		 PGA_FATAL, PGA_INT, (void *) &pop );

    if (p>0 && p<ctx->ga.PopSize)
      ind += p;
    else
      if (p == PGA_TEMP1)
	ind += ctx->ga.PopSize;
      else
	if (p == PGA_TEMP2)
	  ind += ctx->ga.PopSize + 1;
	else
	  PGAError(ctx, "PGAGetIndividual: Invalid value of p:",
		   PGA_FATAL, PGA_INT, (void *) &p );
#endif

    PGADebugExited("PGAGetIndividual");

    return(ind);
}


/*I****************************************************************************
   PGAUpdateAverage - Updates the average fitness statistic for reporting.

   Inputs:
       ctx - context variable
       pop - symbolic constant of the population

   Outputs:

   Example:

**************************************************************************I*/
void PGAUpdateAverage(PGAContext *ctx, int pop)
{
    double ThisGensTotal = 0;
    int p;

    PGADebugEntered("PGAUpdateAverage");
    
    for (p = 0; p < ctx->ga.PopSize; p++)
	if (!PGAGetEvaluationUpToDateFlag(ctx, p, pop))
	    PGAError(ctx, "PGAUpdateOnline: Evaluate function not up to "
		     "date:", PGA_FATAL, PGA_INT, (void *) &p);

    for (p = 0; p < ctx->ga.PopSize; p++)
	ThisGensTotal += PGAGetEvaluation(ctx, p, pop);
    
    ctx->rep.Average = ThisGensTotal / (double)ctx->ga.PopSize;
    
    PGADebugExited("PGAUpdateAverage");
}


/*I****************************************************************************
  PGAUpdateOnline - Updates the online value based on the results in
  the new generation

  Inputs:
     ctx - context variable
     pop - symbolic constant of the population whose statistics to use

  Outputs:
     Updates an internal field in the context variable

  Example:
     PGAContext *ctx;
     :
     PGAUpdateOnline(ctx,PGA_NEWPOP);

**************************************************************************I*/
void PGAUpdateOnline(PGAContext *ctx, int pop)
{
     double ThisGensTotal = 0;
     int p;

    PGADebugEntered("PGAUpdateOnline");

     for (p = 0; p < ctx->ga.PopSize; p++)
          if (!PGAGetEvaluationUpToDateFlag(ctx, p, pop))
               PGAError(ctx, "PGAUpdateOnline: Evaluate function not up to "
                        "date:", PGA_FATAL, PGA_INT, (void *) &p);

     for (p = 0; p < ctx->ga.PopSize; p++)
          ThisGensTotal += PGAGetEvaluation(ctx, p, pop);

     PGADebugPrint(ctx, PGA_DEBUG_PRINTVAR, "PGAUpdateOnline",
                   "ThisGensTotal = ", PGA_DOUBLE, (void *) &ThisGensTotal);

     ctx->rep.Online = (ctx->rep.Online * ctx->ga.PopSize * (ctx->ga.iter - 1)
                        + ThisGensTotal) / ctx->ga.iter / ctx->ga.PopSize;

    PGADebugExited("PGAUpdateOnline");

}

/*I****************************************************************************
  PGAUpdateOffline - Updates the offline value based on the results in
  the new generation

  Inputs:
     ctx - context variable
     pop - symbolic constant of the population whose statistics to use

  Outputs:
     Updates an internal field in the context variable

  Example:
     PGAContext *ctx;
     :
     PGAUpdateOffline(ctx,PGA_NEWPOP);

**************************************************************************I*/
void PGAUpdateOffline(PGAContext *ctx, int pop)
{
     int p;

    PGADebugEntered("PGAUpdateOffline");

     for (p = 0; p < ctx->ga.PopSize; p++)
          if (!PGAGetEvaluationUpToDateFlag(ctx, p, pop))
               PGAError(ctx, "PGAUpdateOffline: Evaluate function not up to "
                        "date:", PGA_FATAL, PGA_INT, (void *) &p);

     p = PGAGetBestIndex(ctx, pop);

     ctx->rep.Offline = ((ctx->ga.iter - 1) * ctx->rep.Offline +
                         PGAGetEvaluation(ctx, p, pop)) / ctx->ga.iter;

    PGADebugExited("PGAUpdateOffline");
}

/*I****************************************************************************
  PGAComputeSimilarity - computes the percentage of the population that have
  the same evaluation function

  Inputs:
     ctx - context variable
     pop - symbolic constant of the population whose statistics to use

  Outputs:
     returns a count of the number of  population members that have the same
     evaluation function value

  Example:
     PGAContext *ctx;
     :
     PGAComputeSimilarity(ctx,PGA_NEWPOP);

**************************************************************************I*/
int PGAComputeSimilarity(PGAContext *ctx, PGAIndividual *pop)
{
     int max = 0, curr = 1, i;
     double prev;

    PGADebugEntered("PGAComputeSimilarity");

     for(i=0; i < ctx->ga.PopSize; i++)
     {
          ctx->scratch.dblscratch[i] = (pop + i)->evalfunc;
          ctx->scratch.intscratch[i] = i;
     }

     PGADblHeapSort(ctx, ctx->scratch.dblscratch, ctx->scratch.intscratch,
                    ctx->ga.PopSize);

     prev = ctx->scratch.dblscratch[0];

     for(i = 1; i < ctx->ga.PopSize; i++)
     {
          if (ctx->scratch.dblscratch[i] == prev)
               curr++;
          else
          {
               if (curr > max)
                    max = curr;
               curr = 1;
          }
     }

    PGADebugExited("PGAComputeSimilarity");

     return(100 * max / ctx->ga.PopSize);
}

/* Used in PGAEvalCompare below */
static inline double CMP (double a, double b)
{
    return (a < b ? -1 : (a > b ? 1 : 0));
}

/*U****************************************************************************
  PGAEvalCompare - Compare two strings for selection.
  This typically simply compares fitness.
  But if auxiliary evaluations are defined, the auxiliary evaluations are
  treated as constraints. This function is a user function and can be
  redefined for other purposes: By redefining this function, e.g.,
  instead of using the aux evaluations for constraint-handling, instead
  (or in addition), multi objective evaluation can be implemented.
  The default handling of auxiliary evaluations is incompatible with
  certain selection schemes, see checks in create.c

  Note that PGAEvalCompare is now used in several contexts, including
  finding the best evaluation. For very badly scaled problems, the
  default fitness computation will degenerate if there are very large
  evaluation values and very small ones. In that case the fitness will
  not reflect the evaluation. Therefore PGAEvalCompare will now always
  sort on evaluation values ignoring the fitness. This improves
  Tournament selection for very badly scaled problems.

  Category: Operators

  Inputs:
    ctx   - context variable
    p1    - first string to compare
    pop1  - symbolic constant of population of first string
    p2    - second string to compare
    pop2  - symbolic constant of population of second string

  Outputs:
    <0 if p2 is "better" than p1
    >0 if p1 is "better" than p2
    0  if both compare equal
    Thinks of this as sorting individuals by decreasing fitness or
    increasing constraint violations.

  Example:
    PGAContext *ctx;
    int result;
    :
    result = PGAEvalCompare(ctx, p1, PGA_OLDPOP, p2, PGA_OLDPOP);

****************************************************************************U*/
int PGAEvalCompare (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    double auxt1 = 0, auxt2 = 0;
    int dir = PGAGetOptDirFlag (ctx);
    PGAIndividual *ind1, *ind2;
    if (!PGAGetEvaluationUpToDateFlag (ctx, p1, pop1)) {
        PGAError
            ( ctx
            , "PGAEvalCompare: first individual not up to date:"
            , PGA_FATAL, PGA_INT, (void *) &p1
            );
    }
    if (!PGAGetEvaluationUpToDateFlag (ctx, p2, pop2)) {
        PGAError
            ( ctx
            , "PGAEvalCompare: second individual not up to date:"
            , PGA_FATAL, PGA_INT, (void *) &p2
            );
    }
    if (ctx->ga.NumAuxEval > 0) {
        auxt1 = PGAGetAuxTotal (ctx, p1, pop1);
        auxt2 = PGAGetAuxTotal (ctx, p2, pop2);
    }
    if (auxt1 || auxt2) {
        return CMP (auxt2, auxt1);
    }
    /* We might use the fitness if both populations are the same
       otherwise fitness values are not comparable. But we now
       use the evaluation in any case.
     */
    ind1 = PGAGetIndividual (ctx, p1, pop1);
    ind2 = PGAGetIndividual (ctx, p2, pop2);
    switch (dir) {
    case PGA_MAXIMIZE:
        return CMP (ind1->evalfunc, ind2->evalfunc);
        break;
    case PGA_MINIMIZE:
        return CMP (ind2->evalfunc, ind1->evalfunc);
        break;
    default:
        PGAError
            (ctx
            , "PGAEvalCompare: Invalid value of PGAGetOptDirFlag:"
            , PGA_FATAL, PGA_INT, (void *) &dir
            );
        break;
    }
    /* notreached */
    return 0;
}

