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

/*!***************************************************************************
* \file
* This file contains routines that perform utility functions.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/
/*!***************************************************************************
 *  \defgroup utility Utility functions
 *  \brief Utility functions for internal and explicit use
 *****************************************************************************/

#include <pgapack.h>

/*!****************************************************************************
    \brief Shuffle a list.
    \ingroup utility
    \param  ctx   Context pointer
    \param  list  array of integers to shuffle
    \param  n     number of elements in array
    \return Shuffled array

    \rst

    Description
    -----------

    We're using Durstenfeld's version of the Fisher-Yates shuffle.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int list [ctx->ga.PopSize];
      ...
      for (i=0; i<ctx->ga.PopSize; i++) {
          list [i] = i;
      }
      PGAShuffle (ctx, list, ctx->ga.PopSize);

    \endrst

******************************************************************************/
void PGAShuffle (PGAContext *ctx, int *list, int n)
{
    int i, j, tmp = 0;

    for (i=0; i<n-1; i++) {
        j = PGARandomInterval (ctx, i, n - 1);
        tmp = list [j];
        list [j] = list [i];
        list [i] = tmp;
    }
}

/*!****************************************************************************
    \brief Shuffle a list of PGAInteger.
    \ingroup utility
    \param  ctx   Context pointer
    \param  list  array of PGAInteger elements to shuffle
    \param  n     number of elements in array
    \return Shuffled array

    \rst

    Description
    -----------

    We're using Durstenfeld's version of the Fisher-Yates shuffle.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      PGAInteger list [ctx->ga.StringLen];
      ...
      for (i=0; i<ctx->ga.StringLen; i++) {
          list [i] = i;
      }
      PGAShuffle (ctx, list, ctx->ga.StringLen);

    \endrst

******************************************************************************/
void PGAShufflePGAInteger (PGAContext *ctx, PGAInteger *list, int n)
{
    int i, j;
    PGAInteger tmp = 0;

    for (i=0; i<n-1; i++) {
        j = PGARandomInterval (ctx, i, n - 1);
        tmp = list [j];
        list [j] = list [i];
        list [i] = tmp;
    }
}

/*!****************************************************************************
    \brief Calculates the mean value of an array of elements
    \ingroup utility
    \param  ctx   context variable
    \param  a     array to take the mean of
    \param  n     number of elements in array a
    \return The mean of the n elements in array a

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      double a [100], mean;

      ...
      mean = PGAMean (ctx, a, 100);

    \endrst

******************************************************************************/
double PGAMean (PGAContext *ctx, double *a, int n)
{
    int i;
    double result;

    PGADebugEntered ("PGAMean");

    result = 0.;
    for (i=n-1; i>=0; i--) {
        result += a [i];
    }

    PGADebugExited ("PGAMean");

    return result / n;
}


/*!****************************************************************************
    \brief Calculate the standard deviation of an array of elements
    \ingroup utility
    \param  ctx   context variable
    \param  a     array to take the standard deviation of
    \param  n     number of elements in array a
    \param  mean  the mean of the elements in array a
    \return The standard deviation of the n elements in array a

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      double a [100], mean, sigma;

      ...
      mean  = PGAMean (ctx, a, 100);
      sigma = PGAStddev (ctx, a, 100, mean);

    \endrst

******************************************************************************/
double PGAStddev (PGAContext *ctx, double *a, int n, double mean)
{
    int i;
    double result;

    PGADebugEntered ("PGAStddev");

    result = 0;
    for (i=n-1; i>=0; i--) {
        result += (a [i] - mean) * (a [i] - mean);
    }
    result = sqrt (result/n);

    PGADebugExited ("PGAStddev");
    return result;
}

/*!****************************************************************************
    \brief Mathematically round a double to an integer, using 0.5 as the
           cutoff value.
    \ingroup utility
    \param   ctx  context variable
    \param   x    the number to be rounded
    \return The rounded number

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int y;

       ...
       y = PGARound (ctx, -78.6);

    \endrst

******************************************************************************/
int PGARound (PGAContext *ctx, double x)
{
    double ipart, frac;

    PGADebugEntered ("PGARound");

    frac = modf (x, &ipart);
    if (frac <= -0.5) {
        ipart--;
    } else if (frac >= 0.5) {
        ipart++;
    }

    PGADebugExited ("PGARound");

    return (int)ipart;
}

/*!****************************************************************************
    \brief Copy individual source to individual dest.
    \ingroup explicit
    \param  src  Individual to copy
    \param  dst  Individual to copy ind1 to
    \return Individual dest is a copy of Individual source

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      PGAIndividual *source, *dest;

      ...
      INDCopyIndividual (ctx, source, dest);

    \endrst

******************************************************************************/
void INDCopyIndividual (PGAIndividual *src, PGAIndividual *dst)
{
    PGAContext *ctx = src->ctx;
    int srcidx, dstidx;
    int spop, dpop;

    PGADebugEntered ("INDCopyIndividual");

    assert (src->ctx != NULL);
    assert (dst->ctx == NULL || dst->ctx == src->ctx);
    dst->ctx              = src->ctx;
    dst->evalue           = src->evalue;
    dst->fitness          = src->fitness;
    dst->evaluptodate     = src->evaluptodate;
    dst->auxtotal         = src->auxtotal;
    dst->auxtotalok       = src->auxtotalok;
    dst->rank             = src->rank;
    if (ctx->ga.NumAuxEval) {
        memcpy
            (dst->auxeval, src->auxeval, ctx->ga.NumAuxEval * sizeof (double));
    } else {
        dst->auxeval = NULL;
    }
    srcidx = src - src->pop;
    dstidx = dst - dst->pop;
    spop   = src->pop == ctx->ga.oldpop ? PGA_OLDPOP : PGA_NEWPOP;
    dpop   = dst->pop == ctx->ga.oldpop ? PGA_OLDPOP : PGA_NEWPOP;

    (*ctx->cops.CopyString)(ctx, srcidx, spop, dstidx, dpop);

    PGADebugExited ("INDCopyIndividual");
}

/*!****************************************************************************
    \brief Copy string p1 in population pop1 to position p2 in population pop2.
    \ingroup explicit
    \param  ctx   context variable
    \param  p1    string to copy
    \param  pop1  symbolic constant of population containing string p1
    \param  p2    string to copy p1 to
    \param  pop2  symbolic constant of population containing string p2
    \return String p2 is an exact copy of string p1

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int i, j;

      ...
      PGACopyIndividual (ctx, i, PGA_OLDPOP, j, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGACopyIndividual (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAIndividual *src, *dst;

    PGADebugEntered ("PGACopyIndividual");

    src = PGAGetIndividual (ctx, p1, pop1);
    dst = PGAGetIndividual (ctx, p2, pop2);

    INDCopyIndividual (src, dst);

    PGADebugExited ("PGACopyIndividual");
}

/*!****************************************************************************
    \brief Map a string to a number to be used as verification check.
    \ingroup utility
    \param  ctx     context variable
    \param  p       string index
    \param  pop     symbolic constant for the population
    \return An integer representing the "value" of the string

    \rst

    Description
    -----------

    The data type :c:macro:`PGA_DATATYPE_USER` is not supported.
    See also the USER function PGAHash :c:macro:`PGA_USERFUNCTION_HASH`
    and the hash implementations for individual data types.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p, sum;

       ...
       sum = PGACheckSum (ctx, p, PGA_NEWPOP);

    \endrst

******************************************************************************/
int PGACheckSum (PGAContext *ctx, int p, int pop)
{
    long stringlen, totalchars, charbits, i, j, checksum, totalbytes, out_bit;
    unsigned char *message, specimen;

    PGADebugEntered ("PGACheckSum");

    stringlen = PGAGetStringLength (ctx);
    switch (ctx->ga.datatype) {
      case PGA_DATATYPE_BINARY:
        totalbytes = ctx->ga.tw * sizeof (PGABinary);
        break;
      case PGA_DATATYPE_INTEGER:
        totalbytes = stringlen * sizeof (PGAInteger);
        break;
      case PGA_DATATYPE_REAL:
        totalbytes = stringlen * sizeof (PGAReal);
        break;
      case PGA_DATATYPE_CHARACTER:
        totalbytes = stringlen * sizeof (PGACharacter);
        break;
      default:
        totalbytes = 0;
        PGAError
            ( ctx, "PGACheckSum: User datatype checksum may be invalid."
            , PGA_WARNING, PGA_VOID, NULL
            );
        break;
    }

    message = (unsigned char *)PGAGetIndividual (ctx, p, pop)->chrom;
    totalchars = totalbytes / sizeof (unsigned char);
    charbits = sizeof (unsigned char) * 8;
    checksum = 0;
    for (i = 0; i < totalchars; i++) {
        specimen = *(message + i);
        for (j = 0; j < charbits; j++) {
            out_bit = (checksum & 0x80000000);
            checksum = (checksum << 1) + ((specimen & 0x80) == 0x80);
            if (out_bit) {
                checksum ^= 0x04c11db7;
            }
            specimen <<= 1;
        }
    }

    PGADebugExited ("PGACheckSum");

    return checksum;
}

/*!***************************************************************************
    \brief Return the index of the string with the worst evaluation
           function value in population pop.
    \ingroup query
    \param  ctx  context variable
    \param  pop  symbolic constant of the population to find the worst string in
    \return Index of the string with the worst evaluation function value

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int worst;

       ...
       worst = PGAGetWorstIndex (ctx, PGA_OLDPOP);

    \endrst

*****************************************************************************/
int PGAGetWorstIndex (PGAContext *ctx, int pop)
{
    int     p, worst_indx = 0;

    PGADebugEntered ("PGAGetWorstIndex");

    for (p=0; p<ctx->ga.PopSize; p++) {
        if (!PGAGetEvaluationUpToDateFlag (ctx, p, pop)) {
            PGAError
                ( ctx, "PGAGetWorstIndex: Evaluate function not up to date:"
                , PGA_FATAL, PGA_INT, (void *) &p
                );
        }
    }

    for (p=1; p<ctx->ga.PopSize; p++) {
        if (PGAEvalCompare (ctx, p, pop, worst_indx, pop) > 0) {
            worst_indx = p;
        }
    }

    PGADebugExited ("PGAGetWorstIndex");

    return (worst_indx);
}

/*!***************************************************************************
    \brief Return the index of the string with the best evaluation
           function value in population pop
    \ingroup query
    \param  ctx    context variable
    \param  popidx symbolic constant of the population to find the best
            string in
    \return Index of the string with the best evaluation function value

    \rst

    Description
    -----------

    Note that in the presence of multiple evaluations, this will return

    - If all strings violate constraints the one with the least constraint
      violation
    - In case there exist strings without constraint violations one
      *randomly* chosen from rank 0.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int best;

       ...
       best = PGAGetBestIndex (ctx, PGA_OLDPOP);

    \endrst

*****************************************************************************/
int PGAGetBestIndex (PGAContext *ctx, int popidx)
{
    int p, Best_indx = 0;
    DECLARE_DYNARRAY (int, bestidxs, ctx->ga.PopSize);
    int nbest = 0;
    int is_multi = ctx->ga.NumAuxEval > ctx->ga.NumConstraint;
    PGAIndividual *pop = PGAGetIndividual (ctx, 0, popidx);

    PGADebugEntered ("PGAGetBestIndex");

    for (p=0; p<ctx->ga.PopSize; p++) {
        if (!PGAGetEvaluationUpToDateFlag (ctx, p, popidx)) {
            PGAError
                ( ctx, "PGAGetBestIndex: Evaluate function not up to date:"
                , PGA_FATAL, PGA_INT, (void *) &p
                );
        }
    }
    if (is_multi && pop->rank == 0 && !INDGetAuxTotal (pop)) {
        bestidxs [nbest++] = 0;
    }

    for (p=1; p<ctx->ga.PopSize; p++) {
        if (PGAEvalCompare (ctx, p, popidx, Best_indx, popidx) < 0) {
            Best_indx = p;
        }
        if (is_multi && (pop + p)->rank == 0 && !INDGetAuxTotal (pop + p)) {
            bestidxs [nbest++] = p;
        }
    }
    /* Shuffle and return first */
    if (is_multi && nbest > 1) {
        PGAShuffle (ctx, bestidxs, nbest);
        Best_indx = bestidxs [0];
        assert ((pop + Best_indx)->rank == 0);
    }

    PGADebugExited ("PGAGetBestIndex");

    return Best_indx;
}


/*!***************************************************************************
    \brief Return the best evaluation value in population pop
           for the given evaluation index.
    \ingroup query
    \param  ctx  context variable
    \param  pop  symbolic constant of the population to find the best
                 eval, only PGA_OLDPOP is allowed
    \param  idx  Index of the evaluation function, 0 for the one return
                 by your evaluation function 1 .. NumAuxEval for the
                 auxiliary evaluations
    \return Index of the string with the best evaluation function value

    \rst

    Description
    -----------

    The evaluation function index must be 0 < ``idx`` <= ``NumAuxEval``.
    So for single evaluation only the index 0 is allowed and the return
    is the evaluation as from :c:func:`PGAGetBestIndex`. Note that this
    accesses the pre-computed statistics and therefore only
    :c:macro:`PGA_OLDPOP` is allowed for the population constant.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double best;

       ...
       best = PGAGetBestReport (ctx, PGA_OLDPOP, 1);

    \endrst

*****************************************************************************/
double PGAGetBestReport (PGAContext *ctx, int pop, int idx)
{
    if (pop != PGA_OLDPOP) {
        PGAErrorPrintf (ctx, PGA_FATAL, "Index must be PGA_OLDPOP");
    }
    if (idx < 0 || idx > ctx->ga.NumAuxEval) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "Index must be 0 <= idx <= NumAuxEval, got %d"
            , idx
            );
    }

    return ctx->rep.Best [idx];
}

/*!***************************************************************************
    \brief Return the index of the string with the best evaluation
           function value in population pop for the given evaluation index.
    \ingroup query
    \param  ctx  context variable
    \param  pop  symbolic constant of the population to find the best
                 string in
    \param  idx  Index of the evaluation function, 0 for the one return
                 by your evaluation function 1 .. NumAuxEval for the
                 auxiliary evaluations
    \return Index of the string with the best evaluation function value

    \rst

    Description
    -----------

    The evaluation function index must be 0 < ``idx`` <= ``NumAuxEval``.
    So for single evaluation only the index 0 is allowed and the return
    is the same as from :c:func:`PGAGetBestIndex`. Note that this
    accesses the pre-computed statistics and therefore only
    :c:macro:`PGA_OLDPOP` is allowed for the population constant.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int best;

       ...
       best = PGAGetBestReportIndex (ctx, PGA_OLDPOP, 1);

    \endrst

*****************************************************************************/
int PGAGetBestReportIndex (PGAContext *ctx, int pop, int idx)
{
    if (pop != PGA_OLDPOP) {
        PGAErrorPrintf (ctx, PGA_FATAL, "Index must be PGA_OLDPOP");
    }
    if (idx < 0 || idx > ctx->ga.NumAuxEval) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "Index must be 0 <= idx <= NumAuxEval, got %d"
            , idx
            );
    }

    return ctx->rep.BestIdx [idx];
}


/*!****************************************************************************
    \brief Return individual defined by population symbolic constant and
           population index.
    \ingroup allele
    \param  ctx  context variable
    \param  p    string index
    \param  pop  symbolic constant of the population the string is in
    \return Address of the PGAIndividual structure for string p in
            population pop

    \rst

    Example
    -------

    .. code-block:: c

      PGAIndividual *source;
      PGAContext *ctx;
      int p;

      ...
      source = PGAGetIndividual (ctx, p, PGA_NEWPOP);

    \endrst

******************************************************************************/
PGAIndividual *PGAGetIndividual (PGAContext *ctx, int p, int pop)
{
    PGAIndividual *ind;

    PGADebugEntered ("PGAGetIndividual");

#ifdef OPTIMIZE
    ind = (pop == PGA_OLDPOP) ? ctx->ga.oldpop : ctx->ga.newpop;

    if (p >= 0) {
        ind += p;
    } else {
        ind += (p == PGA_TEMP1) ? ctx->ga.PopSize : ctx->ga.PopSize+1;
    }
#else
    if (pop == PGA_OLDPOP) {
        ind = ctx->ga.oldpop;
    } else {
        if (pop == PGA_NEWPOP) {
            ind = ctx->ga.newpop;
        } else {
            PGAError
              ( ctx, "PGAGetIndividual: Invalid value of pop:"
              , PGA_FATAL, PGA_INT, (void *) &pop
              );
        }
    }

    if (p>0 && p<ctx->ga.PopSize) {
        ind += p;
    } else {
        if (p == PGA_TEMP1) {
            ind += ctx->ga.PopSize;
        } else {
            if (p == PGA_TEMP2) {
                ind += ctx->ga.PopSize + 1;
            } else {
                PGAError
                    ( ctx, "PGAGetIndividual: Invalid value of p:"
                    , PGA_FATAL, PGA_INT, (void *) &p
                    );
            }
        }
    }
#endif

    PGADebugExited ("PGAGetIndividual");

    return ind;
}


/*!****************************************************************************
    \brief Updates the best fitness statistic for reporting.
    \ingroup internal
    \param  ctx   context variable
    \param  popix symbolic constant of the population
    \return None

    \rst

    Description
    -----------

    Note that in the presence of constraints if no individual without
    constraints is found, the best value (for all functions except
    constraint functions) is ``NAN``.

    \endrst

****************************************************************************/
void PGAUpdateBest (PGAContext *ctx, int popix)
{
    PGAIndividual *ind = popix == PGA_OLDPOP ? ctx->ga.oldpop : ctx->ga.newpop;
    PGAIndividual *pop = ind;
    DECLARE_DYNARRAY (PGAIndividual *, best, ctx->ga.NumAuxEval + 1);
    int numfun = ctx->ga.NumAuxEval - ctx->ga.NumConstraint;
    int validcount = 0;
    int p, k;

    PGADebugEntered ("PGAUpdateBest");

    for (k=0; k<=ctx->ga.NumAuxEval; k++) {
        best [k] = ind;
    }
    if (!INDGetAuxTotal (ind)) {
        validcount++;
    }
    ctx->rep.MinSumConstr = ind->auxtotal;
    ind++;
    for (p=1; p<ctx->ga.PopSize; p++) {
        if (!PGAGetEvaluationUpToDateFlag (ctx, p, popix)) {
            PGAError
                ( ctx, "PGAUpdateBest: Evaluate function not up to date:"
                , PGA_FATAL, PGA_INT, (void *) &p
                );
        }
        if (!INDGetAuxTotal (ind)) {
            validcount++;
        }
        if (ind->auxtotal < ctx->rep.MinSumConstr) {
            ctx->rep.MinSumConstr = ind->auxtotal;
        }
        for (k=0; k<=ctx->ga.NumAuxEval; k++) {
            PGAIndividual *b = best [k];
            double enew;
            double eold;
            if (k == 0) {
                enew = ind->evalue;
                eold = b->evalue;
            } else if (k <= numfun) {
                enew = ind->auxeval [k-1];
                eold = b->auxeval   [k-1];
            } else {
                enew = ind->auxeval [k-1] <= 0 ? 0 : ind->auxeval [k-1];
                eold = b->auxeval [k-1]   <= 0 ? 0 : b->auxeval [k-1];
            }
            if (k > numfun) {
                if (enew < eold) {
                    best [k] = ind;
                }
            } else {
                int funcmp = ctx->ga.optdir == PGA_MINIMIZE
                           ? CMP (enew, eold)
                           : CMP (eold, enew);
                if (  (  INDGetAuxTotal (b)   - ctx->ga.Epsilon  > 0
                      && INDGetAuxTotal (ind) - ctx->ga.Epsilon <= 0
                      )
                   || (  INDGetAuxTotal (ind) - ctx->ga.Epsilon <= 0
                      && funcmp < 0
                      )
                   )
                {
                    best [k] = ind;
                }
            }
        }
        ind++;
    }
    ctx->rep.validcount = validcount;
    for (k=0; k<=ctx->ga.NumAuxEval; k++) {
        if (k<=numfun) {
            double e = (k==0) ? best [k]->evalue : best [k]->auxeval [k-1];
            if (INDGetAuxTotal (best [k]) - ctx->ga.Epsilon > 0) {
                ctx->rep.Best [k] = NAN;
            } else {
                ctx->rep.Best [k] = e;
            }
        } else {
            ctx->rep.Best [k] = best [k]->auxeval [k-1];
        }
        ctx->rep.BestIdx [k] = best [k] - pop;
    }

    PGADebugExited ("PGAUpdateBest");
}

/*!****************************************************************************
    \brief Update the average fitness statistic for reporting.
    \ingroup internal
    \param  ctx  context variable
    \param  pop  symbolic constant of the population
    \return None

    \rst

    Description
    -----------

    Note that in the presence of constraints only the *unconstrained*
    function evaluations are averaged, except for the
    constraint-functions, these are averaged unconditionally.

    \endrst

****************************************************************************/
void PGAUpdateAverage (PGAContext *ctx, int pop)
{
    int validcount = 0;
    PGAIndividual *ind = pop == PGA_OLDPOP ? ctx->ga.oldpop : ctx->ga.newpop;
    int numfun = ctx->ga.NumAuxEval - ctx->ga.NumConstraint;
    int p, k;

    PGADebugEntered ("PGAUpdateAverage");

    /* Intentionally one more than NumAuxEval: evalue + auxeval */
    for (k=0; k<=ctx->ga.NumAuxEval; k++) {
        ctx->rep.Average [k] = 0;
    }
    for (p=0; p<ctx->ga.PopSize; p++) {
        if (!PGAGetEvaluationUpToDateFlag (ctx, p, pop)) {
            PGAError
                ( ctx, "PGAUpdateAverage: Evaluate function not up to date:"
                , PGA_FATAL, PGA_INT, (void *) &p
                );
        }
        if (INDGetAuxTotal (ind) == 0) {
            ctx->rep.Average [0] += ind->evalue;
            for (k=0; k<numfun; k++) {
                ctx->rep.Average [k+1] += ind->auxeval [k];
            }
            validcount++;
        }
        for (k=numfun; k<ctx->ga.NumAuxEval; k++) {
            ctx->rep.Average [k+1] += ind->auxeval [k];
        }
        ind++;
    }
    if (validcount) {
        for (k=0; k<=numfun; k++) {
            ctx->rep.Average [k] /= (double)validcount;
        }
    }
    for (k=numfun+1; k<=ctx->ga.NumAuxEval; k++) {
        ctx->rep.Average [k] /= (double)ctx->ga.PopSize;
    }
    ctx->rep.validcount = validcount;

    PGADebugExited ("PGAUpdateAverage");
}


/*!****************************************************************************
    \brief Update the online value based on the results in the new generation.
    \ingroup internal
    \param  ctx  context variable
    \param  pop  symbolic constant of the population whose statistics to use
    \return Updates an internal field in the context variable

    \rst

    Description
    -----------

    Note that in the presence of constraints only the *unconstrained*
    function evaluations are averaged, except for the
    constraint-functions, these are averaged unconditionally.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGAUpdateOnline (ctx, PGA_NEWPOP);

    \endrst

****************************************************************************/
void PGAUpdateOnline (PGAContext *ctx, int pop)
{
    int validcount = 0;
    PGAIndividual *ind = pop == PGA_OLDPOP ? ctx->ga.oldpop : ctx->ga.newpop;
    int numfun = ctx->ga.NumAuxEval - ctx->ga.NumConstraint;
    DECLARE_DYNARRAY (double, evalsum, ctx->ga.NumAuxEval + 1);
    int p, k;

    PGADebugEntered ("PGAUpdateOnline");

    for (k=0; k<=ctx->ga.NumAuxEval; k++) {
        if (k > numfun) {
            ctx->rep.Online [k] *=
                (double)ctx->ga.PopSize * (double)(ctx->ga.iter - 1);
        } else {
            ctx->rep.Online [k] *= (double)ctx->rep.validonline;
        }
        evalsum [k] = 0;
    }
    for (p=0; p<ctx->ga.PopSize; p++) {
        if (!PGAGetEvaluationUpToDateFlag (ctx, p, pop)) {
            PGAError
                ( ctx, "PGAUpdateOnline: Evaluate function not up to date:"
                , PGA_FATAL, PGA_INT, (void *) &p
                );
        }
        if (INDGetAuxTotal (ind) == 0) {
            evalsum [0] += ind->evalue;
            for (k=0; k<numfun; k++) {
                evalsum [k+1] += ind->auxeval [k];
            }
            validcount++;
        }
        for (k=numfun; k<ctx->ga.NumAuxEval; k++) {
            evalsum [k+1] += ind->auxeval [k];
        }
        ind++;
    }
    ctx->rep.validonline += validcount;
    for (k=0; k<=ctx->ga.NumAuxEval; k++) {
        ctx->rep.Online [k] += evalsum [k];
    }
    if (ctx->rep.validonline) {
        for (k=0; k<=numfun; k++) {
            ctx->rep.Online [k] /= (double)ctx->rep.validonline;
        }
    }
    for (k=numfun+1; k<=ctx->ga.NumAuxEval; k++) {
        ctx->rep.Online [k] /= (double)ctx->ga.PopSize * (double)ctx->ga.iter;
    }

    PGADebugExited ("PGAUpdateOnline");

}

/*!****************************************************************************
    \brief Update the offline value based on the results in the new generation.
    \ingroup internal
    \param  ctx  context variable
    \param  pop  symbolic constant of the population whose statistics to use
    \return Updates an internal field in the context variable

    \rst

    Description
    -----------

    Note that in the presence of constraints only the *unconstrained*
    best function evaluations are averaged, except for the
    constraint-functions, these are averaged unconditionally.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGAUpdateOffline (ctx, PGA_NEWPOP);

    \endrst

****************************************************************************/
void PGAUpdateOffline (PGAContext *ctx, int pop)
{
    int k;
    int numfun = ctx->ga.NumAuxEval - ctx->ga.NumConstraint;

    PGADebugEntered ("PGAUpdateOffline");

    PGAUpdateBest (ctx, pop);

    for (k=0; k<=ctx->ga.NumAuxEval; k++) {
        if (k<=numfun) {
            if (ctx->rep.validcount) {
                assert (!isnan (ctx->rep.Best [k]));
                ctx->rep.Offline [k] =
                    ( ctx->rep.Offline [k] * (double)ctx->rep.validoffline
                    + ctx->rep.Best [k]
                    ) / (double)(ctx->rep.validoffline + 1);
            }
        } else {
            assert (!isnan (ctx->rep.Best [k]));
            ctx->rep.Offline [k] =
                ( ctx->rep.Offline [k] * (double)(ctx->ga.iter - 1)
                + ctx->rep.Best [k]
                ) / (double)ctx->ga.iter;

        }
    }
    if (ctx->rep.validcount) {
        ctx->rep.validoffline += 1;
    }

    PGADebugExited ("PGAUpdateOffline");
}

/*!****************************************************************************
    \brief Compute the percentage of the population that have the same
           evaluation function.
    \ingroup internal
    \param  ctx      context variable
    \param  popindex symbolic constant of the population whose statistics to use
    \return returns a count of the number of  population members that
            have the same evaluation function value

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGAComputeSimilarity (ctx, PGA_NEWPOP);

    \endrst

****************************************************************************/
int PGAComputeSimilarity (PGAContext *ctx, int popindex)
{
    int max = 0, curr = 1, i;
    PGAIndividual *prev;

    PGADebugEntered ("PGAComputeSimilarity");

    /* No need to init the indeces, filled in by PGAEvalSort */
    PGAEvalSort (ctx, popindex, ctx->scratch.intscratch);

    prev = PGAGetIndividual (ctx, 0, popindex);

    for (i=1; i<ctx->ga.PopSize; i++)
    {
        PGAIndividual *ind = PGAGetIndividual (ctx, i, popindex);
        int same = 0;

        if (ind->evalue == prev->evalue) {
            int j;
            same = 1;
            for (j=0; j<ctx->ga.NumAuxEval; j++) {
                if (ind->auxeval [j] != prev->auxeval [j]) {
                    same = 0;
                    break;
                }
            }
        }

        if (same) {
            curr++;
        } else {
            if (curr > max)
                max = curr;
            curr = 1;
        }
        same = 0;
    }

    PGADebugExited ("PGAComputeSimilarity");

    return 100 * max / ctx->ga.PopSize;
}

/*!****************************************************************************
    \brief Compare two individuals by evaluation.
    \ingroup utility
    \param  ind1   Pointer to first individual
    \param  ind2   Pointer to second individual
    \return
    - >0 if p2 is "better" than p1
    - <0 if p1 is "better" than p2
    - 0  if both compare equal

    \rst

    Description
    -----------

    This typically simply compares evaluation taking into account the
    evaluation direction (minimize/maximize). We sort "better" individuals
    first. For more details see :c:func:`PGAEvalCompare`.

    Thinks of this as sorting individuals by decreasing fitness or
    increasing constraint violations. See also :c:func:`PGAEvalCompare`.

    Example
    -------

    .. code-block:: c

      PGAIndividual *ind1, *ind2;
      int result;
      ind1 = PGAGetIndividual (...
      ind2 = PGAGetIndividual (...

      result = INDEvalCompare (ind1, ind2);

    \endrst

******************************************************************************/
int INDEvalCompare (PGAIndividual *ind1, PGAIndividual *ind2)
{
    double auxt1 = 0, auxt2 = 0;
    PGAContext *ctx = ind1->ctx;
    int dir = PGAGetOptDirFlag (ctx);
    assert (ind1->ctx == ind2->ctx);

    if (!ind1->evaluptodate) {
        PGAError
            ( ctx
            , "EvalCompare: first individual not up to date:"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    if (!ind2->evaluptodate) {
        PGAError
            ( ctx
            , "EvalCompare: second individual not up to date:"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    if (ctx->ga.NumConstraint > 0) {
        auxt1 = INDGetAuxTotal (ind1);
        auxt2 = INDGetAuxTotal (ind2);
    }
    if (ctx->ga.Epsilon) {
        auxt1 -= ctx->ga.Epsilon;
        auxt2 -= ctx->ga.Epsilon;
        if (auxt1 < 0) {
            auxt1 = 0;
        }
        if (auxt2 < 0) {
            auxt2 = 0;
        }
    }
    if (auxt1 || auxt2) {
        return CMP (auxt1, auxt2);
    }
    /* Multiobjective? */
    if (ctx->ga.NumAuxEval - ctx->ga.NumConstraint > 0) {
        return CMP (ind1->rank, ind2->rank);
    }
    /* We might use the fitness if both populations are the same
       otherwise fitness values are not comparable. But we now
       use the evaluation in any case. This avoids overflows on fitness
       computation, see documentation.
     */
    switch (dir) {
    case PGA_MAXIMIZE:
        return CMP (ind2->evalue, ind1->evalue);
        break;
    case PGA_MINIMIZE:
        return CMP (ind1->evalue, ind2->evalue);
        break;
    default:
        PGAError
            (ctx
            , "EvalCompare: Invalid value of PGAGetOptDirFlag:"
            , PGA_FATAL, PGA_INT, (void *) &dir
            );
        break;
    }
    /* notreached */
    return 0;
}

/*!****************************************************************************
    \brief Compare two strings for selection.
    \ingroup utility
    \param  ctx    context variable
    \param  p1     first string to compare
    \param  pop1   symbolic constant of population of first string
    \param  p2     second string to compare
    \param  pop2   symbolic constant of population of second string
    \return
    - >0 if p2 is "better" than p1
    - <0 if p1 is "better" than p2
    - 0  if both compare equal

    \rst

    Description
    -----------

    Thinks of this as sorting individuals by decreasing fitness or
    increasing constraint violations.

    This typically simply compares evaluation taking into account the
    evaluation direction (minimize/maximize).
    But if auxiliary evaluations are defined, the auxiliary evaluations are
    treated as constraints or for multi-objective optimization.
    The default handling of auxiliary evaluations is incompatible with
    certain selection schemes, see checks in create.c
    We handle constraints to compare first: If two constrained individuals
    are compared, the one with less constraint violations wins. If a
    constrained individual is compared to an unconstrained one, the latter
    wins. If two unconstrained individuals are compared, the (single)
    evaluation is compared depending on the direction of optimization
    (minimization or maximization).

    For multi-objective optimization we do not compare the evaluations but
    only the rank (as computed by the NSGA-II algorithm). Note that many
    individuals may have the same rank.

    Note that :c:func:`PGAEvalCompare` is now used in several contexts,
    including finding the best evaluation. For very badly scaled
    problems, the default fitness computation will degenerate if there
    are very large evaluation values and very small ones. In that case
    the fitness will not reflect the evaluation. Therefore
    :c:func:`PGAEvalCompare` will now always sort on evaluation values
    ignoring the fitness. This improves Tournament selection for very
    badly scaled problems.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int result;

      ...
      result = PGAEvalCompare (ctx, p1, PGA_OLDPOP, p2, PGA_OLDPOP);

    \endrst

******************************************************************************/
int PGAEvalCompare (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAIndividual *ind1 = PGAGetIndividual (ctx, p1, pop1);
    PGAIndividual *ind2 = PGAGetIndividual (ctx, p2, pop2);
    return INDEvalCompare (ind1, ind2);
}

/*!****************************************************************************
    \brief Compare two PGAIndividual pointers for qsort
    \ingroup internal
    \param a1 void pointer to first individual
    \param a2 void pointer to second individual
    \return -1 if less, 0 if equal, 1 if greater

    \rst

    Description
    -----------

    Used as the comparison function in ``qsort``.

    \endrst

******************************************************************************/
int PGAEvalSortHelper (const void *a1, const void *a2)
{
    /* Note: We cast away the const because INDGetAuxTotal does caching
     * and writes to the individual
     */
    PGAIndividual **i1 = (void *)a1;
    PGAIndividual **i2 = (void *)a2;
    return INDEvalCompare (*i1, *i2);
}

/*!****************************************************************************
    \brief Sort list of population indeces by evaluation.
    \ingroup utility
    \param  ctx    context variable
    \param  pop    symbolic constant of population to sort
    \param  idx    Array of integer indeces into population pop
    \return Array idx sorted by evalation

    \rst

    Description
    -----------

    Note that the given list ``idx`` of indeces need not be initialized.
    It just has to have the correct size of course.
    Think of this as sorting individuals by "better" evaluation or
    increasing constraint violations. The best individals are sorted
    first.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int indexes [ctx->ga.PopSize];

      ...
      PGAEvalSort (ctx, PGA_OLDPOP, indexes);

    \endrst

******************************************************************************/
void PGAEvalSort (PGAContext *ctx, int pop, int *idx)
{
    int i;
    DECLARE_DYNARRAY (PGAIndividual *, sorttmp, ctx->ga.PopSize);
    /* No need to check validity of pop, done by PGAIndividual */
    PGAIndividual *first = PGAGetIndividual (ctx, 0, pop);

    for (i=0; i<ctx->ga.PopSize; i++) {
        sorttmp [i] = first + i;
    }
    qsort (sorttmp, ctx->ga.PopSize, sizeof (sorttmp [0]), PGAEvalSortHelper);
    for (i=0; i<ctx->ga.PopSize; i++) {
        idx [i] = sorttmp [i] - first;
    }
}

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

#define HASH_JEN_MIX(a,b,c)              \
do {                                     \
  a -= b; a -= c; a ^= ( c >> 13 );      \
  b -= c; b -= a; b ^= ( a << 8 );       \
  c -= a; c -= b; c ^= ( b >> 13 );      \
  a -= b; a -= c; a ^= ( c >> 12 );      \
  b -= c; b -= a; b ^= ( a << 16 );      \
  c -= a; c -= b; c ^= ( b >> 5 );       \
  a -= b; a -= c; a ^= ( c >> 3 );       \
  b -= c; b -= a; b ^= ( a << 10 );      \
  c -= a; c -= b; c ^= ( b >> 15 );      \
} while (0)

#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/*!****************************************************************************
    \brief Hashing utility function.
    \ingroup standard-api
    \param data  pointer to data to hash
    \param len   length of data
    \param hashv previous hash, use PGA_INITIAL_HASH if first
    \return hash over data

    \rst

    Description
    -----------

    This is Bob Jenkins' hash function as implemented in the
    `uthash project`_.
    The last parameter is the previous hash (if hashing different
    pieces), use :c:macro:`PGA_INITIAL_HASH` for the first hash.

    .. _`uthash project`:
        https://github.com/troydhanson/uthash/blob/master/src/uthash.h

    \endrst

******************************************************************************/

PGAHash PGAUtilHash (const void *data, size_t len, PGAHash hashv)
{
    PGAHash _hj_i, _hj_j, _hj_k;
    unsigned const char *_hj_key = data;
    _hj_i = _hj_j = 0x9e3779b9u;
    _hj_k = (PGAHash)(len);
    while (_hj_k >= 12U) {
        _hj_i +=    (_hj_key[0] + ( (PGAHash)_hj_key[1] << 8 )
            + ( (PGAHash)_hj_key[2] << 16 )
            + ( (PGAHash)_hj_key[3] << 24 ) );
        _hj_j +=    (_hj_key[4] + ( (PGAHash)_hj_key[5] << 8 )
            + ( (PGAHash)_hj_key[6] << 16 )
            + ( (PGAHash)_hj_key[7] << 24 ) );
        hashv += (_hj_key[8] + ( (PGAHash)_hj_key[9] << 8 )
            + ( (PGAHash)_hj_key[10] << 16 )
            + ( (PGAHash)_hj_key[11] << 24 ) );

        HASH_JEN_MIX(_hj_i, _hj_j, hashv);

        _hj_key += 12;
        _hj_k -= 12U;
    }
    hashv += (PGAHash)(len);
    switch ( _hj_k ) {
      case 11: hashv += ( (PGAHash)_hj_key[10] << 24 ); /* FALLTHROUGH */
      case 10: hashv += ( (PGAHash)_hj_key[9] << 16 );  /* FALLTHROUGH */
      case 9:  hashv += ( (PGAHash)_hj_key[8] << 8 );   /* FALLTHROUGH */
      case 8:  _hj_j += ( (PGAHash)_hj_key[7] << 24 );  /* FALLTHROUGH */
      case 7:  _hj_j += ( (PGAHash)_hj_key[6] << 16 );  /* FALLTHROUGH */
      case 6:  _hj_j += ( (PGAHash)_hj_key[5] << 8 );   /* FALLTHROUGH */
      case 5:  _hj_j += _hj_key[4];                      /* FALLTHROUGH */
      case 4:  _hj_i += ( (PGAHash)_hj_key[3] << 24 );  /* FALLTHROUGH */
      case 3:  _hj_i += ( (PGAHash)_hj_key[2] << 16 );  /* FALLTHROUGH */
      case 2:  _hj_i += ( (PGAHash)_hj_key[1] << 8 );   /* FALLTHROUGH */
      case 1:  _hj_i += _hj_key[0];                      /* FALLTHROUGH */
      default: ;
    }
    HASH_JEN_MIX(_hj_i, _hj_j, hashv);
    return hashv;
}

/*!****************************************************************************
    \brief Compute Hash for individual in given population
    \ingroup internal
    \param   ctx   context variable
    \param   p     string index
    \param   pop   symbolic constant of the population containing string p
    \return  Hash of given individual

    \rst

    Description
    -----------

    Call the hash user function, by default :c:func:`PGACreate` makes
    sure the correct one is set for the data type.

    \endrst

******************************************************************************/
size_t PGAIndividualHashIndex (PGAContext *ctx, int p, int pop)
{
    PGABinary hash = 0;
    assert (ctx->ga.NoDuplicates);

    if (ctx->fops.Hash) {
        int fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p + 1;
        hash = (*ctx->fops.Hash)(&ctx, &fp, &pop);
    } else {
        hash = (*ctx->cops.Hash)(ctx, p, pop);
    }
    return (size_t)(hash % ctx->ga.PopSize);
}

/*!****************************************************************************
    \brief Calculate the mean genetic distance for a population.
    \ingroup explicit
    \param  ctx  context variable
    \param  pop  symbolic constant of the population for which the
                 genetic distance is to be calculated
    \return The mean genetic distance in the population

    \rst

    Description
    -----------

    This function has effort :math:`O(n^2)` in the population size
    :math:`n`. It is very useful for detecting premature convergence and
    is used when genetic distance reporting has been turned on with
    :c:macro:`PGA_REPORT_GENE_DISTANCE`.

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        double gd;

        ...
        gd = PGAGeneDistance (ctx, PGA_NEWPOP);

    \endrst

******************************************************************************/

double PGAGeneDistance (PGAContext *ctx, int pop)
{
    int i, j, count = 0;
    double sum_dist = 0;
    for (i=0; i<ctx->ga.PopSize - 1; i++) {
        for (j=i + 1; j<ctx->ga.PopSize; j++) {
            sum_dist += PGAUserFunctionGeneDistance (ctx, i, pop, j, pop);
            count++;
        }
    }
    return sum_dist / count;
}
