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
* This file contains the routines specific to the floating point data
* structure.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include <pgapack.h>

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* Helper for bounds/bounce check */
static void bouncheck
    ( PGAContext *ctx, int idx, int boundflag, int bounceflag
    , PGAReal *child, PGAReal minp, PGAReal maxp
    )
{
    if (boundflag || bounceflag) {
        if (child [idx] < ctx->init.RealMin [idx]) {
            if (bounceflag) {
                child [idx] = PGARandomUniform
                    (ctx, ctx->init.RealMin [idx], minp);
            } else {
                child [idx] = ctx->init.RealMin [idx];
            }
        }
        if (child [idx] > ctx->init.RealMax [idx]) {
            if (bounceflag) {
                child [idx] = PGARandomUniform
                    (ctx, maxp, ctx->init.RealMax [idx]);
            } else {
                child [idx] = ctx->init.RealMax [idx];
            }
        }
    }
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Set the value of real-valued allele i in string p in
           population pop
    \ingroup allele
    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \param   i    allele index
    \param   val  real value to set the allele to
    \return  The specified allele in p is modified by side-effect.

    \rst

    Example
    -------

    Sets the value of allele ``i`` of string ``p`` in population
    :c:macro:`PGA_NEWPOP` to 1.57

    .. code-block:: c

       PGAContext *ctx;
       int i, p;

       ...
       PGASetRealAllele (ctx, p, PGA_NEWPOP, i, 1.57);

    \endrst

******************************************************************************/
void PGASetRealAllele (PGAContext *ctx, int p, int pop, int i, double val)
{
    PGAIndividual *ind;
    PGAReal      *chrom;

    PGADebugEntered  ("PGASetRealAllele");
    PGACheckDataType ("PGASetRealAllele", PGA_DATATYPE_REAL);

    ind = PGAGetIndividual (ctx, p, pop);
    chrom = (PGAReal *)ind->chrom;
    chrom [i] = val;

    PGADebugExited ("PGASetRealAllele");
}

/*!****************************************************************************
    \brief Return the value of real-valued allele i in string p in
           population pop
    \ingroup allele

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \param   i    allele index
    \return  The value of allele i

    \rst

    Example
    -------

    Returns the value of real-valued allele ``i`` of string ``p``
    in population :c:macro:`PGA_NEWPOP`

    .. code-block:: c

       PGAContext *ctx;
       int p, i
       double d;

       ...
       d = PGAGetRealAllele (ctx, p, PGA_NEWPOP, i);

    \endrst

******************************************************************************/
double PGAGetRealAllele (PGAContext *ctx, int p, int pop, int i)
{
    PGAIndividual *ind;
    PGAReal      *chrom;

    PGADebugEntered  ("PGAGetRealAllele");
    PGACheckDataType ("PGAGetRealAllele", PGA_DATATYPE_REAL);

    ind = PGAGetIndividual (ctx, p, pop);
    chrom = (PGAReal *)ind->chrom;

    PGADebugExited ("PGAGetRealAllele");

    return (double) chrom [i];
}

/*!****************************************************************************
    \brief Set the upper and lower bounds for randomly initializing
           real-valued genes.
    \ingroup init
    \param   ctx      context variable
    \param   median   an array containing the mean value of the interval
    \param   frac     an array containing the fraction of median to add and
                      subtract to/from the median to define the interval
    \return  None

    \rst

    Description
    -----------

    For each gene these bounds define an interval from which the initial
    allele value is selected uniformly randomly.  With this routine the
    user specifies a median value and a fraction of the median for each allele.


    Example
    -------

    Set the initialization routines to select a value for each real-valued
    gene ``i`` uniformly randomly from the interval :math:`[i-v,i+v]`, where
    :math:`v = i/2`.
    Assumes all strings are the same length.

    .. code-block:: c

       PGAContext *ctx;
       double *median, *frac;
       int i, stringlen;

       ...
       stringlen = PGAGetStringLength (ctx);
       median = malloc (stringlen * sizeof(double));
       frac   = malloc (stringlen * sizeof(double));
       for (i=0; i<stringlen; i++) {
          median [i] = (double) i;
          frac   [i] = 0.5;
       }
       PGASetRealInitPercent (ctx, median, frac);

    \endrst

******************************************************************************/
void PGASetRealInitFraction (PGAContext *ctx, double *median, double *frac)
{
    int i;
    int stringlen;
    double offset;

    PGADebugEntered  ("PGASetRealInitPercent");
    PGAFailIfSetUp   ("PGASetRealInitPercent");
    PGACheckDataType ("PGASetRealInitPercent", PGA_DATATYPE_REAL);

    stringlen = PGAGetStringLength (ctx);
    for (i=0; i<stringlen; i++) {
         offset = fabs (median [i] * frac [i]);
         ctx->init.RealMin [i] = median [i] - offset;
         ctx->init.RealMax [i] = median [i] + offset;
    }
    ctx->init.RealType = PGA_RINIT_PERCENT;

    PGADebugExited ("PGASetRealInitPercent");
}

/*!****************************************************************************
    \brief Set the upper and lower bounds for randomly initializing
           real-valued genes.
    \ingroup init

    \param   ctx context variable
    \param   min array containing the lower bound of the interval for each gene
    \param   max array containing the upper bound of the interval for each gene
    \return  None

    \rst

    Description
    -----------

    For each gene these bounds define an interval from which the initial
    allele value is selected uniformly randomly.  The user specifies two
    arrays containing lower and upper bound for each gene to define the
    interval.  This is the default strategy for initializing real-valued
    strings.  The default interval is :math:`[0,1.0]` for each gene.

    Example
    -------

    Set the initialization routines to select a value for each real-valued
    gene :math:`i` uniformly randomly from the interval :math:`[-10.,i]`
    Assumes all strings are of the same length.

    .. code-block:: c

       PGAContext *ctx;
       double *low, *high;
       int i, stringlen;

       ...
       stringlen = PGAGetStringLength (ctx);
       low  = malloc (stringlen * sizeof(double));
       high = malloc (stringlen * sizeof(double));
       for (i=0; i<stringlen; i++) {
          low  [i] = -10.0;
          high [i] = i;
       }
       PGASetRealInitRange (ctx, low, high);

    \endrst

******************************************************************************/
void PGASetRealInitRange (PGAContext *ctx, const double *min, const double *max)
{
    int i;
    PGADebugEntered  ("PGASetRealInitRange");
    PGAFailIfSetUp   ("PGASetRealInitRange");
    PGACheckDataType ("PGASetRealInitRange", PGA_DATATYPE_REAL);

    for (i=ctx->ga.StringLen-1; i>=0; i--) {
        if (max [i] < min [i]) {
            PGAError
                ( ctx
                , "PGASetRealInitRange: Lower bound exceeds upper bound "
                  "for allele #"
                , PGA_FATAL, PGA_INT, (void *) &i
                );
        } else {
             ctx->init.RealMin [i] = min [i];
             ctx->init.RealMax [i] = max [i];
        }
    }
    ctx->init.RealType = PGA_RINIT_RANGE;

    PGADebugExited ("PGASetRealInitRange");
}

/*!***************************************************************************
    \brief Returns the minimum value used to randomly initialize allele
           i in a real string.
    \ingroup query
    \param   ctx  context variable
    \param   i    an allele position
    \return  The minimum value used to randomly initialize allele i

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        int min;

        ...
        min = PGAGetMinRealInitValue (ctx, 0);
    \endrst

*****************************************************************************/
double PGAGetMinRealInitValue (PGAContext *ctx, int i)
{
    PGADebugEntered   ("PGAGetMinRealInitValue");
    PGAFailIfNotSetUp ("PGAGetMinRealInitValue");
    PGACheckDataType  ("PGAGetMinRealInitValue", PGA_DATATYPE_REAL);

    if (i < 0 || i >= ctx->ga.StringLen) {
        PGAError
            ( ctx, "PGAGetMinRealInitValue: Index out of range:"
            , PGA_FATAL, PGA_INT, (int *) &i
            );
    }

    PGADebugExited ("PGAGetMinRealInitValue");

    return ctx->init.RealMin [i];
}

/*!***************************************************************************
    \brief Return the maximum value used to randomly initialize allele i
           in a real string.
    \ingroup query
    \param   ctx  context variable
    \param   i    an allele position
    \return  The maximum value used to randomly initialize allele i

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        int max;

        ...
        max = PGAGetMaxRealInitValue (ctx, 0);
    \endrst

*****************************************************************************/
double PGAGetMaxRealInitValue (PGAContext *ctx, int i)
{
    PGADebugEntered   ("PGAGetMaxRealInitValue");
    PGAFailIfNotSetUp ("PGAGetMaxRealInitValue");
    PGACheckDataType  ("PGAGetMaxRealInitValue", PGA_DATATYPE_REAL);

    if (i < 0 || i >= ctx->ga.StringLen) {
        PGAError
            ( ctx, "PGAGetMaxRealInitValue: Index out of range:"
            , PGA_FATAL, PGA_INT, (int *) &i
            );
    }

    PGADebugExited ("PGAGetMaxRealInitValue");

    return ctx->init.RealMax [i];
}

/*!***************************************************************************
    \brief Return the type of scheme used to randomly initialize strings
           of data type real.
    \ingroup query
    \param  ctx  context variable
    \return Returns the integer corresponding to the symbolic constant
            used to specify the scheme used to initialize real strings

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        int inittype;

        ...
        inittype = PGAGetRealInitType (ctx);
        switch (inittype) {
        case PGA_RINIT_PERCENT:
            printf ("Data Type = PGA_RINIT_PERCENT\n");
            break;
        case PGA_RINIT_RANGE:
            printf ("Data Type = PGA_RINIT_RANGE\n");
            break;
        }
    \endrst

*****************************************************************************/
int PGAGetRealInitType (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetRealInitType");
    PGAFailIfNotSetUp ("PGAGetRealInitType");
    PGACheckDataType  ("PGAGetRealInitType", PGA_DATATYPE_REAL);

    PGADebugExited ("PGAGetRealInitType");

    return ctx->init.RealType;
}

/*!****************************************************************************
    \brief Allocate memory for a string of type real.
    \ingroup internal
    \param   ctx       context variable
    \param   p         string index
    \param   pop       symbolic constant of the population string p is in
    \param   initflag  A boolean flag to indicate random initialization
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the create
    string user function for the real datatype by default.
    Parameter ``initflag`` is used in conjunction with
    ``ctx->ga.RandomInit`` to initialize the string either randomly or
    set to zero.

    Example
    -------

    Allocates memory and assigns the address of the allocated memory to
    the real string field ``ind->chrom`` of the individual.  Also, clears
    the string.

    .. code-block:: c

       PGAContext *ctx;
       int p;

       ...
       PGARealCreateString (ctx, p, PGA_NEWPOP, PGA_FALSE);
    \endrst

******************************************************************************/
void PGARealCreateString (PGAContext *ctx, int p, int pop, int initflag)
{
    PGAIndividual *new = PGAGetIndividual (ctx, p, pop);
    int i, fp;
    PGAReal *c;

    PGADebugEntered ("PGARealCreateString");

    new->chrom = (void *) malloc (ctx->ga.StringLen * sizeof(PGAReal));
    if (new->chrom == NULL) {
        PGAError
            ( ctx, "PGARealCreateString: No room to allocate new->chrom"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    c = (PGAReal *)new->chrom;
    if (initflag) {
        if (ctx->fops.InitString) {
            fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
            (*ctx->fops.InitString)(&ctx, &fp, &pop);
        } else {
            (*ctx->cops.InitString)(ctx, p, pop);
        }
    } else {
        for (i=ctx->ga.StringLen-1; i>=0; i--) {
            c[i] = 0.0;
        }
    }

    PGADebugExited ("PGARealCreateString");
}

/*!****************************************************************************
    \brief Randomly mutates a floating point string with probability mr.
    \ingroup internal
    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population string p is in
    \param   mr   probability of mutating a real-valued gene
    \return  The number of mutations performed

    \rst

    Description
    -----------

    Three of the four mutation operators are of the form
    :math:`v = v +- p \cdot v`.
    That is, the new value of :math:`v` (allele :math:`i`) is the old
    value + or - a percentage, :math:`p`, of the old value. There are
    three possibilities for choosing :math:`p`: (1) constant value (0.01
    by default), (2) selected uniformly on :math:`(0,\text{UB})`
    (UB is .1 by default), and (3) selected from a Gaussian distribution
    (with mean 0 and standard deviation .1 be default).  The change to
    an allele, :math:`p \cdot v`, is added or subtracted to the old value with
    a probability of .5. The fourth option is to replace :math:`v` with a value
    selected uniformly random from the initialization range of that
    gene. Alleles to mutate are randomly selected.  The value set by the
    routine :c:func:`PGASetMutationRealValue` is used as :math:`p`, UB,
    and sigma in cases 1, 2, and 3, respectively.

    Note that this function is set in :c:func:`PGASetUp` as the mutation
    user function for the real datatype by default.

    Example
    -------

    Mutate string ``p`` in population :c:macro:`PGA_NEWPOP` with
    probability 0.001.

    .. code-block:: c

       PGAContext *ctx;
       int NumMutations, p;

       ...
       NumMutations = PGARealMutation (ctx, p, PGA_NEWPOP, .001);

    \endrst

******************************************************************************/
int PGARealMutation (PGAContext *ctx, int p, int pop, double mr)
{
    PGAReal *c;
    int i, j;
    int count = 0;
    double val = 0.0;
    /* The following are used for DE variants */
    int midx = 0;
    int do_best = (ctx->ga.DEVariant == PGA_DE_VARIANT_BEST);
    int nrand  = 2 * ctx->ga.DENumDiffs + (!do_best);
    int maxidx = 2 * ctx->ga.DENumDiffs + 1;
    DECLARE_DYNARRAY (PGAReal *, indivs, maxidx);
    int do_crossover = 1;
    static double de_dither = 0.0;
    static int last_iter = -1;


    PGADebugEntered ("PGARealMutation");
    if (ctx->ga.MutationType == PGA_MUTATION_DE) {
        DECLARE_DYNARRAY (int, idx, maxidx);
        PGASampleState sstate;
        int best = 0;
        int avoid = 1;
        if (do_best) {
            best = PGAGetBestIndex (ctx, PGA_OLDPOP);
        }

        if (ctx->ga.DEDither > 0) {
            if (ctx->ga.DEDitherPerIndividual || last_iter != ctx->ga.iter) {
                de_dither = ctx->ga.DEDither * (PGARandom01 (ctx, 0) - 0.5);
                last_iter = ctx->ga.iter;
            }
        }

        /* We rely on the fact that we operate on an individual of
         * PGA_NEWPOP and take data from PGA_OLDPOP. Assert this is the
         * case.
         */
        if (pop != PGA_NEWPOP) {
            PGAError
                ( ctx, "PGARealMutation: Invalid value of pop:"
                , PGA_FATAL, PGA_INT, (void *) &(pop)
                );
        }

        /* for BEST strategy we have to avoid collision with running
         * index p *and* the best individual
         */
        if (do_best && best != p) {
            avoid = 2;
        }
        /* Use (PopSize - avoid) to avoid collision with running index
         * and with the best index (unless this is the same)
         * We may use this form of selection of the indeces needed for
         * DE for linear selection, for other selection schemes we use
         * the selected individuals and re-sample if we get collisions.
         */
        if (ctx->ga.SelectType == PGA_SELECT_LINEAR) {
            PGARandomSampleInit (ctx, &sstate, nrand, ctx->ga.PopSize - avoid);
            for (i=0; i<nrand; i++) {
                int rawidx = PGARandomNextSample (&sstate);
                /* Avoid collision with p and optionally best, samples
                 * are drawn with reduced upper bound, see
                 * PGARandomSampleInit above
                 */
                if (rawidx >= p) {
                    rawidx += 1;
                }
                /* The best individual can be in the part of the old
                 * population that is copied verbatim to the new
                 * population. In that case never increment the index
                 */
                if (do_best && best != p && rawidx >= best) {
                    rawidx += 1;
                }
                idx [i] = rawidx;
            }
        } else {
            for (i=0; i<nrand; i++) {
                do {
                    idx [i] = PGASelectNextIndex (ctx, PGA_OLDPOP);
                } while
                    (  idx [i] == p
                    || (do_best && idx [i] == best)
                    || (i > 0 && idx [i] == idx [i-1])
                    || (i > 1 && idx [i] == idx [i-2])
                    || (i > 2 && idx [i] == idx [i-3])
                    || (i > 3 && idx [i] == idx [i-4])
                    );
            }
        }
        /* Since indices from PGARandomNextSample are
         * returned in order we need to shuffle
         */
        PGAShuffle (ctx, idx, nrand);
        /* Now we have a list of shuffled indexes that do not
         * collide with one-another or with p or best
         * Add best index as last to the list if do_best
         */
        if (do_best) {
            idx [maxidx - 1] = best;
        }
        for (i=0; i<maxidx; i++) {
            indivs [i] = (PGAReal *)
                PGAGetIndividual (ctx, idx [i], PGA_OLDPOP)->chrom;
        }
        /* Index of allele that is mutated in any case */
        midx = PGARandomInterval (ctx, 0, ctx->ga.StringLen - 1);
    }

    c = (PGAReal *)PGAGetIndividual (ctx, p, pop)->chrom;
    for (i=0; i<ctx->ga.StringLen; i++) {
        double old_value = c [i];
        int idx = i;

        switch (ctx->ga.MutationType) {
            case PGA_MUTATION_RANGE:
            case PGA_MUTATION_CONSTANT:
            case PGA_MUTATION_UNIFORM:
            case PGA_MUTATION_GAUSSIAN:
            case PGA_MUTATION_POLY:
                /* randomly choose an allele   */
                if ( PGARandomFlip(ctx, mr) ) {
                    /* generate on range, or calculate multiplier */
                    switch (ctx->ga.MutationType) {
                      case PGA_MUTATION_RANGE:
                        c [i] = PGARandomUniform
                            (ctx, ctx->init.RealMin [i], ctx->init.RealMax [i]);
                        break;
                      case PGA_MUTATION_CONSTANT:
                        val = ctx->ga.MutateRealValue;
                        break;
                      case PGA_MUTATION_UNIFORM:
                        val = PGARandomUniform
                            (ctx, 0.0, ctx->ga.MutateRealValue);
                        break;
                      case PGA_MUTATION_GAUSSIAN:
                        val = PGARandomGaussian
                            (ctx, 0.0, ctx->ga.MutateRealValue);
                        break;
                      case PGA_MUTATION_POLY:
                      {
                        double u = PGARandom01 (ctx, 0);
                        double eta = PGAGetMutationPolyEta (ctx) + 1;
                        double delta;
                        if (u < 0.5) {
                            delta = pow (2 * u, 1.0 / eta) - 1.0;
                        } else {
                            delta = 1.0 - pow (2 * (1 - u), 1.0 / eta);
                        }
                        if (ctx->ga.MutatePolyValue >= 0) {
                            c [i] += delta * ctx->ga.MutatePolyValue;
                        } else {
                            if (delta < 0) {
                                val = fabs (c [i] - ctx->init.RealMin [i]);
                            } else {
                                val = fabs (ctx->init.RealMax [i] - c [i]);
                            }
                            c [i] += delta * val;
                        }
                        break;
                      }
                    }
                    /* apply multiplier calculated in switch above */
                    if ( (ctx->ga.MutationType == PGA_MUTATION_CONSTANT) ||
                         (ctx->ga.MutationType == PGA_MUTATION_UNIFORM)  ||
                         (ctx->ga.MutationType == PGA_MUTATION_GAUSSIAN)
                       )
                    {
                         /* add/subtract from allele */
                        if ( PGARandomFlip(ctx, .5) )
                            c [i] += val * c [i];
                        else
                            c [i] -= val * c [i];
                    }
                    /* increment mutation count */
                    count++;
                }
                break;
            case PGA_MUTATION_DE:
                switch (ctx->ga.DECrossoverType) {
                case PGA_DE_CROSSOVER_BIN:
                    do_crossover =
                        (  idx == midx
                        || PGARandomFlip (ctx, ctx->ga.DECrossoverProb)
                        );
                    break;
                case PGA_DE_CROSSOVER_EXP:
                    /* The first index copied is midx, then all indices
                     * are copied while the coin flip is valid
                     */
                    if (do_crossover) {
                        idx = (midx + i) % ctx->ga.StringLen;
                        if (i > 0) {
                            do_crossover =
                                (PGARandomFlip (ctx, ctx->ga.DECrossoverProb));
                        }
                    }
                    break;
                default:
                    PGAError
                        ( ctx, "PGARealMutation: Invalid DE crossover type:"
                        , PGA_FATAL, PGA_INT
                        , (void *) &(ctx->ga.DECrossoverType)
                        );
                    break;
                }
                if (do_crossover){
                    double f = ctx->ga.DEScaleFactor + de_dither;
                    if (ctx->ga.DEJitter > 0) {
                        f += ctx->ga.DEJitter * (PGARandom01 (ctx, 0) - 0.5);
                    }
                    switch (ctx->ga.DEVariant) {
                    case PGA_DE_VARIANT_RAND:
                    case PGA_DE_VARIANT_BEST:
                        /* the last element is either a random individual
                         * or the best depending on variant
                         */
                        c [idx] = indivs [maxidx - 1][idx];
                        /* Add difference vectors */
                        for (j=0; j < (maxidx - 1); j+=2) {
                            c [idx] +=
                                f * (indivs [j][idx] - indivs [j+1][idx]);
                        }
                        break;
                    case PGA_DE_VARIANT_EITHER_OR:
                        /* We use only 1 difference and ignore DENumDiffs */
                        if (PGARandom01 (ctx, 0) < ctx->ga.DEProbabilityEO) {
                            c [idx] = indivs [0][idx]
                                  + f * (indivs [1][idx] - indivs [2][idx]);
                        } else {
                            double k = ctx->ga.DEAuxFactor;
                            c [idx] = indivs [0][idx]
                                  + k * ( indivs [1][idx] + indivs [2][idx]
                                        - 2 * indivs [0][idx]
                                        );
                        }
                        break;
                    default:
                        PGAError(ctx, "PGARealMutation: Invalid value of "
                                 "ga.DEVariant:", PGA_FATAL, PGA_INT,
                                 (void *) &(ctx->ga.DEVariant));
                        break;
                    }
                    count++;
                }
                break;
            default:
                PGAError(ctx, "PGARealMutation: Invalid value of "
                         "ga.MutationType:", PGA_FATAL, PGA_INT,
                         (void *) &(ctx->ga.MutationType));
                break;
        }

        /* reset to min/max or bounce if outside range */
        bouncheck
            ( ctx, idx, ctx->ga.MutateBoundedFlag, ctx->ga.MutateBounceFlag
            , c, old_value, old_value
            );
    }

    PGADebugExited ("PGARealMutation");

    return count;
}

/*!****************************************************************************
    \brief This routine performs one point crossover on two parent
           strings, producing (via side effect) the crossed children
           child1 and child2.
    \ingroup internal
    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing
                   string p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population to contain
                   string c1 and c2
    \return  c1 and c2 in population pop2 are modified by side-effect

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the real datatype when selecting
    one-point crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGARealOneptCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);
    \endrst

******************************************************************************/
void PGARealOneptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAReal *parent1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *parent2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    PGAReal *child1  = (PGAReal *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    PGAReal *child2  = (PGAReal *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    int i, xsite;

    PGADebugEntered ("PGARealOneptCrossover");

    xsite = PGARandomInterval (ctx, 1, ctx->ga.StringLen - 1);

    for (i=0; i<xsite; i++) {
        child1 [i] = parent1 [i];
        child2 [i] = parent2 [i];
    }

    for (i=xsite; i<ctx->ga.StringLen; i++) {
        child1 [i] = parent2 [i];
        child2 [i] = parent1 [i];
    }

    PGADebugExited ("PGARealOneptCrossover");
}


/*!****************************************************************************
    \brief Perform two-point crossover on two parent strings producing
           two children via side-effect.
    \ingroup internal
    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing
                   string p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population to contain
                   string c1 and c2
    \return  c1 and c2 in population pop2 are modified by side-effect

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the real datatype when selecting
    two-point crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGARealTwoptCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);
    \endrst

******************************************************************************/
void PGARealTwoptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAReal *parent1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *parent2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    PGAReal *child1  = (PGAReal *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    PGAReal *child2  = (PGAReal *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    int i, temp, xsite1, xsite2;

    PGADebugEntered ("PGARealTwoptCrossover");

    /* pick two cross sites such that xsite2 > xsite1 */
    xsite1 = PGARandomInterval (ctx, 1, ctx->ga.StringLen - 1);
    xsite2 = xsite1;
    while (xsite2 == xsite1) {
        xsite2 = PGARandomInterval (ctx, 1, ctx->ga.StringLen - 1);
    }
    if (xsite1 > xsite2) {
        temp   = xsite1;
        xsite1 = xsite2;
        xsite2 = temp;
    }

    for (i=0; i<xsite1; i++) {
        child1 [i] = parent1 [i];
        child2 [i] = parent2 [i];
    }

    for (i=xsite1; i<xsite2; i++) {
        child1 [i] = parent2 [i];
        child2 [i] = parent1 [i];
    }

    for (i=xsite2; i<ctx->ga.StringLen; i++) {
        child1 [i] = parent1 [i];
        child2 [i] = parent2 [i];
    }

    PGADebugExited ("PGARealTwoptCrossover");
}


/*!****************************************************************************
    \brief Perform uniform crossover on two parent strings producing two
           children via side-effect
    \ingroup internal
    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing
                   string p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population to contain
                   string c1 and c2
    \return  c1 and c2 in population pop2 are modified by side-effect

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the real datatype when selecting
    uniform crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGARealUniformCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);
    \endrst

******************************************************************************/
void PGARealUniformCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAReal *parent1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *parent2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    PGAReal *child1  = (PGAReal *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    PGAReal *child2  = (PGAReal *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    int i;

    PGADebugEntered ("PGARealUniformCrossover");

    for (i=0; i<ctx->ga.StringLen; i++) {
        if (parent1 [i] == parent2 [i]) {
            child1 [i] = parent1 [i];
            child2 [i] = parent2 [i];
        } else {
            if (PGARandomFlip(ctx, ctx->ga.UniformCrossProb)) {
                child1 [i] = parent1 [i];
                child2 [i] = parent2 [i];
            } else {
                child1 [i] = parent2 [i];
                child2 [i] = parent1 [i];
            }
        }
    }

    PGADebugExited ("PGARealUniformCrossover");
}

/*!****************************************************************************
    \brief Perform simulated binary crossover (SBX) on two parent
           strings producing two children via side-effect.
    \ingroup internal
    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing
                   string p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population to contain
                   string c1 and c2
    \return  c1 and c2 in population pop2 are modified by side-effect

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the real datatype when selecting
    simulated binary crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGARealSBXCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);
    \endrst

******************************************************************************/
void PGARealSBXCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAReal *parent1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *parent2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    PGAReal *child1  = (PGAReal *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    PGAReal *child2  = (PGAReal *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    int i;
    double u = 0;

    if (ctx->ga.CrossSBXOnce) {
        u = PGARandom01 (ctx, 0);
    }

    for (i=0; i<ctx->ga.StringLen; i++) {
        if (  parent1 [i] == parent2 [i]
           || !PGARandomFlip (ctx, ctx->ga.UniformCrossProb)
           )
        {
            child1 [i] = parent1 [i];
            child2 [i] = parent2 [i];
        } else {
            int j;
            PGAReal minp =
                parent1 [i] < parent2 [i] ? parent1 [i] : parent2 [i];
            PGAReal maxp =
                parent1 [i] > parent2 [i] ? parent1 [i] : parent2 [i];
            if (!ctx->ga.CrossSBXOnce) {
                u = PGARandom01 (ctx, 0);
            }
            PGACrossoverSBX
                (ctx, parent1 [i], parent2 [i], u, child1 + i, child2 + i);
            for (j=0; j<2; j++) {
                bouncheck
                    ( ctx, i, ctx->ga.CrossBoundedFlag, ctx->ga.CrossBounceFlag
                    , j ? child2 : child1, minp, maxp
                    );
            }
        }
    }
}

/*!****************************************************************************
    \brief Write a real-valued string to a file.
    \ingroup internal
    \param   ctx  context variable
    \param   fp   file pointer to file to write the string to
    \param   p    index of the string to write out
    \param   pop  symbolic constant of the population string p is in
    \return  None

    \rst

    Example
    -------

    Write string ``s`` to ``stdout``.

    .. code-block:: c

       PGAContext *ctx;
       int s;

       ...
       PGARealPrintString (ctx, stdout, s, PGA_NEWPOP);
    \endrst

******************************************************************************/
void PGARealPrintString (PGAContext *ctx, FILE *fp, int p, int pop)
{
    PGAReal *c = (PGAReal *)PGAGetIndividual (ctx, p, pop)->chrom;
    int i;

    PGADebugEntered ("PGARealPrintString");

    for (i = 0; i < ctx->ga.StringLen; i++) {
        switch (i % 5) {
        case 0:
            fprintf (fp, "#%4d: [%11.7g]",i,c[i]);
            break;
        case 1:
        case 2:
        case 3:
            fprintf (fp, ", [%11.7g]",c[i]);
            break;
        case 4:
            fprintf (fp, ", [%11.7g]",c[i]);
            if (i+1 < ctx->ga.StringLen) {
                fprintf (fp, "\n");
            }
            break;
        }
    }
    fprintf (fp, "\n");

    PGADebugExited ("PGARealPrintString");
}


/*!****************************************************************************
    \brief Copy one real-valued string string to another
    \ingroup internal
    \param   ctx   context variable
    \param   p1    string to copy
    \param   pop1  symbolic constant of population containing string p1
    \param   p2    string to copy p1 to
    \param   pop2  symbolic constant of population containing string p2
    \return  String p2 in population pop2 is modified to be a copy of
             string p1 in population pop1.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the copy
    string user function for the real datatype by default.

    Example
    -------

    Copy string ``x`` in old population to ``y`` in new population.

    .. code-block:: c

       PGAContext *ctx;
       int x, y;

       ...
       PGARealCopyString (ctx, x, PGA_OLDPOP, y, PGA_NEWPOP);
    \endrst

******************************************************************************/
void PGARealCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAReal *source = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *dest   = (PGAReal *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int i;

    PGADebugEntered ("PGARealCopyString");

    for (i=ctx->ga.StringLen-1; i>=0; i--) {
        *(dest++) = *(source++);
    }

    PGADebugExited ("PGARealCopyString");
}

/*!****************************************************************************
    \brief Returns true if real-valued string a is a duplicate of
           real-valued string b, else returns false.
    \ingroup internal
    \param   ctx   context variable
    \param   p1    string index of the first string to compare
    \param   pop1  symbolic constant of the population string p1 is in
    \param   p2    string index of the second string to compare
    \param   pop2  symbolic constant of the population string p2 is in
    \return  Returns true/false if strings are duplicates

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    duplicate checking user function for the real datatype by default.

    Example
    -------

    Compare strings ``x`` with ``y`` to see if they are duplicates.

    .. code-block:: c

       PGAContext *ctx;
       int x, y;

       ...
       if (PGARealDuplicate (ctx, x, PGA_OLDPOP, y, PGA_OLDPOP)) {
           printf ("strings are duplicates\n");
       }
    \endrst

******************************************************************************/
int PGARealDuplicate (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAReal *a = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *b = (PGAReal *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int count = 0;

    PGADebugEntered ("PGARealDuplicate");

    for (count=0; count<ctx->ga.StringLen; count++) {
        if (a [count] != b [count]) {
            break;
        }
    }

    PGADebugExited ("PGARealDuplicate");

    return count == ctx->ga.StringLen ? PGA_TRUE : PGA_FALSE;
}

/*!****************************************************************************
    \brief Return hash value of given gene.
    \ingroup internal
    \param   ctx   context variable
    \param   p     string index of the string to hash
    \param   pop   symbolic constant of the population string p is in
    \return  Hash value for string

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    hash user function for the real datatype by default.

    \endrst

******************************************************************************/
PGAHash PGARealHash (PGAContext *ctx, int p, int pop)
{
    void *a = PGAGetIndividual (ctx, p, pop)->chrom;
    PGAHash hash = PGAUtilHash
        (a, sizeof (PGAReal) * ctx->ga.StringLen, PGA_INITIAL_HASH);
    return hash;
}

/*!****************************************************************************
    \brief Randomly initialize a string of type real.
    \ingroup internal
    \param   ctx  context variable
    \param   p    index of string to randomly initialize
    \param   pop  symbolic constant of the population string p is in
    \return  String p in population pop is randomly initialized by side-effect

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    init string user function for the real datatype by default.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p;

       ...
       PGARealInitString (ctx, p, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGARealInitString (PGAContext *ctx, int p, int pop)
{
    int i;
    PGAReal *c = (PGAReal *)PGAGetIndividual (ctx, p, pop)->chrom;

    PGADebugEntered ("PGARealInitString");

    for (i = 0; i < ctx->ga.StringLen; i++) {
        c[i] = PGARandomUniform
            (ctx, ctx->init.RealMin[i], ctx->init.RealMax[i]);
    }

    PGADebugExited ("PGARealInitString");
}

/*!****************************************************************************
    \brief Build an MPI datatype for a string.
    \ingroup internal
    \param   ctx    context variable
    \param   p      index of string
    \param   pop    symbolic constant of population string p is in
    \return  An MPI_Datatype

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext   *ctx;
       int           p;
       MPI_Datatype  dt;

       ...
       dt = PGARealBuildDatatype(ctx, p, pop);
    \endrst

******************************************************************************/
MPI_Datatype PGARealBuildDatatype (PGAContext *ctx, int p, int pop)
{
    int            idx = 0;
    /* Number of elements in each block (array of integer) */
    int            counts [PGA_MPI_HEADER_ELEMENTS + 1];
    /* byte displacement of each block (array of integer) */
    MPI_Aint       displs [PGA_MPI_HEADER_ELEMENTS + 1];
    /* type of elements in each block (array of handles to datatype objects) */
    MPI_Datatype   types  [PGA_MPI_HEADER_ELEMENTS + 1];
    MPI_Datatype   individualtype; /* new datatype (handle) */
    PGAIndividual *traveller;      /* address of individual in question */

    PGADebugEntered  ("PGARealBuildDatatype");

    traveller = PGAGetIndividual  (ctx, p, pop);
    idx = PGABuildDatatypeHeader (ctx, p, pop, counts, displs, types);

    MPI_Get_address (traveller->chrom, &displs [idx]);
    counts [idx] = ctx->ga.StringLen;
    types  [idx] = MPI_DOUBLE;
    idx++;

    MPI_Type_create_struct (idx, counts, displs, types, &individualtype);
    MPI_Type_commit (&individualtype);

    PGADebugExited ("PGARealBuildDatatype");

    return individualtype;
}

/*!****************************************************************************
    \brief Compute genetic difference of two strings.
    \ingroup internal

    \param   ctx    context variable
    \param   p1     first string index
    \param   pop1   symbolic constant of the population the first string is in
    \param   p2     second string index
    \param   pop2   symbolic constant of the population the second string is in
    \return  genetic distance of the two strings

    \rst

    Description
    -----------

    Sum of the absolute values of the differences of each allele.
    So this is a Manhattan distance (mainly for performance reasons).
    Internal function.  Use :c:func:`PGAUserFunctionGeneDistance`.

    \endrst

******************************************************************************/
double PGARealGeneDistance (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAReal *c1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *c2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    double ret = 0.0;
    int i;

    PGADebugEntered ("PGARealGeneDistance");
    for (i=0; i<ctx->ga.StringLen; i++) {
        ret += fabs (c1 [i] - c2 [i]);
    }
    PGADebugExited ("PGARealGeneDistance");
    return ret;
}

/*!****************************************************************************
    \brief Compute genetic difference of two strings.
    \ingroup standard-api
    \param   ctx    context variable
    \param   p1     first string index
    \param   pop1   symbolic constant of the population the first string is in
    \param   p2     second string index
    \param   pop2   symbolic constant of the population the second string is in
    \return  Genetic euclidian distance of the two strings

    \rst

    Description
    -----------

    This uses the Euclidian distance metric, the square-root of the sum
    of all squared differences of each allele. It can be used to
    override the default real genetic distance function (which uses a
    manhattan distance metric) using :c:func:`PGASetUserFunction` with
    the :c:macro:`PGA_USERFUNCTION_GEN_DISTANCE` setting.

    Example
    -------

    Override genetic distance function:

    .. code-block:: c

       PGAContext *ctx;

       ...
       assert (PGAGetDataType (ctx) == PGA_DATATYPE_REAL);
       PGASetUserFunction
        (ctx, PGA_USERFUNCTION_GEN_DISTANCE, PGARealEuclidianDistance);

    \endrst

******************************************************************************/
double PGARealEuclidianDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAReal *c1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *c2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop2)->chrom;

    PGADebugEntered ("PGARealGeneDistance");
    PGADebugExited  ("PGARealGeneDistance");
    return LIN_euclidian_distance (ctx->ga.StringLen, c1, c2);
}
