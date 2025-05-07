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
* This file contains the routines specific to the integer data structure.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include <stdint.h>
#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* Helper for bounds/bounce check */
static void bouncheck
    ( PGAContext *ctx, int idx, int boundflag, int bounceflag
    , PGAInteger *child, PGAInteger minp, PGAInteger maxp
    )
{
    if (boundflag || bounceflag) {
        if (child [idx] < ctx->init.IntegerMin [idx]) {
            if (bounceflag) {
                child [idx] = PGARandomInterval
                    (ctx, ctx->init.IntegerMin [idx], minp);
            } else {
                child [idx] = ctx->init.IntegerMin [idx];
            }
        }
        if (child [idx] > ctx->init.IntegerMax [idx]) {
            if (bounceflag) {
                child [idx] = PGARandomInterval
                    (ctx, maxp, ctx->init.IntegerMax [idx]);
            } else {
                child [idx] = ctx->init.IntegerMax [idx];
            }
        }
    }
}
/* Helper for sorting / searching */
static int intcmp (const void *v1, const void *v2)
{
    const PGAInteger *i1 = v1;
    const PGAInteger *i2 = v2;
    if (*i1 < *i2) {
        return -1;
    }
    if (*i1 > *i2) {
        return 1;
    }
    return 0;
}
/* Helpers for checking an allele against fixed edges */

/* Search for node in rev edges and return pointer or NULL if not found */
static PGAInteger (*rev_edge (PGAContext *ctx, PGAInteger node))[2]
{
    if (ctx->ga.n_edges == 0) {
        return NULL;
    }
    return bsearch
        ( &node
        , ctx->ga.r_edge
        , ctx->ga.n_edges
        , 2 * sizeof (PGAInteger)
        , intcmp
        );
}

/* Search for node in forward edges and return pointer or NULL if not found */
static PGAFixedEdge *get_edge (PGAContext *ctx, PGAInteger node)
{
    if (ctx->ga.n_edges == 0) {
        return NULL;
    }
    return bsearch
        ( &node
        , ctx->ga.edges
        , ctx->ga.n_edges
        , sizeof (PGAFixedEdge)
        , intcmp
        );
}

#ifdef DEBUG
static void assert_has_edges (PGAContext *ctx, PGAInteger *a)
{
    int i;
    int l = ctx->ga.StringLen;
    if (!ctx->ga.n_edges) {
        return;
    }
    for (i=0; i<l; i++) {
        PGAFixedEdge *p = get_edge (ctx, a [i]);
        if (p && !p->prev) {
            int j = 0, inc = 0;
            while (p) {
                /* Construct to easily set breakpoint on failing assert */
                if (a [(i + j + l) % l] != p->lhs) {
                    assert (a [(i + j + l) % l] == p->lhs);
                }
                /* First item */
                if (inc == 0) {
                    if (a [(i + 1) % l] == p->rhs) {
                        inc = 1;
                    } else if (a [(i - 1 + l) % l] == p->rhs) {
                        inc = -1;
                    } else {
                        assert (0);
                    }
                } else { /* Subsequent items */
                    assert (a [(i + j + inc + l) % l] == p->rhs);
                }
                p = p->next;
                j += inc;
            }
            if (j > 0) {
                i += j;
            }
        }
    }
}
#endif /* DEBUG */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/*!****************************************************************************
    \brief Set the value of a (integer) allele.
    \ingroup allele

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \param   i    allele index
    \param   val  integer value to set the allele to
    \return  None

    \rst

    Example
    -------

    Set the value of allele ``i`` of string ``p`` in population
    :c:macro:`PGA_NEWPOP` to 64.

    .. code-block:: c

       PGAContext *ctx;
       int p, i;

       ...
       PGASetIntegerAllele (ctx, p, PGA_NEWPOP, i, 64)

    \endrst

******************************************************************************/
void PGASetIntegerAllele (PGAContext *ctx, int p, int pop, int i, int val)
{
    PGAIndividual *ind;
    PGAInteger     *chrom;

    PGADebugEntered("PGASetIntegerAllele");
    PGACheckDataType("PGASetIntegerAllele", PGA_DATATYPE_INTEGER);

    ind = PGAGetIndividual ( ctx, p, pop );
    chrom = (PGAInteger *)ind->chrom;
    chrom[i] = val;

    PGADebugExited("PGASetIntegerAllele");
}

/*!****************************************************************************
    \brief Return the value of allele i of member p in population pop.
    \ingroup allele

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \param   i    allele index
    \return  The value of allele i in string p

    \rst

    Description
    -----------

    Assumes the data type is :c:macro:`PGA_DATATYPE_INTEGER`.

    Example
    -------

    Returns the value of integer allele ``i`` of string ``p`` in
    population :c:macro:`PGA_NEWPOP`.

    .. code-block:: c

       PGAContext *ctx;
       int p, i, k;

       ...
       k =  PGAGetIntegerAllele (ctx, p, PGA_NEWPOP, i)

    \endrst

******************************************************************************/
int PGAGetIntegerAllele (PGAContext *ctx, int p, int pop, int i)
{
    PGAIndividual *ind;
    PGAInteger     *chrom;

    PGADebugEntered("PGAGetIntegerAllele");

    PGACheckDataType("PGAGetIntegerAllele", PGA_DATATYPE_INTEGER);

    ind = PGAGetIndividual ( ctx, p, pop );
    chrom = (PGAInteger *)ind->chrom;

    PGADebugExited("PGAGetIntegerAllele");

    return( (int) chrom[i] );
}

/*!****************************************************************************
    \brief Set a flag to tell the initialization routines to set each
           integer-valued gene to a random permutation of the values given
           by an upper and lower bound.
    \ingroup init

    \param   ctx  context variable
    \param   min  the lower bound of numbers used in the permutation
    \param   max  the upper bound of numbers used in the permutation
    \return  None

    \rst

    Description
    -----------

    The length of the interval must be the same
    as the string length.  This is the default strategy for initializing
    integer-valued strings. The default interval is :math:`[0,L-1]`
    where :math:`L` is the string length.  No string initialization is
    done by this call.

    Example
    -------

    Set the initialization routines to set each gene to a random and
    unique value from the interval [500,599].

    .. code-block:: c

        PGAContext *ctx;

        ...
        PGASetIntegerInitPermute (ctx, 500, 599)}

    \endrst
******************************************************************************/
void PGASetIntegerInitPermute (PGAContext *ctx, int min, int max)
{
    int i, range;

    PGADebugEntered  ("PGASetIntegerInitPermute");
    PGAFailIfSetUp   ("PGASetIntegerInitPermute");
    PGACheckDataType ("PGASetIntegerInitPermute", PGA_DATATYPE_INTEGER);

    range = max - min + 1;
    if (max <= min) {
        PGAError
            ( ctx, "PGASetIntegerInitPermute: max does not exceed min:"
            , PGA_FATAL, PGA_INT, (void *) &max
            );
    } else if (range != ctx->ga.StringLen) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetIntegerInitPermute: range of %d does not equal "
              "StringLen %d"
            , range, ctx->ga.StringLen
            );
    } else {
        ctx->init.IntegerType = PGA_IINIT_PERMUTE;
        for (i = 0; i < ctx->ga.StringLen; i++) {
            ctx->init.IntegerMin [i] = min;
            ctx->init.IntegerMax [i] = max;
        }
    }

    PGADebugExited ("PGASetIntegerInitPermute");
}

/*!****************************************************************************
    \brief Set a flag to tell the initialization routines to set each
           integer-valued gene to a value chosen randomly from the
           interval given by an upper and lower bound.
    \ingroup init
    \param   ctx  context variable
    \param   min  array of lower bounds that define the interval the gene is
                  initialized from
    \param   max  array of upper bounds that define the interval the gene is
                  initialized from
    \return  None

    \rst

    Example
    -------

    Set the initialization routines to select a value for gene ``i``
    uniformly randomly from the interval :math:`[0,i]`.  Assumes all
    strings are of the same length.

    .. code-block:: c

        PGAContext *ctx;
        int *low, *high, stringlen, i;

        ...
        stringlen = PGAGetStringLength (ctx);
        low  = malloc (stringlen * sizeof (int));
        high = malloc (stringlen * sizeof (int));
        for (i=0; i<stringlen; i++) {
            low  [i] = 0;
            high [i] = i;
        }
        PGASetIntegerInitRange (ctx, low, high);

    \endrst

******************************************************************************/
void PGASetIntegerInitRange (PGAContext *ctx, const int *min, const int *max)
{
    int i;

    PGADebugEntered  ("PGASetIntegerInitRange");
    PGAFailIfSetUp   ("PGASetIntegerInitRange");
    PGACheckDataType ("PGASetIntegerInitRange", PGA_DATATYPE_INTEGER);

    for (i = 0; i < ctx->ga.StringLen; i++)
    {
        if (max[i] < min[i]) {
            PGAError
                ( ctx
                , "PGASetIntegerInitRange: Lower bound exceeds upper "
                  "bound for allele #"
                , PGA_FATAL, PGA_INT, (void *) &i
                );
        } else {
            ctx->init.IntegerMin [i] = min [i];
            ctx->init.IntegerMax [i] = max [i];
        }
    }
    ctx->init.IntegerType = PGA_IINIT_RANGE;

    PGADebugExited ("PGASetIntegerInitRange");
}

/*!***************************************************************************
    \brief Return the type of scheme used to randomly initialize strings
           of data type integer.
    \ingroup query
    \param    ctx  context variable
    \return Returns the integer corresponding to the symbolic constant
            used to specify the scheme used to initialize integer strings.

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        int inittype;

        ...
        inittype = PGAGetIntegerInitType (ctx);
        switch (inittype) {
        case PGA_IINIT_PERMUTE:
            printf ("Data Type = PGA_IINIT_PERMUTE\n");
            break;
        case PGA_IINIT_RANGE:
            printf ("Data Type = PGA_IINIT_RANGE\n");
            break;
        }

    \endrst

*****************************************************************************/
int PGAGetIntegerInitType (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetIntegerInitType");
    PGAFailIfNotSetUp ("PGAGetIntegerInitType");
    PGACheckDataType  ("PGAGetIntegerInitType", PGA_DATATYPE_INTEGER);

    PGADebugExited ("PGAGetIntegerInitType");

    return ctx->init.IntegerType;
}

/*!***************************************************************************
    \brief Return the minimum of the range of integers used to randomly
           initialize integer strings.
    \ingroup query

    \param   ctx  context variable
    \param   i    Index of the initialization range
    \return The minimum of the range of integers used to randomly initialize
            integer strings

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int min;
       int idx;

       ...
       min = PGAGetMinIntegerInitValue (ctx, idx);

    \endrst

*****************************************************************************/
int PGAGetMinIntegerInitValue (PGAContext *ctx, int i)
{
    PGADebugEntered   ("PGAGetMinIntegerInitValue");
    PGAFailIfNotSetUp ("PGAGetMinIntegerInitValue");
    PGACheckDataType  ("PGASetIntegerAllele", PGA_DATATYPE_INTEGER);

    if (i < 0 || i >= ctx->ga.StringLen) {
        PGAError
            ( ctx, "PGAGetMinIntegerInitValue: Index out of range:"
            , PGA_FATAL, PGA_INT, (int *) &i
            );
    }

    PGADebugExited ("PGAGetMinIntegerInitValue");

    return ctx->init.IntegerMin [i];
}

/*!***************************************************************************
    \brief Return the maximum of the range of integers used to randomly
           initialize integer strings.
    \ingroup query

    \param   ctx  context variable
    \param   i    Index of the initialization range
    \return  The maximum of the range of integers used to randomly initialize
             integer strings

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int idx;
       int max;

       ...
       max = PGAGetMaxIntegerInitValue (ctx, idx);

    \endrst

*****************************************************************************/
int PGAGetMaxIntegerInitValue (PGAContext *ctx, int i)
{
    PGADebugEntered   ("PGAGetMaxIntegerInitValue");
    PGAFailIfNotSetUp ("PGAGetMaxIntegerInitValue");
    PGACheckDataType  ("PGAGetMaxIntegerInitValue", PGA_DATATYPE_INTEGER);

    if (i < 0 || i >= ctx->ga.StringLen) {
        PGAError
            ( ctx, "PGAGetMaxIntegerInitValue: Index out of range:"
            , PGA_FATAL, PGA_INT, (int *) &i
            );
    }

    PGADebugExited ("PGAGetMaxIntegerInitValue");

    return ctx->init.IntegerMax [i];
}


/*!****************************************************************************
    \brief Allocate memory for a string of type PGAInteger, and
           initializes or clears the string according to initflag.
    \ingroup internal
    \param   ctx       context variable
    \param   p         string index
    \param   pop       symbolic constant of the population string p is in
    \param   initflag  A true/false flag used in conjunction with
                       ctx->ga.RandomInit to initialize the string
                       either randomly or set to zero
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    create string user function for the integer datatype.

    Example
    -------

    Allocates and clears memory and assigns the address of the allocated
    memory to the string field ``ind->chrom`` of the individual.

    .. code-block:: c

       PGAContext *ctx;
       PGAIndividual *ind;

       ...
       PGAIntegerCreateString (ctx, ind, PGA_FALSE);

    \endrst

******************************************************************************/
void PGAIntegerCreateString (PGAContext *ctx, int p, int pop, int
initflag)
{
    int i, fp;
    PGAInteger *c;
    PGAIndividual *new = PGAGetIndividual(ctx, p, pop);

    PGADebugEntered ("PGAIntegerCreateString");

    new->chrom = (void *)malloc(ctx->ga.StringLen * sizeof(PGAInteger));
    if (new->chrom == NULL) {
        PGAError
            ( ctx
            , "PGAIntegerCreateString: No room to allocate new->chrom"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    c = (PGAInteger *)new->chrom;
    if (initflag) {
        if (ctx->fops.InitString) {
            fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
            (*ctx->fops.InitString)(&ctx, &fp, &pop);
        } else {
            (*ctx->cops.InitString)(ctx, p, pop);
        }
    } else {
        for (i=0; i<ctx->ga.StringLen; i++) {
            c[i] = 0;
        }
    }

    PGADebugExited ("PGAIntegerCreateString");
}

/*!****************************************************************************
    \brief Randomly mutates an integer-valued gene with a specified
           probability.
    \ingroup internal

    \param   ctx       context variable
    \param   p         string index
    \param   pop       symbolic constant of the population string p is in
    \param   mr        probability of mutating an integer-valued gene
    \return  Returns the number of mutations

    \rst

    Description
    -----------

    This routine is called from :c:func:`PGAMutate`.

    Note that this function is set in :c:func:`PGASetUp` as the mutation
    user function for the integer datatype by default.

    \endrst

******************************************************************************/
int PGAIntegerMutation (PGAContext *ctx, int p, int pop, double mr)
{
    PGAInteger *c;
    int i, j, temp;
    int count = 0;
    /* The following are used for DE variants */
    int midx = 0;
    int maxidx = 2 * ctx->ga.DENumDiffs + 1;
    DECLARE_DYNARRAY (PGAInteger *, indivs, maxidx);
    int do_crossover = 1;
    static double de_dither = 0.0;

    PGADebugEntered ("PGAIntegerMutation");

    c = (PGAInteger *)PGAGetIndividual (ctx, p, pop)->chrom;

    switch (ctx->ga.MutationType) {
    case PGA_MUTATION_DE:
      {
        DECLARE_DYNARRAY (int, idx, maxidx);
        de_dither = PGASetupDE (ctx, p, pop, maxidx, idx);
        for (i=0; i<maxidx; i++) {
            indivs [i] = (PGAInteger *)
                PGAGetIndividual (ctx, idx [i], PGA_OLDPOP)->chrom;
        }
        /* Index of allele that is mutated in any case */
        midx = PGARandomInterval (ctx, 0, ctx->ga.StringLen - 1);
        break;
      }
    case PGA_MUTATION_SCRAMBLE:
      {
        int l = 0;
        int pos;
        if (PGARandomFlip (ctx, mr)) {
            l   = PGARandomInterval (ctx, 2, ctx->ga.MutateScrambleMax);
            pos = PGARandomInterval (ctx, 0, ctx->ga.StringLen - l - 1);
            PGAShufflePGAInteger (ctx, c + pos, l);
        }
        return l;
        break;
      }
    case PGA_MUTATION_POSITION:
      {
        int l = ctx->ga.StringLen;
        int pos1, pos2;
        if (!PGARandomFlip (ctx, mr)) {
            return 0;
        }
        pos1 = PGARandomInterval (ctx, 0, l - 1);
        pos2 = PGARandomInterval (ctx, 0, l - 1);
        PGAInteger tmp = c [pos1];
        if (pos1 < pos2) {
            count = pos2 - pos1;
            for (i=pos1+1; i<=pos2; i++) {
                c [i-1] = c [i];
            }
        } else {
            count = pos1 - pos2;
            for (i=pos1; i>pos2; i--) {
                c [i] = c [i-1];
            }
        }
        c [pos2] = tmp;
        return count;
        break;
      }
    }

    for (i=0; i<ctx->ga.StringLen; i++) {
        int old_value = c [i];
        int idx = i;
        /* Do not permute fixed edges */
        if (  ctx->ga.n_edges
           && (get_edge (ctx, c [i]) || rev_edge (ctx, c [i]))
           )
        {
            continue;
        }

        /* randomly choose an allele   */
        if (ctx->ga.MutationType == PGA_MUTATION_DE || PGARandomFlip (ctx, mr)){
            /* apply appropriate mutation operator */
            switch (ctx->ga.MutationType) {
            case PGA_MUTATION_CONSTANT:
                /* add or subtract from allele */
                if (PGARandomFlip(ctx, .5)) {
                    c [i] += ctx->ga.MutateIntegerValue;
                } else {
                    c [i] -= ctx->ga.MutateIntegerValue;
                }
                break;
            case PGA_MUTATION_PERMUTE:
              {
                /* could check for j == i if we were noble */
                /* edd: 16 Jun 2007  applying patch from Debian bug
                 * report #333381 correcting an 'off-by-one' here
                 * bu reducing StringLen by 1
                 */
                j = PGARandomInterval (ctx, 0, ctx->ga.StringLen - 1);
                if (ctx->ga.n_edges) {
                    while (get_edge (ctx, c [j]) || rev_edge (ctx, c [j])) {
                        j = PGARandomInterval (ctx, 0, ctx->ga.StringLen - 1);
                    }
                }
                temp  = c [i];
                c [i] = c [j];
                c [j] = temp;
                break;
              }
            case PGA_MUTATION_RANGE:
                c [i] = PGARandomInterval
                    (ctx, ctx->init.IntegerMin [i], ctx->init.IntegerMax [i]);
                break;
            case PGA_MUTATION_POLY:
              {
                double u = PGARandom01 (ctx, 0);
                double eta = PGAGetMutationPolyEta (ctx) + 1;
                double delta, val;
                if (u < 0.5) {
                    delta = pow (2 * u, 1.0 / eta) - 1.0;
                } else {
                    delta = 1.0 - pow (2 * (1 - u), 1.0 / eta);
                }
                if (ctx->ga.MutatePolyValue >= 0) {
                    c [i] += (int)round (delta * ctx->ga.MutatePolyValue);
                } else {
                    if (delta < 0) {
                        val = fabs (c [i] - ctx->init.IntegerMin [i] + 0.4999);
                    } else {
                        val = fabs (ctx->init.IntegerMax [i] - c [i] + 0.4999);
                    }
                    c [i] += (int)round (delta * val);
                }
                break;
              }
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
                        ( ctx, "PGAIntegerMutation: Invalid DE crossover type:"
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
                            c [idx] += (int)round
                                (f * (indivs [j][idx] - indivs [j+1][idx]));
                        }
                        break;
                    case PGA_DE_VARIANT_EITHER_OR:
                        /* We use only 1 difference and ignore DENumDiffs */
                        if (PGARandom01 (ctx, 0) < ctx->ga.DEProbabilityEO) {
                            c [idx] = indivs [0][idx] + (int)round
                                (f * (indivs [1][idx] - indivs [2][idx]));
                        } else {
                            double k = ctx->ga.DEAuxFactor;
                            c [idx] = indivs [0][idx] + (int)round
                                (k * ( indivs [1][idx] + indivs [2][idx]
                                     - 2 * indivs [0][idx]
                                     )
                                );
                        }
                        break;
                    default:
                        PGAError(ctx, "PGAIntegerMutation: Invalid value of "
                                 "ga.DEVariant:", PGA_FATAL, PGA_INT,
                                 (void *) &(ctx->ga.DEVariant));
                        break;
                    }
                    count++;
                }
                break;
            default:
                PGAError
                    ( ctx
                    , "PGAIntegerMutation: Invalid value of ga.MutationType:"
                    , PGA_FATAL, PGA_INT, (void *) &(ctx->ga.MutationType)
                    );
                break;
            }

            /* reset to min/max or bounce if outside range */
            bouncheck
                ( ctx, i, ctx->ga.MutateBoundedFlag, ctx->ga.MutateBounceFlag
                , c, old_value, old_value
                );
            count++;
        }
    }
    PGADebugExited ("PGAIntegerMutation");
    return count;
}

/*!****************************************************************************
    \brief Perform one-point crossover on two parent strings producing
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
    crossover user function for the integer datatype when selecting
    one-point crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``,
    producing children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerOneptCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGAIntegerOneptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
     PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual
        (ctx, p1, pop1)->chrom;
     PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual
        (ctx, p2, pop1)->chrom;
     PGAInteger *child1  = (PGAInteger *)PGAGetIndividual
        (ctx, c1, pop2)->chrom;
     PGAInteger *child2  = (PGAInteger *)PGAGetIndividual
        (ctx, c2, pop2)->chrom;
     int i, xsite;

    PGADebugEntered ("PGAIntegerOneptCrossover");

    xsite = PGARandomInterval (ctx, 1,ctx->ga.StringLen-1);

    for (i=0; i<xsite; i++) {
        child1 [i] = parent1 [i];
        child2 [i] = parent2 [i];
    }

    for (i=xsite; i<ctx->ga.StringLen; i++) {
        child1 [i] = parent2 [i];
        child2 [i] = parent1 [i];
    }

    PGADebugExited ("PGAIntegerOneptCrossover");
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
    crossover user function for the integer datatype when selecting
    two-point crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerTwoptCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGAIntegerTwoptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
     PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual
        (ctx, p1, pop1)->chrom;
     PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual
        (ctx, p2, pop1)->chrom;
     PGAInteger *child1  = (PGAInteger *)PGAGetIndividual
        (ctx, c1, pop2)->chrom;
     PGAInteger *child2  = (PGAInteger *)PGAGetIndividual
        (ctx, c2, pop2)->chrom;
     int i, temp, xsite1, xsite2;

    PGADebugEntered ("PGAIntegerTwoptCrossover");

    /* pick two cross sites such that xsite2 > xsite1 */
    xsite1 = PGARandomInterval (ctx, 1,ctx->ga.StringLen-1);
    xsite2 = xsite1;
    while (xsite2 == xsite1) {
        xsite2 = PGARandomInterval (ctx, 1,ctx->ga.StringLen-1);
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

    PGADebugExited ("PGAIntegerTwoptCrossover");
}


/*!****************************************************************************
    \brief Perform uniform crossover on two parent strings producing two
           children via side-effect.
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
    crossover user function for the integer datatype when selecting
    uniform crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerUniformCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGAIntegerUniformCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual
        (ctx, p1, pop1)->chrom;
    PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual
        (ctx, p2, pop1)->chrom;
    PGAInteger *child1  = (PGAInteger *)PGAGetIndividual
        (ctx, c1, pop2)->chrom;
    PGAInteger *child2  = (PGAInteger *)PGAGetIndividual
        (ctx, c2, pop2)->chrom;
    int i;

    PGADebugEntered ("PGAIntegerUniformCrossover");

    for (i=0; i<ctx->ga.StringLen; i++) {
        if (parent1 [i] == parent2 [i]) {
            child1 [i] = parent1 [i];
            child2 [i] = parent2 [i];
        } else {
            if(PGARandomFlip (ctx, ctx->ga.UniformCrossProb)) {
                child1 [i] = parent1 [i];
                child2 [i] = parent2 [i];
            } else {
                child1 [i] = parent2 [i];
                child2 [i] = parent1 [i];
            }
        }
    }
    PGADebugExited ("PGAIntegerUniformCrossover");
}

/*!****************************************************************************
    \brief Performs simulated binary crossover (SBX) on two parent
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    simulated binary crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerSBXCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGAIntegerSBXCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual
        (ctx, p1, pop1)->chrom;
    PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual
        (ctx, p2, pop1)->chrom;
    PGAInteger *child1  = (PGAInteger *)PGAGetIndividual
        (ctx, c1, pop2)->chrom;
    PGAInteger *child2  = (PGAInteger *)PGAGetIndividual
        (ctx, c2, pop2)->chrom;
    int i;
    double u = 0;

    if (ctx->ga.CrossSBXOnce) {
        u = PGARandom01 (ctx, 0);
    }

    for (i=0; i<ctx->ga.StringLen; i++) {
        if (  parent1[i] == parent2[i]
           || !PGARandomFlip (ctx, ctx->ga.UniformCrossProb)
           )
        {
            child1 [i] = parent1 [i];
            child2 [i] = parent2 [i];
        } else {
            int j;
            double c1, c2;
            PGAInteger minp =
                parent1 [i] < parent2 [i] ? parent1 [i] : parent2 [i];
            PGAInteger maxp =
                parent1 [i] > parent2 [i] ? parent1 [i] : parent2 [i];
            if (!ctx->ga.CrossSBXOnce) {
                u = PGARandom01 (ctx, 0);
            }
            PGACrossoverSBX (ctx, parent1 [i], parent2 [i], u, &c1, &c2);
            child1 [i] = (int)round (c1);
            child2 [i] = (int)round (c2);
            for (j=0; j<2; j++) {
                bouncheck
                    ( ctx, i, ctx->ga.CrossBoundedFlag, ctx->ga.CrossBounceFlag
                    , j ? child2 : child1, minp, maxp
                    );
            }
        }
    }
}

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

static void append_edge (PGAContext *ctx, PGAInteger n1, PGAInteger n2)
{
    int i;
    PGAInteger *em = ctx->scratch.edgemap [n1];
    for (i=0; i<4; i++) {
        if (em [i] == 0) {
            em [i] = n2 + 1;
            return;
        } else if (em [i] == n2 + 1) {
            em [i] = -em [i];
            return;
        } else {
            assert (-em [i] != n2 + 1);
        }
    }
}

static void build_edge_map (PGAContext *ctx, PGAInteger **parent)
{
    int j;
    PGAInteger i, l = ctx->ga.StringLen;
    memset (ctx->scratch.edgemap, 0, sizeof (PGAInteger) * 4 * l);
    unsigned long long s [2] = {0, 0};
    for (i=0; i<l; i++) {
        for (j=0; j<2; j++) {
            PGAInteger n1 = parent [j][i];
            PGAInteger n2 = parent [j][(i + 1) % l];
            if (n1 < 0 || n1 >= l || n2 < 0 || n2 >= l) {
                PGAErrorPrintf
                    (ctx, PGA_FATAL, "Parent gene is no permutation");
            }
            append_edge (ctx, n1, n2);
            append_edge (ctx, n2, n1);
            s [j] += parent [j][i];
        }
    }
    for (j=0; j<2; j++) {
        if (s [j] != (unsigned long long)(l) * (unsigned long long)(l-1) / 2) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Parent gene is no permutation");
        }
    }
}

/* Remove new edge from all right sides, note that only the edges in
 * the edge table for the new edge have this edge on the right side
 */
static void remove_edge_from_right (PGAContext *ctx, PGAInteger cidx)
{
    int i, j;
    for (j=0; j<4; j++) {
        PGAInteger v = labs (ctx->scratch.edgemap [cidx][j]);
        if (v == 0) {
            continue;
        }
        v -= 1;
        for (i=0; i<4; i++) {
            if (labs (ctx->scratch.edgemap [v][i]) - 1 == cidx) {
                ctx->scratch.edgemap [v][i] = 0;
                break;
            }
        }
    }
}

static void fix_edge_map (PGAContext *ctx)
{
    size_t i, j, k;
    PGAFixedEdge *e;
    /* We only need to do something if we have fixed edges */
    if (ctx->ga.n_edges == 0) {
        return;
    }
    for (i=0; i<ctx->ga.n_edges; i++) {
        e = ctx->ga.edges + i;
        if (ctx->ga.symmetric) {
            /* If intermediate edge we don't need to do anything */
            if (e->prev && e->next) {
                continue;
            }
            /* Ensure that fixed edge is always preferred */
            for (j=0; j<2; j++) {
                PGAInteger node  = j == 0 ? e->lhs : e->rhs;
                PGAInteger other = j == 0 ? e->rhs : e->lhs;
                PGAInteger *em = ctx->scratch.edgemap [node];
                /* Ignore warning if compiling without assertions */
                #pragma GCC diagnostic push
                #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
                int found = 0;
                #pragma GCC diagnostic pop
                if (j == 0 && e->prev) {
                    continue;
                }
                if (j == 1 && e->next) {
                    continue;
                }
                for (k=0; k<4; k++) {
                    if (abs (em [k]) - 1 == other) {
                        /* This assertion fails if not both of the
                         * parents have the fixed edge
                         */
                        assert (em [k] < 0);
                        found = 1;
                    } else {
                        em [k] = abs (em [k]);
                    }
                }
                assert (found);
            }
        } else {
            /* Asymmetric case, edge must have given orientation */
            /* right side may never occur in the right side of edge table */
            remove_edge_from_right (ctx, e->rhs);
            /* Entry in edgemap can only have fixed right side */
            ctx->scratch.edgemap [e->lhs][0] = e->rhs;
            ctx->scratch.edgemap [e->lhs][1]
                = ctx->scratch.edgemap [e->lhs][2]
                = ctx->scratch.edgemap [e->lhs][3] = 0;
        }
    }
}

static int count_edges (PGAContext *ctx, PGAInteger idx)
{
    int j;
    int c = 0;
    PGAInteger *em = ctx->scratch.edgemap [idx];
    for (j=0; j<4; j++) {
        if (em [j]) {
            c++;
        }
    }
    return c;
}

static void next_edge (PGAContext *ctx, PGAInteger *child, PGAInteger idx)
{
    PGAInteger i, j;
    PGAInteger *em = ctx->scratch.edgemap [child [idx]];
    assert (idx < ctx->ga.StringLen - 1);
    /* Prefer common edges */
    if (em [0] < 0 && em [1] >= 0) {
        child [idx + 1] = -em [0] - 1;
    } else if (em [1] < 0 && em [0] >= 0) {
        child [idx + 1] = -em [1] - 1;
    } else if (em [2] < 0) {
        child [idx + 1] = -em [2] - 1;
    } else {
        PGAInteger idxm = 0;
        PGAInteger emin [4];
        int minv = -1;
        for (j=0; j<4; j++) {
            int v;
            if (em [j] == 0) {
                continue;
            }
            v = count_edges (ctx, labs (em [j]) - 1);
            if (idxm == 0 || v < minv) {
                minv = v;
                emin [0] = labs (em [j]) - 1;
                idxm = 1;
            } else if (v == minv) {
                emin [idxm++] = labs (em [j]) - 1;
            }
        }
        if (idxm == 1) {
            child [idx + 1] = emin [0];
        } else if (idxm > 1) {
            child [idx + 1] = emin [PGARandomInterval (ctx, 0, idxm - 1)];
        } else {
            PGAInteger used [idx + 1];
            PGAInteger mini = -1;
            PGAInteger lastu = 0;
            memcpy (used, child, sizeof (PGAInteger) * (idx + 1));
            qsort (used, idx + 1, sizeof (PGAInteger), intcmp);
            for (j=0; j<idx+1; j++) {
                for (i=lastu; i<used [j]; i++) {
                    int ec = count_edges (ctx, i);
                    PGAFixedEdge *p = get_edge (ctx, i);
                    /* May not start in the middle of a run of fixed edges */
                    if (p && p->prev) {
                        continue;
                    }
                    /* If we have asymmetric (directed) fixed edges and
                     * this is the last of a run (intermediate cases
                     * above) we may not start with that edge
                     */
                    if (!ctx->ga.symmetric && NULL != rev_edge (ctx, i)) {
                        continue;
                    }
                    if (minv < 0 || ec < minv) {
                        minv = ec;
                        mini = i;
                    }
                }
                lastu = used [j] + 1;
            }
            /* don't randomize different indexes with same count */
            if (mini < 0) {
                /* The current (sorted) list is tight
                 * Or did contain only invalid starting points
                 * (right sides of fixed edges)
                 */
                for (j=idx+1; j<ctx->ga.StringLen; j++) {
                    PGAFixedEdge *p = get_edge (ctx, j);
                    /* Don't use intermediate fixed edge */
                    if (p && p->prev) {
                        continue;
                    }
                    if (ctx->ga.symmetric || NULL == rev_edge (ctx, j)) {
                        child [idx + 1] = j;
                        break;
                    }
                }
            } else {
                child [idx + 1] = mini;
            }
        }
    }
    remove_edge_from_right (ctx, child [idx + 1]);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Perform Edge Recombination on two parent strings producing
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    edge crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerEdgeCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerEdgeCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    int j, oc;
    PGAInteger i, ci;
    PGAInteger l = ctx->ga.StringLen;
    int p_ok [] = {1, 1};
    PGAInteger ok [2];

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    /* Build Edge-Map, 0 is unallocated, a node is represented by index+1 */
    build_edge_map (ctx, parent);
    fix_edge_map   (ctx);
    /* Find a node with edge-count < 4 for each child
     * If not possible (e.g. for the 2nd child) we use first node of one
     * of the parents.
     */
    ci = 0;
    oc = 0;
    for (i=0; i<l; i++) {
        PGAFixedEdge *p = get_edge (ctx, i);
        /* May not start with intermediate edge of fixed-edge run */
        if (p && p->prev) {
            continue;
        }
        /* May not start with end node if asymmetric fixed edges */
        if (!ctx->ga.symmetric && rev_edge (ctx, i)) {
            continue;
        }
        for (j=0; j<4; j++) {
            if (!ctx->scratch.edgemap [i][j]) {
                child [ci++][0] = i;
                break;
            }
        }
        if (oc < 2 && j == 4) {
            ok [oc++] = i;
        }
        if (ci >= 2) {
            break;
        }
    }
    if (ci < 2) {
        for (i=0; i<2; i++) {
            PGAFixedEdge *p = get_edge (ctx, parent [i][0]);
            /* May not start with intermediate edge of fixed-edge run */
            if (p && p->prev) {
                p_ok [i] = 0;
                continue;
            }
            /* May not start with end node if asymmetric fixed edges */
            if (!ctx->ga.symmetric && rev_edge (ctx, i)) {
                p_ok [i] = 0;
                continue;
            }
        }
        assert (oc == 2);
        if (ci == 0) {
            if (p_ok [0]) {
                child [0][0] = parent [0][0];
            } else {
                child [0][0] = ok [0];
            }
            if (p_ok [1]) {
                child [1][0] = parent [1][0];
            } else {
                child [1][0] = ok [1];
            }
        } else {
            if (child [0][0] == parent [0][0] && p_ok [1]) {
                child [1][0] = parent [1][0];
            } else if (child [0][0] == parent [1][0] && p_ok [0]) {
                child [1][0] = parent [0][0];
            } else if (!(p_ok [0] && p_ok [1])) {
                child [1][0] = ok [0];
            } else {
                child [1][0] = parent [PGARandomFlip (ctx, 0.5)][0];
            }
        }
    }
    remove_edge_from_right (ctx, child [0][0]);
    /* Now do the actual crossover */
    for (ci=0; ci<2; ci++) {
        for (i=0; i<l-1; i++) {
            next_edge (ctx, child [ci], i);
        }
        #ifdef DEBUG
        assert_has_edges (ctx, child [ci]);
        #endif /* DEBUG */

        /* Need to rebuild, edge-map is consumed above */
        build_edge_map (ctx, parent);
        fix_edge_map   (ctx);
        remove_edge_from_right (ctx, child [1][0]);
    }
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

/*!****************************************************************************
    \brief Perform Partially Mapped Crossover on two parent strings producing
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerPartiallyMappedCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerPartiallyMappedCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger i;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger pos1, pos2;
    PGAInteger *c0_seen = ctx->scratch.pgaintscratch [0];
    PGAInteger *c1_seen = ctx->scratch.pgaintscratch [1];
    PGAInteger imin = ctx->init.IntegerMin [0];
    for (i=0; i<l; i++) {
        c0_seen [i] = c1_seen [i] = -1;
    }

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    /* Chose two positions */
    pos1 = PGARandomInterval (ctx, 0, l - 1);
    pos2 = PGARandomInterval (ctx, 0, l - 1);
    for (i=pos1; i != (pos2 + 1) % l; i=(i+1) % l) {
        if (parent [1][i] - imin < 0 || parent [1][i] - imin >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: Second Parent is no permutation");
        }
        child [0][i] = parent [1][i];
        c0_seen [parent [1][i] - imin] = parent [0][i];
        if (parent [0][i] - imin < 0 || parent [0][i] - imin >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: First Parent is no permutation");
        }
        child [1][i] = parent [0][i];
        c1_seen [parent [0][i] - imin] = parent [1][i];
    }
    for (i = (pos2 + 1) % l; i != pos1; i=(i+1) % l) {
        int j;
        PGAInteger v;
        v = parent [0][i];
        for (j=0; j<l + 5 && c0_seen [v - imin] >= 0; j++) {
            v = c0_seen [v - imin];
        }
        /* Should have been checked above, possibly a bug */
        if (i > l) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Crossover: no permutation");
        }
        child [0][i] = v;
        v = parent [1][i];
        for (j=0; j<l + 5 && c1_seen [v - imin] >= 0; j++) {
            v = c1_seen [v - imin];
        }
        /* Should have been checked above, possibly a bug */
        if (i > l) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Crossover: no permutation");
        }
        child [1][i] = v;
    }
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

static void copy_middle_part
    ( PGAContext *ctx
    , PGAInteger *parent, PGAInteger *child
    , PGAInteger flag, PGAInteger pos1, PGAInteger pos2
    )
{
    PGAInteger i;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger *seen = ctx->scratch.pgaintscratch [0];
    PGAInteger imin = ctx->init.IntegerMin [0];

    for (i=pos1; i != pos2; i=(i+1) % l) {
        child [i] = parent [i];
        if (parent [i] - imin < 0 || parent [i] - imin >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: Parent is no permutation");
        }
        seen [parent [i] - imin] |= flag;
    }
}

static PGAInteger copy_rest
    ( PGAContext *ctx
    , PGAInteger *parent, PGAInteger *child
    , PGAInteger flag, PGAInteger pidx, PGAInteger sc, PGAInteger ec
    )
{
    PGAInteger i;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger *seen = ctx->scratch.pgaintscratch [0];
    PGAInteger imin = ctx->init.IntegerMin [0];
    /* Copy rest from other parent, retaining order in other parent */
    for (i=sc; i<ec; i++) {
        while (seen [parent [pidx] - imin] & flag) {
            pidx = (pidx + 1) % l;
        }
        child [i] = parent [pidx];
        pidx = (pidx + 1) % l;
    }
    return pidx;
}

/*!****************************************************************************
    \brief Perform Modified Crossover on two parent strings producing
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerModifiedCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerModifiedCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger pos;
    PGAInteger *seen = ctx->scratch.pgaintscratch [0];
    memset (seen, 0, sizeof (PGAInteger) * l);

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);


    /* Chose one position and copy first part */
    pos = PGARandomInterval (ctx, 0, l - 1);
    copy_middle_part (ctx, parent [0], child [0], 1, 0, pos);
    copy_middle_part (ctx, parent [1], child [1], 2, 0, pos);
    /* Copy rest from other parent, retaining order in other parent */
    (void)copy_rest (ctx, parent [1], child [0], 1, 0, pos, l);
    (void)copy_rest (ctx, parent [0], child [1], 2, 0, pos, l);
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

/*!****************************************************************************
    \brief Perform Order Crossover on two parent strings producing
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerOrderCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerOrderCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger j;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger pos1, pos2;
    PGAInteger *seen = ctx->scratch.pgaintscratch [0];
    memset (seen, 0, sizeof (PGAInteger) * l);

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    /* Chose two positions and copy middle part, note that 'middle' part
     * can wrap
     */
    pos1 = PGARandomInterval (ctx, 0, l - 1);
    pos2 = PGARandomInterval (ctx, 0, l - 1);
    copy_middle_part (ctx, parent [0], child [0], 1, pos1, pos2);
    copy_middle_part (ctx, parent [1], child [1], 2, pos1, pos2);
    /* Copy rest from other parent, retaining order in other parent */
    if (pos2 >= pos1) {
        j =   copy_rest (ctx, parent [1], child [0], 1, pos2, pos2, l);
        (void)copy_rest (ctx, parent [1], child [0], 1, j,    0, pos1);
        j =   copy_rest (ctx, parent [0], child [1], 2, pos2, pos2, l);
        (void)copy_rest (ctx, parent [0], child [1], 2, j,    0, pos1);
    } else {
        (void)copy_rest (ctx, parent [1], child [0], 1, pos2, pos2, pos1);
        (void)copy_rest (ctx, parent [0], child [1], 2, pos2, pos2, pos1);
    }
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

/*!****************************************************************************
    \brief Perform Cycle Crossover on two parent strings producing
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerCycleCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerCycleCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger i;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger pos;
    PGAInteger imin = ctx->init.IntegerMin [0];
    PGAInteger *seen = ctx->scratch.pgaintscratch [0];
    PGAInteger *pidx = ctx->scratch.pgaintscratch [1];
    memset (seen, 0, sizeof (PGAInteger) * l);

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    /* Build parent index */
    for (i=0; i<l; i++) {
        if (parent [0][i] - imin < 0 || parent [0][i] - imin >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: First Parent is no permutation");
        }
        pidx [parent [0][i] - imin] = i;
    }
    pos = PGARandomInterval (ctx, 0, l - 1);
    for (i=pos; !seen [i]; i=pidx [parent [1][i] - imin]) {
        if (i < 0 || i >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: Second Parent is no permutation");
        }
        seen [i] = 1;
        child [0][i] = parent [0][i];
        child [1][i] = parent [1][i];
    }
    for (i=0; i<l; i++) {
        if (seen [i]) {
            continue;
        }
        child [0][i] = parent [1][i];
        child [1][i] = parent [0][i];
    }
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

/*!****************************************************************************
    \brief Perform Order Based Crossover on two parent strings producing
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerOrderBasedCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerOrderBasedCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger i, j0, j1;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger imin = ctx->init.IntegerMin [0];
    PGAInteger *idx0 = ctx->scratch.pgaintscratch [0];
    PGAInteger *idx1 = ctx->scratch.pgaintscratch [1];
    PGAInteger *val0 = ctx->scratch.pgaintscratch [2];
    PGAInteger *val1 = ctx->scratch.pgaintscratch [3];

    for (i=0; i<l; i++) {
        idx0 [i] = idx1 [i] = val0 [i] = val1 [i] = -1;
    }

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    j0 = 0;
    for (i=0; i<l; i++) {
        if (PGARandomFlip (ctx, 0.5)) {
            if (parent [0][i] - imin < 0 || parent [0][i] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: First Parent is no permutation"
                    );
            }
            idx0 [parent [0][i] - imin] = j0;
            val0 [j0] = parent [0][i];
            if (parent [1][i] - imin < 0 || parent [1][i] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: Second Parent is no permutation"
                    );
            }
            idx1 [parent [1][i] - imin] = j0;
            val1 [j0] = parent [1][i];
            j0++;
        }
    }
    j0 = j1 = 0;
    for (i=0; i<l; i++) {
        if (parent [0][i] - imin < 0 || parent [0][i] - imin >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: First Parent is no permutation");
        }
        if (parent [1][i] - imin < 0 || parent [1][i] - imin >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: Second Parent is no permutation");
        }
        if (idx1 [parent [0][i] - imin] < 0) {
            child [0][i] = parent [0][i];
        } else {
            child [0][i] = val1 [j0++];
        }
        if (idx0 [parent [1][i] - imin] < 0) {
            child [1][i] = parent [1][i];
        } else {
            child [1][i] = val0 [j1++];
        }
    }
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

/*!****************************************************************************
    \brief Perform Position Based Crossover on two parent strings producing
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerPositionBasedCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerPositionBasedCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger *seen0 = ctx->scratch.pgaintscratch [0];
    PGAInteger *seen1 = ctx->scratch.pgaintscratch [1];
    PGAInteger i, j0, j1;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger imin = ctx->init.IntegerMin [0];
    memset (seen0, 0, sizeof (PGAInteger) * l);
    memset (seen1, 0, sizeof (PGAInteger) * l);

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    for (i=0; i<l; i++) {
        ctx->scratch.intscratch [i] = 0;
        if (PGARandomFlip (ctx, 0.5)) {
            child [0][i] = parent [0][i];
            child [1][i] = parent [1][i];
            if (parent [0][i] - imin < 0 || parent [0][i] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: First Parent is no permutation"
                    );
            }
            if (parent [1][i] - imin < 0 || parent [1][i] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: Second Parent is no permutation"
                    );
            }
            assert (!seen0 [child [0][i] - imin]);
            assert (!seen1 [child [1][i] - imin]);
            seen0 [child [0][i] - imin] = 1;
            seen1 [child [1][i] - imin] = 1;
            ctx->scratch.intscratch [i] = 1;
        }
    }
    j0 = j1 = 0;
    for (i=0; i<l; i++) {
        if (!ctx->scratch.intscratch [i]) {
            for ( ; seen0 [parent [1][j0] - imin]; j0++)
                ;
            for ( ; seen1 [parent [0][j1] - imin]; j1++)
                ;
            assert (j0 < l && j1 < l);
            child [0][i] = parent [1][j0];
            child [1][i] = parent [0][j1];
            j0++;
            j1++;
        }
    }
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

/*!****************************************************************************
    \brief Perform Uniform Order Based Crossover on two parent strings
           producing two children via side-effect.
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Note that this crossover variant is identical to
    PGAIntegerPositionBasedCrossover for the first child.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerUniformOrderBasedCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerUniformOrderBasedCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger *seen0 = ctx->scratch.pgaintscratch [0];
    PGAInteger *seen1 = ctx->scratch.pgaintscratch [1];
    PGAInteger i, j0, j1;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger imin = ctx->init.IntegerMin [0];
    memset (seen0, 0, sizeof (PGAInteger) * l);
    memset (seen1, 0, sizeof (PGAInteger) * l);

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    for (i=0; i<l; i++) {
        ctx->scratch.intscratch [i] = 0;
        if (PGARandomFlip (ctx, 0.5)) {
            child [0][i] = parent [0][i];
            if (child [0][i] - imin < 0 || child [0][i] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: First Parent is no permutation"
                    );
            }
            assert (!seen0 [child [0][i] - imin]);
            seen0 [child [0][i] - imin] = 1;
            ctx->scratch.intscratch [i] = 1;
        } else {
            child [1][i] = parent [1][i];
            if (child [1][i] - imin < 0 || child [1][i] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: Second Parent is no permutation"
                    );
            }
            assert (!seen1 [child [1][i] - imin]);
            seen1 [child [1][i] - imin] = 1;
        }
    }
    j0 = j1 = 0;
    for (i=0; i<l; i++) {
        if (!ctx->scratch.intscratch [i]) {
            for ( ; seen0 [parent [1][j0] - imin]; j0++)
                ;
            assert (j0 < l);
            child [0][i] = parent [1][j0];
            j0++;
        } else {
            for ( ; seen1 [parent [0][j1] - imin]; j1++)
                ;
            assert (j1 < l);
            child [1][i] = parent [0][j1];
            j1++;
        }
    }
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

static void ae_fill_child
    ( PGAContext *ctx
    , PGAInteger *parent [2], PGAInteger *child
    , PGAInteger pidx, PGAInteger off
    )
{
    PGAInteger *idx [2];
    PGAInteger *val  = ctx->scratch.pgaintscratch [2];
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger i;
    PGAInteger imin = ctx->init.IntegerMin [0];
    idx [0] = ctx->scratch.pgaintscratch [0];
    idx [1] = ctx->scratch.pgaintscratch [1];
    memset (val, 0, sizeof (PGAInteger) * l);

    child [off] = parent [pidx][off];
    if (child [off] - imin < 0 || child [off] - imin >= l) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "Crossover: Parent %d is no permutation"
            , (pidx) + 1
            );
    }
    val [child [off] - imin] = 1;
    off = (off + 1) % l;
    child [off] = parent [pidx][off];
    if (child [off] - imin < 0 || child [off] - imin >= l) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "Crossover: Parent %d is no permutation"
            , (pidx) + 1
            );
    }
    val [child [off] - imin] = 1;
    off = (off + 1) % l;
    for (i=2; i<l; i++) {
        PGAInteger v = child [(off + l - 1) % l];
        PGAInteger pos = idx [!pidx][v - imin];
        if (!val [parent [!pidx][(pos + 1) % l] - imin]) {
            child [off] = parent [!pidx][(pos + 1) % l];
            if (child [off] - imin < 0 || child [off] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: Parent %d is no permutation"
                    , (!pidx) + 1
                    );
            }
            assert (val [child [off] - imin] == 0);
            val [child [off] - imin] = 1;
            pidx = !pidx;
            off  = (off + 1) % l;
            continue;
        }
        if (!val [parent [!pidx][(pos + l - 1) % l] - imin]) {
            child [off] = parent [!pidx][(pos + l - 1) % l];
            if (child [off] - imin < 0 || child [off] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: Parent %d is no permutation"
                    , (!pidx) + 1
                    );
            }
            assert (val [child [off] - imin] == 0);
            val [child [off] - imin] = 1;
            pidx = !pidx;
            off  = (off + 1) % l;
            continue;
        }
        /* If not found in other parent try current one */
        pos = idx [pidx][v - imin];
        if (!val [parent [pidx][(pos + 1) % l] - imin]) {
            child [off] = parent [pidx][(pos + 1) % l];
            if (child [off] - imin < 0 || child [off] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: Parent %d is no permutation"
                    , pidx + 1
                    );
            }
            assert (val [child [off] - imin] == 0);
            val [child [off] - imin] = 1;
            off  = (off + 1) % l;
            continue;
        }
        if (!val [parent [pidx][(pos + l - 1) % l] - imin]) {
            child [off] = parent [pidx][(pos + l - 1) % l];
            if (child [off] - imin < 0 || child [off] - imin >= l) {
                PGAErrorPrintf
                    ( ctx, PGA_FATAL
                    , "Crossover: Parent %d is no permutation"
                    , pidx + 1
                    );
            }
            assert (val [child [off] - imin] == 0);
            val [child [off] - imin] = 1;
            off  = (off + 1) % l;
            continue;
        }
        /* Select rand index and search from there for free position */
        pos = PGARandomInterval (ctx, 0, l - 1);
        for (; val [pos]; pos = (pos + 1) % l)
            ;
        child [off] = pos + imin;
        val [pos] = 1;
        off  = (off + 1) % l;
    }
}

/*!****************************************************************************
    \brief Perform Alternating Edge Crossover on two parent strings producing
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Note that the original paper [GGRG85]_ mandates that the search
    starts with a random position. We start with a random position but
    keep the absolute position in the child (and don't copy from the
    middle of the parent to the start of the child).
    This may make the crossover work for other problems than just TSP.
    From the paper it is unclear if edge reversals are possible, we
    allow them but prefer non-reversed egdes.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerAlternatingEdgeCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerAlternatingEdgeCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger *idx [2];
    PGAInteger i, off;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger imin = ctx->init.IntegerMin [0];
    idx [0] = ctx->scratch.pgaintscratch [0];
    idx [1] = ctx->scratch.pgaintscratch [1];
    for (i=0; i<l; i++) {
        idx [0][i] = idx[1][i] = -1;
    }

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    for (i=0; i<l; i++) {
        if (parent [0][i] - imin < 0 || parent [0][i] - imin >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: First Parent is no permutation");
        }
        if (parent [1][i] - imin < 0 || parent [1][i] - imin >= l) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Crossover: Second Parent is no permutation");
        }
        idx [0][parent [0][i] - imin] = i;
        idx [1][parent [1][i] - imin] = i;
    }
    off = PGARandomInterval (ctx, 0, l - 1);
    ae_fill_child (ctx, parent, child [0], 0, off);
    ae_fill_child (ctx, parent, child [1], 1, off);
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

/*!****************************************************************************
    \brief Perform Non-wrapping Order Crossover on two parent strings
           producing two children via side-effect.
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
    \return  c1 and c2 in population pop2 are modified by side-effect.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the integer datatype when selecting
    partially mapped crossover.

    The operation produces permutations of the integer genes of both
    parents. The result is a permutation again for both children.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGAIntegerNonWrappingOrderCrossover
           (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/

void PGAIntegerNonWrappingOrderCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    PGAInteger j;
    PGAInteger l = ctx->ga.StringLen;
    PGAInteger pos1, pos2;
    PGAInteger *seen = ctx->scratch.pgaintscratch [0];
    memset (seen, 0, sizeof (PGAInteger) * l);

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    DEBUG_CHECK_PERMUTE (ctx, parent [0]);
    DEBUG_CHECK_PERMUTE (ctx, parent [1]);

    /* Chose two positions and copy middle part, note that 'middle' part
     * can wrap
     */
    pos1 = PGARandomInterval (ctx, 0, l - 1);
    pos2 = PGARandomInterval (ctx, 0, l - 1);
    copy_middle_part (ctx, parent [0], child [0], 1, pos1, pos2);
    copy_middle_part (ctx, parent [1], child [1], 2, pos1, pos2);
    /* Copy rest from other parent, retaining order in other parent */
    if (pos2 >= pos1) {
        j =   copy_rest (ctx, parent [1], child [0], 1, 0, 0, pos1);
        (void)copy_rest (ctx, parent [1], child [0], 1, j, pos2, l);
        j =   copy_rest (ctx, parent [0], child [1], 2, 0, 0, pos1);
        (void)copy_rest (ctx, parent [0], child [1], 2, j, pos2, l);
    } else {
        (void)copy_rest (ctx, parent [1], child [0], 1, 0, pos2, pos1);
        (void)copy_rest (ctx, parent [0], child [1], 2, 0, pos2, pos1);
    }
    DEBUG_CHECK_PERMUTE (ctx, child [0]);
    DEBUG_CHECK_PERMUTE (ctx, child [1]);
}

/*!****************************************************************************
    \brief Set edges that have to be present.
    \ingroup standard-api

    \param   ctx        context variable
    \param   n          Number of edges
    \param   edge       Pointer to edges, each edge consists of two indices
    \param   symmetric  Flag that indicates if edges are allowed in reverse
                        direction
    \return  None

    \rst

    Description
    -----------

    This is used only in Edge Crossover.
    Note: The edges data structure passed into the function is copied
    and ownership remains with the caller.
    It is admissible that the edges data is in automatic
    variables allocated on the stack.

    Example
    -------

    Set edges (0, 213) and (7, 11) as fixed edges.

    .. code-block:: c

       PGAContext *ctx;
       PGAInteger edge[][2] = {(PGAInteger []){0, 213}, (PGAInteger []){7, 11}};
       int  n = sizeof (edge) / (2 * sizeof (PGAInteger));

       ...
       PGAIntegerSetFixedEdges (ctx, n, edge, PGA_TRUE);

    \endrst

******************************************************************************/
void PGAIntegerSetFixedEdges
    (PGAContext *ctx, size_t n, PGAInteger (*edge)[2], int symmetric)
{
    size_t i, j;
    PGAInteger prev;
    /* Do not allocate twice */
    assert (ctx->ga.edges == NULL);
    ctx->ga.edges = malloc (sizeof (PGAFixedEdge) * n);
    if (ctx->ga.edges == NULL) {
        PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate edges");
    }
    ctx->ga.r_edge = malloc (2 * sizeof (PGAInteger) * n);
    if (ctx->ga.r_edge == NULL) {
        PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate rev-edges");
    }
    for (i=0; i<n; i++) {
        ctx->ga.edges  [i].lhs  = edge [i][0];
        ctx->ga.edges  [i].rhs  = edge [i][1];
        ctx->ga.edges  [i].next = ctx->ga.edges [i].prev = NULL;
    }
    qsort (ctx->ga.edges,  n, sizeof (PGAFixedEdge),   intcmp);
    for (i=0; i<n; i++) {
        ctx->ga.r_edge [i][0] = ctx->ga.edges [i].rhs;
        ctx->ga.r_edge [i][1] = i;
    }
    qsort (ctx->ga.r_edge, n, 2 * sizeof (PGAInteger), intcmp);
    ctx->ga.n_edges = n;
    ctx->ga.symmetric = symmetric ? PGA_TRUE : PGA_FALSE;
    /* predecessors, cycle check */
    prev = -1;
    for (i=0; i<n; i++) {
        PGAFixedEdge *p;
        if (ctx->ga.edges [i].lhs == prev) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Fixed edges may not share nodes");
        }
        prev = ctx->ga.edges [i].lhs;
        /* Try finding predecessor */
        p = get_edge (ctx, ctx->ga.edges [i].rhs);
        if (p != NULL) {
            if (p->prev) {
                PGAErrorPrintf
                    (ctx, PGA_FATAL, "Fixed edges may not share nodes");
            } else {
                p->prev = ctx->ga.edges + i;
            }
        }
    }
    /* Now loop again over all nodes and perform cycle check:
     * We follow the predecessor indexes until we encounter a stop (-1)
     * or we have exhausted n.
     */
    for (i=0; i<n; i++) {
        if (ctx->ga.edges [i].prev) {
            PGAFixedEdge *p = ctx->ga.edges + i;
            for (j=0; j<n; j++) {
                if (p->prev) {
                    PGAFixedEdge *np = p->prev;
                    if (np->next) {
                        assert (np->next == p);
                    } else {
                        np->next = p;
                    }
                    p = np;
                } else {
                    break;
                }
            }
            if (j == n) {
                PGAErrorPrintf (ctx, PGA_FATAL, "Fixed edges contain a cycle");
            }
        }
    }
}

/*!****************************************************************************
    \brief Write an integer-valued string to a file.
    \ingroup internal

    \param   ctx  context variable
    \param   fp   file pointer to file to write the string to
    \param   p    index of the string to write out
    \param   pop  symbolic constant of the population string p is in
    \return  None

    \rst

    Example
    -------

    Write member ``p`` in population :c:macro:`PGA_NEWPOP` to ``stdout``.

    .. code-block:: c

       PGAContext *ctx;
       int  p;

       ...
       PGAIntegerPrintString (ctx, stdout, p, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGAIntegerPrintString ( PGAContext *ctx, FILE *fp, int p, int pop)
{
    PGAInteger *c = (PGAInteger *)PGAGetIndividual (ctx, p, pop)->chrom;
    int i;

    PGADebugEntered ("PGAIntegerPrintString");

    for (i = 0; i < ctx->ga.StringLen; i++) {
        switch (i % 6) {
        case 0:
            fprintf (fp, "#%5d: [%8ld]",i,c[i]);
            break;
        case 1:
        case 2:
        case 3:
        case 4:
            fprintf (fp, ", [%8ld]",c[i]);
            break;
        case 5:
            fprintf (fp, ", [%8ld]",c[i]);
            if (i+1 < ctx->ga.StringLen) {
                fprintf ( fp, "\n");
            }
            break;
        }
    }
    fprintf (fp, "\n");

    PGADebugExited ("PGAIntegerPrintString");
}

/*!****************************************************************************
    \brief Copy one integer-valued string to another.
    \ingroup internal

    \param   ctx   context variable
    \param   p1    string to copy
    \param   pop1  symbolic constant of population containing string p1
    \param   p2    string to copy p1 to
    \param   pop2  symbolic constant of population containing string p2
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the copy
    string user function for the integer datatype by default.

    \endrst

******************************************************************************/
void PGAIntegerCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAInteger *source = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAInteger *dest   = (PGAInteger *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int i;

    PGADebugEntered ("PGAIntegerCopyString");

    for (i = 0; i < ctx->ga.StringLen; i++) {
        dest [i] = source [i];
    }

    PGADebugExited ("PGAIntegerCopyString");
}

/*!****************************************************************************
    \brief Return true if string a is a duplicate of string b, else
           returns false.
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

    Note that this function is set in :c:func:`PGASetUp` as the duplicate
    checking user function for the integer datatype by default.

    \endrst

******************************************************************************/
int PGAIntegerDuplicate (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAInteger *a = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAInteger *b = (PGAInteger *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int count = 0;

    PGADebugEntered ("PGAIntegerDuplicate");

    for (count=0; count<ctx->ga.StringLen; count++) {
        if (a [count] != b [count]) {
            break;
        }
    }

    PGADebugExited ("PGAIntegerDuplicate");

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
    hash user function for the integer datatype by default.

    \endrst

******************************************************************************/
PGAHash PGAIntegerHash (PGAContext *ctx, int p, int pop)
{
    void *a = PGAGetIndividual (ctx, p, pop)->chrom;
    PGAHash hash = PGAUtilHash
        (a, sizeof (PGAInteger) * ctx->ga.StringLen, PGA_INITIAL_HASH);
    return hash;
}

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* Helper function for computing index into gene for given value if it
 * is part of a fixed edge
 */
static void compute_idx
    (PGAContext *ctx, size_t (*x)[2], size_t k, PGAInteger a)
{
    PGAInteger (*r)[2];
    PGAFixedEdge *p;

    p = get_edge (ctx, a);
    if (p != NULL) {
        size_t off = p - ctx->ga.edges;
        x [off][0] = k;
        if (p->prev) {
            assert (p->prev->rhs == a);
            off = p->prev - ctx->ga.edges;
            x [off][1] = k;
        }
    } else if (NULL != (r = rev_edge (ctx, a))) {
        size_t off = (*r) [1];
        assert ((ctx->ga.edges + off)->rhs == a);
        x [off][1] = k;
    }
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/*!****************************************************************************
    \brief Randomly initialize a string of data type integer.
    \ingroup internal

    \param   ctx  context variable
    \param   p    index of string to randomly initialize
    \param   pop  symbolic constant of the population string p is in
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    init string user function for the integer datatype by default.

    \endrst

******************************************************************************/

void PGAIntegerInitString (PGAContext *ctx, int p, int pop)
{
    int len = ctx->ga.StringLen;
    PGAInteger *c = (PGAInteger *)PGAGetIndividual (ctx, p, pop)->chrom;

    PGADebugEntered ("PGAIntegerInitString");

    switch (ctx->init.IntegerType) {
    case PGA_IINIT_PERMUTE:
      {
        PGAInteger ii;
        /* We need to initialize the array in reverse order to emulate
         * the old shuffling below
         */
        for (ii=0; ii<len; ii++) {
             c [len - ii - 1] = ii + ctx->init.IntegerMin [0];
        }
        /* We should really use PGAShufflePGAInteger here but this would
         * destroy all test cases that rely on random number sequence
         * and initial population. So we emulate the original
         * initialization but with only a single tmp variable not with a
         * whole array (which was allocated dynamically).
         */
        for (ii=0; ii<len; ii++) {
            int j = PGARandomInterval (ctx, 0, len - ii - 1);
            PGAInteger tmp = c [ii];
            c [ii] = c [len - j - 1];
            c [len - j - 1] = tmp;
        }
        /* Ensure fixed edges if configured, all fixed edges are used in
         * forward direction regardless of the symmetric flag
         */
        if (ctx->ga.n_edges) {
            size_t (*x)[2] = malloc (2 * ctx->ga.n_edges * sizeof (size_t));
            size_t k, j = 0;
            if (x == NULL) {
                PGAErrorPrintf
                    (ctx, PGA_FATAL, "PGAIntegerInitString: Cannot allocate");
            }
            for (k=0; k<ctx->ga.n_edges; k++) {
                x [k][0] = x [k][1] = SIZE_MAX;
            }
            for (k=0; k<(size_t)len; k++) {
                compute_idx (ctx, x, k, c [k]);
            }
            for (k=0; k<ctx->ga.n_edges; k++) {
                assert (x [k][0] != SIZE_MAX && x [k][1] != SIZE_MAX);
            }
            for (k=0; k<ctx->ga.n_edges; k++) {
                PGAInteger v;
                PGAFixedEdge *p = ctx->ga.edges + k;
                size_t ix;
                if (p->prev) {
                    continue;
                }
                /* We start at 0, otherwise sequences may overlap and
                 * destroy each other
                 */
                ix = k;
                /* Relocate first item of sequence */
                v = c [j];
                c [j] = c [x [ix][0]];
                c [x [ix][0]] = v;
                compute_idx (ctx, x, x [ix][0], v);
                j++;
                do {
                    size_t tmp;
                    v = c [j];
                    c [j] = c [x [ix][1]];
                    c [x [ix][1]] = v;
                    compute_idx (ctx, x, x [ix][1], v);
                    j++;
                    p = p->next;
                    if (p) {
                        tmp = p - ctx->ga.edges;
                        assert  (x [tmp][0] == x [ix][1]);
                        x [tmp][0] = j;
                        ix = tmp;
                    }
                } while (p);
            }
            free (x);
            #ifdef DEBUG
            assert_has_edges (ctx, c);
            # endif /* DEBUG */
        }
        break;
      }
    case PGA_IINIT_RANGE:
      {
        int i;
        for (i=0; i<len; i++)
            c [i] = PGARandomInterval
                (ctx, ctx->init.IntegerMin [i], ctx->init.IntegerMax [i]);
        break;
      }
    }

    PGADebugExited ("PGAIntegerInitString");
}

/*!****************************************************************************
    \brief Build an MPI datatype for a string of data type integer.
    \ingroup internal

    \param    ctx  context variable
    \param    p    index of string
    \param    pop  symbolic constant of the population string p is in
    \return   MPI_Datatype

    \rst

    Description
    -----------

    Called only by MPI routines.  Not for user consumption.

    \endrst

******************************************************************************/
MPI_Datatype PGAIntegerBuildDatatype (PGAContext *ctx, int p, int pop)
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

    PGADebugEntered ("PGAIntegerBuildDatatype");

    traveller = PGAGetIndividual (ctx, p, pop);
    idx = PGABuildDatatypeHeader (ctx, p, pop, counts, displs, types);

    MPI_Get_address (traveller->chrom, &displs [idx]);
    counts [idx] = ctx->ga.StringLen;
    types  [idx] = MPI_LONG;
    idx++;

    MPI_Type_create_struct (idx, counts, displs, types, &individualtype);
    MPI_Type_commit (&individualtype);

    PGADebugExited ("PGAIntegerBuildDatatype");

    return individualtype;
}

/*!****************************************************************************
    \brief Compute genetic distance of two strings.
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
    Internal function.  Use :c:func:`PGAUserFunctionGeneDistance`.

    \endrst

******************************************************************************/
double PGAIntegerGeneDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAInteger *c1 = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAInteger *c2 = (PGAInteger *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int ret = 0;
    int i;

    PGADebugEntered ("PGAIntegerGeneDistance");
    for (i=0; i<ctx->ga.StringLen; i++) {
        ret += labs (c1 [i] - c2 [i]);
    }
    PGADebugExited ("PGAIntegerGeneDistance");
    return ret;
}

/*!****************************************************************************
    \brief Compute genetic distance of two strings.
    \ingroup standard-api

    \param   ctx    context variable
    \param   p1     first string index
    \param   pop1   symbolic constant of the population the first string is in
    \param   p2     second string index
    \param   pop2   symbolic constant of the population the second string is in
    \return  genetic euclidian distance of the two strings

    \rst

    Description
    -----------

    This uses the Euclidian distance metric, the square-root of the sum
    of all squared differences of each allele. It can be used to
    override the default integer genetic distance function (which uses a
    manhattan distance metric) using :c:func:`PGASetUserFunction` with
    the :c:macro:`PGA_USERFUNCTION_GEN_DISTANCE` setting.

    Example
    -------

    Override genetic distance function:

    .. code-block:: c

       PGAContext *ctx;

       ...
       assert (PGAGetDataType (ctx) == PGA_DATATYPE_INTEGER);
       PGASetUserFunction
        (ctx, PGA_USERFUNCTION_GEN_DISTANCE, PGAIntegerEuclidianDistance);

    \endrst

******************************************************************************/
double PGAIntegerEuclidianDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAInteger *c1 = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAInteger *c2 = (PGAInteger *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    double ret = 0.0;
    int i;

    PGADebugEntered ("PGAIntegerGeneDistance");
    for (i=0; i<ctx->ga.StringLen; i++) {
        ret += pow (c1 [i] - c2 [i], 2);
    }
    PGADebugExited ("PGAIntegerGeneDistance");
    return sqrt (ret);
}
