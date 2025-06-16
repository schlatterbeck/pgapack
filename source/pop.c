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
* This file contains routines that act on entire populations.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

#define OPT_DIR_CMP(ctx, e1, e2) \
    (ctx->ga.optdir == PGA_MAXIMIZE ? CMP ((e2), (e1)) : CMP ((e1), (e2)))
#define OPT_DIR_CMP_EV(ctx, e1, e2) \
    (((ctx)->nsga.is_ev) ? OPT_DIR_CMP((ctx), (e1), (e2)) : CMP ((e1), (e2)))
#define NORMALIZE(ctx, e, u) \
    (ctx->ga.optdir == PGA_MAXIMIZE ? ((u) - (e)) : ((e) - (u)))
#define DENORMALIZE(ctx, e, u) \
    (ctx->ga.optdir == PGA_MAXIMIZE ? ((u) - (e)) : ((e) + (u)))

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Creates an (internal) array of indices according to one of
           three criteria.
    \ingroup explicit
    \param   ctx  context variable
    \param   pop  symbolic constant of the population from which to
                  create the sorted array
    \return  An internal array of indices sorted according to one of
             three criteria is created

    \rst

    Description
    -----------

    If :c:macro:`PGA_POPREPL_BEST` is used (the default) the array is
    sorted from most fit to least fit. If :c:macro:`PGA_POPREPL_RANDOM_REP`
    is used the indices in the array are selected randomly with replacement.
    If :c:macro:`PGA_POPREPL_RANDOM_NOREP` is used the indices in the
    array are selected randomly without replacement. The function
    :c:func:`PGASetPopReplaceType` is used to specify which strategy is
    used.  The indices of the sorted population members may then be
    accessed from the internal array via :c:func:`PGAGetSortedPopIndex`.
    This routine is typically used during population replacement.

    Example
    -------

    Copy the five best strings from the old population into the new
    population.  The rest of the new population will be created by
    recombination, and is not shown.

    .. code-block:: c

       PGAContext *ctx;
       int i, j;

       ...
       PGASetPopReplaceType (ctx, PGA_POPREPL_BEST)
       ...
       PGASortPop (ctx, PGA_OLDPOP);
       for (i=0; i<5; i++) {
           j = PGAGetSortedPopIndex (ctx, i);
           PGACopyIndividual (ctx, j, PGA_OLDPOP, i, PGA_NEWPOP);
       }

    \endrst

******************************************************************************/
void PGASortPop (PGAContext *ctx, int pop)
{
    int i,j;

    PGADebugEntered ("PGASortPop");
    if (pop != PGA_OLDPOP && pop != PGA_NEWPOP) {
        PGAFatalPrintf (ctx, "PGASort: Invalid value of pop: %d", pop);
    }
    switch (ctx->ga.PopReplace) {
    case PGA_POPREPL_BEST:
        /* No need to init ga.sorted, done by PGAEvalSort */
        PGAEvalSort (ctx, pop, ctx->ga.sorted);
        break;
    case PGA_POPREPL_RANDOM_REP:
        for (i = 0; i < ctx->ga.PopSize; i++) {
            ctx->scratch.intscratch [i] = i;
        };
        for (i = 0; i < ctx->ga.PopSize; i++) {
            j = PGARandomInterval (ctx, 0, ctx->ga.PopSize-1);
            ctx->ga.sorted [i] = ctx->scratch.intscratch [j];
        };
        break;
    case PGA_POPREPL_RANDOM_NOREP:
        for (i = 0; i < ctx->ga.PopSize; i++) {
            ctx->scratch.intscratch [i] = i;
        };
        for (i = 0; i < ctx->ga.PopSize; i++) {
            j = PGARandomInterval (ctx, 0, ctx->ga.PopSize-i-1);
            ctx->ga.sorted [i] = ctx->scratch.intscratch [j];
            ctx->scratch.intscratch [j] =
                ctx->scratch.intscratch [ctx->ga.PopSize-i-1];
        };
        break;
    }
    PGADebugExited ("PGASortPop");
}


/*!***************************************************************************
    \brief Return the population size
    \ingroup query
    \param   ctx  context variable
    \return  The population size

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int popsize;

       ...
       popsize = PGAGetPopSize (ctx);

    \endrst

*****************************************************************************/
int PGAGetPopSize (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetPopSize");
    PGAFailIfNotSetUp ("PGAGetPopSize");

    PGADebugExited ("PGAGetPopSize");

    return ctx->ga.PopSize;
}

/*!***************************************************************************
    \brief Return the maximum number of strings to replace in each generation.
    \ingroup query
    \param   ctx  context variable
    \return The maximum number number of strings to replace each generation

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int numreplace;

       ...
       numreplace = PGAGetNumReplaceValue (ctx);

    \endrst

*****************************************************************************/
int PGAGetNumReplaceValue (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetNumReplaceValue");
    PGAFailIfNotSetUp ("PGAGetNumReplaceValue");

    PGADebugExited ("PGAGetNumReplaceValue");

    return ctx->ga.NumReplace;
}

/*!***************************************************************************
    \brief Return the symbolic constant used to determine which strings
           to copy from the old population to the new population.
    \ingroup query
    \param   ctx  context variable
    \return  The symbolic constant of the replacement strategy

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int popreplace;

       ...
       popreplace = PGAGetPopReplaceType (ctx);
       switch (popreplace) {
       case PGA_POPREPL_BEST:
           printf ("Replacement Strategy = PGA_POPREPL_BEST\n");
           break;
       case PGA_POPREPL_RANDOM_REP:
           printf ("Replacement Strategy = PGA_POPREPL_RANDOM_REP\n");
           break;
       case PGA_POPREPL_RANDOM_NOREP:
           printf ("Replacement Strategy = PGA_POPREPL_RANDOM_NOREP\n");
           break;
       default:
           printf ("Another Replacement Strategy\n");
           break;
       }

    \endrst

******************************************************************************/
int PGAGetPopReplaceType (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetPopReplaceType");
    PGAFailIfNotSetUp ("PGAGetPopRelaceType");

    PGADebugExited ("PGAGetPopReplaceType");

    return ctx->ga.PopReplace;
}

/*!****************************************************************************
    \brief Set window size used for restricted tournament replacement.
    \ingroup init
    \param   ctx         context variable
    \param   windowsize  size of the window for restricted tournament
                         replacement
    \return  None

    \rst

    Description
    -----------

    This function will have no effect unless :c:macro:`PGA_POPREPL_RTR`
    was specified as  the population replacement strategy with
    :c:func:`PGASetPopReplaceType`.
    The window size must be smaller than the population size.
    The default is :math:`\min (n, N/20)` where :math:`n` is the string length
    and :math:`N` is the population size.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetRTRWindowSize (ctx, windowsize);

    \endrst

******************************************************************************/
void PGASetRTRWindowSize (PGAContext *ctx, int windowsize)
{
    PGADebugEntered ("PGASetRTRWindowSize");

    ctx->ga.RTRWindowSize = windowsize;

    PGADebugExited ("PGASetRTRWindowSize");
}

/*!***************************************************************************
    \brief Return the window size for restricted tournament replacement.
    \ingroup query
    \param   ctx  context variable
    \return  The size of the window for restricted tournament selection

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int windowsize;

       ...
       windowsize = PGAGetRTRWindowSize (ctx);

    \endrst

*****************************************************************************/
int PGAGetRTRWindowSize (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetRTRWindowSize");

    PGADebugExited ("PGAGetRTRWindowSize");

    return ctx->ga.RTRWindowSize;
}

/*!****************************************************************************
    \brief Return a population string index from the array created by
           sorting of the population.
    \ingroup query
    \param   ctx       context variable
    \param   n         specified which index element is to be returned.
    \return  A population string index from the array created by
             \ref PGASortPop

    \rst

    Example
    -------

    Copy the five best strings from the old population into the new
    population.  The rest of the new population will be created by
    recombination, and is not shown.

    .. code-block:: c

       PGAContext *ctx;
       int i, j;

       ...
       PGASetPopReplaceType (ctx,PGA_POPREPL_BEST)
       PGASortPop (ctx, PGA_OLDPOP);
       for (i=0; i<5; i++) {
           j = PGAGetSortedPopIndex (ctx, i);
           PGACopyIndividual (ctx, j, PGA_OLDPOP, i, PGA_NEWPOP);
       }

    \endrst

******************************************************************************/
int PGAGetSortedPopIndex (PGAContext *ctx, int n)
{
    int temp = 0;

    PGADebugEntered ("PGAGetSortedPopIndex");
    if (n >= 0 && n < ctx->ga.PopSize) {
        temp = ctx->ga.sorted [n];
    } else {
        PGAFatalPrintf (ctx, "PGAGetSorted: Invalid value of n: %d", n);
    }

    PGADebugExited ("PGAGetSortedPopIndex");

    return temp;
}

/*!****************************************************************************
    \brief Specify the size of the genetic algorithm population.
    \ingroup init

    \param   ctx      context variable
    \param   popsize  the genetic algorithm population size to use
    \return  None

    \rst

    Description
    -----------

    The default population size is 100, unless reference directions or
    reference points have been specified for NSGA-III replacement in
    which case the default is the number of reference points plus the
    number of reference directions.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetPopSize (ctx, 200);

    \endrst

******************************************************************************/
void PGASetPopSize (PGAContext *ctx, int popsize)
{

    PGADebugEntered ("PGASetPopSize");
    PGAFailIfSetUp  ("PGASetPopSize");

    if (popsize < 1 || popsize % 2) {
        PGAFatalPrintf
            (ctx, "PGASetPopSize: Invalid value of popsize: %d", popsize);
    } else {
        ctx->ga.PopSize = popsize;
    }

    PGADebugExited ("PGASetPopSize");
}



/*!****************************************************************************
    \brief Specify the number of new strings to create each generation.
    \ingroup init
    \param   ctx          context variable
    \param   pop_replace  the number of population members to create
                          each generation
    \return  None

    \rst

    Description
    -----------

    The default is ten percent of the population size.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetNumReplaceValue (ctx, 35);

    \endrst

******************************************************************************/
void PGASetNumReplaceValue (PGAContext *ctx, int pop_replace)
{
    PGADebugEntered ("PGASetNumReplaceValue");

    if (pop_replace < 0) {
        PGAFatalPrintf
            ( ctx
            , "PGASetNumReplaceValue: Invalid value of pop_replace: %d"
            , pop_replace
            );
    } else {
        ctx->ga.NumReplace = pop_replace;
    }

    PGADebugExited ("PGASetNumReplaceValue");
}




/*!****************************************************************************
    \brief Choose method of replacing strings in the new population.
    \ingroup init
    \param   ctx          context variable
    \param   pop_replace  symbolic constant to specify the population
                          replacement strategy
    \return  None

    \rst

    Description
    -----------

    Valid choices are :c:macro:`PGA_POPREPL_BEST`,
    :c:macro:`PGA_POPREPL_RANDOM_NOREP`, or
    :c:macro:`PGA_POPREPL_RANDOM_REP` for copying the best
    strings, or  random string, with or without replacement, respectively,
    from the old population into the new population. Additional
    replacement types are :c:macro:`PGA_POPREPL_RTR` for restricted
    tournament replacement, :c:macro:`PGA_POPREPL_PAIRWISE_BEST` for
    pairwise comparison of each individual in the old/new population,
    and :c:macro:`PGA_POPREPL_NSGA_II` and :c:macro:`PGA_POPREPL_NSGA_III`
    for multiobjective optimization using the Nondominated Sorting
    Genetic Algorithm (NSGA-II or NSGA-III).
    The default is :c:macro:`PGA_POPREPL_BEST`.
    See :ref:`group:const-poprep` for the constants and section
    :ref:`sec:population-replacement` in the user guide for details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetPopReplaceType (ctx, PGA_POPREPL_RANDOM_NOREP);

    \endrst

******************************************************************************/
void PGASetPopReplaceType (PGAContext *ctx, int pop_replace)
{
    PGADebugEntered ("PGASetPopReplaceType");

    switch (pop_replace) {
    case PGA_POPREPL_BEST:
    case PGA_POPREPL_RANDOM_NOREP:
    case PGA_POPREPL_RANDOM_REP:
    case PGA_POPREPL_RTR:
    case PGA_POPREPL_PAIRWISE_BEST:
    case PGA_POPREPL_NSGA_II:
    case PGA_POPREPL_NSGA_III:
        ctx->ga.PopReplace = pop_replace;
        break;
    default:
        PGAFatalPrintf
            ( ctx
            , "PGASetPopReplaceType: Invalid value of pop_replace: %d"
            , pop_replace
            );
        break;
    }

    PGADebugExited ("PGASetPopReplaceType");
}

/*!****************************************************************************
    \brief Set reference points on reference hyperplane for NSGA-III.
    \ingroup init
    \param   ctx      context variable
    \param   npoints  Number of points
    \param   points   Pointer to points
    \return  None

    \rst

    Description
    -----------

    Reference points should be allocated by :c:func:`LIN_dasdennis`.
    They are freed when calling :c:func:`PGADestroy`. If no explicit
    reference points are set a default will be computed by
    :c:func:`PGASetUp`.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       ...
       int dim = PGAGetNumAuxEval (ctx) - PGAGetNumConstraint (ctx) + 1;
       void *p = NULL;
       int np  = 0;

       ...
       np = LIN_dasdennis (dim, 3, &p, 0, 1, NULL);
       PGASetReferencePoints (ctx, np, p);

    \endrst

******************************************************************************/
void PGASetReferencePoints (PGAContext *ctx, size_t npoints, void *points)
{
    if (ctx->ga.nrefpoints) {
        PGAFatalPrintf (ctx, "Can't set reference points twice");
    }
    if (points == NULL) {
        PGAFatalPrintf (ctx, "Need non-NULL points");
    }
    ctx->ga.nrefpoints = npoints;
    ctx->ga.refpoints  = points;
}

/*!****************************************************************************
    \brief Set reference directions for NSGA-III.
    \ingroup init
    \param   ctx    context variable
    \param   ndirs  Number of directions
    \param   dirs   Pointer to directions
    \param   npart  Number of Das / Dennis partitions
    \param   scale  Scale factor for constructed Das / Dennis points,
                    must be 0 < scale <= 1 but will typically be < 0.5
    \return  None

    \rst

    Description
    -----------

    A direction is a point in objective space and can be seen as a vector
    from the origin to that point. During optimization the reference
    directions are mapped to the reference hyperplane and a scaled
    Das/Dennis hyperplane is constructed around that point.
    Each direction consists of ``dimension`` ``double`` variables.

    Example
    -------

    Asume 3 dimensions, i.e. 3 evaluation functions

    .. code-block:: c


       PGAContext *ctx;
       double dirs [][3] = {{1, 2, 3}, {4, 5, 6}};
       ...
       int dim = PGAGetNumAuxEval (ctx) - PGAGetNumConstraint (ctx) + 1;

       ...
       PGASetReferenceDirections (ctx, 2, dirs, 5, 0.1);

    \endrst

******************************************************************************/
void PGASetReferenceDirections
    (PGAContext *ctx, size_t ndirs, void *dirs, int npart, double scale)
{
    if (ctx->ga.nrefdirs) {
        PGAFatalPrintf (ctx, "Can't set reference directions twice");
    }
    if (dirs == NULL) {
        PGAFatalPrintf (ctx, "Need non-NULL directions");
    }
    ctx->ga.nrefdirs   = ndirs;
    ctx->ga.refdirs    = dirs;
    ctx->ga.ndir_npart = npart;
    ctx->ga.dirscale   = scale;
}

/*!****************************************************************************
   \brief Perform restricted tournament replacement.
   \ingroup explicit
   \param   ctx          context variable
   \return  None

   \rst

   Description
   -----------

   For each individual in :c:macro:`PGA_NEWPOP` we select a window of
   individuals from :c:macro:`PGA_OLDPOP`, find the one genetically most
   like the new candidate and replace the individual if the new
   candidate has better evalutation. Note that we may not use the
   fitness here: Fitness from two different populations are
   uncompareable! After this populations are swapped (exchange of
   :c:macro:`PGA_NEWPOP` and :c:macro:`PGA_OLDPOP`) for further processing.

   Example
   -------

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGARestrictedTournamentReplacement (ctx);

    \endrst

******************************************************************************/
void PGARestrictedTournamentReplacement (PGAContext *ctx)
{
    int i, j;
    int popsize = PGAGetPopSize (ctx);
    int numreplace = PGAGetNumReplaceValue (ctx);
    PGASampleState state;
    PGAIndividual *temp;
    int oldpop = PGA_OLDPOP;
    int newpop = PGA_NEWPOP;

    PGADebugEntered ("PGARestrictedTournamentReplacement");
    for (i=popsize - numreplace; i<popsize; i++) {
        double dist = -1.0;
        int closest = -1;
        PGARandomSampleInit
            (ctx, &state, ctx->ga.RTRWindowSize, ctx->ga.PopSize);
        /* Avoid duplicates */
        if (PGADuplicate (ctx, i, PGA_NEWPOP, PGA_OLDPOP)) {
            continue;
        }
        for (j=0; j<ctx->ga.RTRWindowSize; j++) {
            double d;
            int idx = PGARandomNextSample (&state);
            if (ctx->fops.GeneDistance) {
                d = (*ctx->fops.GeneDistance)
                    (&ctx, &idx, &oldpop, &i, &newpop);
            } else {
                d = (*ctx->cops.GeneDistance)
                    (ctx, idx, PGA_OLDPOP, i, PGA_NEWPOP);
            }
            if (dist < 0 || d < dist) {
                dist = d;
                closest = idx;
            }
        }
        assert (closest >= 0);

        /* If new population individual is better */
        if (PGAEvalCompare (ctx, i, PGA_NEWPOP, closest, PGA_OLDPOP) <= 0) {
            /* Remove old individual from hash */
            PGAUnHashIndividual (ctx, closest, PGA_OLDPOP);
            /* Copy i in PGA_NEWPOP to closest in PGA_OLDPOP */
            PGACopyIndividual (ctx, i, PGA_NEWPOP, closest, PGA_OLDPOP);
            /* Add new individual to hash */
            PGAHashIndividual (ctx, closest, PGA_OLDPOP);
        }
    }
    /* Exchange old/newpop, will be done again by PGAUpdateGeneration,
     * we just make sure that the result of above replacements is now
     * newpop.
     */
    temp           = ctx->ga.oldpop;
    ctx->ga.oldpop = ctx->ga.newpop;
    ctx->ga.newpop = temp;
    PGADebugExited ("PGARestrictedTournamentReplacement");
}

/*!****************************************************************************

    \brief Perform pairwise best replacement.
    \ingroup explicit
    \param   ctx          context variable
    \return  None

    \rst

    Description
    -----------

    Compare individuals with same index in :c:macro:`PGA_OLDPOP` and
    :c:macro:`PGA_NEWPOP` and select the one with better evalutation.
    Note that we may not use the fitness here: Fitness from two
    different populations are uncompareable!
    This replacement strategy is used in evolutionary algorithms that
    modify a single individual and replace the parent if the offspring is
    better. A popular example is Differential Evolution (DE).
    After this populations are swapped (exchange of
    :c:macro:`PGA_NEWPOP` and :c:macro:`PGA_OLDPOP`) for further processing.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGAPairwiseBestReplacement (ctx);

    \endrst

******************************************************************************/
void PGAPairwiseBestReplacement (PGAContext *ctx)
{
    int i;
    int popsize = PGAGetPopSize (ctx);
    int numreplace = PGAGetNumReplaceValue (ctx);
    PGAIndividual *temp;

    PGADebugEntered ("PGAPairwiseBestReplacement");
    for (i=popsize - numreplace; i<popsize; i++) {
        /* Avoid duplicates */
        if (PGADuplicate (ctx, i, PGA_NEWPOP, PGA_OLDPOP)) {
            continue;
        }
        /* Note the '<=' comparison, differential evolution can walk across
         * areas with equal evaluation this way
         */
        if (PGAEvalCompare (ctx, i, PGA_NEWPOP, i, PGA_OLDPOP) <= 0) {
            /* Remove old individual from hash */
            PGAUnHashIndividual (ctx, i, PGA_OLDPOP);
            /* Copy i in PGA_NEWPOP to i in PGA_OLDPOP */
            PGACopyIndividual (ctx, i, PGA_NEWPOP, i, PGA_OLDPOP);
            /* Add new individual to hash */
            PGAHashIndividual (ctx, i, PGA_OLDPOP);
        }
    }
    /* Exchange old/newpop, will be done again by PGAUpdateGeneration,
     * we just make sure that the result of above replacements is now
     * newpop.
     */
    temp           = ctx->ga.oldpop;
    ctx->ga.oldpop = ctx->ga.newpop;
    ctx->ga.newpop = temp;
    PGADebugExited ("PGAPairwiseBestReplacement");
}

/* Helper functions for PGA_NSGA_II_Replacement */
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

#define NONNEGEVAL(v, is_ev) ((is_ev) ? (v) : (((v) < 0) ? 0 : (v)))
#define GETEVAL(ind, fidx, is_ev)                       \
    ( ((fidx) == 0)                                     \
    ? (ind)->evalue                                     \
    : NONNEGEVAL ((ind)->auxeval [(fidx) - 1], (is_ev)) \
    )
#define GETEVAL_EV(ind, fidx) \
    GETEVAL((ind), (fidx) + (ind)->ctx->nsga.base, (ind)->ctx->nsga.is_ev)

static int crowdsort_cmp (const void *a1, const void *a2)
{
    PGAIndividual **i1 = (void *)a1;
    PGAIndividual **i2 = (void *)a2;
    double e1 = GETEVAL_EV (*i1, (*i1)->funcidx);
    double e2 = GETEVAL_EV (*i2, (*i2)->funcidx);
    return CMP (e1, e2);
}

static int nondom_cmp (const void *a1, const void *a2)
{
    const PGAIndividual * const *i1 = a1;
    const PGAIndividual * const *i2 = a2;
    if ((*i1)->rank < (*i2)->rank) {
        return -1;
    }
    if ((*i1)->rank > (*i2)->rank) {
        return 1;
    }
    return CMP ((*i2)->crowding, (*i1)->crowding);
}

/* typedef to make it easier to pass crowding functions as parameter */
typedef void (* crowding_t)
    (PGAContext *, PGAIndividual **, size_t, unsigned int);

/* Compute crowding distance over the given individuals
 * This is specific to NSGA-II.
 */
STATIC void crowding
    (PGAContext *ctx, PGAIndividual **start, size_t n, unsigned int rank)
{
    size_t i;
    int k;
    size_t nrank = 0;
    PGAIndividual **crowd = ctx->scratch.nsga_tmp.ind_tmp;
    /* declare values of functions as arrays on the stack, this should
     * never exceed memory limits
     */
    DECLARE_DYNARRAY (double, f_min, ctx->nsga.nfun);
    DECLARE_DYNARRAY (double, f_max, ctx->nsga.nfun);

    for (i=0; i<n; i++) {
        PGAIndividual *ind = start [i];
        ind->crowding = 0;
        if (ind->rank == rank) {
            for (k=0; k<ctx->nsga.nfun; k++) {
                double e = GETEVAL_EV (ind, k);
                if (nrank == 0 || f_min [k] > e) {
                    f_min [k] = e;
                }
                if (nrank == 0 || f_max [k] < e) {
                    f_max [k] = e;
                }
            }
            crowd [nrank++] = ind;
        }
    }
    assert (nrank > 0);
    for (k=0; k<ctx->nsga.nfun; k++) {
        double norm = f_max [k] - f_min [k];
        for (i=0; i<nrank; i++) {
            (crowd [i])->funcidx = k;
        }
        qsort (crowd, nrank, sizeof (crowd [0]), crowdsort_cmp);
        (crowd [0])->crowding = DBL_MAX;
        (crowd [nrank-1])->crowding = DBL_MAX;
        for (i=1; i<nrank-1; i++) {
            if ((crowd [i])->crowding != DBL_MAX) {
                (crowd [i])->crowding +=
                    ( GETEVAL_EV (crowd [i+1], k)
                    - GETEVAL_EV (crowd [i-1], k)
                    ) / norm;
            }
        }
    }
}

/* Compute utopian point as the minimum over all solutions */
STATIC void compute_utopian (PGAContext *ctx, PGAIndividual **start, int n)
{
    int i, j;
    int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;

    if (!ctx->ga.utopian_valid) {
        for (j=0; j<dim; j++) {
            ctx->ga.utopian [j] = GETEVAL (start [0], j, 1);
        }
        ctx->ga.utopian_valid = PGA_TRUE;
    }
    for (i=0; i<n; i++) {
        for (j=0; j<dim; j++) {
            double e = GETEVAL (start [i], j, 1);
            if (OPT_DIR_CMP (ctx, ctx->ga.utopian [j], e) > 0) {
                ctx->ga.utopian [j] = e;
            }
        }
    }
}

/* When computing ASF set values in vector < EPS_VAL to 0
 * Taken from pymoo code, not documented in any paper.
 */
#define EPS_VAL 1e-3
#define EPS_ASF 1e-6
#define EPS_NAD 1e-6
#define WEIGHT(a, j) ((a) == (j) ? 1 : EPS_ASF)

/* Compute ASF (see nsga-iii paper), axis is <= dim - 1 */
STATIC double compute_asf (PGAContext *ctx, const double *point, int axis)
{
    int j;
    int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
    double *utop = ctx->ga.utopian;
    double asf = NORMALIZE (ctx, point [0], utop [0]) / WEIGHT (axis, 0);

    for (j=1; j<dim; j++) {
        double a = NORMALIZE (ctx, point [j], utop [j]);
        if (a < EPS_VAL) {
            a = 0;
        }
        a /= WEIGHT (axis, j);
        if (a > asf) {
            asf = a;
        }
    }
    return asf;
}

/* Compute extreme points */
STATIC void compute_extreme (PGAContext *ctx, PGAIndividual **start, int n)
{
    int i, j;
    int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
    DECLARE_DYNPTR (double, extreme, dim) = ctx->ga.extreme;
    DECLARE_DYNARRAY (double, asfmin, dim);

    if (!ctx->ga.extreme_valid) {
        for (i=0; i<dim; i++) {
            for (j=0; j<dim; j++) {
                DEREF2_DYNPTR (extreme, dim, i, j) = GETEVAL (start [0], j, 1);
            }
        }
    }
    ctx->ga.extreme_valid = PGA_TRUE;
    for (j=0; j<dim; j++) {
        asfmin [j] = compute_asf (ctx, DEREF1_DYNPTR (extreme, dim, j), j);
    }

    for (i=0; i<n; i++) {
        DECLARE_DYNARRAY (double, e, dim);
        for (j=0; j<dim; j++) {
            e [j] = GETEVAL (start [i], j, 1);
        }
        for (j=0; j<dim; j++) {
            double asf = compute_asf (ctx, e, j);
            if (asf < asfmin [j]) {
                asfmin [j] = asf;
                memcpy
                    (DEREF1_DYNPTR (extreme, dim, j), e, sizeof (double) * dim);
            }
        }
    }
}

/* Preferred nadir estimate via extreme points and axes intersect */
STATIC int compute_intersect (PGAContext *ctx, PGAIndividual **start, int n)
{
    int i, j, d;
    int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
    DECLARE_DYNPTR (double, extreme, dim) = ctx->ga.extreme;
    DECLARE_DYNARRAY2 (double, m, dim, dim);
    DECLARE_DYNARRAY (double, x, dim);
    DECLARE_DYNARRAY (double, v, dim);

    for (d=0; d<dim; d++) {
        int result;
        for (i=0; i<dim-1; i++) {
            for (j=0; j<dim; j++) {
                /* No need to normalize with utopian here, cancels */
                DEREF2_DYNPTR (m, dim, j, i)
                    = DEREF2_DYNPTR (extreme, dim, i+1, j)
                    - DEREF2_DYNPTR (extreme, dim, 0, j);
            }
        }
        for (j=0; j<dim; j++) {
            DEREF2_DYNPTR (m, dim, j, dim-1) = j == d ? -1 : 0;
            v [j] = -NORMALIZE
                (ctx, DEREF2_DYNPTR (extreme, dim, 0, j), ctx->ga.utopian [j]);
        }
        result = LIN_solve (dim, m, v);
        /* Matrix singular? */
        if (result != 0) {
            break;
        }
        x [d] = v [dim - 1];
        /* Intercept too small or negative? */
        if (x [d] <= EPS_NAD) {
            break;
        }
    }
    if (d >= dim) {
        /* Success: Use nadir estimate from hyper-plane */
        for (j=0; j<dim; j++) {
            ctx->ga.nadir [j] = DENORMALIZE (ctx, x [j], ctx->ga.utopian [j]);
        }
        return 0;
    }
    /* Fail: No nadir estimate via extreme points, fall back to wof0 */
    PGAErrorPrintf
        ( ctx
        , PGA_WARNING
        , "Intercept computation failed in Generation %d\n"
        , ctx->ga.iter
        );
    return 1;
}

/* Compute worst of population and nadir estimate (worst of front 0) */
STATIC void compute_worst
    (PGAContext *ctx, PGAIndividual **start, int n, double *wpop, double *wof0)
{
    int i, j;
    int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
    int wv = 0;
    for (i=0; i<n; i++) {
        for (j=0; j<dim; j++) {
            double e = GETEVAL (start [i], j, 1);
            if (  !ctx->ga.worst_valid
               || OPT_DIR_CMP (ctx, ctx->ga.worst [j], e) < 0
               )
            {
                ctx->ga.worst [j] = e;
            }
            if (i==0 || OPT_DIR_CMP (ctx, wpop [j], e) < 0) {
                wpop [j] = e;
            }
            if (  start [i]->rank == 0
               && (wv <= j || OPT_DIR_CMP (ctx, wof0 [j], e) < 0)
               )
            {
                wof0 [j] = e;
                if (wv <= j) {
                    wv++;
                }
            }
        }
        ctx->ga.worst_valid = PGA_TRUE;
    }
    assert (wv == dim);
}

STATIC void compute_nadir (PGAContext *ctx, PGAIndividual **start, int n)
{
    int j;
    int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
    DECLARE_DYNARRAY (double, wof0, dim); /* Worst of front 0 */
    DECLARE_DYNARRAY (double, wpop, dim); /* Worst of population */
    int ret;

    ret = compute_intersect (ctx, start, n);
    compute_worst (ctx, start, n, wpop, wof0);
    /* If a component in estimated nadir is *worse* than corresponding
     * component of the worst point ever, replace this component.
     * Taken from pymoo, not documented in any paper.
     */
    if (ret == 0) {
        for (j=0; j<dim; j++) {
            if (OPT_DIR_CMP (ctx, ctx->ga.worst [j], ctx->ga.nadir [j]) < 0) {
                ctx->ga.nadir [j] = ctx->ga.worst [j];
            }
        }
    } else {
        /* Matrix inversion failed, need to fall back on wof0/wpop */
        for (j=0; j<dim; j++) {
            if (fabs (wof0 [j] - ctx->ga.utopian [j]) < EPS_NAD) {
                ctx->ga.nadir [j] = wpop [j];
            } else {
                ctx->ga.nadir [j] = wof0 [j];
            }
        }
    }
}

/* Get points in refpoints and refdir-cloud */
static double *get_point (PGAContext *ctx, size_t idx)
{
    int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
    DECLARE_DYNPTR (double, normdirs,  dim) = ctx->ga.normdirs;
    DECLARE_DYNPTR (double, refpoints, dim) = ctx->ga.refpoints;
    assert (idx < ctx->ga.ndpoints * ctx->ga.nrefdirs + ctx->ga.nrefpoints);
    if (idx >= ctx->ga.nrefpoints) {
        return DEREF1_DYNPTR (normdirs, dim, idx - ctx->ga.nrefpoints);
    }
    return DEREF1_DYNPTR (refpoints, dim, idx);
}

static int assoc_cmp (const void *a1, const void *a2)
{
    const PGAIndividual * const *i1 = a1;
    const PGAIndividual * const *i2 = a2;
    if ((*i1)->point_idx < (*i2)->point_idx) {
        return -1;
    }
    if ((*i1)->point_idx > (*i2)->point_idx) {
        return 1;
    }
    if ((*i1)->rank < (*i2)->rank) {
        return -1;
    }
    if ((*i1)->rank > (*i2)->rank) {
        return 1;
    }
    return CMP ((*i1)->distance, (*i2)->distance);
}

/* Compute niche preservation algorithm over the given individuals
 * This is specific to NSGA-III.
 */
static void niching
    (PGAContext *ctx, PGAIndividual **start, size_t n, unsigned int rank)
{
    size_t i, j;
    size_t dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
    size_t npoints = ctx->ga.ndpoints * ctx->ga.nrefdirs + ctx->ga.nrefpoints;
    DECLARE_DYNARRAY (double, point, dim);

    compute_utopian (ctx, start, n);
    compute_extreme (ctx, start, n);
    compute_nadir   (ctx, start, n);
    /* Normalize points to hyperplane */
    for (i=0; i<n; i++) {
        PGAIndividual *ind = start [i];
        /* Init crowding metric for all individuals */
        ind->crowding = 0;
        for (j=0; j<dim; j++) {
            double e = GETEVAL (start [i], j, 1);
            e  = NORMALIZE (ctx, e, ctx->ga.utopian [j]);
            e /= NORMALIZE (ctx, ctx->ga.nadir [j], ctx->ga.utopian [j]);
            ind->normalized [j] = e;
        }
        LIN_normalize_to_refplane (dim, ind->normalized);
    }
    /* map reference directions to hyperplane
     * and compute the resulting points
     */
    for (i=0; i<ctx->ga.nrefdirs; i++) {
        size_t sz = ctx->ga.ndpoints * sizeof (double) * dim;
        DECLARE_DYNPTR (double, refdirs, dim) = ctx->ga.refdirs;
        memcpy (point, DEREF1_DYNPTR (refdirs, dim, i), sizeof (double) * dim);
        for (j=0; j<dim; j++) {
            point [j]  = NORMALIZE (ctx, point [j], ctx->ga.utopian [j]);
            point [j] /= NORMALIZE
                (ctx, ctx->ga.nadir [j], ctx->ga.utopian [j]);
        }
        LIN_dasdennis_allocated
            ( dim, ctx->ga.ndir_npart
            , ctx->ga.dirscale, point
            , ctx->ga.ndpoints, ((char *)ctx->ga.normdirs) + i * sz
            );
    }
    for (j=0; j<n; j++) {
        PGAIndividual *ind = start [j];
        double mindist = LIN_euclidian_distance
            (dim, ind->normalized, get_point (ctx, 0));
        int minidx = 0;
        for (i=1; i<npoints; i++) {
            double *point = get_point (ctx, i);
            double d = LIN_euclidian_distance
                (dim, ind->normalized, point);
            if (d < mindist) {
                mindist = d;
                minidx  = i;
            }
        }
        ind->distance = mindist;
        ind->point_idx = minidx;
    }
    /* Sort individuals by associated point index, rank, distance */
    qsort (start, n, sizeof (PGAIndividual *), assoc_cmp);
    /* Iterate over individuals and set crowding metric to
     * the negative of the count of the current point
     * We don't care that points with lower rank have now a crowding
     * metric, too.
     */
    {
        int last_pointidx = -1;
        int pointcount = 0;
        for (j=0; j<n; j++) {
            PGAIndividual *ind = start [j];
            if (last_pointidx != ind->point_idx) {
                last_pointidx = ind->point_idx;
                pointcount = 0;
            } else {
                pointcount++;
            }
            ind->crowding = -pointcount;
        }
    }
}

/* Comparison function for sorting individuals by objectives */
static int obj_sort_cmp (const void *a1, const void *a2)
{
    PGAIndividual * const *i1 = (PGAIndividual * const *)a1;
    PGAIndividual * const *i2 = (PGAIndividual * const *)a2;
    PGAContext *ctx = (*i1)->ctx;

    double e1_i1 = GETEVAL_EV (*i1, 0);
    double e1_i2 = GETEVAL_EV (*i2, 0);
    int cmp = OPT_DIR_CMP_EV (ctx, e1_i1, e1_i2);

    if (cmp == 0) {
        /* If first objectives are equal, sort by second objective */
        double e2_i1 = GETEVAL_EV (*i1, 1);
        double e2_i2 = GETEVAL_EV (*i2, 1);
        cmp = OPT_DIR_CMP_EV (ctx, e2_i1, e2_i2);
    }

    return cmp;
}

/* Same as above but in addition reverse by rank */
static int obj_sort_cmp_rank (const void *a1, const void *a2)
{
    PGAIndividual * const *i1 = (PGAIndividual * const *)a1;
    PGAIndividual * const *i2 = (PGAIndividual * const *)a2;
    int r = 0;

    if ((r = obj_sort_cmp (a1, a2))) {
        return r;
    }
    return CMP ((*i2)->rank, (*i1)->rank);
}

/* Helper function to compare two solutions for domination
 * Special case when check_eq is 1 we return -1 if s1 == s2 in all
 * objectives. This case happens when nd_helper_b has partitioned
 * l and h into two sets (which indicates that l dominates h in a higher
 * dimension). Now if these are equal in all lower dimensions we still
 * need to return -1.
 */
static int compare_solutions
    (PGAContext *ctx, PGAIndividual *s1, PGAIndividual *s2, int m, int check_eq)
{
    int k;
    int cmp = 0;

    for (k = 0; k < m; k++) {
        double e1, e2;
        int ncmp;
        e1 = GETEVAL_EV (s1, k);
        e2 = GETEVAL_EV (s2, k);
        ncmp = OPT_DIR_CMP_EV (ctx, e1, e2);
        if (cmp && ncmp && ncmp != cmp) {
            return 0; /* Non-dominated */
        }
        /* Don't allow cmp to get zero again */
        if (ncmp) {
            cmp = ncmp;
        }
    }
    if (check_eq && cmp == 0) {
        return -1;
    }
    return cmp;
}

/* For sorting an array of double */
static int double_cmp (const void *a1, const void *a2)
{
    const double *d1 = a1;
    const double *d2 = a2;
    return CMP (*d1, *d2);
}

/*
 * Compute median by sorting and taking the middle elements.
 * Used for small sizes of values array.
 */
static double nlogn_select (double *values, size_t lo, size_t hi, size_t k)
{
    qsort (values + lo, hi - lo, sizeof (*values), double_cmp);
    return values [lo + k];
}

/* Need forward declaration for recursive call from pick_pivot */
static double quickselect (double *values, size_t lo, size_t hi, size_t k);

/* Find median of medians to use as pivot */
static double pick_pivot (double *values, size_t lo, size_t hi)
{
    size_t i;
    double *medians = NULL;
    double pivot = 0;

    assert (hi > lo);
    medians = malloc (sizeof (double) * (hi - lo) / 5);
    if (medians == NULL) {
        /* Cannot use PGAFatalPrintf here, no context variable */
        fprintf (stderr, "Out of memory in pick_pivot\nPGAError: Fatal\n");
        exit (-1);
    }

    /* Small case: Just return median */
    if (hi - lo <= 5) {
        return nlogn_select (values, lo, hi, (lo + hi) / 2);
    }

    /* This skips the last non-full chunk if (hi - lo) is not divisible
     * by 5. This doesn't affect the pivot computation too much.
     */
    for (i=lo; i<hi; i+=5) {
        if (i + 5 > hi) {
            break;
        }
        medians [(i - lo) / 5] = nlogn_select (values, i, i + 5, 2);
    }
    /* Should be same with truncation */
    assert ((i - lo) / 5 == (hi - lo) / 5);
    /* Recursive call to find median of medians */
    pivot = quickselect (medians, 0, (i - lo) / 5, (i - lo) / 10);
    free (medians);
    return pivot;
}

/* Linear-time index value selection algorithm
 * (quickselect with median-of-medians pivot)
 * We separate values into three sections: values < pivot (at the start),
 * values equal to pivot (at the end) and values > pivot in the middle.
 */
static double quickselect (double *values, size_t lo, size_t hi, size_t k)
{
    size_t i, j, pidx;
    double pivot, tmp;
    int n_pivots = 0;

    assert (hi > lo);
    assert (k >= 0 && lo + k < hi);
    /* For small arrays, sort and return median */
    if (hi - lo < 10) {
        return nlogn_select (values, lo, hi, k);
    }

    pivot = pick_pivot (values, lo, hi);

    /* Partition values into <= pivot and > pivot */
    for (pidx = lo; pidx < hi; pidx++) {
        if (values [pidx] == pivot) {
            for (j=hi-1; hi > 0 && j>=pidx && values [j] == pivot; j--) {
                n_pivots++;
                hi--;
            }
            if (pidx >= hi) {
                break;
            }
            tmp = values [pidx];
            values [pidx] = values [hi - 1];
            values [hi - 1] = tmp;
            hi--;
            n_pivots++;
        }
        assert (values [pidx] != pivot);
        if (values [pidx] >= pivot) {
            break;
        }
    }
    for (i = pidx + 1; i < hi; i++) {
        if (values [i] < pivot) {
            tmp = values [i];
            values [i] = values [pidx];
            values [pidx] = tmp;
            pidx++;
        } else if (values [i] == pivot) {
            for (j=hi-1; j>=i && values [j] == pivot; j--) {
                n_pivots++;
                hi--;
            }
            if (i >= hi) {
                break;
            }
            tmp = values [i];
            values [i] = values [hi - 1];
            values [hi - 1] = tmp;
            hi--;
            n_pivots++;
            /* Element from above might have been < pivot */
            if (values [i] < pivot) {
                tmp = values [i];
                values [i] = values [pidx];
                values [pidx] = tmp;
                pidx++;
            }
        }
    }
    /* pick_pivot always returns an *existing* value in values */
    assert (n_pivots);

    if (k < pidx - lo) {
        assert (pidx <= hi);
        return quickselect (values, lo, pidx, k);
    } else if (k < (pidx - lo) + n_pivots) {
        return pivot;
    } else {
        assert ((pidx - lo) + n_pivots <= k);
        return quickselect (values, pidx, hi, k - (pidx - lo) - n_pivots);
    }
}

/* Find median value for objective m in the set using linear-time
 * selection algorithm.
 * We do not find the "real" median, when the number of elements is
 * even. We only use the median for partitioning, so it makes no sense
 * to call quickselect twice just to get the average of the two middle
 * values.
 */
STATIC double find_median (PGAIndividual **s, size_t n, int m)
{
    PGAContext *ctx = (*s)->ctx;
    double *values = ctx->scratch.nsga_tmp.medval;
    size_t i;

    assert (n > 0);
    assert (m >= 0 && m <= ctx->ga.NumAuxEval);
    for (i = 0; i < n; i++) {
        values [i] = GETEVAL_EV (s [i], m);
    }
    return quickselect (values, 0, n, n / 2);
}

/* Split function to divide a set into two based on the median value of
 * objective m. We're doing this in-place and return two pointers.
 * We split so that the pivot ends up in the *upper* half. This is
 * because the median implementation when called with n/2 for 2 elements
 * makes n/2 = 1 the pivot. So this ensures that for this median call we
 * end up with to equal-sized halves.
 * The variable split_lower should be 0 by default indicates if the
 * pivot belongs to the lower set (default is the upper set)
 */
STATIC void split_set
    ( PGAIndividual **s, size_t n
    , PGAIndividual ***l, size_t *nl
    , PGAIndividual ***h, size_t *nh
    , int m, double pivot, int split_lower
    )
{
    size_t i, pidx = 0;
    PGAContext *ctx = (*s)->ctx;

    *nl = 0;
    *nh = 0;

    for (pidx=0; pidx<n; pidx++) {
        double val = GETEVAL_EV (s [pidx], m);
        if (OPT_DIR_CMP_EV (ctx, val, pivot) > 0) {
            break;
        }
        if (!split_lower && val == pivot) {
            break;
        }
    }
    for (i = pidx + 1; i<n; i++) {
        double val = GETEVAL_EV (s [i], m);
        if (  OPT_DIR_CMP_EV (ctx, val, pivot) < 0
           || (split_lower && val == pivot)
           )
        {
            PGAIndividual *tmp = s [i];
            s [i] = s [pidx];
            s [pidx] = tmp;
            pidx++;
        }
    }
    *l  = s;
    *nl = pidx;
    *h  = s + pidx;
    *nh = n - pidx;
}

/* 2D ranking of H according to L
 * Precondition is that ranks in L are already correctly computed.
 * We may have gaps in the ranks.
 */
STATIC void rank_2d_b
    ( PGAContext *ctx
    , PGAIndividual **l, size_t nl
    , PGAIndividual **h, size_t nh
    )
{
    size_t i;
    size_t l_idx = 0;
    size_t num_fronts = 0;
    /* Note: we're keeping only the last in each front */
    PGAIndividual **fronts = NULL;
    PGAIndividual *cur_l = NULL;
    PGAIndividual *last = NULL;
    size_t max_front = 0;
    double e1_l, e2_l;

    fronts = malloc (sizeof (*fronts) * ctx->ga.PopSize * 2);
    if (fronts == NULL) {
        PGAFatalPrintf (ctx, "Out of memory in rank_2d_b");
    }

    /* Initialize to zero */
    memset (fronts, 0, ctx->ga.PopSize * sizeof (PGAIndividual *));

    /* nl must be > 0 */
    assert (nl > 0);
    /* And it doesn't make sense to call this with h empty */
    assert (nh > 0);

    /* Sort L by first/second objective (ascending) */
    qsort (l, nl, sizeof (*l), obj_sort_cmp);

    /* Sort H by first/second objective (ascending) */
    qsort (h, nh, sizeof (*h), obj_sort_cmp);

    cur_l = l [l_idx];
    e1_l = GETEVAL_EV (cur_l, 0);
    e2_l = GETEVAL_EV (cur_l, 1);

    /* Build fronts, we sweep over L and H simultaneously */
    for (i = 0; i < nh; i++) {
        PGAIndividual *cur_h = h [i];
        double e1_h = GETEVAL_EV (cur_h, 0);
        double e2_h = GETEVAL_EV (cur_h, 1);
        int e1c = OPT_DIR_CMP_EV (ctx, e1_l, e1_h);
        int e2c = OPT_DIR_CMP_EV (ctx, e2_l, e2_h);
        double e1_last, e2_last;
        int e1cmp, e2cmp;

        /* This *continues* to sweep l */
        while (l_idx < nl && (e1c < 0 || (e1c == 0 && e2c <= 0))) {
            int j;
            size_t front, k;
            front = cur_l->rank;
            /* Get *next* front that is non-zero */
            for (k=front; k<=max_front; k++) {
                if (fronts [k] != NULL) {
                    break;
                }
            }
            if (k <= max_front && fronts [k] != NULL) {
                double ef2 = GETEVAL_EV (fronts [k], 1);
                int e2fc = OPT_DIR_CMP_EV (ctx, e2_l, ef2);
                /* Dominated by this front, happens due to higher dimension */
                if (e2fc >= 0) {
                    l_idx++;
                    if (l_idx < nl) {
                        cur_l = l [l_idx];
                        e1_l = GETEVAL_EV (cur_l, 0);
                        e2_l = GETEVAL_EV (cur_l, 1);
                        e1c  = OPT_DIR_CMP_EV (ctx, e1_l, e1_h);
                        e2c  = OPT_DIR_CMP_EV (ctx, e2_l, e2_h);
                    }
                    continue;
                }
            }
            if (fronts [front] != NULL) {
                /* Going to overwrite existing */
                num_fronts--;
            }
            fronts [front] = cur_l;
            num_fronts++;
            if (max_front < front) {
                max_front = front;
            }
            /* Need to invalidate lower fronts with larger eval
             * This is the case where stairs are discontinued in the
             * published algorithm
             */
            for (j = front - 1; j >= 0; j--) {
                double f_e1, f_e2;
                int f1_c, f2_c;
                if (fronts [j] == NULL) {
                    continue;
                }
                f_e1 = GETEVAL_EV (fronts [j], 0);
                f_e2 = GETEVAL_EV (fronts [j], 1);
                f1_c = OPT_DIR_CMP_EV (ctx, f_e1, e1_l);
                f2_c = OPT_DIR_CMP_EV (ctx, f_e2, e2_l);
                /* This should hold by construction: */
                assert (f1_c <= 0);
                if (f2_c > 0 || (f2_c == 0 && f1_c <= 0)) {
                    fronts [j] = NULL;
                    num_fronts--;
                } else {
                    break;
                }
            }
            l_idx++;
            if (l_idx < nl) {
                cur_l = l [l_idx];
                e1_l = GETEVAL_EV (cur_l, 0);
                e2_l = GETEVAL_EV (cur_l, 1);
                e1c  = OPT_DIR_CMP_EV (ctx, e1_l, e1_h);
                e2c  = OPT_DIR_CMP_EV (ctx, e2_l, e2_h);
            }
        }
        /* If cur_h is stricly smaller than first cur_l do nothing */
        if (!num_fronts) {
            continue;
        }
        last = fronts [max_front];
        e1_last = GETEVAL_EV (last, 0);
        e2_last = GETEVAL_EV (last, 1);
        e1cmp = OPT_DIR_CMP_EV (ctx, e1_last, e1_h);
        e2cmp = OPT_DIR_CMP_EV (ctx, e2_last, e2_h);

        /* When both evals are equal the higher dimensions already
         * established dominance, so that case is also a dominance case
         */
        if (e2cmp < 0 || (e2cmp == 0 && e1cmp <= 0)) {
            /* Current individual dominated by last front */
            if (cur_h->rank <= last->rank + 1) {
                cur_h->rank = last->rank + 1;
            }
        } else {
            /* Binary search to find right front
             * Note that we have to take gaps in the fronts into account
             */
            int low = 0;
            int high = max_front;
            PGAIndividual *lo_ind;
            double e1_lo, e2_lo;

            for (low = 0; fronts [low] == NULL; low++)
                ;
            while (low < high) {
                size_t mid = (low + high) / 2;
                PGAIndividual *i_mid;
                double e1_mid, e2_mid;

                for (; fronts [mid] == NULL; mid--)
                    ;
                i_mid = fronts [mid];
                e1_mid = GETEVAL_EV (i_mid, 0);
                e2_mid = GETEVAL_EV (i_mid, 1);
                e2cmp = OPT_DIR_CMP_EV (ctx, e2_mid, e2_h);
                e1cmp = OPT_DIR_CMP_EV (ctx, e1_mid, e1_h);

                if (e2cmp > 0) {
                    high = mid;
                } else {
                    for (low = mid + 1; fronts [low] == NULL; low++)
                        ;
                }
            }
            /* Skip down to next which should dominate
             */
            for (low = low - 1; low >= 0 && fronts [low] == NULL; low--)
                ;
            if (low >= 0) {
                lo_ind = fronts [low];
                e1_lo = GETEVAL_EV (lo_ind, 0);
                e2_lo = GETEVAL_EV (lo_ind, 1);
                e2cmp = OPT_DIR_CMP_EV (ctx, e2_lo, e2_h);
                e1cmp = OPT_DIR_CMP_EV (ctx, e1_lo, e1_h);
                assert (e1cmp <= 0 && e2cmp <= 0);
                /* Dominated even when equal */
                if (e1cmp <= 0 && e2cmp <= 0) {
                    if (cur_h->rank < fronts [low]->rank + 1) {
                        cur_h->rank = fronts [low]->rank + 1;
                    }
                }
            }
        }
    }
    free (fronts);
}

/* Assign the front numbers to the solutions in h according to solutions
 * in l. The solutions in l are assumed to have the correct front
 * numbers in l [k]->rank already.
 */
static void nd_helper_b
    ( PGAContext *ctx
    , PGAIndividual **l, size_t nl
    , PGAIndividual **h, size_t nh, int m
    )
{
    size_t i;

    assert (m > 0);
    assert (m < ctx->nsga.nfun);
    /* Handle special cases */
    if (nl == 0 || nh == 0) {
        return;
    } else if (nl == 1) {
        /* Compare single L solution to all H solutions */
        /* Note: If l and h got separated into two sets
         * l already dominated h in one objective. So if we find that
         * now l == h for all other objectives, l still dominates h.
         * So we also need to check for equalness if cmp == 0.
         */
        PGAIndividual *l1 = l [0];
        for (i = 0; i < nh; i++) {
            int cmp = compare_solutions (ctx, l1, h [i], m, 1);
            if (cmp < 0) { /* l1 dominates h [i] */
                if (h [i]->rank < l1->rank + 1) {
                    h [i]->rank = l1->rank + 1;
                }
            }
        }
    } else if (nh == 1) {
        /* Compare single H solution to all L solutions
         * See note above.
         */
        PGAIndividual *h1 = h [0];
        for (i = 0; i < nl; i++) {
            int cmp = compare_solutions (ctx, l [i], h1, m, 1);
            if (cmp < 0) { /* l [i] dominates h1 */
                if (h1->rank < l [i]->rank + 1) {
                    h1->rank = l [i]->rank + 1;
                }
            }
        }
    } else if (m == 2) {
        /* 2D ranking case */
        rank_2d_b (ctx, l, nl, h, nh);
    } else {
        double l_max_m, l_min_m, h_max_m, h_min_m;

        l_max_m = l_min_m = GETEVAL_EV (l [0], m-1);
        h_max_m = h_min_m = GETEVAL_EV (h [0], m-1);
        for (i = 1; i < nl; i++) {
            double val = GETEVAL_EV (l [i], m-1);
            int maxcmp = OPT_DIR_CMP_EV (ctx, val, l_max_m);
            int mincmp = OPT_DIR_CMP_EV (ctx, val, l_min_m);
            if (maxcmp > 0) l_max_m = val;
            if (mincmp < 0) l_min_m = val;
        }
        for (i = 1; i < nh; i++) {
            double val = GETEVAL_EV (h [i], m-1);
            int maxcmp = OPT_DIR_CMP_EV (ctx, val, h_max_m);
            int mincmp = OPT_DIR_CMP_EV (ctx, val, h_min_m);
            if (maxcmp > 0) h_max_m = val;
            if (mincmp < 0) h_min_m = val;
        }

        if (OPT_DIR_CMP_EV (ctx, l_max_m, h_min_m) <= 0) {
            /* Objective m can be ignored, recursive call with m reduced */
            nd_helper_b (ctx, l, nl, h, nh, m-1);
        } else {
            /* Check if L and H overlap */
            if (OPT_DIR_CMP_EV (ctx, l_min_m, h_max_m) <= 0) {
                PGAIndividual **l1, **l2, **h1, **h2;
                size_t nl1 = 0, nl2 = 0, nh1 = 0, nh2 = 0;
                /* L and H overlap, split and recurse */
                double x_m_split;
                int alternative;

                if (nl > nh) {
                    x_m_split = find_median (l, nl, m - 1);
                } else {
                    x_m_split = find_median (h, nh, m - 1);
                }
                for (alternative=0;alternative<4;alternative++) {
                    int maxi = ctx->ga.optdir == PGA_MAXIMIZE;
                    /* Split L and H and assert parts are smaller */
                    split_set
                        (l, nl, &l1, &nl1, &l2, &nl2, m - 1, x_m_split, maxi);
                    split_set
                        (h, nh, &h1, &nh1, &h2, &nh2, m - 1, x_m_split, maxi);
                    /* At least one set got smaller */
                    if ((nl1 > 0 && nl2 > 0) || (nh1 > 0 && nh2 > 0)) {
                        break;
                    }
                    switch (alternative) {
                    case 0:
                        /* Split smaller */
                        if (nl > nh) {
                            x_m_split = find_median (h, nh, m - 1);
                        } else {
                            x_m_split = find_median (l, nl, m - 1);
                        }
                        break;
                    case 1:
                        /* Split left max */
                        x_m_split = l_max_m;;
                        break;
                    case 2:
                        /* Split right min */
                        x_m_split = h_min_m;
                        break;
                    case 3:
                        /* Split left min */
                        x_m_split = l_min_m;
                        break;
                    case 4:
                        /* Split right max */
                        x_m_split = h_max_m;
                        break;
                    default:
                        assert (0);
                        break;
                    }
                }

                /* Recursive calls */
                nd_helper_b (ctx, l1, nl1, h1, nh1, m);
                nd_helper_b (ctx, l1, nl1, h2, nh2, m-1);
                nd_helper_b (ctx, l2, nl2, h2, nh2, m);
            }
        }
    }
}

STATIC void rank_2d_a (PGAContext *ctx, PGAIndividual **s, size_t n)
{
    const size_t max_front = ctx->ga.PopSize * 2;
    size_t last_front = 0;
    PGAIndividual **fronts = NULL;
    unsigned int rank = UINT_MAX;
    int i, j;

    assert (n >= 1);
    /* Initialize to zero */
    fronts = malloc (sizeof (*fronts) * max_front);
    if (fronts == NULL) {
        PGAFatalPrintf (ctx, "Out of memory in rank_2d_a");
    }
    memset (fronts, 0, sizeof (*fronts) * max_front);
    /* Sort s by first/second objective (ascending) and rank (descending) */
    qsort (s, n, sizeof (*s), obj_sort_cmp_rank);
    rank = s [0]->rank;
    assert (rank < max_front);
    fronts [rank] = s [0];
    last_front = rank;

    for (i=1; i<(int)n; i++) {
        PGAIndividual *last = fronts [last_front];
        PGAIndividual *cur  = s [i];
        double e1 = GETEVAL_EV (cur, 0);
        double e2 = GETEVAL_EV (cur, 1);
        double e1_last = GETEVAL_EV (last, 0);
        double e2_last = GETEVAL_EV (last, 1);
        int e1cmp = OPT_DIR_CMP_EV (ctx, e1_last, e1);
        int e2cmp = OPT_DIR_CMP_EV (ctx, e2_last, e2);
        /* Current individual dominated by last front? */
        if (e2cmp <= 0 && e2cmp <= 0) {
            /* When both evals are equal the individual is *not* dominated
             * but in that case it needs to have the same eval as the other
             */
            int add = (e2cmp < 0 || e1cmp < 0);
            if (cur->rank <= last->rank + add) {
                cur->rank = last->rank + add;
            }
        } else {
            /* Binary search to find right front
             * Note that we have to take gaps in the fronts into account
             */
            size_t low = 0;
            size_t high = last_front;
            PGAIndividual *lo_ind;
            double e1_lo, e2_lo;

            for (low = 0; fronts [low] == NULL; low++)
                ;
            while (low < high) {
                size_t mid = (low + high) / 2;
                PGAIndividual *i_mid;
                double e1_mid, e2_mid;

                for (; fronts [mid] == NULL; mid--)
                    ;
                i_mid = fronts [mid];
                e1_mid = GETEVAL_EV (i_mid, 0);
                e2_mid = GETEVAL_EV (i_mid, 1);
                e1cmp = OPT_DIR_CMP_EV (ctx, e1_mid, e1);
                e2cmp = OPT_DIR_CMP_EV (ctx, e2_mid, e2);

                if (e2cmp > 0 || (e2cmp == 0 && e1cmp > 0)) {
                    high = mid;
                } else {
                    for (low = mid + 1; fronts [low] == NULL; low++)
                        ;
                }
            }
            /* Skip down to next which should dominate
             * Note unsigned comparison for <= due to wrap
             */
            for (low = low - 1; low <= max_front && fronts [low] == NULL; low--)
                ;
            if (low <= max_front) {
                int add = 0;
                lo_ind = fronts [low];
                e1_lo = GETEVAL_EV (lo_ind, 0);
                e2_lo = GETEVAL_EV (lo_ind, 1);
                e2cmp = OPT_DIR_CMP_EV (ctx, e2_lo, e2);
                e1cmp = OPT_DIR_CMP_EV (ctx, e1_lo, e1);
                assert (e1cmp <= 0 && e2cmp <= 0);
                add = 0;
                /* Strictly lower */
                if (e2cmp < 0 || e1cmp < 0) {
                    add = 1;
                }
                if (cur->rank < fronts [low]->rank + add) {
                    cur->rank = fronts [low]->rank + add;
                }
            }
        }
        rank = cur->rank;
        assert (rank < max_front);
        fronts [rank] = cur;
        if (rank > last_front) {
            last_front = rank;
        }
        /* Need to invalidate lower fronts with larger eval
         * This is the case where stairs are discontinued in the
         * published algorithm
         */
        for (j = rank - 1; j >= 0; j--) {
            double f_e1, f_e2;
            if (fronts [j] == NULL) {
                continue;
            }
            f_e1 = GETEVAL_EV (fronts [j], 0);
            f_e2 = GETEVAL_EV (fronts [j], 1);
            e1cmp = OPT_DIR_CMP_EV (ctx, f_e1, e1);
            e2cmp = OPT_DIR_CMP_EV (ctx, f_e2, e2);
            if (e2cmp > 0 || (e2cmp == 0 && e1cmp <= 0)) {
                fronts [j] = NULL;
            } else {
                break;
            }
        }
    }
    free (fronts);
}

/* Create a non-dominated sorting of s on the first m objectives.
 * The front numbers in s [k]->rank are taken as basis for sorting.
 */
static void nd_helper_a (PGAContext *ctx, PGAIndividual **s, size_t n, int m)
{
    size_t i;

    assert (m >= 2);
    /* Handle base case */
    if (n <= 1) {
        return;
    } else if (n == 2) {
        /* Compare the two solutions */
        int cmp = compare_solutions (ctx, s [0], s [1], m, 0);
        if (cmp < 0) { /* s [0] dominates s [1] */
            if (s [1]->rank < s [0]->rank + 1) {
                s [1]->rank = s [0]->rank + 1;
            }
        } else if (cmp > 0) { /* s [1] dominates s [0] */
            if (s [0]->rank < s [1]->rank + 1) {
                s [0]->rank = s [1]->rank + 1;
            }
        }
    } else if (m == 2) {
        rank_2d_a (ctx, s, n);
    } else {
        /* Check if all values for objective m-1 are identical */
        int all_identical = 1;
        double first_val = GETEVAL_EV (s [0], m-1);

        for (i = 1; i < n && all_identical; i++) {
            if (GETEVAL_EV (s [i], m-1) != first_val) {
                all_identical = 0;
                break;
            }
        }
        if (all_identical) {
            assert (m > 2);
            /* If all values are identical, skip this dimension */
            nd_helper_a (ctx, s, n, m-1);
        } else {
            int maxi = ctx->ga.optdir == PGA_MAXIMIZE;
            PGAIndividual **l, **h;
            size_t nl = 0, nh = 0;

            /* Split the set and recurse */
            double x_m_split = find_median (s, n, m - 1);

            /* Split S into L and H
             * We know that not all values are identical. So if we're
             * getting all elements in one of the split sets we know
             * that more than half the elements are identical to the
             * pivot and that the pivot is the lowest element.
             * This is reversed when we do maximization, so the maxi
             * flags needs to be checked.
             */
            split_set (s, n, &l, &nl, &h, &nh, m - 1, x_m_split, maxi);
            if (nl == 0 || nh == 0) {
                assert ((!maxi && nl == 0) || (maxi && nh == 0));
                split_set (s, n, &l, &nl, &h, &nh, m - 1, x_m_split, !maxi);
            }

            /* Split should have split off something, we asserted above
             * that we have >2 elements
             */
            assert (nl > 0 && nh > 0);

            /* Recursive calls */
            nd_helper_a (ctx, l, nl, m);
            nd_helper_b (ctx, l, nl, h, nh, m-1);
            nd_helper_a (ctx, h, nh, m);
        }
    }
}

static unsigned int max_rank
    (PGAContext *ctx, PGAIndividual **start, size_t n, int goal)
{
    size_t i;
    unsigned int max_rank = 0;
    size_t rankcount = 0;
    const size_t max_front = ctx->ga.PopSize * 2;
    size_t *front_sizes = ctx->scratch.nsga_tmp.front_sizes;

    /* Count individuals in each front and determine max_rank */
    memset (front_sizes, 0, sizeof (*front_sizes) * max_front);
    for (i = 0; i < n; i++) {
        unsigned int rank = start [i]->rank;
        assert (rank < n);
        front_sizes [rank]++;
        if (rank > max_rank) {
            max_rank = rank;
        }
    }
    for (i=0; i<max_rank+1; i++) {
        rankcount += front_sizes [i];
        if (rankcount >= (size_t)goal) {
            max_rank = i;
            break;
        }
    }
    assert (rankcount <= n);

    /* All ranks beyond max_rank need to be set to UINT_MAX */
    for (i = 0; i < n; i++) {
        if (start [i]->rank > max_rank) {
            start [i]->rank = UINT_MAX;
        }
    }

    /* No need for crowding computation if we hit goal exactly */
    if (rankcount == (size_t)goal) {
        return UINT_MAX;
    }

    return max_rank;
}

/* Specialized ranking function for the two-objective case */
static unsigned int ranking_2_objectives
    (PGAContext *ctx, PGAIndividual **start, size_t n, int goal)
{
    /* Call 2d special case */
    rank_2d_a (ctx, start, n);
    /* Compute max_rank and return it */
    return max_rank (ctx, start, n, goal);
}

/* Main non-dominated sorting function */
static unsigned int ranking_3_plus_objectives
    (PGAContext *ctx, PGAIndividual **start, size_t n, int goal)
{
    /* Call the recursive helper function */
    nd_helper_a (ctx, start, n, ctx->nsga.nfun);

    /* Compute max_rank and return it */
    return max_rank (ctx, start, n, goal);
}

/*
 * The is_ev flag decides if we're ranking the evaluation functions or
 * if we're ranking constraint violations, it is 1 for eval functions.
 * The base specifies which is the first aux evaluation. This is 0 if
 * ranking eval function (not constraint violations) and is the index of
 * the first constraint otherwise. The nfun variable is the number of
 * functions to test in each case.
 */
STATIC void set_nsga_state (PGAContext *ctx, int is_ev)
{
    int nc = ctx->ga.NumConstraint;
    int na = ctx->ga.NumAuxEval;
    ctx->nsga.is_ev = is_ev;
    ctx->nsga.nfun = is_ev ? (na - nc + 1) : nc;
    ctx->nsga.base = is_ev ? 0 : (na - nc + 1);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Perform nondominated sorting with O(n**2) algorithm
    \ingroup internal
    \param   ctx             context variable
    \param   start           pointer to individuals
    \param   n               number of individuals
    \param   goal            number of individuals needed in next generation
    \return  Maximum rank given or UINT_MAX if goal was reached exactly
             (in which case no crowding is necessary)

    \rst

    Description
    -----------

    Perform dominance computation known as nondominated sorting or
    ranking. This is the old algorithm which is O(n**2).
    First compute a dominance matrix of N x N bits. The rows are the
    dominated-by relation. We loop over all n^2 pairs of individuals and
    fill the matrix. Init all ranks with -1.  Then starting with rank0:

    - Get all rows of the matrix which are 0 and where the individual has
      no rank yet: These are the currently non-dominated individuals,
      assign the current rank
    - Loop over all individuals with the current rank and remove their
      bits from the dominance matrix
    - Increment the rank counter
    - Stop early when goal is reached

    \endrst
******************************************************************************/
unsigned int PGASortND_NSquare
    (PGAContext *ctx, PGAIndividual **start, size_t n, int goal)
{
    size_t i, j;
    unsigned int rank;
    int nranked = 0;
    size_t intsforn = (n + WL - 1) / WL;
    DECLARE_DYNPTR (PGABinary, dominance, intsforn) =
        (void *)(ctx->scratch.dominance);

    /* Initialize rank, dominance, crowding
     * We initialize crowding here, too because later we compute
     * crowding metric only for the last dominance rank but we sort
     * *all* individuals by crowding.
     */
    for (i=0; i<n; i++) {
        (start [i])->rank     = UINT_MAX;
        (start [i])->crowding = 0;
        for (j=0; j<intsforn; j++) {
            DEREF2_DYNPTR (dominance, intsforn, i, j) = 0;
        }
    }
    for (i=0; i<n; i++) {
        for (j=i+1; j<n; j++) {
            int cmp = compare_solutions
                (ctx, start [i], start [j], ctx->nsga.nfun, 0);
            /* Non-dominated if cmp == 0 */
            /* j dominated by i */
            if (cmp < 0) {
                SET_BIT (DEREF1_DYNPTR (dominance, intsforn, j), i);
            /* i dominated by j */
            } else if (cmp > 0) {
                SET_BIT (DEREF1_DYNPTR (dominance, intsforn, i), j);
            }
        }
    }
    /* Now loop over individuals and establish rank */
    nranked = 0;
    for (rank=0; rank<n; rank++) {
        for (i=0; i<n; i++) {
            if ((*(start+i))->rank != UINT_MAX) {
                continue;
            }
            for (j=0; j<intsforn; j++) {
                if (DEREF2_DYNPTR (dominance, intsforn, i, j)) {
                    break;
                }
            }
            /* Non-dominated in this rank */
            if (j == intsforn) {
                (*(start+i))->rank = rank;
                nranked ++;
            }
        }
        /* Need to rank only goal individuals */
        if (nranked >= goal || (size_t)nranked >= n) {
            break;
        }
        /* Remove dominance bits for this rank */
        for (i=0; i<n; i++) {
            if ((*(start+i))->rank != rank) {
                continue;
            }
            for (j=0; j<n; j++) {
                CLEAR_BIT (DEREF1_DYNPTR (dominance, intsforn, j), i);
            }
        }
    }
    /* No need for crowding computation if we hit goal exactly */
    if (nranked == goal) {
        return UINT_MAX;
    }
    return rank;
}

/*!****************************************************************************
    \brief Perform nondominated sorting, compare old and new implementation
    \ingroup internal
    \param   ctx             context variable
    \param   start           pointer to individuals
    \param   n               number of individuals
    \param   goal            number of individuals needed in next generation
    \return  Maximum rank given or UINT_MAX if goal was reached exactly
             (in which case no crowding is necessary)

    \rst

    Description
    -----------

    Perform dominance computation known as nondominated sorting or
    ranking. This version compares the old O(n**2) implementation against
    the new O(n*(log(n))**m) algorithm (where n is the population size
    and m is the number of objectives). It asserts that the two produce
    identical results.

    \endrst
******************************************************************************/
unsigned int PGASortND_Both
    (PGAContext *ctx, PGAIndividual **start, size_t n, int goal)
{
    DECLARE_DYNARRAY (unsigned int, rank1, n);
    DECLARE_DYNARRAY (unsigned int, rank2, n);
    const size_t max_front = ctx->ga.PopSize * 2;
    size_t *front_sizes = ctx->scratch.nsga_tmp.front_sizes;
    PGAIndividual **cpy_start = ctx->scratch.nsga_tmp.ind_tmp;
    unsigned int max_rank = 0, mr_new = 0, mr_old = 0;
    size_t i;

    if (!n) {
        return 0;
    }
    /* Initialize all ranks to 0 and crowding to 0 */
    for (i=0; i<n; i++) {
        start [i]->rank = 0;
        start [i]->crowding = 0;
    }
    memcpy (cpy_start, start, sizeof (*start) * n);
    mr_old = max_rank = PGASortND_NSquare (ctx, start, n, goal);
    /* Copy ranks and re-initialize to 0 */
    for (i=0; i<n; i++) {
        rank1 [i] = start [i]->rank;
        start [i]->rank = 0;
    }
    assert (ctx->nsga.nfun >= 2);
    if (ctx->nsga.nfun == 2) {
        mr_new = ranking_2_objectives (ctx, cpy_start, n, goal);
    } else {
        mr_new = ranking_3_plus_objectives (ctx, cpy_start, n, goal);
    }
    for (i=0; i<n; i++) {
        rank2 [i] = start [i]->rank;
    }
    if (max_rank == UINT_MAX) {
        int s = 0;
        memset (front_sizes, 0, sizeof (*front_sizes) * max_front);
        for (i=0; i<n; i++) {
            if (start [i]->rank < 2 * (unsigned int)ctx->ga.PopSize) {
                front_sizes [start [i]->rank] += 1;
            }
        }
        for (i=0; i<2u*ctx->ga.PopSize; i++) {
            s += front_sizes [i];
            max_rank = i;
            if (s >= goal) {
                break;
            }
        }
    }
    for (i=0; i<n; i++) {
        /* Dump evaluations in error case */
        if (rank1 [i] != rank2 [i]) {
#define DEBUG_RANKING_DUMP
#ifdef DEBUG_RANKING_DUMP
            int j, k;
            for (j=0; j<(int)n; j++) {
                int n_ev = ctx->ga.NumAuxEval + 1;
                PGAIndividual *ind = start [j];
                printf ("%c{", j == 0 ? '{' : ',');
                for (k=0; k<n_ev; k++) {
                    double ev = (k == 0) ? ind->evalue : ind->auxeval [k - 1];
                    if (n_ev > 3 && k == 0) {
                        printf (" ");
                    }
                    printf ("%.18g", ev);
                    if (k == n_ev - 1) {
                        if (n_ev > 3) {
                            printf ("\n }\n");
                        } else {
                            printf ("}\n");
                        }
                    } else if (n_ev > 3 && ((k + 1) % 3) == 0) {
                        printf ("\n , ");
                    } else {
                        printf (", ");
                    }
                }
            }
            printf ("};\n\n\n");
            fflush (stdout);
#endif /* DEBUG_RANKING_DUMP */
            PGAFatalPrintf (ctx, "Non-dominated sorting of old/new differs");
        }
    }
    if (mr_old != mr_new) {
        PGAFatalPrintf (ctx, "Non-dominated sorting: max_rank differs");
    }
    return mr_new;
}

/*!****************************************************************************
    \brief Perform nondominated sorting using Jensen's algorithm
    \ingroup internal
    \param   ctx             context variable
    \param   start           pointer to individuals
    \param   n               number of individuals
    \param   goal            number of individuals needed in next generation
    \return  Maximum rank given or UINT_MAX if goal was reached exactly
             (in which case no crowding is necessary)

    \rst

    Description
    -----------

    Perform dominance computation known as nondominated sorting or
    ranking. This is the new algorithm originally by Jensen [Jen03]_,
    modified to correctly handle duplicates by Fortin et. al. [FGP13]_
    and (for a slightly modified version) shown to be O(N*log(N)**(M-1))
    by Buzdalov and Shalyto [BS14]_.

    The algorithm selects the appropriate algorithm based on the number
    of objectives (there is a special case for m=2 objectives).

    \endrst
******************************************************************************/
unsigned int PGASortND_Jensen
    (PGAContext *ctx, PGAIndividual **start, size_t n, int goal)
{
    size_t i;
    PGAIndividual **cpy_start = ctx->scratch.nsga_tmp.ind_tmp;
    memcpy (cpy_start, start, sizeof (*start) * n);

    if (!n) {
        return 0;
    }
    /* Initialize all ranks to 0 and crowding to 0 */
    for (i = 0; i < n; i++) {
        start [i]->rank = 0;
        start [i]->crowding = 0;
    }
    if (ctx->nsga.nfun == 2) {
        return ranking_2_objectives (ctx, cpy_start, n, goal);
    } else {
        return ranking_3_plus_objectives (ctx, cpy_start, n, goal);
    }
}

/*!****************************************************************************
    \brief Perform NSGA Replacement
    \ingroup internal
    \param   ctx             context variable
    \param   crowding_method Method used for crowding sort
    \return  None

    \rst

    Description
    -----------

    - Perform dominance computation (ranking)
    - Perform crowding computation specific to the NSGA-Variant given as
      the parameter crowding_method
    - Sort individuals and replace into next generation

    Note that the crowding_method makes the difference between NSGA-II
    and NSGA-III.

    \endrst
******************************************************************************/
static void PGA_NSGA_Replacement (PGAContext *ctx, crowding_t crowding_method)
{
    int i;
    int n_unc_ind, n_con_ind;
    int popsize = PGAGetPopSize (ctx);
    int numreplace = PGAGetNumReplaceValue (ctx);
    PGAIndividual **all_individuals = ctx->scratch.nsga_tmp.ind_all;
    PGAIndividual **constrained = all_individuals + popsize + numreplace;
    PGAIndividual **unconstrained = all_individuals;
    PGAIndividual *oldpop = ctx->ga.oldpop;
    PGAIndividual *newpop = ctx->ga.newpop;
    PGAIndividual *temp;
    size_t         n_dupes = 0;

    PGADebugEntered ("PGA_NSGA_Replacement");

    /* We keep two pointers into the all_individuals array. One with
     * constrained individuals starts from the end. The other with
     * unconstrained individuals starts from the start.
     * Note that with the NSGA algorithms it does not make much sense to
     * have numreplace < popsize: The algorithm is elitist and will keep
     * the good individuals anyway. If numreplace < popsize we do not
     * use the individuals that were copied to the new population,
     * otherwise we would produce duplicates.
     */

    /* Before we put constrained and unconstrained individuals into the
     * start array, if the NoDuplicates flag is set, we sort all the
     * duplicates to the end of the array.
     * A precondition is that neither the individuals in oldpop nor the
     * individuals in newpop contain duplicates. And all the individuals
     * in newpop are already in the hash. So we just loop over oldpop
     * and put all duplicates at the end of the array.
     * if numreplace < popsize we rely on the fact that the first
     * popsize-numreplace individuals are the same in oldpop and newpop.
     * The individuals from newpop are not used in that case, so we do
     * *not* mark the first popsize-numreplace individuals as
     * duplicates.
     */
    if (ctx->ga.NoDuplicates) {
        for (i=popsize - numreplace; i<popsize; i++) {
            if (PGADuplicate (ctx, i, PGA_OLDPOP, PGA_NEWPOP)) {
                constrained--;
                *constrained = oldpop + i;
                n_dupes++;
            }
        }
    }

    /* First loop over all old individuals and put them into the
     * all_individuals array. Note that we compare using the current
     * Epsilon, so this uses Epsilon-Constrained optimization if
     * enabled.
     * Dupes are already handled above, skip those.
     */
    for (i=0; i<popsize; i++) {
        /* Dupes are handled above */
        if (ctx->ga.NoDuplicates && i >= popsize - numreplace) {
            if (PGADuplicate (ctx, i, PGA_OLDPOP, PGA_NEWPOP)) {
                continue;
            }
        }
        if (INDGetAuxTotal (oldpop + i) - ctx->ga.Epsilon > 0) {
            constrained--;
            *constrained = oldpop + i;
        } else {
            *unconstrained = oldpop + i;
            unconstrained++;
        }
    }
    /* Now put all the new individuals into the same array */
    for (i=popsize - numreplace; i<popsize; i++) {
        if (INDGetAuxTotal (newpop + i) - ctx->ga.Epsilon > 0) {
            constrained--;
            *constrained = newpop + i;
        } else {
            *unconstrained = newpop + i;
            unconstrained++;
        }
    }
    /* When array has been filled both pointers point to the same
     * individual somewhere in the middle
     */
    assert (constrained == unconstrained);
    n_unc_ind = unconstrained - all_individuals;
    n_con_ind = popsize + numreplace - n_unc_ind - n_dupes;

    /* First perform non-dominated sorting on unconstrained individuals
     * Normal sorting if only one eval function
     */
    if (ctx->ga.NumConstraint == ctx->ga.NumAuxEval) {
        qsort
            ( all_individuals
            , n_unc_ind
            , sizeof (all_individuals [0])
            , PGAEvalSortHelper
            );
    } else {
        unsigned int rank;
        set_nsga_state (ctx, 1);
        rank = ctx->cops.SortND (ctx, all_individuals, n_unc_ind, popsize);
        if (n_unc_ind >= popsize && rank != UINT_MAX) {
            crowding_method (ctx, all_individuals, n_unc_ind, rank);
        }
        qsort \
            ( all_individuals
            , n_unc_ind
            , sizeof (all_individuals [0])
            , nondom_cmp
            );
    }

    /* Only sort constrained individuals if some of them go into next
     * generation, i.e. if number of unconstrained individuals is below
     * popsize.
     */
    if (n_unc_ind < popsize) {
        assert (ctx->ga.NumConstraint);
        /* Normal sorting if only one constraint function or
         * SumConstraints is PGA_TRUE
         */
        if (ctx->ga.SumConstraints || ctx->ga.NumConstraint < 2) {
            qsort
                ( all_individuals + n_unc_ind
                , n_con_ind
                , sizeof (all_individuals [0])
                , PGAEvalSortHelper
                );
        } else {
            unsigned int rank;
            set_nsga_state (ctx, 0);
            rank = ctx->cops.SortND
                ( ctx
                , all_individuals + n_unc_ind
                , n_con_ind
                , popsize - n_unc_ind
                );
            if (rank != UINT_MAX) {
                crowding (ctx, all_individuals + n_unc_ind, n_con_ind, rank);
            }
            qsort \
                ( all_individuals + n_unc_ind
                , n_con_ind
                , sizeof (all_individuals [0])
                , nondom_cmp
                );
        }
    }
    /* Now we can put all better individuals sorted first into the
     * places of the individuals sorted last, if the source is in
     * PGA_NEWPOP and the target is in PGA_OLDPOP.
     */

    constrained = all_individuals;
    unconstrained = all_individuals + popsize + numreplace;
    while (constrained - all_individuals < popsize) {
        while ((*constrained)->pop == ctx->ga.oldpop) {
            constrained++;
        }
        while ((*(unconstrained-1))->pop == ctx->ga.newpop) {
            unconstrained--;
        }
        if (constrained - all_individuals < popsize) {
            assert (constrained < unconstrained);
            INDCopyIndividual (*constrained, *(unconstrained-1));
            constrained++;
            unconstrained--;
        }
    }
    assert (constrained - all_individuals <= popsize + numreplace);
    assert (constrained - all_individuals >= popsize);

    /* Exchange old/newpop, will be done again by PGAUpdateGeneration,
     * we just make sure that the result of above replacements is now
     * newpop.
     */
    temp           = ctx->ga.oldpop;
    ctx->ga.oldpop = ctx->ga.newpop;
    ctx->ga.newpop = temp;
    PGADebugExited ("PGA_NSGA_Replacement");
}

/*!****************************************************************************
    \brief Perform NSGA-II Replacement.
    \ingroup explicit
    \param   ctx          context variable
    \return  None

    \rst

    Description
    -----------

    - Perform dominance computation (ranking)
    - Perform crowding computation specific to NSGA-II
    - Sort individuals and replace into next generation

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGA_NSGA_II_Replacement (ctx);

    \endrst

******************************************************************************/

void PGA_NSGA_II_Replacement (PGAContext *ctx)
{
    PGA_NSGA_Replacement (ctx, crowding);
}

/*!****************************************************************************
    \brief Perform NSGA-III Replacement.
    \ingroup explicit
    \param   ctx          context variable
    \return  None

    \rst

    Description
    -----------

    - Perform dominance computation (ranking)
    - Perform crowding computation specific to NSGA-III
    - Sort individuals and replace into next generation

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGA_NSGA_II_Replacement (ctx);

    \endrst

******************************************************************************/

void PGA_NSGA_III_Replacement (PGAContext *ctx)
{
    PGA_NSGA_Replacement (ctx, niching);
}
