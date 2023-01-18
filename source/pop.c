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
    (ctx->ga.optdir == PGA_MAXIMIZE ? CMP ((e1), (e2)) : CMP ((e2), (e1)))
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
    \return  An inteneral array of indices sorted according to one of
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
        PGAError
            (ctx, "PGASort: Invalid value of pop:"
            , PGA_FATAL, PGA_INT, (void *) &pop
            );
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
            ctx->ga.sorted[i] = ctx->scratch.intscratch[j];
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
        PGAError
            ( ctx, "PGAGetSorted: Invalid value of n:"
            , PGA_FATAL, PGA_INT, (void *) &n
            );
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
        PGAError
            ( ctx, "PGASetPopSize: Invalid value of popsize:"
            , PGA_FATAL, PGA_INT, (void *) &popsize
            );
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
        PGAError
            ( ctx, "PGASetNumReplaceValue: Invalid value of pop_replace:"
            , PGA_FATAL, PGA_INT, (void *) &pop_replace
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
        PGAError
            ( ctx, "PGASetPopReplaceType: Invalid value of pop_replace:"
            , PGA_FATAL, PGA_INT, (void *) &pop_replace
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
        PGAErrorPrintf (ctx, PGA_FATAL, "Can't set reference points twice");
    }
    if (points == NULL) {
        PGAErrorPrintf (ctx, PGA_FATAL, "Need non-NULL points");
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
        PGAErrorPrintf (ctx, PGA_FATAL, "Can't set reference directions twice");
    }
    if (dirs == NULL) {
        PGAErrorPrintf (ctx, PGA_FATAL, "Need non-NULL directions");
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
    int popsize = PGAGetPopSize(ctx);
    int numreplace = PGAGetNumReplaceValue(ctx);
    PGASampleState state;
    PGAIndividual *temp;
    int oldpop = PGA_OLDPOP;
    int newpop = PGA_NEWPOP;

    PGADebugEntered ("PGARestrictedTournamentReplacement");
    for (i=popsize - numreplace; i<popsize; i++) {
        double dist = -1.0;
        int closest = 0;
        PGARandomSampleInit
            (ctx, &state, ctx->ga.RTRWindowSize, ctx->ga.PopSize);
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

        /* If new population individual is better */
        if (PGAEvalCompare (ctx, i, PGA_NEWPOP, closest, PGA_OLDPOP) <= 0) {
            /* Copy i in PGA_NEWPOP to closest in PGA_OLDPOP */
            PGACopyIndividual (ctx, i, PGA_NEWPOP, closest, PGA_OLDPOP);
        }
    }
    /* Exchange old/newpop, will be done again by PGAUpdateGeneration,
     * we just make sure that the result of above replacements is now
     * newpop.
     */
    temp           = ctx->ga.oldpop;
    ctx->ga.oldpop = ctx->ga.newpop;
    ctx->ga.newpop = temp;
    PGADebugExited("PGARestrictedTournamentReplacement");
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
    int popsize = PGAGetPopSize(ctx);
    int numreplace = PGAGetNumReplaceValue(ctx);
    PGAIndividual *temp;

    PGADebugEntered ("PGAPairwiseBestReplacement");
    for (i=popsize - numreplace; i<popsize; i++) {
        /* Note the '<=' comparison, differential evolution can walk across
         * areas with equal evaluation this way
         */
        if (PGAEvalCompare (ctx, i, PGA_NEWPOP, i, PGA_OLDPOP) <= 0) {
            /* Copy i in PGA_NEWPOP to i in PGA_OLDPOP */
            PGACopyIndividual (ctx, i, PGA_NEWPOP, i, PGA_OLDPOP);
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

static int crowdsort_cmp (const void *a1, const void *a2)
{
    PGAIndividual **i1 = (void *)a1;
    PGAIndividual **i2 = (void *)a2;
    int is_ev = INDGetAuxTotal (*i1) ? 0 : 1;
    double e1 = GETEVAL(*i1, (*i1)->funcidx, is_ev);
    double e2 = GETEVAL(*i2, (*i2)->funcidx, is_ev);
    return CMP(e1, e2);
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
    return CMP((*i2)->crowding, (*i1)->crowding);
}

/* typedef to make it easier to pass crowding functions as parameter */
typedef void (* crowding_t)
    (PGAContext *, PGAIndividual **, size_t, unsigned int);

/* Compute crowding distance over the given individuals
 * This is specific to NSGA-II.
 * For explanation of is_ev see ranking
 */
STATIC void crowding
    (PGAContext *ctx, PGAIndividual **start, size_t n, unsigned int rank)
{
    size_t i;
    int k;
    size_t nrank = 0;
    int is_ev = INDGetAuxTotal (*start) ? 0 : 1;
    DECLARE_DYNARRAY (PGAIndividual *, crowd, n);
    int nc = ctx->ga.NumConstraint;
    int na = ctx->ga.NumAuxEval;
    int base = is_ev ? 0 : (na - nc);
    int nfun = is_ev ? (na - nc + 1) : nc;
    DECLARE_DYNARRAY (double, f_min, nfun);
    DECLARE_DYNARRAY (double, f_max, nfun);

    for (i=0; i<n; i++) {
        PGAIndividual *ind = start [i];
        ind->crowding = 0;
        if (ind->rank == rank) {
            for (k=0; k<nfun; k++) {
                double e = GETEVAL (ind, k + base, is_ev);
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
    for (k=0; k<nfun; k++) {
        double norm = f_max [k] - f_min [k];
        for (i=0; i<nrank; i++) {
            (crowd [i])->funcidx = k + base;
        }
        qsort (crowd, nrank, sizeof (crowd [0]), crowdsort_cmp);
        (crowd [0])->crowding = DBL_MAX;
        (crowd [nrank-1])->crowding = DBL_MAX;
        for (i=1; i<nrank-1; i++) {
            if ((crowd [i])->crowding != DBL_MAX) {
                (crowd [i])->crowding +=
                    ( GETEVAL(crowd [i+1], k + base, is_ev)
                    - GETEVAL(crowd [i-1], k + base, is_ev)
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
            if (OPT_DIR_CMP (ctx, e, ctx->ga.utopian [j]) > 0) {
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
STATIC double compute_asf (PGAContext *ctx, double *point, int axis)
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
        ( ctx, PGA_WARNING
        , "Intercept computation failed in Generation %d\n", ctx->ga.iter
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
               || OPT_DIR_CMP (ctx, e, ctx->ga.worst [j]) < 0
               )
            {
                ctx->ga.worst [j] = e;
            }
            if (i==0 || OPT_DIR_CMP (ctx, e, wpop [j]) < 0) {
                wpop [j] = e;
            }
            if (  start [i]->rank == 0
               && (wv <= j || OPT_DIR_CMP (ctx, e, wof0 [j]) < 0)
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
            if (OPT_DIR_CMP (ctx, ctx->ga.nadir [j], ctx->ga.worst [j]) < 0) {
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
    return CMP((*i1)->distance, (*i2)->distance);
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

/* Dominance computation, return the maximum rank given or UINT_MAX if
 * goal was reached exactly (in which case no crowding is necessary)
 * First compute a dominance matrix of N x N bits. The rows are the
 * dominated-by relation. We loop over all n^2 pairs of individuals and
 * fill the matrix. Initit all ranks with -1.
 * Then starting with rank0:
 * - Get all rows of the matrix which are 0 and where the individual has
 *   no rank yet: These are the currently non-dominated individuals,
 *   assign the current rank
 * - Loop over all individuals with the current rank and remove their
 *   bits from the dominance matrix
 * - Increment the rank counter
 * The is_ev flag decides if we're ranking the evaluation functions or
 * if we're ranking constraint violations, it is 0 for eval functions.
 * The base specifies which is the first aux evaluation. This is 0 if
 * ranking eval function (not constraint violations) and is the index of
 * the first constraint otherwise. The nfun variable is the number of
 * functions to test in each case.
 */
STATIC unsigned int ranking
    (PGAContext *ctx, PGAIndividual **start, size_t n, int goal)
{
    size_t i, j;
    int k;
    int is_ev = INDGetAuxTotal (*start) - ctx->ga.Epsilon <= 0;
    unsigned int rank;
    int nranked = 0;
    int nc = ctx->ga.NumConstraint;
    int na = ctx->ga.NumAuxEval;
    int base = is_ev ? 0 : (na - nc);
    int nfun = is_ev ? (na - nc + 1) : nc;
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
            int cmp = 0;
            for (k=0; k<nfun; k++) {
                double e1, e2;
                int ncmp;
                e1 = GETEVAL (start [i], k + base, is_ev);
                e2 = GETEVAL (start [j], k + base, is_ev);
                if (is_ev && ctx->ga.optdir == PGA_MAXIMIZE) {
                    ncmp = CMP (e2, e1);
                } else {
                    ncmp = CMP (e1, e2);
                }
                if (cmp && ncmp && ncmp != cmp) {
                    break;
                }
                /* Don't allow cmp to get zero again */
                if (ncmp) {
                    cmp = ncmp;
                }
            }
            /* Non-dominated? */
            if (!cmp || k<nfun) {
                continue;
            }
            /* j dominated by i */
            if (cmp < 0) {
                SET_BIT (DEREF1_DYNPTR (dominance, intsforn, j), i);
            /* i dominated by j */
            } else {
                SET_BIT (DEREF1_DYNPTR (dominance, intsforn, i), j);
            }
        }
    }
    /* Now repeatedly loop over individuals and establish rank */
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
        if (nranked >= goal) {
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

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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
    DECLARE_DYNARRAY (PGAIndividual *, all_individuals, popsize + numreplace);
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
        rank = ranking (ctx, all_individuals, n_unc_ind, popsize);
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
        if (ctx->ga.SumConstraints || ctx->ga.NumConstraint == 1) {
            qsort
                ( all_individuals + n_unc_ind
                , n_con_ind
                , sizeof (all_individuals [0])
                , PGAEvalSortHelper
                );
        } else {
            unsigned int rank;
            rank = ranking
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
