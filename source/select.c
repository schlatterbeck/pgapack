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
* This file contains the routines that have to do with selection.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

/*!****************************************************************************
    \brief Compute total value over all constraint violations.
    \ingroup standard-api
    \param  ind    Pointer to Individual
    \return Computed or cached total value over all constraint violations

    \rst

    Description
    -----------

    This returns the sum of all *positive* individual aux evaluations
    that are used for constraints.
    The semantics is a total value of all constraint violations.

    Example
    -------

    .. code-block:: c

      PGAIndividual *ind = PGAGetIndividual (ctx, p, PGA_OLDPOP);
      double result;

      ...
      result = INDGetAuxTotal (ind);

    \endrst

******************************************************************************/
double INDGetAuxTotal (PGAIndividual *ind)
{
    PGAContext *ctx = ind->ctx;
    if (!ind->auxtotalok) {
        int i;
        double s = 0;
        int numaux = ctx->ga.NumAuxEval;
        int numcon = ctx->ga.NumConstraint;
        for (i=numaux - numcon; i<numaux; i++) {
            if (ind->auxeval [i] > 0) {
                s += ind->auxeval [i];
            }
        }
        ind->auxtotal   = s;
        ind->auxtotalok = PGA_TRUE;
    }
    return ind->auxtotal;
}

/*!****************************************************************************
    \brief Compute total value over all constraint violations.
    \ingroup query
    \param  ctx    context variable
    \param  p      index of individual
    \param  pop    population
    \return Computed or cached total value over all constraint violations

    \rst

    Description
    -----------

    This returns the sum of all *positive* individual aux evaluations
    that are used for constraints.
    The semantics is a total value of all constraint violations.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      double result;

      ...
      result = PGAGetAuxTotal (ctx, p, PGA_OLDPOP);

    \endrst

******************************************************************************/
double PGAGetAuxTotal (PGAContext *ctx, int p, int pop)
{
    PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
    return INDGetAuxTotal (ind);
}

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* Forward declarations of static functions */
static int  PGASelectLinear       (PGAContext *ctx, PGAIndividual *pop);
static int  PGASelectProportional (PGAContext *ctx, PGAIndividual *pop);
static void PGASelectSUS          (PGAContext *ctx, PGAIndividual *pop);
static int  PGASelectTournament   (PGAContext *ctx, int pop);
static int  PGASelectPTournament  (PGAContext *ctx, int pop);
static int  PGASelectTruncation   (PGAContext *ctx, int pop);

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Perform genetic algorithm selection using the defined
           selection scheme.
    \ingroup explicit
    \param  ctx    context variable
    \param  popix  symbolic constant of population to select from
    \return An array used by \ref PGASelectNextIndex is
            created which contains the population indices of the
            selected individuals

    \rst

    Description
    -----------

    The selection scheme used is either the default selection scheme or
    that specified with :c:func:`PGASetSelectType`.

    Valid selection methods are proportional, stochastic universal,
    tournament, probabilistic tournament selection, truncation
    selection, or linear selection with macros
    :c:macro:`PGA_SELECT_PROPORTIONAL`, :c:macro:`PGA_SELECT_SUS`,
    :c:macro:`PGA_SELECT_TOURNAMENT`, :c:macro:`PGA_SELECT_PTOURNAMENT`,
    :c:macro:`PGA_SELECT_TRUNCATION`, :c:macro:`PGA_SELECT_LINEAR`
    respectively. This function updates an internal array with the
    indices of members of ``popix`` selected for recombination.  These
    indices may be accessed with :c:func:`PGASelectNextIndex`.
    See :ref:`group:const-selection` for the constants and section
    :ref:`sec:selection` in the user guide for details.


    Example
    -------

    .. code-block:: c

      PGAContext *ctx,

      ...
      PGASelect (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
void PGASelect (PGAContext *ctx, int popix)
{
    int i;                   /* not to intefere with dummy argument        */
    int j;                   /* random number                              */
    int temp;                /* for shuffling selected indices US          */
    PGAIndividual *pop;      /* pointer to appropriate population          */

    PGADebugEntered ("PGASelect");

    pop = PGAGetIndividual (ctx, 0, popix);

    switch (ctx->ga.SelectType) {

    case PGA_SELECT_PROPORTIONAL:  /* proportional selection             */
        for (i=0; i<ctx->ga.PopSize; i++) {
            ctx->ga.selected [i] = PGASelectProportional (ctx, pop);
        }
        break;
    case PGA_SELECT_SUS:           /* stochastic universal selection     */
        PGASelectSUS (ctx, pop);
        break;
    case PGA_SELECT_TOURNAMENT:    /* tournament selection               */
        for (i=0; i<ctx->ga.PopSize; i++) {
            ctx->ga.selected [i] = PGASelectTournament (ctx, popix);
        }
        break;
    case PGA_SELECT_PTOURNAMENT:   /* probabilistic tournament selection */
        for (i=0; i<ctx->ga.PopSize; i++) {
            ctx->ga.selected [i] = PGASelectPTournament (ctx, popix);
        }
        break;
    case PGA_SELECT_TRUNCATION:   /* truncation selection */
        for (i=0; i<ctx->ga.PopSize; i++) {
            ctx->ga.selected [i] = PGASelectTruncation (ctx, popix);
        }
        break;
    case PGA_SELECT_LINEAR:      /* linear selection */
        for (i=0; i<ctx->ga.PopSize; i++) {
            ctx->ga.selected [i] = PGASelectLinear (ctx, pop);
        }
        break;
    default:
        PGAError
            ( ctx, "PGASelect: Invalid value of SelectType:"
            , PGA_FATAL, PGA_INT, (void *) &(ctx->ga.SelectType)
            );
        break;
    }

    /* randomize selected string locations
     * Note that for all selection schemes above *except* SUS the items
     * are already randomized. So we randomize only if the selection
     * scheme is SUS *or* we have the backward-compatibility flag
     * ctx->ga.RandomizeSelect set. Note that the point of linear
     * selection is that the sequence is *not* randomized. So in this
     * case we never randomize the sequence.
     */
    if  (  ctx->ga.SelectType == PGA_SELECT_SUS
        || (  ctx->ga.RandomizeSelect
           && ctx->ga.SelectType != PGA_SELECT_LINEAR
           )
        )
    {
        for (i=0; i<ctx->ga.PopSize; i++) {
            j          = PGARandomInterval (ctx, 0, ctx->ga.PopSize-1);
            temp       = ctx->ga.selected [j];
            ctx->ga.selected [j] = ctx->ga.selected [i];
            ctx->ga.selected [i] = temp;
        }
    }

    PGADebugExited ("PGASelect");
}

/*!****************************************************************************
    \brief Return the index of next individual in internal array.
    \ingroup explicit
    \param  ctx     context variable
    \param  popix   the population index, typically PGA_OLDPOP
    \return A population index for the next selected creature

    \rst

    Description
    -----------

    The internal array used contains the indices determined by
    :c:func:`PGASelect`.


    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int l;

      ...
      l = PGASelectNextIndex (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
int PGASelectNextIndex (PGAContext *ctx, int popix)
{
    PGADebugEntered ("PGASelectNextIndex");

    /* We allow the select to consume more than ga.PopSize items
     * If more are consumed we need to prepare the next batch.
     */
    if (ctx->ga.SelectIndex >= ctx->ga.PopSize) {
        PGASelect (ctx, popix);
        ctx->ga.SelectIndex = 0;
    }

    PGADebugExited ("PGASelectNextIndex");
    return ctx->ga.selected [ctx->ga.SelectIndex++];
}

/*!****************************************************************************
    \brief Specify the type of selection to use.
    \ingroup init

    \param   ctx          context variable
    \param   select_type  symbolic constant to specify selection type
    \return  None

    \rst

    Description
    -----------

    Valid choices are :c:macro:`PGA_SELECT_PROPORTIONAL`,
    :c:macro:`PGA_SELECT_SUS`, :c:macro:`PGA_SELECT_TOURNAMENT`,
    :c:macro:`PGA_SELECT_PTOURNAMENT`, :c:macro:`PGA_SELECT_TRUNCATION`,
    and :c:macro:`PGA_SELECT_LINEAR` for proportional, stochastic
    universal selection, tournament, probabilistic tournament selection,
    truncation selection and linear selection, respectively.  The
    default is :c:macro:`PGA_SELECT_TOURNAMENT`. See
    :ref:`group:const-selection` for the constants and section
    :ref:`sec:selection` in the user guide for details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetSelectType (ctx, PGA_SELECT_SUS);

    \endrst

******************************************************************************/
void PGASetSelectType (PGAContext *ctx, int select_type)
{

    PGADebugEntered ("PGASetSelectType");

    switch (select_type) {
        case PGA_SELECT_PROPORTIONAL:
        case PGA_SELECT_SUS:
        case PGA_SELECT_TOURNAMENT:
        case PGA_SELECT_PTOURNAMENT:
        case PGA_SELECT_TRUNCATION:
        case PGA_SELECT_LINEAR:
            ctx->ga.SelectType = select_type;
            break;
        default:
            PGAError
                ( ctx, "PGASetSelectType: Invalid value of select_type:"
                , PGA_FATAL, PGA_INT, (void *) &select_type
                );
        break;
    }

    PGADebugExited ("PGASetSelectType");
}

/*!***************************************************************************
    \brief Return the type of selection selected.
    \ingroup query
    \param   ctx  context variable
    \return  Return the integer corresponding to the symbolic constant
             used to specify the type of selection specified

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int selecttype;

       ...
       selecttype = PGAGetSelectType (ctx);
       switch (selecttype) {
       case PGA_SELECT_PROPORTIONAL:
           printf ("Selection Type = PGA_SELECT_PROPORTIONAL\n");
           break;
       case PGA_SELECT_SUS:
           printf ("Selection Type = PGA_SELECT_SUS\n");
           break;
       case PGA_SELECT_TOURNAMENT:
           printf ("Selection Type = PGA_SELECT_TOURNAMENT\n");
           break;
       case PGA_SELECT_PTOURNAMENT:
           printf ("Selection Type = PGA_SELECT_PTOURNAMENT\n");
           break;
       case PGA_SELECT_TRUNCATION:
           printf ("Selection Type = PGA_SELECT_TRUNCATION\n");
           break;
       }

    \endrst

*****************************************************************************/
int PGAGetSelectType (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetSelectType");
    PGAFailIfNotSetUp ("PGAGetSelectType");

    PGADebugExited ("PGAGetSelectType");

    return ctx->ga.SelectType;
}


/*!****************************************************************************
    \brief Specify the probability that the string that wins a binary
           tournament will be selected.
    \ingroup init
    \param   ctx              context variable
    \param   ptournament_prob the probability of selecting the better string
    \return  None

    \rst

    Description
    -----------

    This function will have no effect unless :c:macro:`PGA_SELECT_PTOURNAMENT`
    was specified as the type of selection to use with
    :c:func:`PGASetSelectType`.  The default value is 0.6.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetPTournamentProb (ctx, 0.8);

    \endrst

******************************************************************************/
void PGASetPTournamentProb (PGAContext *ctx, double ptournament_prob)
{
    PGADebugEntered ("PGASetPTournamentProb");

    ctx->ga.PTournamentProb = ptournament_prob;

    PGADebugExited ("PGASetPTournamentProb");
}

/*!***************************************************************************
    \brief Return the probability of selecting the best string in a
           probabilistic binary tournament.
    \ingroup query
    \param   ctx  context variable
    \return  The probabilistic binary tournament selection probability

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double pt;

       ...
       pt = PGAGetPTournamentProb (ctx);

    \endrst

*****************************************************************************/
double PGAGetPTournamentProb (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetPTournamentProb");
    PGAFailIfNotSetUp ("PGAGetPTournamentProb");

    PGADebugExited ("PGAGetPTournamentProb");

    return ctx->ga.PTournamentProb;
}

/*!****************************************************************************
    \brief Specify the number of participants in a non-probabilistic Tournament.
    \ingroup init
    \param   ctx               context variable
    \param   tournament_size   the size of the tournament
    \return  None

    \rst

    Description
    -----------

    This function will have no effect unless
    :c:macro:`PGA_SELECT_TOURNAMENT` was specified as the type of
    selection to use with :c:func:`PGASetSelectType`.  The default value
    is 2.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetTournamentSize (ctx, 3);

    \endrst

******************************************************************************/
void PGASetTournamentSize (PGAContext *ctx, double tournament_size)
{
    PGADebugEntered ("PGASetTournamentSize");

    ctx->ga.TournamentSize = tournament_size;

    PGADebugExited ("PGASetTournamentSize");
}

/*!***************************************************************************
    \brief Return the number of participants in a tournament
    \ingroup query
    \param   ctx  context variable
    \return  The number of participants in a non probabilistic tournament

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double sz;

       ...
       sz = PGAGetTournamentSize (ctx);

    \endrst

*****************************************************************************/
double PGAGetTournamentSize (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetTournamentSize");
    PGAFailIfNotSetUp ("PGAGetTournamentSize");

    PGADebugExited ("PGAGetTournamentSize");

    return ctx->ga.TournamentSize;
}
/*!****************************************************************************
    \brief Specify if tournament is with or without replacement.
    \ingroup init
    \param   ctx  context variable
    \param   v    flag indicating if replacement is used
    \return  None

    \rst

    Description
    -----------

    This function will have no effect unless
    :c:macro:`PGA_SELECT_TOURNAMENT` was specified as the type of
    selection to use with :c:func:`PGASetSelectType`. The default value
    is :c:macro:`PGA_TRUE` indicating tournament with replacement.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetTournamentWithReplacement (ctx, PGA_FALSE);

    \endrst

******************************************************************************/
void PGASetTournamentWithReplacement (PGAContext *ctx, int v)
{
    ctx->ga.TournamentWithRepl = v;
}

/*!***************************************************************************
    \brief Return the setting for tournament sampling, true if with
           replacement.
    \ingroup query
    \param   ctx  context variable
    \return  The setting of tournament with/without replacement, true
             if with replacement

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int v;

       ...
       v = PGAGetTournamentWithReplacement (ctx);

    \endrst

*****************************************************************************/
int PGAGetTournamentWithReplacement (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetTournamentWithReplacement");
    return ctx->ga.TournamentWithRepl;
}
/*!****************************************************************************
    \brief Specify the proportion of selected individuals for truncation
           selection.
    \ingroup init
    \param   ctx         context variable
    \param   proportion  The value, 0 < proportion <= 1
    \return  None

    \rst

    Description
    -----------

    This function will have no effect unless
    :c:macro:`PGA_SELECT_TRUNCATION` was specified as the type of
    selection to use with :c:func:`PGASetSelectType`. The default value
    is 0.5.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetTruncationProportion (ctx, 0.7);

    \endrst

******************************************************************************/
void PGASetTruncationProportion (PGAContext *ctx, double proportion)
{
    PGAFailIfSetUp ("PGASetTruncationProportion");
    ctx->ga.TruncProportion = proportion;
}

/*!***************************************************************************
    \brief Return the proportion of best individuals selected in
           truncation selection.
    \ingroup query
    \param   ctx  context variable
    \return  Proportion of best individuals selected in truncation selection

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double v;

       ...
       v = PGAGetTruncationProportion (ctx);

    \endrst

*****************************************************************************/
double PGAGetTruncationProportion (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetTruncationProportion");
    return ctx->ga.TruncProportion;
}

/*!****************************************************************************
    \brief Specify if during PGASelect the chosen individuals should be
           randomized again.
    \ingroup init
    \param   ctx  context variable
    \param   v    flag, true if randomized again
    \return  None

    \rst

    Description
    -----------

    All selection schemes except :c:macro:`PGA_SELECT_SUS` already
    return the individuals in randomized order, previously this was
    randomized again. With this method you can re-enable the
    randomization for selection schemes other than
    :c:macro:`PGA_SELECT_SUS` (for which a randomization step is always
    performed).

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetRandomizeSelect (ctx, PGA_TRUE);

    \endrst

******************************************************************************/
void PGASetRandomizeSelect (PGAContext *ctx, int v)
{
    ctx->ga.RandomizeSelect = v;
}

/*!***************************************************************************
    \brief Return the setting for additional select randomization.
    \ingroup query
    \param   ctx  context variable
    \return  The setting of select randomization, true if on

    \rst

    Description
    -----------

    This function will return :c:macro:`PGA_TRUE` if a second
    randomization step after selection is performed. All selection
    schemes except :c:macro:`PGA_SELECT_SUS` already return the
    individuals in randomized order, previously  this  was randomized
    again.  With this method you can find out if the randomization for
    selection schemes other than :c:macro:`PGA_SELECT_SUS` (for which a
    randomization step is always performed)


    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int v;

       ...
       v = PGAGetRandomizeSelect (ctx);

    \endrst

*****************************************************************************/
int PGAGetRandomizeSelect (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetTruncationProportion");
    return ctx->ga.RandomizeSelect;
}



/*!****************************************************************************
    \brief Select a parent for the next generation using a linear search
           through a (fitness) weighted 'roulette wheel'.
    \ingroup internal
    \param  ctx    context variable
    \param  pop    pointer to first individual of population
    \return index of the selected string

    \rst
    .. |_| unicode:: U+00A0 .. Non-breaking space
       :trim:

    Description
    -----------

    The probability of selection of individual :math:`i` with fitness
    :math:`f_i` is given by Goldberg [Gol89]_, p. |_| 11 as
    :math:`p_i = \frac{f_i}{\sum_{i} f_i}`

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      int l;

      ...
      l = PGASelectProportional (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
static int PGASelectProportional (PGAContext *ctx, PGAIndividual *pop)
{
    double sum, sumfitness, r;
    int i;

    PGADebugEntered ("PGASelectProportional");

    sumfitness = 0.0;
    for (i=0; i<ctx->ga.PopSize; i++) {
        sumfitness += (pop+i)->fitness;
    }

    i = 0;
    sum = (pop+i)->fitness;

    r = sumfitness * PGARandom01 (ctx, 0);
    while (r > sum || i==ctx->ga.PopSize) {
        i++;
        sum += (pop+i)->fitness;
    }

    PGADebugExited ("PGASelectProportional");

    return i;
}

/*!****************************************************************************
    \brief A select routine using stochastic universal sampling
    \ingroup internal
    \param  ctx    context variable
    \param  pop    pointer to first individual of population
    \return the array ga.selected [] created via side effect.

    \rst
    .. |_| unicode:: U+00A0 .. Non-breaking space
       :trim:

    Description
    -----------

    Perform stochastic universal sampling selection [Bak87]_, p. |_| 16.
    This routine creates the entire selected population with one call.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,

      ...
      PGASelectSUS (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
static void PGASelectSUS (PGAContext *ctx, PGAIndividual *pop)
{
    int i;
    int k;                          /* index to fill samples array    */
    double davg;                    /* population average fitness     */
    double sum;                     /* running sum of expected values */
    double r;                       /* random number                  */

    PGADebugEntered ("PGASelectSUS");

    /* fill the expected value array */
    davg = 0.0;
    for (i=0; i<ctx->ga.PopSize; i++) {
        davg += (pop+i)->fitness;
    }
    davg /=  (double) ctx->ga.PopSize;
    for (i=0; i<ctx->ga.PopSize; i++) {
        ctx->scratch.dblscratch [i] = (pop+i)->fitness / davg;
    }

    /* select ctx->ga.PopSize as follows */
    sum = 0;
    k   = 0;
    r   = PGARandom01 (ctx, 0);
    for (i=0; i<ctx->ga.PopSize; i++) {
        for (sum+=ctx->scratch.dblscratch [i]; sum>r; r++) {
            ctx->ga.selected [k++] = i;
        }
    }

    PGADebugExited ("PGASelectSUS");
}


/*!****************************************************************************
    \brief Choose N strings randomly and return the one with best evaluation.
    \ingroup internal
    \param  ctx  context variable
    \param  pop  symbolic constant of population to select from
    \return index of the selected string

    \rst

    Description
    -----------

    The configuration parameter N is the value set with
    :c:func:`PGASetTournamentSize`, the default is 2.
    The selection happens *with* replacement.
    This is a generalization of Goldbergs description [Gol89]_,
    for the generalization see, e.g. in Goldberg and Deb [GD91]_.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      int l;

      ...
      l = PGASelectTournamentWithReplacement (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
static
int PGASelectTournamentWithReplacement (PGAContext *ctx, int pop)
{
    int i;
    int m = PGARandomInterval (ctx, 0, ctx->ga.PopSize-1);
    int t = (int)ctx->ga.TournamentSize;
    if (ctx->ga.TournamentSize - t) {
        t += PGARandomFlip (ctx, ctx->ga.TournamentSize - t);
    }
    for (i=1; i<t; i++) {
        int mn = PGARandomInterval (ctx, 0, ctx->ga.PopSize-1);
        /* use '<=' for backwards-compat with prev. binary tournament */
        if (PGAEvalCompare (ctx, mn, pop, m, pop) <= 0) {
            m = mn;
        }
    }
    return m;
}

/*!****************************************************************************
    \brief Choose N strings randomly and return the one with best evaluation.
    \ingroup internal
    \param  ctx  context variable
    \param  pop  symbolic constant of population to select from
    \return index of the selected string

    \rst

    Description
    -----------

    The configuration parameter N is the value set with
    :c:func:`PGASetTournamentSize`, the default is 2.
    The selection happens *without* replacement.
    This means if we select N individuals with a tournament
    size of 2, each individual is participating in exactly two
    tournaments. This does *not* mean that a single individual cannot be
    returned more than once. For implementation notes on the algorithm
    see [GKD89]_, p. 504.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      int l;

      ...
      l = PGASelectTournamentWithoutReplacement (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* Helper function to compute permuted list */
static void _shuffle (PGAContext *ctx, int k, int init)
{
    if (init) {
        int i = 0;
        for (i=0; i<k; i++) {
            ctx->scratch.permute [i] = i;
        }
    }
    PGAShuffle (ctx, ctx->scratch.permute, k);
}
#define NEXT_IDX(ctx, k, init)                        \
    ((ctx)->ga.perm_idx >= (k))                       \
    ? ( _shuffle((ctx), (k), (init))                  \
      , (ctx)->ga.perm_idx = 1                        \
      , (ctx)->scratch.permute [0]                    \
      )                                               \
    : ctx->scratch.permute [ctx->ga.perm_idx++]

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

static
int PGASelectTournamentWithoutReplacement (PGAContext *ctx, int pop)
{
    int m;
    int i;
    int t = (int)ctx->ga.TournamentSize;

    if (ctx->ga.TournamentSize - t) {
        t += PGARandomFlip (ctx, ctx->ga.TournamentSize - t);
    }

    m = NEXT_IDX (ctx, ctx->ga.PopSize, 1);
    assert (0 <= m && m < ctx->ga.PopSize);
    for (i=1; i<t; i++) {
        int mn = NEXT_IDX (ctx, ctx->ga.PopSize, 1);
        assert (0 <= mn && mn < ctx->ga.PopSize);
        if (PGAEvalCompare (ctx, mn, pop, m, pop) <= 0) {
            m = mn;
        }
    }
    return m;
}

/*!****************************************************************************
    \brief Choose all strings that are not already copied to the next
           generation due to elitist strategies.
    \ingroup internal
    \param  ctx    context variable
    \param  pop    pointer to first individual of population
    \return index of the selected string

    \rst

    Description
    -----------

    Note that this 'selection' scheme isn't a selection scheme in the
    genetic sense, it has no selection pressure. Note that the indeces
    are *not* randomized.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      int l;

      ...
      l = PGASelectLinear (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
static int PGASelectLinear (PGAContext *ctx, PGAIndividual *pop)
{
    int numreplace = PGAGetNumReplaceValue (ctx);
    int popsize = PGAGetPopSize (ctx);

    if (ctx->ga.last_iter != ctx->ga.iter || ctx->ga.perm_idx >= popsize) {
        ctx->ga.perm_idx  = popsize - numreplace;
        ctx->ga.last_iter = ctx->ga.iter;
    }
    return ctx->ga.perm_idx++;
}

/*!****************************************************************************
    \brief Choose the best k strings and return them in random order.
    \ingroup internal
    \param  ctx  context variable
    \param  pop  symbolic constant of population to select from
    \return index of the selected string

    \rst

    Description
    -----------

    The value k is (N * TruncationProportion) rounded to the
    next integer.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      int l;

      ...
      l = PGASelectTruncation (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
int PGASelectTruncation (PGAContext *ctx, int pop)
{
    int m = -1;
    int k = (int)(ctx->ga.PopSize * ctx->ga.TruncProportion + 0.5);

    if (k < 1) {
        k = 1;
    }
    if (k > ctx->ga.PopSize) {
        k = ctx->ga.PopSize;
    }

    if (ctx->ga.last_iter != ctx->ga.iter) {
        int i;
        DECLARE_DYNARRAY (int, bestidx, ctx->ga.PopSize);
        /* Returns sorted list of indeces in bestidx, no need to initialize */
        PGAEvalSort (ctx, pop, bestidx);
        /* Copy the k first indeces */
        for (i=0; i<k; i++) {
            ctx->scratch.permute [i] = bestidx [i];
        }
        ctx->ga.last_iter = ctx->ga.iter;
        ctx->ga.perm_idx  = k;
    }

    m = NEXT_IDX (ctx, k, 0);
    assert (0 <= m && m < ctx->ga.PopSize);
    return m;
}

/*!****************************************************************************
    \brief Choose N strings randomly and return the one with best evaluation.
    \ingroup internal
    \param  ctx  context variable
    \param  pop  symbolic constant of population to select from
    \return index of the selected string

    \rst

    Description
    -----------

    The configuration parameter N is the value set with
    :c:func:`PGASetTournamentSize`, the default is 2.
    Depending on the setting of :c:func:`PGASetTournamentWithReplacement`
    calls one of two local functions to use the right sampling.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      int l;

      ...
      l = PGASelectTournament (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
static int PGASelectTournament (PGAContext *ctx, int pop)
{
    PGADebugEntered ("PGASelectTournament");
    if (ctx->ga.TournamentWithRepl) {
        return PGASelectTournamentWithReplacement (ctx, pop);
    } else {
        return PGASelectTournamentWithoutReplacement (ctx, pop);
    }
    PGADebugExited ("PGASelectTournament");
}

/*!****************************************************************************
    \brief Choose two strings randomly and return the one with better
           evaluation with a specified probability.
    \ingroup internal
    \param  ctx  context variable
    \param  pop  symbolic constant of population to select from
    \return index of the selected string

    \rst
    .. |_| unicode:: U+00A0 .. Non-breaking space
       :trim:

    Description
    -----------

    See description in Goldbergs classic book [Gol89]_, p. |_| 121.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      int l;

      ...
      l = PGASelectPTournament (ctx, PGA_OLDPOP);

    \endrst

******************************************************************************/
int PGASelectPTournament (PGAContext *ctx, int pop)
{
    int m1, m2;
    int RetVal;
    double drand;

    PGADebugEntered ("PGASelectPTournament");

    m1 = PGARandomInterval (ctx, 0, ctx->ga.PopSize-1);
    m2 = PGARandomInterval (ctx, 0, ctx->ga.PopSize-1);
    drand = PGARandom01 (ctx, 0);

    if (PGAEvalCompare (ctx, m1, pop, m2, pop) < 0) {
        RetVal = drand < ctx->ga.PTournamentProb ? m1 : m2;
    } else {
        RetVal = drand < ctx->ga.PTournamentProb ? m2 : m1;
    }

    PGADebugExited ("PGASelectPTournament");
    return RetVal;
}
