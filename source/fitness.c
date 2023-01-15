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
* This file contains the routines that have to do with fitness calculations.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

static void PGAFitnessLinearRank (PGAContext *ctx, int popindex);
static void PGAFitnessLinearNormal (PGAContext *ctx, int popindex);
static void PGAFitnessMinCmax (PGAContext *ctx, PGAIndividual *pop);
static void PGAFitnessMinReciprocal (PGAContext *ctx, PGAIndividual *pop);

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* Helper function for fitness computation */
static void remap_to_positive (PGAContext *ctx, PGAIndividual *pop)
{
    int i;
    double mineval;
    /* put raw fitness into fitness field */

    for (i=0; i<ctx->ga.PopSize; i++) {
        (pop+i)->fitness = (pop+i)->evalue;
    }

    /* translate to all positive sequence (if necessary) */

    mineval = ctx->sys.PGAMaxDouble;
    for (i=0; i<ctx->ga.PopSize; i++) {
        if ((pop+i)->fitness < mineval) {
            mineval = (pop+i)->fitness;
        }
    }
    if (mineval < 0.0) {
        mineval = (-1.01) * mineval;
        for (i=0; i<ctx->ga.PopSize; i++) {
            double fitness = (pop+i)->fitness;
            if (  (mineval != 0 && fitness + mineval == fitness)
               || (fitness != 0 && fitness + mineval == mineval)
               )
            {
                PGAErrorPrintf(ctx, PGA_FATAL
                    , "PGAFitness: Overflow computing positive fitness\n"
                      "Values: mineval=%e, fitness=%e\n"
                      "You can use a ranking variant of fitness computation\n"
                    , mineval, fitness
                    );
            }
            (pop+i)->fitness = fitness + mineval;
        }
    }
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Map the user's evaluation function value to a fitness value.
    \ingroup evaluation
    \param  ctx      context variable
    \param  popindex symbolic constant of the population to calculate
                     fitness for
    \return Calculate the fitness for each string in the population via
            side effect

    \rst

    Description
    -----------

    First, the user's evaluation function value is translated to all
    positive values if any are negative.  Next, this positive sequence
    is translated to a maximization problem if the user's optimization
    direction was minimization.  This positive sequence is then mapped
    to a fitness value using linear ranking, linear normalization
    fitness, or the identity (i.e., the evaluation function value).
    See :ref:`group:const-fitness` in the user guide for allowed
    values. This routine is usually used after :c:func:`PGAEvaluate` is
    called.

    Example
    -------

    Calculate the fitness of all strings in population
    :c:macro:`PGA_NEWPOP` after calling :c:func:`PGAEvaluate` to
    calculate the strings evaluation value.

    .. code-block:: c

       double energy (PGAContext *ctx, int p, int pop, double *aux);
       PGAContext *ctx;
       MPI_Comm comm;

       ...
       PGAEvaluate (ctx, PGA_NEWPOP, energy, comm);
       PGAFitness  (ctx, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGAFitness (PGAContext *ctx, int popindex)
{
    int i;
    /* set pointer to appropriate population, validity check below */
    PGAIndividual *pop =
        (popindex == PGA_OLDPOP) ? ctx->ga.oldpop : ctx->ga.newpop;

    PGADebugEntered ("PGAFitness");

    if (popindex != PGA_OLDPOP && popindex != PGA_NEWPOP) {
        PGAError
            ( ctx, "PGAFitness: Invalid value of popindex:"
            , PGA_FATAL, PGA_INT, (void *) &popindex
            );
    }

    /* make sure all evaluation function values are up-to-date */

    for (i=0; i<ctx->ga.PopSize; i++) {
        if ((pop+i)->evaluptodate != PGA_TRUE) {
            PGAError
                ( ctx, "PGAFitness: evaluptodate not PGA_TRUE for:"
                , PGA_FATAL, PGA_INT, (void *) &i
                );
        }
    }

    /* translate to maximization problem  (if necessary)
     * Only necessary for raw fitness
     */
    if (  ctx->ga.optdir == PGA_MINIMIZE
       && ctx->ga.FitnessType != PGA_FITNESS_RANKING
       )
    {
        switch (ctx->ga.FitnessMinType) {
        case PGA_FITNESSMIN_RECIPROCAL:
            remap_to_positive (ctx, pop);
            PGAFitnessMinReciprocal (ctx, pop);
            break;
        case PGA_FITNESSMIN_CMAX:
            /* Needs no remapping to positive value */
            PGAFitnessMinCmax (ctx, pop);
            break;
        default:
            PGAError( ctx,
                     "PGAFitness: Invalid FitnessMinType:",
                      PGA_FATAL,
                      PGA_INT,
                      (void *) &(ctx->ga.FitnessMinType) );
            break;
        }
    }

    /* last step in fitness calculation */

    switch (ctx->ga.FitnessType) {
    case PGA_FITNESS_RAW:
        if (ctx->ga.optdir != PGA_MINIMIZE) {
            /* minimization already did the remapping */
            remap_to_positive (ctx, pop);
        }
        break;
    case PGA_FITNESS_NORMAL:
        if (ctx->ga.optdir != PGA_MINIMIZE) {
            /* minimization already did the remapping */
            remap_to_positive (ctx, pop);
        }
        PGAFitnessLinearNormal (ctx, popindex);
        break;
    case PGA_FITNESS_RANKING:
        /* Needs no remapping */
        PGAFitnessLinearRank (ctx, popindex);
        break;
    default:
        PGAError( ctx,
                 "PGAFitness: Invalid FitnessType:",
                  PGA_FATAL,
                  PGA_INT,
                  (void *) &(ctx->ga.FitnessType) );
        break;
    }

    PGADebugExited("PGAFitness");
}


/*!****************************************************************************
    \brief Return the rank of a string in a population.
    \ingroup evaluation

    \param  ctx    context variable
    \param  p      the index of the string whose rank is desired
    \param  order  an array containing a unique rank for each string
    \param  n      the size of the array order
    \return The rank of string p

    \rst

    Description
    -----------

    This is a value between 1,...,N (the population size).  The most fit
    string has rank 1, the least fit string has rank N.

    Example
    -------

    Determine the rank of string ``p``.

    .. code-block:: c

      PGAContext *ctx;
      int i, popsize, rank;
      int popsize = PGAGetPopsize (ctx)
      int order [popsize];

      ...
      PGAEvalSort (ctx, pop, order);
      rank = PGARank (ctx, p, order, popsize)

    \endrst

******************************************************************************/
int PGARank (PGAContext *ctx, int p, int *order, int n)
{
    int i;

    PGADebugEntered ("PGARank");

    /*  If the user gives us PGA_TEMP1 or PGA_TEMP2 (or, gasp, some random
     *  number that is not in the population), fail.
     */
    if ((p<0) || (p >= PGAGetPopSize (ctx))) {
        PGAError
            ( ctx, "PGARank: Not a valid population member, p = "
            , PGA_FATAL, PGA_INT, (void *)&p
            );
    }

    /*  Search through all the orderings until we find the one that
     *  matches the given string.  Return the index number.  If we do not
     *  find one, something is _very_ bad; terminate with a fatal error.
     */
    for (i=0; i<n; i++) {
        if (order[i] == p) {
            PGADebugExited ("PGARank");
            return i+1;
        }
    }

    /*  Ideally, we should print out the order array, but, well, ideally,
     *  we should never get here anyway...Also, to make some compilers
     *  shut up, return(0) is here, even though PGAError doesn't return.
     */
    PGAError
        ( ctx, "PGARank: Bottom of loop in rank, p = ", PGA_FATAL
        , PGA_INT, (void *) &p
        );
    return 0;
}

/*!***************************************************************************
    \brief Return the fitness value for a string.
    \ingroup evaluation

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \return  The fitness value for string p in population pop

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p;
       double fit;

       ...
       fit = PGAGetFitness (ctx, p, PGA_NEWPOP);

    \endrst

*****************************************************************************/
double PGAGetFitness (PGAContext *ctx, int p, int pop)
{
    PGAIndividual *ind;

    PGADebugEntered ("PGAGetFitness");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR,"PGAGetFitness", "p = "
        , PGA_INT, (void *) &p
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR,"PGAGetFitness", "pop = "
        , PGA_INT, (void *) &pop
        );

    ind = PGAGetIndividual (ctx, p, pop);

    PGADebugExited ("PGAGetFitness");

    return ind->fitness;
}

/*!***************************************************************************
    \brief Return the type of fitness transformation used.
    \ingroup query

    \param   ctx  context variable
    \return  Returns the integer corresponding to the symbolic constant
             to specify the type of fitness transformation used

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int fittype;

       ...
       fittype = PGAGetFitnessType (ctx);
       switch (fittype) {
       case PGA_FITNESS_RAW:
           printf ("Fitness Type = PGA_FITNESS_RAW\n");
           break;
       case PGA_FITNESS_NORMAL:
           printf ("Fitness Type = PGA_FITNESS_NORMAL\n");
           break;
       case PGA_FITNESS_RANKING:
           printf ("Fitness Type = PGA_FITNESS_RANKING\n");
           break;
       }

    \endrst

*****************************************************************************/
int PGAGetFitnessType (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetFitnessType");
    PGAFailIfNotSetUp ("PGAGetFitnessType");

    PGADebugExited ("PGAGetFitnessType");

    return(ctx->ga.FitnessType);
}

/*!***************************************************************************
    \brief Return the type of fitness transformation used for
           minimization problems.
    \ingroup query

    \param   ctx  context variable
    \return  Returns the integer corresponding to the symbolic constant
             used to specify the type of fitness transformation used for
             minimization problems

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int fitmintype;

       ...
       fitmintype = PGAGetFitnessMinType (ctx);
       switch (fitmintype) {
       case PGA_FITNESSMIN_RECIPROCAL:
           printf ("Fitness Minimization Type = PGA_FITNESSMIN_RECIPROCAL\n");
           break;
       case PGA_FITNESSMIN_CMAX:
           printf ("Fitness Minimization Type = PGA_FITNESSMIN_CMAX\n");
           break;
       }

    \endrst

*****************************************************************************/
int PGAGetFitnessMinType (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetFitnessMinType");
    PGAFailIfNotSetUp ("PGAGetFitnessType");

    PGADebugExited ("PGAGetFitnessMinType");

    return(ctx->ga.FitnessMinType);
}

/*!***************************************************************************
    \brief Return the maximum value used in rank-based fitness.
    \ingroup query

    \param   ctx  context variable
    \return  The value of MAX used in rank-based fitness

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double max;

       ...
       max = PGAGetMaxFitnessRank (ctx);

    \endrst

*****************************************************************************/
double PGAGetMaxFitnessRank (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetMaxFitnessRank");
    PGAFailIfNotSetUp ("PGAGetMaxFitnessRank");

    PGADebugExited ("PGAGetMaxFitnessRank");

    return ctx->ga.FitnessRankMax;
}

/*!****************************************************************************
    \brief Set the type of fitness algorithm to use.
    \ingroup init

    \param   ctx           context variable
    \param   fitness_type  symbolic constant to specify fitness type
    \return  None

    \rst

    Description
    -----------

    Valid choices are :c:macro:`PGA_FITNESS_RAW`,
    :c:macro:`PGA_FITNESS_NORMAL`, or :c:macro:`PGA_FITNESS_RANKING` for
    raw fitness (the evaluation function value), linear normalization,
    or linear ranking, respectively.  The default is
    :c:macro:`PGA_FITNESS_RAW`. See :ref:`group:const-fitness` for the
    constants and section :ref:`sec:evaluation` in the user guide for
    details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetFitnessType (ctx, PGA_FITNESS_RANKING);

    \endrst

******************************************************************************/
void PGASetFitnessType (PGAContext *ctx, int fitness_type)
{

    PGADebugEntered ("PGASetFitnessType");

    switch (fitness_type) {
        case PGA_FITNESS_RAW:
        case PGA_FITNESS_NORMAL:
        case PGA_FITNESS_RANKING:
            ctx->ga.FitnessType = fitness_type;
            break;
        default:
            PGAError
                ( ctx, "PGASetFitnessType: Invalid value of fitness_type:"
                , PGA_FATAL, PGA_INT, (void *) &fitness_type
                );
            break;
    }

    PGADebugExited ("PGASetFitnessType");
}

/*!****************************************************************************
    \brief Set the type of algorithm used if a minimization problem is
           specified to determine how values are remapped for maximization.
    \ingroup init

    \param ctx          context variable
    \param fitness_type symbolic constant to specify fitness minimization type
    \return None

    \rst

    Description
    -----------

    Valid choices are :c:macro:`PGA_FITNESSMIN_RECIPROCAL` and
    :c:macro:`PGA_FITNESSMIN_CMAX` to do the mapping using the
    reciprocal of the evaluation function, or by subtracting the worst
    evaluation function value from each evaluation function value,
    respectively.  The default is :c:macro:`PGA_FITNESSMIN_CMAX`.
    See :ref:`group:const-fitness-min` for the constants and section
    :ref:`sec:evaluation` in the user guide for details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetFitnessMinType (ctx, PGA_FITNESSMIN_CMAX);

   \endrst

******************************************************************************/
void PGASetFitnessMinType (PGAContext *ctx, int fitness_type)
{

    PGADebugEntered ("PGASetFitnessMinType");

    switch (fitness_type) {
        case PGA_FITNESSMIN_RECIPROCAL:
        case PGA_FITNESSMIN_CMAX:
            ctx->ga.FitnessMinType = fitness_type;
            break;
        default:
            PGAError
                ( ctx, "PGASetFitnessMinType: Invalid value of fitness_type:"
                , PGA_FATAL, PGA_INT, (void *) &fitness_type
                );
        break;
    }

    PGADebugExited ("PGASetFitnessMinType");
}

/*!****************************************************************************
    \brief The value of the parameter Max when using linear ranking for
           fitness determination.
    \ingroup init

    \param   ctx  context variable
    \param   max  the value of the parameter Max when using linear ranking
    \return None

    \rst

    Description
    -----------

    The default value is 1.2.  The value must be from the interval
    :math:`[1.0, 2.0]`. The fitness type must have been set to
    :c:macro:`PGA_FITNESS_RANKING` with :c:func:`PGASetFitnessType` for
    this function call to have any effect.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMaxFitnessRank (ctx, 1.1);

    \endrst

******************************************************************************/
void PGASetMaxFitnessRank (PGAContext *ctx, double max)
{
    PGADebugEntered ("PGASetMaxFitnessRank");

    if ((max < 1.0) || (max > 2.0)) {
        PGAError
            ( ctx, "PGASetMaxFitnessRank: Invalid value of max:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &max
            );
    } else {
        ctx->ga.FitnessRankMax = max;
    }

    PGADebugExited ("PGASetMaxFitnessRank");
}

/*!****************************************************************************
    \brief Calculates fitness using a ranking method and linear ordering.
    \ingroup internal

    \param  ctx      context variable
    \param  popindex population index to calculate fitness for
    \return Calculates the fitness for each string in the population via
            side effect

    \rst

    Description
    -----------

    The fitness function is of the form
    :math:`u(x) = K - (\text{rank} \cdot \sigma)` with the constant
    :math:`K` equal to the mean of the evaluation functions, and the
    decrement :math:`\sigma` equal to the standard deviation of the same
    as defined by Davis [Dav91]_, p.\ 33.

    \endrst


******************************************************************************/
static void PGAFitnessLinearNormal (PGAContext *ctx, int popindex)
{

    int i;
    double K, sigma, mean;
    PGAIndividual *pop =
        (popindex == PGA_OLDPOP) ? ctx->ga.oldpop : ctx->ga.newpop;

    PGADebugEntered ("PGAFitnessLinearNormal");

    /* Sort by *eval* (not fitness), no need to init array */
    PGAEvalSort (ctx, popindex, ctx->scratch.intscratch);

    /* calculate parameters for linear normalization */

    for (i=0; i<ctx->ga.PopSize; i++) {
        ctx->scratch.dblscratch [i] = (pop+i)->fitness;
    }
    mean  = PGAMean   (ctx, ctx->scratch.dblscratch, ctx->ga.PopSize);
    sigma = PGAStddev (ctx, ctx->scratch.dblscratch, ctx->ga.PopSize, mean);
    if (sigma == 0) {
         sigma = 1;
    }
    K = sigma * (double) ctx->ga.PopSize;

    for (i=0; i<ctx->ga.PopSize; i++) {
        (pop+i)->fitness =
            K - ( sigma
                * (double)PGARank
                    (ctx,i,ctx->scratch.intscratch,ctx->ga.PopSize)
                );
    }

    PGADebugExited ("PGAFitnessLinearNormal");
}

/*!****************************************************************************
    \brief Calculate fitness using linear ranking.
    \ingroup internal

    \param  ctx      context variable
    \param  popindex population index to calculate fitness for
    \return Calculates the fitness for each string in the population via
            side effect

    \rst
    .. |_| unicode:: U+00A0 .. Non-breaking space
       :trim:


    Description
    -----------

    The fitness function is of the form

    .. math::

        \frac{1}{N} \cdot \left(\max - (\max-\min) * \frac{i-1}{N-1}\right)

    where :math:`\min = 2-\max` and :math:`1 \le \max \le 2`.
    See Baker [Bak87]_, Bäck and Hoffmeister [BH91]_, Grefenstette and
    Baker [GB89]_ and Whitley's ``linear`` function on p. |_| 121 in [Whi89]_.

    \endrst

******************************************************************************/
static void PGAFitnessLinearRank (PGAContext *ctx, int popindex)
{
    double max, min, popsize, rpopsize;
    int *scratch = ctx->scratch.intscratch;
    int i;
    PGAIndividual *pop =
        (popindex == PGA_OLDPOP) ? ctx->ga.oldpop : ctx->ga.newpop;

    PGADebugEntered ("PGAFitnessLinearRank");

    max      = ctx->ga.FitnessRankMax;
    min      = 2. - max;
    popsize  = (double) ctx->ga.PopSize;
    rpopsize = 1.0/popsize;

    /* Sort by *eval* (not fitness), no need to init array */
    PGAEvalSort (ctx, popindex, scratch);

    for(i=0;i<ctx->ga.PopSize;i++) {
        (pop+i)->fitness =
            ( rpopsize
            * ( max
              - ( (max - min)
                * ( ((double)PGARank (ctx,i,scratch,ctx->ga.PopSize) - 1.)
                  / (popsize - 1.)
                  )
                )
              )
            );
    }

    PGADebugExited ("PGAFitnessLinearRank");
}


/*!****************************************************************************
    \brief Calculate fitness in the case of a minimization problem using
           the reciprocal of the evaluation function.
    \ingroup internal
    \param  ctx   context variable
    \param  pop   population pointer to calculate fitness for
    \return Calculates the fitness for each string in the population via
            side effect

    \rst

    Description
    -----------

    This is a power law
    :math:`u(x) = (a f(x) + b)^k` with :math:`a=1, b=0, k=-1`.

    \endrst

******************************************************************************/
static void PGAFitnessMinReciprocal (PGAContext *ctx, PGAIndividual *pop)
{
    int i;

    PGADebugEntered ("PGAFitnessMinReciprocal");

    for (i=0; i<ctx->ga.PopSize; i++) {
        if ((pop+i)->fitness != 0.) {
            (pop+i)->fitness = 1. / (pop+i)->fitness;
        } else {
            PGAError
                ( ctx, "PGAFitnessReciprocal: Value 0.0 for fitness member:"
                , PGA_FATAL, PGA_INT, (void *) &i
                );
        }
    }

    PGADebugExited ("PGAFitnessMinReciprocal");
}


/*!****************************************************************************
    \brief Calculate fitness in the case of a minimization problem by
           subtracting the worst evaluation function value from each
           evaluation function.
    \ingroup internal

    \param  ctx   context variable
    \param  pop   population pointer to calculate fitness for
    \return Calculates the fitness for each string in the population via
            side effect

    \rst

    Description
    -----------

    This is a dynamic linear fitness function
    :math:`u(x) = a f(x) + b(t)` with :math:`a=-1, b(t) = 1.1 * \max f(x)`

    \endrst
******************************************************************************/
static void PGAFitnessMinCmax (PGAContext *ctx, PGAIndividual *pop)
{
    int i;
    double cmax;

    PGADebugEntered ("PGAFitnessMinCmax");

    cmax = 0.;

    for (i=0; i<ctx->ga.PopSize; i++) {
        if ((pop+i)->evalue > cmax) {
            cmax = (pop+i)->evalue;
        }
    }

    cmax *= ctx->ga.FitnessCmaxValue; /* so worst string has nonzero fitness */

    for (i=0; i<ctx->ga.PopSize; i++) {
        /* Check that we're not mapping distinct evaluations to the same
         * fitness
         */
        double evalue = (pop+i)->evalue;
        if (cmax == cmax - evalue && evalue != 0) {
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "PGAFitness: Overflow computing cmax fitness\n"
                  "Values: cmax=%e, evalue=%e\n"
                  "You can use a ranking variant of fitness computation\n"
                , cmax, evalue
                );
        }
        (pop+i)->fitness = cmax - evalue;
    }

    PGADebugExited ("PGAFitnessMinCmax");
}


/*!****************************************************************************
    \brief Set value of the multiplier used by MinCmax fitness algorithm
           so that the worst string has a nonzero fitness.
    \ingroup init

    \param   ctx  context variable
    \param   val  the value of the multiplier
    \return  None

    \rst

    Description
    -----------

    The default value is 1.01.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetFitnessCmaxValue (ctx, 1.2);

    \endrst

******************************************************************************/
void PGASetFitnessCmaxValue (PGAContext *ctx, double val)
{
    PGADebugEntered ("PGASetFitnessCmaxValue");
    ctx->ga.FitnessCmaxValue = val;
    PGADebugExited ("PGASetFitnessCmaxValue");
}

/*!***************************************************************************
    \brief Return the value of the multiplier used in the MinCmax
           fitness algorithm.
    \ingroup query

    \param  ctx  context variable
    \return The value of Cmax used

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double cmax;

       ...
       cmax = PGAGetFitnessCmaxValue (ctx);

    \endrst

*****************************************************************************/
double PGAGetFitnessCmaxValue (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetFitnessCmaxValue");
    PGAFailIfNotSetUp ("PGAGetFitnessCmaxValue");
    PGADebugExited ("PGAGetFitnessCmaxValue");
    return ctx->ga.FitnessCmaxValue;
}
