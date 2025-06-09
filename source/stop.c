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
* This file contains routines related to the stopping conditions for the GA.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

/*!****************************************************************************
    \brief Return true if the stopping conditions have been met.
    \ingroup explicit
    \param   ctx   context variable
    \param   comm  an MPI communicator
    \return  return true if at least one of the termination
             conditions has been met

    \rst

    Description
    -----------

    Calls exactly one of the user defined C or Fortran or system
    :c:func:`PGACheckStoppingConditions` stopping condition functions.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGADone (ctx, comm);

    \endrst

******************************************************************************/
int PGADone (PGAContext *ctx, MPI_Comm comm)
{
    int rank, size, done;

    PGADebugEntered ("PGADone");

    rank = PGAGetRank (ctx, comm);
    size = PGAGetNumProcs (ctx, comm);

    if (rank == 0) {
        if (ctx->fops.StopCond) {
            done = (*ctx->fops.StopCond)(&ctx);
        } else if (ctx->cops.StopCond) {
            done = (*ctx->cops.StopCond)(ctx);
        } else {
            done = PGACheckStoppingConditions (ctx);
        }
    }

    if (size > 1) {
        MPI_Bcast (&done, 1, MPI_INT, 0, comm);
    }

    PGADebugExited ("PGADone");

    return done;
}

/*!****************************************************************************
    \brief Return boolean to indicate if the configured PGAPack
           termination conditions have been met.
    \ingroup standard-api
    \param   ctx   context variable
    \return  return true if at least one of the termination
             conditions has been met

    \rst

    Description
    -----------

    The default termination conditions are given in
    :ref:`group:const-stop` and more details are found in section
    :ref:`sec:stopping-criteria` of the user guide.

    Example
    -------

    Useful in a user-defined function that is registered as a stopping
    condition function. We can use this to keep the builtin stopping
    conditions in addition to a user-defined condition.

    .. code-block:: c

      int StopCond (PGAContext *ctx)
      {
          if (my_stop_check (ctx)) {
              return PGA_TRUE;
          }
          return PGACheckStoppingConditions (ctx);
      }

    \endrst

******************************************************************************/
int PGACheckStoppingConditions (PGAContext *ctx)
{
    int done = PGA_FALSE;

    PGADebugEntered ("PGACheckStoppingConditions");

    /* Since the check happens *after* the generation, test for >= not > */
    if (  ((ctx->ga.StoppingRule & PGA_STOP_MAXITER) == PGA_STOP_MAXITER)
       && (ctx->ga.iter >= ctx->ga.MaxIter)
       )
    {
        done = PGA_TRUE;
    }

    if (  ((ctx->ga.StoppingRule & PGA_STOP_NOCHANGE) == PGA_STOP_NOCHANGE)
       && (ctx->ga.ItersOfSame >= ctx->ga.MaxNoChange)
       )
    {
        done = PGA_TRUE;
    }

    if (  ((ctx->ga.StoppingRule & PGA_STOP_TOOSIMILAR) == PGA_STOP_TOOSIMILAR)
       && (ctx->ga.PercentSame >= ctx->ga.MaxSimilarity)
       )
    {
        done = PGA_TRUE;
    }

    if (ctx->ga.Epsilon > 0) {
        done = PGA_FALSE;
    }

    PGADebugExited ("PGACheckStoppingConditions");
    return done;
}

/*!****************************************************************************
    \brief Specify a stopping criterion.
    \ingroup init

    \param   ctx       context variable
    \param   stoprule  symbolic constant to specify stopping rule
    \return  None

    \rst

    Description
    -----------

    If called more than once the different stopping criterion are ORed
    together.  Valid choices are :c:macro:`PGA_STOP_MAXITER`,
    :c:macro:`PGA_STOP_TOOSIMILAR`, or :c:macro:`PGA_STOP_NOCHANGE` to
    specify iteration limit reached, population too similar, or no change in
    the best solution found in a given number of iterations, respectively.
    The default is to stop when a maximum iteration limit is reached (by
    default, 1000 iterations). The constants can be found in
    :ref:`group:const-stop` and more details are in section
    :ref:`sec:stopping-criteria` of the user guide.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetStoppingRuleType (ctx, PGA_STOP_TOOSIMILAR);

    \endrst

******************************************************************************/
void PGASetStoppingRuleType (PGAContext *ctx, int stoprule)
{

    PGADebugEntered ("PGASetStoppingRuleType");
    PGAFailIfSetUp  ("PGASetStoppingRuleType");

    switch (stoprule) {
        case PGA_STOP_MAXITER  :
        case PGA_STOP_NOCHANGE :
        case PGA_STOP_TOOSIMILAR :
            ctx->ga.StoppingRule |= stoprule;
            break;
        default:
            PGAFatalPrintf
                ( ctx
                , "PGASetStoppingRuleType: Invalid value of stoprule: %d"
                , stoprule
                );
    }

    PGADebugExited ("PGASetStoppingRuleType");
}

/*!***************************************************************************
    \brief Return a symbolic constant that defines the termination criteria.
    \ingroup query
    \param   ctx  context variable
    \return  Return an integer which is an ORed mask of the symbolic constants
             used to specify the stopping rule(s).

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int stop;

       ...
       stop = PGAGetStoppingRuleType (ctx);
       if (stop & PGA_STOP_MAXITER) {
           printf ("Stopping Rule = PGA_STOP_MAXITER\n");
       }
       if (stop & PGA_STOP_NOCHANGE) {
           printf ("Stopping Rule = PGA_STOP_NOCHANGE\n");
       }
       if (stop & PGA_STOP_TOOSIMILAR) {
           printf ("Stopping Rule = PGA_STOP_TOOSIMILAR\n");
       }

    \endrst

*****************************************************************************/
int PGAGetStoppingRuleType (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetStoppingRuleType");
    PGAFailIfNotSetUp ("PGAGetStoppingRuleType");

    PGADebugExited ("PGAGetStoppingRuleType");

    return ctx->ga.StoppingRule;
}

/*!****************************************************************************
    \brief Specify the maximum number of iterations.
    \ingroup init

    \param   ctx     context variable
    \param   maxiter the maximum number of GA iterations to run before stopping
    \return  None

    \rst

    Description
    -----------

    The stopping rule :c:macro:`PGA_STOP_MAXITER` is the default
    stopping rule and is always in effect.
    The default value is 1000 iterations.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMaxGAIterValue (ctx, 5000);

    \endrst

******************************************************************************/
void PGASetMaxGAIterValue (PGAContext *ctx, int maxiter)
{

    PGADebugEntered ("PGASetMaxGAIterValue");
    PGAFailIfSetUp  ("PGASetMaxGAIterValue");

    if (maxiter < 1) {
        PGAFatalPrintf
            ( ctx
            , "PGASetMaxGAIterValue: Invalid value of maxiter: %d"
            , maxiter
            );
    } else {
        ctx->ga.MaxIter = maxiter;
    }

    PGADebugExited ("PGASetMaxGAIterValue");
}

/*!***************************************************************************
    \brief Return the maximum number of iterations to run.
    \ingroup query
    \param   ctx  context variable
    \return  The maximum number of iterations to run

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int maxiter;

       ...
       maxiter = PGAGetMaxGAIterValue (ctx);

    \endrst

*****************************************************************************/
int PGAGetMaxGAIterValue (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetMaxGAIterValue");
    PGAFailIfNotSetUp ("PGAGetMaxGAIterValue");

    PGADebugExited ("PGAGetMaxGAIterValue");

    return ctx->ga.MaxIter;
}

/*!****************************************************************************
    \brief Specify maximum number of iterations of no change in the
           evaluation function value of the best string before stopping.
    \ingroup init
    \param   ctx           context variable
    \param   max_no_change the maximum number of GA iterations allowed
                           with no change in the best evaluation
                           function value
    \return  None

    \rst

    Description
    -----------

    The default value is 100.  The stopping rule
    :c:macro:`PGA_STOP_NOCHANGE` must have been set by
    :c:func:`PGASetStoppingRuleType` for this function
    call to have any effect.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMaxNoChangeValue (ctx, 100);

    \endrst

******************************************************************************/
void PGASetMaxNoChangeValue (PGAContext *ctx, int max_no_change)
{
    PGADebugEntered ("PGASetMaxNoChangeValue");
    PGAFailIfSetUp  ("PGASetMaxNoChangeValue");

    if (max_no_change <= 0) {
        PGAFatalPrintf
            ( ctx
            , "PGASetMaxNoChangeValue: max_no_change invalid: %d"
            , max_no_change
            );
    }

    ctx->ga.MaxNoChange = max_no_change;

    PGADebugExited ("PGASetMaxNoChangeValue");
}

/*!****************************************************************************
    \brief Specify the maximum percent of homogeneity of the population
           before stopping.
    \ingroup init
    \param   ctx             context variable
    \param   max_similarity  the maximum percent of the population that can
                             share the same evaluation function value
    \return  None

    \rst

    Description
    -----------

    The similarity measure is the same evaluation function value. The
    default value is that 95 percent of the population have the same
    evaluation function value. The stopping rule
    :c:macro:`PGA_STOP_TOOSIMILAR` must have been set by
    :c:func:`PGASetStoppingRuleType` for this function call to have any
    effect.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMaxSimilarityValue (ctx, 99);

    \endrst

******************************************************************************/
void PGASetMaxSimilarityValue (PGAContext *ctx, int max_similarity)
{
    PGADebugEntered ("PGASetMaxSimilarityValue");
    PGAFailIfSetUp  ("PGASetMaxSimilarityValue");

    if ((max_similarity <= 0) || (max_similarity > 100)) {
        PGAFatalPrintf
            ( ctx
            , "PGASetMaxSimilarityValue: max_similarity invalid: %d"
            , max_similarity
            );
    }

    ctx->ga.MaxSimilarity = max_similarity;

    PGADebugExited ("PGASetMaxSimilarityValue");
}

/*!****************************************************************************
    \brief Get the maximum percent of homogeneity of the population
           before stopping.
    \ingroup query
    \param   ctx             context variable
    \return  the maximum percent of the population that can share the
             same evaluation function value

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int max_similarity;

       ...
       max_similarity = PGAGetMaxSimilarityValue (ctx);

    \endrst

******************************************************************************/
int PGAGetMaxSimilarityValue (PGAContext *ctx)
{
    return ctx->ga.MaxSimilarity;
}
