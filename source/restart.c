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
* This file contains the routines needed to handle the restart operator,
* and restarting the GA.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

/*!****************************************************************************
    \brief Reseed a population from the best string
    \ingroup explicit
    \param   ctx          context variable
    \param   source_pop   symbolic constant of the source population
    \param   dest_pop     symbolic constant of the destination population
    \return  dest_pop is modified by side-effect

    \rst

    Description
    -----------

    Perform mutation on the best string.
    For integers and reals, the amount by which to change is set with
    :c:func:`PGASetMutationIntegerValue` and
    :c:func:`PGASetMutationRealValue`,
    respectively.  For binary strings, the bits are complemented.

    Example
    -------

    Perform an unspecified test to determine if the current evolution is
    not evolving fast enough, and if so, restart the evolution.

    .. code-block:: c

       PGAContext *ctx;
       ...
       PGAEvaluate (ctx, PGA_OLDPOP, f, comm);
       PGAFitness  (ctx, PGA_OLDPOP);

       ...
       if (StagnantEvolution ()) {
           PGARestart  (ctx, PGA_OLDPOP, PGA_NEWPOP);
           PGAEvaluate (ctx, PGA_NEWPOP, EvalFunc);
           PGAUpdateGeneration (ctx);
       }

    \endrst

******************************************************************************/
void PGARestart (PGAContext *ctx, int source_pop, int dest_pop)
{
    int dest_p, old_mut_type, source_p;
    double val;

    PGADebugEntered ("PGARestart");

    fprintf (ctx->ga.OutputFile, "Restarting the algorithm . . . \n");
    fflush (ctx->ga.OutputFile);
    source_p = PGAGetBestIndex (ctx, source_pop);
    if (source_p != 0 || source_pop != dest_pop) {
        PGACopyIndividual (ctx, source_p, source_pop, 0, dest_pop);
    }
    PGASetEvaluationUpToDateFlag (ctx, 0, dest_pop, PGA_FALSE);
    old_mut_type = PGAGetMutationType (ctx);
    ctx->ga.MutationType = PGA_MUTATION_CONSTANT;
    val = ctx->ga.restartAlleleProb;

    if (ctx->fops.Mutation) {
        for (dest_p = 2; dest_p <= ctx->ga.PopSize; dest_p++) {
            PGACopyIndividual (ctx, 0, dest_pop, dest_p-1, dest_pop);
            (*ctx->fops.Mutation)(&ctx, &dest_p, &dest_pop, &val);
            PGASetEvaluationUpToDateFlag (ctx, dest_p-1, dest_pop, PGA_FALSE);
        }
    } else {
        for (dest_p = 1; dest_p < ctx->ga.PopSize; dest_p++) {
            PGACopyIndividual (ctx, 0, dest_pop, dest_p, dest_pop);
            (*ctx->cops.Mutation)(ctx, dest_p, dest_pop, val);
            PGASetEvaluationUpToDateFlag (ctx, dest_p, dest_pop, PGA_FALSE);
        }
    }
    ctx->ga.MutationType = old_mut_type;

    PGADebugExited ("PGARestart");
}

/*!****************************************************************************
    \brief Specify whether the algorithm should employ the restart operator.
    \ingroup init
    \param   ctx  context variable
    \param   val  boolean variable
    \return  None

    \rst

    Description
    -----------

    By default no restart is performed.

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;

        ...
        PGASetRestartFlag (ctx, PGA_TRUE);

    \endrst

******************************************************************************/
void PGASetRestartFlag (PGAContext *ctx, int val)
{
    PGADebugEntered ("PGASetRestartFlag");

    switch (val) {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.restart = val;
         break;
    default:
         PGAError
            ( ctx, "PGASetRestartFlag: Invalid value for restart:"
            , PGA_FATAL, PGA_INT, (void *) &val
            );
         break;
    }

    PGADebugExited ("PGASetRestartFlag");
}

/*!****************************************************************************
    \brief Return whether the algorithm should employ the restart operator.
    \ingroup query
    \param   ctx  context variable
    \return  true if restarting is enabled

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int val;

       ...
       val = PGAGetRestartFlag (ctx);

    \endrst

******************************************************************************/
int PGAGetRestartFlag (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetRestartFlag");
    PGAFailIfNotSetUp ("PGAGetRestartFlag");

    PGADebugExited ("PGAGetRestartFlag");

    return ctx->ga.restart;
}

/*!****************************************************************************
    \brief Specify the number of iterations of no change in the best
           string after which the algorithm should restart.
    \ingroup init
    \param    ctx      context variable
    \param    numiter  number of changeless iterations
    \return  None

    \rst

    Description
    -----------

    By default no restarts are performed, see
    :c:func:`PGASetRestartFlag`. If restarts are performed, the default
    is after 50 iterations of no change.

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;

        ...
        PGASetRestartFrequencyValue (ctx, 100);

    \endrst

******************************************************************************/
void PGASetRestartFrequencyValue (PGAContext *ctx, int numiter)
{
    PGADebugEntered ("PGASetRestartFrequencyValue");

    if (numiter > 0) {
         ctx->ga.restartFreq = numiter;
    } else {
         PGAError
            ( ctx
            , "PGASetRestartFrequencyValue: Invalid value for restart freqency:"
            , PGA_FATAL, PGA_INT, (void *) &numiter
            );
    }

    PGADebugExited ("PGASetRestartFrequencyValue");
}

/*!****************************************************************************
    \brief Return the number of iterations of no change in the best
           string after which the algorithm should restart.
    \ingroup query
    \param    ctx      context variable
    \return  The number of iteration of no change required for a restart

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        int frq;

        ...
        frq = PGAGetRestartFrequencyValue (ctx);

    \endrst

******************************************************************************/
int PGAGetRestartFrequencyValue (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetRestartFrequencyValue");
    PGAFailIfNotSetUp ("PGAGetRestartFrequencyValue");

    PGADebugExited ("PGAGetRestartFrequencyValue");

    return ctx->ga.restartFreq;
}

/*!****************************************************************************
    \brief Specify the probability with which an allele will be mutated
           during a restart.
    \ingroup init
    \param   ctx   context variable
    \param   prob  probability of mutation
    \return  None

    \rst

    Description
    -----------

    By default the change probability for allele mutations during
    restart is 0.5.

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;

        ...
        PGASetRestartAlleleChangeProb (ctx, 0.7);

    \endrst

******************************************************************************/
void PGASetRestartAlleleChangeProb (PGAContext *ctx, double prob)
{
    PGADebugEntered ("PGASetRestartAlleleChangeProb");

    if (prob >= 0.0 && prob <= 1.0) {
         ctx->ga.restartAlleleProb = prob;
    } else {
         PGAError
            ( ctx, "PGASetRestartAlleleChangeProb: Invalid probability:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &prob
            );
    }

    PGADebugExited ("PGASetRestartAlleleChangeProb");
}

/*!****************************************************************************
    \brief Return the probability with which an allele will be mutated
           during a restart.
    \ingroup query
    \param   ctx  context variable
    \return  The probability of mutating an allele during a restart

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        double prob;

        ...
        prob = PGAGetRestartAlleleChangeProb (ctx);

    \endrst

******************************************************************************/
double PGAGetRestartAlleleChangeProb (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetRestartAlleleChangeProb");
    PGAFailIfNotSetUp ("PGAGetRestartAlleleChangeProb");

    PGADebugExited ("PGAGetRestartAlleleChangeProb");

    return ctx->ga.restartAlleleProb;
}
