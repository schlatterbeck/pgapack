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
* This file contains all the routines that are data structure neutral.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* Utility function to reset the hash for all individuals */
static void reset_hash (PGAContext *ctx, int pop)
{
    int i = 0;
    if (ctx->ga.NoDuplicates) {
        PGAIndividual *ind = PGAGetIndividual (ctx, 0, pop);
        size_t hashsize = sizeof (PGAIndividual *) * ctx->ga.PopSize;
        memset (ctx->scratch.hashed, 0, hashsize);
        for (i=0; i<ctx->ga.PopSize; i++, ind++) {
            ind->next_hash = NULL;
        }
    }
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Highest level routine to execute the genetic algorithm.
    \ingroup standard-api
    \param  ctx       context variable
    \param  evaluate  a pointer to the user's evaluation function, which
                      must have the calling sequence shown in the
                      example
    \return None

    \rst

    Description
    -----------

    It is called after :c:func:`PGACreate` and :c:func:`PGASetup` have
    been called.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      double f (PGAContext *ctx, int p, int pop, double *aux);

      ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
      PGASetUp (ctx);
      PGARun (ctx, f);
      PGADestroy (ctx);

    \endrst

******************************************************************************/
void PGARun
    ( PGAContext *ctx
    , double (*evaluate)(PGAContext *c, int p, int pop, double *)
    )
{
     MPI_Comm comm;                  /* value of default communicator */
     int nprocs;                     /* number of processes in above  */
     int npops;                      /* number of populations         */
     int ndemes;                     /* number of demes               */


     PGADebugEntered   ("PGARun");
     PGAFailIfNotSetUp ("PGARun");

     comm   = PGAGetCommunicator (ctx);
     nprocs = PGAGetNumProcs     (ctx, comm);
     npops  = PGAGetNumIslands   (ctx);
     ndemes = PGAGetNumDemes     (ctx);

     /**********************************************************************/
     /*              Global model, one island, one deme                    */
     /**********************************************************************/
     if ((npops == 1) && (ndemes == 1)) {
         PGARunGM (ctx, evaluate, comm);
     }

     /**********************************************************************/
     /*              Island model, > one island, one deme                  */
     /**********************************************************************/
     else if ((npops > 1) && (ndemes == 1)) {
         if (nprocs == 1) {
             PGAError
                ( ctx, "PGARun: island model with one process"
                , PGA_FATAL, PGA_VOID, (void *) &nprocs
                );
         } else if (nprocs != npops) {
             PGAError
                ( ctx, "PGARun: island model no. processes != no. pops"
                , PGA_FATAL, PGA_VOID, (void *) &nprocs
                );
         }
         PGARunIM (ctx, evaluate, comm);
     }

     /**********************************************************************/
     /*              Neighborhood model, one island, > one deme            */
     /**********************************************************************/
     else if ((npops == 1) && (ndemes > 1)) {
         if (nprocs == 1) {
             PGAError
                ( ctx, "PGARun: neighborhood model with one process"
                , PGA_FATAL, PGA_VOID, (void *) &nprocs
                );
         } else if (nprocs != ndemes) {
             PGAError
                ( ctx, "PGARun: neighborhood model no. processes != no. demes"
                , PGA_FATAL, PGA_VOID, (void *) &nprocs
                );
         }
         PGARunNM (ctx, evaluate, comm);
     }

     /**********************************************************************/
     /*              Mixed model, > one island, > one deme                 */
     /**********************************************************************/
     else if ((npops > 1) && (ndemes > 1)) {
         PGAError
            ( ctx, "PGARun: Cannot execute mixed models"
            , PGA_FATAL, PGA_VOID, (void *) &nprocs
            );
     }

     /**********************************************************************/
     /*                        E R R O R                                   */
     /**********************************************************************/
     else {
         PGAError
            ( ctx
            , "PGARun: Invalid combination of numislands, ndemes, and nprocs"
            , PGA_FATAL, PGA_VOID, (void *) &nprocs
            );
     }

     /**********************************************************************/
     /*                         E X I T                                    */
     /**********************************************************************/
     PGADebugExited ("PGARun");
}


/*!****************************************************************************
    \brief Perform crossover and mutation from one population to create
           the next.
    \ingroup explicit
    \param  ctx     context variable
    \param  oldpop  symbolic constant of old population
    \param  newpop  symbolic constant of new population
    \return newpop is modified by side-effect

    \rst

    Description
    -----------

    Assumes :c:func:`PGASelect` has been called.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGARunMutationAndCrossover (ctx, PGA_OLDPOP, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGARunMutationAndCrossover (PGAContext *ctx, int oldpop, int newpop)
{
    int i, j, n, m1, m2;
    int popsize, numreplace;
    double pc;

    PGADebugEntered ("PGARunMutationAndCrossover");

    popsize = PGAGetPopSize (ctx);
    numreplace = PGAGetNumReplaceValue (ctx);
    reset_hash (ctx, newpop);
    /*** first, copy n best strings to new pop ***/
    /*** Note that we do not need to do this for PGA_POPREPL_RTR   ***/
    /*** And neither for PGA_POPREPL_PAIRWISE_BEST                 ***/
    /*** And neither for PGA_POPREPL_NSGA_II                       ***/
    n = popsize - numreplace;
    if (  ctx->ga.PopReplace != PGA_POPREPL_RTR
       && ctx->ga.PopReplace != PGA_POPREPL_PAIRWISE_BEST
       && ctx->ga.PopReplace != PGA_POPREPL_NSGA_II
       )
    {
        PGASortPop (ctx, oldpop);
        for (i=0; i < n; i++) {
            j = PGAGetSortedPopIndex (ctx, i);
            PGACopyIndividual (ctx, j, oldpop, i, newpop);
            PGAHashIndividual (ctx, i, newpop);
        }
    }
    pc = PGAGetCrossoverProb (ctx);
    /*** reproduce to create the rest of the new population ***/
    while (n < popsize) {
        m1 = PGASelectNextIndex (ctx, oldpop);
        m2 = PGASelectNextIndex (ctx, oldpop);
        if (PGARandomFlip (ctx, pc)) {
            PGACrossover (ctx, m1, m2, oldpop, PGA_TEMP1, PGA_TEMP2, newpop);

            /*** mutate and copy first string to new population ***/
            PGAMutate (ctx, PGA_TEMP1, newpop);
            while (PGADuplicate (ctx, PGA_TEMP1, newpop, newpop)) {
                PGAChange (ctx, PGA_TEMP1, newpop);
            }
            PGACopyIndividual (ctx, PGA_TEMP1, newpop, n, newpop);
            PGAHashIndividual (ctx, n, newpop);
            n++;

            if (n < popsize) {
                /*** mutate and copy second string to new population ***/
                PGAMutate (ctx, PGA_TEMP2, newpop);
                while (PGADuplicate (ctx, PGA_TEMP2, newpop, newpop)) {
                    PGAChange (ctx, PGA_TEMP2, newpop);
                }
                PGACopyIndividual (ctx, PGA_TEMP2, newpop, n, newpop);
                PGAHashIndividual (ctx, n, newpop);
                n++;
            }
        } else {
            PGACopyIndividual (ctx, m1, oldpop, n, newpop);
            if (ctx->ga.MixingType == PGA_MIX_TRADITIONAL) {
                PGAMutate (ctx, n, newpop);
            }
            while (PGADuplicate (ctx, n, newpop, newpop)) {
                PGAChange (ctx, n, newpop);
            }
            PGAHashIndividual (ctx, n, newpop);
            n++;
            if (n < ctx->ga.PopSize) {
                PGACopyIndividual (ctx, m2, oldpop, n, newpop);
                if (ctx->ga.MixingType == PGA_MIX_TRADITIONAL) {
                    PGAMutate (ctx, n, newpop);
                }
                while (PGADuplicate (ctx, n, newpop, newpop)) {
                    PGAChange (ctx, n, newpop);
                }
                PGAHashIndividual (ctx, n, newpop);
                n++;
            }
        }
    }

    PGADebugExited ("PGARunMutationAndCrossover");
}


/*!****************************************************************************
    \brief Perform crossover or mutation (but not both) from one
           population to create the next.
    \ingroup explicit
    \param  ctx     context variable
    \param  oldpop  symbolic constant of old population
    \param  newpop  symbolic constant of new population
    \return newpop is modified by side-effect

    \rst

    Description
    -----------

    Assumes :c:func:`PGASelect` has been called.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGARunMutationOrCrossover (ctx, PGA_OLDPOP, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGARunMutationOrCrossover (PGAContext *ctx, int oldpop, int newpop)
{
    int i, j, n, m1, m2;
    int popsize, numreplace;
    double pc;

    PGADebugEntered ("PGARunMutationOrCrossover");

    popsize = PGAGetPopSize (ctx);
    numreplace = PGAGetNumReplaceValue (ctx);
    reset_hash (ctx, newpop);
    /*** first, copy n best strings to new pop ***/
    /*** Note that we do not need to do this for PGA_POPREPL_RTR   ***/
    /*** And neither for PGA_POPREPL_PAIRWISE_BEST                 ***/
    /*** And neither for PGA_POPREPL_NSGA_II                       ***/
    n = popsize - numreplace;
    if (  ctx->ga.PopReplace != PGA_POPREPL_RTR
       && ctx->ga.PopReplace != PGA_POPREPL_PAIRWISE_BEST
       && ctx->ga.PopReplace != PGA_POPREPL_NSGA_II
       )
    {
        PGASortPop (ctx, oldpop);
        for (i=0; i < n; i++) {
            j = PGAGetSortedPopIndex (ctx, i);
            PGACopyIndividual (ctx, j, oldpop, i, newpop);
            PGAHashIndividual (ctx, i, newpop);
        }
    }
    pc = PGAGetCrossoverProb (ctx);
    /*** reproduce to create the rest of the new population ***/
    while (n < popsize) {
        m1 = PGASelectNextIndex (ctx, oldpop);
        m2 = PGASelectNextIndex (ctx, oldpop);
        if (PGARandomFlip (ctx, pc)) {
            PGACrossover (ctx, m1, m2, oldpop, PGA_TEMP1, PGA_TEMP2, newpop);

            /*** copy first string to new population ***/
            while (PGADuplicate(ctx, PGA_TEMP1, newpop,  newpop)) {
                PGAChange (ctx, PGA_TEMP1, newpop);
            }
            PGACopyIndividual (ctx, PGA_TEMP1, newpop, n, newpop);
            PGAHashIndividual (ctx, n, newpop);
            n++;

            if (n < popsize) {
                 /*** copy second string to new population ***/
                 while (PGADuplicate(ctx, PGA_TEMP2, newpop,  newpop)) {
                     PGAChange (ctx, PGA_TEMP2, newpop);
                 }
                 PGACopyIndividual (ctx, PGA_TEMP2, newpop, n, newpop);
                 PGAHashIndividual (ctx, n, newpop);
                 n++;
            }
        } else {
             PGACopyIndividual (ctx, m1, oldpop, PGA_TEMP1, newpop);
             PGAMutate (ctx, PGA_TEMP1, newpop);
             while (PGADuplicate (ctx, PGA_TEMP1, newpop, newpop)) {
                 PGAChange (ctx, PGA_TEMP1, newpop);
             }
             PGACopyIndividual (ctx, PGA_TEMP1, newpop, n, newpop);
             PGAHashIndividual (ctx, n, newpop);
             n++;

             if (n < popsize) {
                 PGACopyIndividual (ctx, m2, oldpop, PGA_TEMP2, newpop);
                 PGAMutate (ctx, PGA_TEMP2, newpop);
                 while (PGADuplicate (ctx, PGA_TEMP2, newpop, newpop)) {
                     PGAChange (ctx, PGA_TEMP2, newpop);
                 }
                 PGACopyIndividual (ctx, PGA_TEMP2, newpop, n, newpop);
                 PGAHashIndividual (ctx, n, newpop);
                 n++;
             }
        }
    }

    PGADebugExited ("PGARunMutationOrCrossover");
}


/*!****************************************************************************
    \brief Perform only mutation
    \ingroup explicit
    \param  ctx     context variable
    \param  oldpop  symbolic constant of old population
    \param  newpop  symbolic constant of new population
    \return newpop is modified by side-effect

    \rst

    Description
    -----------

    Assumes :c:func:`PGASelect` has been called.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGARunMutationOnly (ctx, PGA_OLDPOP, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGARunMutationOnly (PGAContext *ctx, int oldpop, int newpop)
{
    int i, j, n, m;
    int popsize, numreplace;

    PGADebugEntered ("PGARunMutationOnly");

    popsize = PGAGetPopSize (ctx);
    numreplace = PGAGetNumReplaceValue (ctx);
    reset_hash (ctx, newpop);
    /*** first, copy n best strings to new pop ***/
    /*** Note that we do not need to do this for PGA_POPREPL_RTR   ***/
    /*** And neither for PGA_POPREPL_PAIRWISE_BEST                 ***/
    /*** And neither for PGA_POPREPL_NSGA_II                       ***/
    n = popsize - numreplace;
    if (  ctx->ga.PopReplace != PGA_POPREPL_RTR
       && ctx->ga.PopReplace != PGA_POPREPL_PAIRWISE_BEST
       && ctx->ga.PopReplace != PGA_POPREPL_NSGA_II
       )
    {
        PGASortPop (ctx, oldpop);
        for (i=0; i < n; i++) {
            j = PGAGetSortedPopIndex (ctx, i);
            PGACopyIndividual (ctx, j, oldpop, i, newpop);
            PGAHashIndividual (ctx, i, newpop);
        }
    }
    /*** reproduce with mutation only to create rest of the new population ***/
    while (n < popsize) {
        m = PGASelectNextIndex (ctx, oldpop);
        PGACopyIndividual (ctx, m, oldpop, n, newpop);
        PGAMutate (ctx, n, newpop);
        while (PGADuplicate (ctx, n, newpop, newpop)) {
            PGAChange (ctx, n, newpop);
        }
        PGAHashIndividual (ctx, n, newpop);
        n++;
    }

    PGADebugExited ("PGARunMutationOnly");
}


/*!****************************************************************************
    \brief Update internal data structures for the next genetic
           algorithm iteration, and check if the termination
           conditions, both user and PGAPack, have been met.
    \ingroup explicit
    \param   ctx   context variable
    \param   comm  an MPI communicator
    \return  PGA_TRUE if the genetic algorithm has terminated, otherwise
             PGA_FALSE

    \rst

    Description
    -----------

    This routine must be called by both rank-0 and worker processes at
    the end of each GA generation.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGAUpdateGeneration (ctx, MPI_COMM_WORLD);

    \endrst

******************************************************************************/
void PGAUpdateGeneration (PGAContext *ctx, MPI_Comm comm)
{
    PGAIndividual *temp;
    int rank;

    PGADebugEntered ("PGAUpdateGeneration");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR,"PGAUpdateGeneration"
        , "ga.iter = ", PGA_INT, (void *) &(ctx->ga.iter)
        );

    rank = PGAGetRank (ctx, comm);

    ctx->ga.iter++;

    if (rank == 0) {
        DECLARE_DYNARRAY (double, oldbest, ctx->ga.NumAuxEval + 1);
        /* The three replacement schemes PGA_POPREPL_PAIRWISE_BEST,
         * PGA_POPREPL_RTR, and PGA_POPREPL_NSGA_II all replace some new
         * individuals from PGA_NEWPOP into PGA_OLDPOP. Then
         * PGA_NEWPOP/PGA_OLDPOP are switched (resulting in the current
         * population to be PGA_NEWPOP). They are switched *again* at
         * the end of this function to be ready for the next generation.
         * Note that these functions may not use the fitness because
         * this is not comparable across populations.
         */
        if (ctx->ga.PopReplace == PGA_POPREPL_RTR) {
            /* This replaces the selected individuals into OLDPOP.
             * Then OLDPOP/NEWPOP are exchanged
             */
            PGARestrictedTournamentReplacement (ctx);
        }
        else if (ctx->ga.PopReplace == PGA_POPREPL_PAIRWISE_BEST) {
            /* This compares individual in OLDPOP with the individual
             * with same index in NEWPOP and puts the better one into
             * OLDPOP.
             * Then OLDPOP/NEWPOP are exchanged
             */
            PGAPairwiseBestReplacement (ctx);
        }
        else if (ctx->ga.PopReplace == PGA_POPREPL_NSGA_II) {
            /* This performs nondominated sorting and replaces the best
             * individuals over both populations into OLDPOP.
             * Then OLDPOP/NEWPOP are exchanged
             */
            PGA_NSGA_II_Replacement (ctx);
        } else if (ctx->ga.PopReplace == PGA_POPREPL_NSGA_III) {
            /* This performs nondominated sorting and replaces the best
             * individuals over both populations into OLDPOP.
             * Then OLDPOP/NEWPOP are exchanged
             */
            PGA_NSGA_III_Replacement (ctx);
        }

        if (ctx->rep.PrintOptions & PGA_REPORT_AVERAGE) {
            PGAUpdateAverage(ctx, PGA_NEWPOP);
        }

        if (ctx->rep.PrintOptions & PGA_REPORT_ONLINE) {
            PGAUpdateOnline(ctx, PGA_NEWPOP);
        }

        if (ctx->rep.PrintOptions & PGA_REPORT_OFFLINE) {
            PGAUpdateOffline(ctx, PGA_NEWPOP);
        }


        memcpy (oldbest, ctx->rep.Best, sizeof (oldbest));
        PGAUpdateBest (ctx, PGA_NEWPOP);
        if ((ctx->ga.StoppingRule & PGA_STOP_NOCHANGE) || ctx->ga.restart) {
            double *best = ctx->rep.Best;
            int k;
            int equal = 1;
            for (k=0; k<ctx->ga.NumAuxEval+1;k++) {
                if (! (  (isnan (oldbest [k]) && isnan (best [k]))
                      || oldbest [k] == best [k]
                      )
                   )
                {
                    equal = 0;
                    break;
                }
            }
            if (equal) {
                ctx->ga.ItersOfSame++;
            } else {
                ctx->ga.ItersOfSame = 1;
            }
        }

        if (ctx->ga.StoppingRule & PGA_STOP_TOOSIMILAR) {
            ctx->ga.PercentSame = PGAComputeSimilarity (ctx, PGA_NEWPOP);
        }

        /*  Clear this twice in case the user EOG calls PGASelect.  */
        ctx->ga.SelectIndex = 0;

        if (ctx->fops.EndOfGen) {
            (*ctx->fops.EndOfGen)(&ctx);
        }
        if (ctx->cops.EndOfGen) {
            (*ctx->cops.EndOfGen)(ctx);
        }

        ctx->ga.SelectIndex = 0;
        temp           = ctx->ga.oldpop;
        ctx->ga.oldpop = ctx->ga.newpop;
        ctx->ga.newpop = temp;
    }

    PGADebugExited ("PGAUpdateGeneration");
}


/*!***************************************************************************
    \brief Return the data type used by the given context.
    \ingroup query
    \param   ctx  context variable
    \return  The integer corresponding to the symbolic constant used to
             specify the data type

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int datatype;

       ...
       datatype = PGAGetDataType (ctx);
       switch (datatype) {
       case PGA_DATATYPE_BINARY:
           printf ("Data Type = PGA_DATATYPE_BINARY\n");
           break;
       case PGA_DATATYPE_CHARACTER:
           printf ("Data Type = PGA_DATATYPE_CHARACTER\n");
           break;
       case PGA_DATATYPE_INTEGER:
           printf ("Data Type = PGA_DATATYPE_INTEGER\n");
           break;
       case PGA_DATATYPE_REAL:
           printf ("Data Type = PGA_DATATYPE_REAL\n");
           break;
       case PGA_DATATYPE_USER:
           printf ("Data Type = PGA_DATATYPE_USER\n");
           break;
       }

    \endrst

*****************************************************************************/
int PGAGetDataType (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetDataType");

    PGADebugExited  ("PGAGetDataType");

    return ctx->ga.datatype;
}

/*!***************************************************************************
    \brief Return a symbolic constant that represents the direction of
           optimization.
    \ingroup query
    \param   ctx  context variable
    \return  The integer corresponding to the symbolic constant used to
             specify the  direction of optimization

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int optdir;

       ...
       optdir = PGAGetOptDirFlag (ctx);
       switch (optdir) {
       case PGA_MAXIMIZE:
           printf ("Optimization direction = PGA_MAXIMIZE\n");
           break;
       case PGA_MINIMIZE:
           printf ("Optimization direction = PGA_MINIMIZE\n");
           break;
       }

    \endrst

*****************************************************************************/
int PGAGetOptDirFlag (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetOptDirFlag");

    PGADebugExited  ("PGAGetOptDirFlag");

    return ctx->ga.optdir;
}

/*!***************************************************************************
    \brief Return the string length
    \ingroup query
    \param   ctx  context variable
    \return  The string length

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int stringlen;

       ...
       stringlen = PGAGetStringLength (ctx);

    \endrst

*****************************************************************************/
int PGAGetStringLength (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetStringLength");

    PGADebugExited  ("PGAGetStringLength");

    return ctx->ga.StringLen;
}

/*!***************************************************************************
    \brief Return the number of the current genetic algorithm generation.
    \ingroup query
    \param   ctx  context variable
    \return  The genetic algorithm generation number

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        int g;

        ...
        g = PGAGetGAIterValue (ctx);

    \endrst

*****************************************************************************/
int PGAGetGAIterValue (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetGAIterValue");
    PGAFailIfNotSetUp ("PGAGetGAIterValue");

    PGADebugExited ("PGAGetGAIterValue");

    return ctx->ga.iter;
}

/*!***************************************************************************
    \brief Return the number of function evaluations so far.
    \ingroup query
    \param   ctx  context variable
    \return  The number of function evaluations

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        int g;

        ...
        g = PGAGetEvalCount (ctx);

    \endrst

*****************************************************************************/
int PGAGetEvalCount (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetEvalCount");

    return ctx->rep.nevals;
}

/*!****************************************************************************
    \brief Set a boolean flag to indicate if recombination uses exactly
           one of crossover or mutation on selected strings.
    \ingroup deprecated
    \param   ctx   context variable
    \param   flag  to indicate if mutation uses exactly one of crossover
                   or mutation
    \return  None

    \rst

    Description
    -----------

    Note: This is a legacy interface, use :c:func:`PGASetMixingType` instead.
    Do not use this for new code.

    \endrst

******************************************************************************/
void PGASetMutationOrCrossoverFlag (PGAContext *ctx, int flag)
{
    PGADebugEntered ("PGASetMutationOrCrossoverFlag");

    if (flag) {
        ctx->ga.MixingType = PGA_MIX_MUTATE_OR_CROSS;
    } else {
        ctx->ga.MixingType = PGA_MIX_MUTATE_AND_CROSS;
    }

    PGADebugExited ("PGASetMutationOrCrossoverFlag");
}

/*!****************************************************************************
    \brief Set a boolean flag to indicate if recombination uses both
           crossover and mutation on selected strings.
    \ingroup deprecated
    \param   ctx   context variable
    \param   flag  PGA_TRUE (default) or PGA_FALSE
    \return  None

    \rst

    Description
    -----------

    Note: This is a legacy interface, use :c:func:`PGASetMixingType` instead.
    Do not use this for new code.

    \endrst

******************************************************************************/
void PGASetMutationAndCrossoverFlag (PGAContext *ctx, int flag)
{
    PGADebugEntered ("PGASetMutationAndCrossoverFlag");

    if (flag) {
        ctx->ga.MixingType = PGA_MIX_MUTATE_AND_CROSS;
    } else {
        ctx->ga.MixingType = PGA_MIX_MUTATE_OR_CROSS;
    }

    PGADebugExited ("PGASetMutationAndCrossoverFlag");
}
/*!***************************************************************************
    \brief Return true if mutation only occurs when crossover does not.
    \ingroup deprecated

    \param   ctx  context variable
    \return  Return true if mutation only occurs when crossover does not

    \rst

    Description
    -----------

    Note: This is a legacy interface. If mixing types other than
    :c:macro:`PGA_MIX_MUTATE_OR_CROSS` and
    :c:macro:`PGA_MIX_MUTATE_AND_CROSS` have been set this might return
    wrong values, use :c:func:`PGAGetMixingType` instead.
    Do not use this for new code.

    \endrst

*****************************************************************************/
int PGAGetMutationOrCrossoverFlag (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetMutationOrCrossoverFlag");
    PGAFailIfNotSetUp ("PGAGetMutationOrCrossoverFlag");

    PGADebugExited ("PGAGetMutationOrCrossoverFlag");

    return ( ctx->ga.MixingType == PGA_MIX_MUTATE_OR_CROSS
           ? PGA_TRUE : PGA_FALSE
           );
}

/*!***************************************************************************
    \brief Return true if mutation occurs only when crossover does.
    \ingroup deprecated
    \param   ctx  context variable
    \return Return true if mutation is applied to crossed-over strings

    \rst

    Description
    -----------

    Note: This is a legacy interface. If mixing types other than
    :c:macro:`PGA_MIX_MUTATE_OR_CROSS` and
    :c:macro:`PGA_MIX_MUTATE_AND_CROSS` have been set this might return
    wrong values, use :c:func:`PGAGetMixingType` instead.
    Do not use this for new code.

    \endrst

*****************************************************************************/
int PGAGetMutationAndCrossoverFlag (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetMutationAndCrossoverFlag");
    PGAFailIfNotSetUp ("PGAGetMutationAndCrossoverFlag");

    PGADebugExited ("PGAGetMutationAndCrossoverFlag");

    return ( ctx->ga.MixingType == PGA_MIX_MUTATE_AND_CROSS
           ? PGA_TRUE : PGA_FALSE
           );
}

/*!****************************************************************************
    \brief Set a boolean flag to indicate that recombination uses mutation only.
    \ingroup deprecated

    \param   ctx   context variable
    \param   flag  to indicate if only mutation is used
    \return  None

    \rst

    Description
    -----------

    Note: This is a legacy interface, use :c:func:`PGASetMixingType` instead.
    Note: This will override settings of
    :c:func:`PGASetMutationOrCrossoverFlag` and
    :c:func:`PGASetMutationAndCrossoverFlag` and will set the default
    (:c:func:`PGASetMutationOrCrossoverFlag`) when using
    :c:macro:`PGA_FALSE` as the flag. Do not use this for new code.

    \endrst

******************************************************************************/
void PGASetMutationOnlyFlag (PGAContext *ctx, int flag)
{
    if (flag) {
        ctx->ga.MixingType = PGA_MIX_MUTATE_ONLY;
    } else {
        ctx->ga.MixingType = PGA_MIX_MUTATE_OR_CROSS;
    }
}

/*!***************************************************************************
    \brief Return true if only mutation is used.
    \ingroup deprecated

    \param   ctx  context variable
    \return flag to indicate if only mutation is applied

    \rst

    Description
    -----------

    Note: This is a legacy interface, use :c:func:`PGAGetMixingType` instead.
    Do not use this for new code.

    \endrst

*****************************************************************************/
int PGAGetMutationOnlyFlag (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetMutationOnlyFlag");
    return (ctx->ga.MixingType == PGA_MIX_MUTATE_ONLY ? PGA_TRUE : PGA_FALSE);
}

/*!****************************************************************************
    \brief Set strategy for combining mutation and crossover.
    \ingroup init

    \param   ctx   context variable
    \param   type  Type of Mutation/Crossover combination
    \return  None

    \rst

    Example
    -------

    Set the genetic algorithm to use mutation only.

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMixingType (ctx, PGA_MIX_MUTATE_ONLY);

    \endrst

******************************************************************************/
void PGASetMixingType (PGAContext *ctx, int type)
{
    switch (type) {
    case PGA_MIX_MUTATE_OR_CROSS:
    case PGA_MIX_MUTATE_AND_CROSS:
    case PGA_MIX_MUTATE_ONLY:
    case PGA_MIX_TRADITIONAL:
        ctx->ga.MixingType = type;
        break;
    default:
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetMixingType: Invalid value of type: %d"
            , type
            );
        break;
    }
}

/*!***************************************************************************
    \brief Return the strategy setting for combination of mutation and
           crossover.
    \ingroup query
    \param   ctx  context variable
    \return  Return the mixing type

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int mixtype;

       ...
       mixtype = PGAGetMixingType (ctx);

    \endrst

*****************************************************************************/
int PGAGetMixingType (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetMixingType");
    return ctx->ga.MixingType;
}
