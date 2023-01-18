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

/*!****************************************************************************
* \file
* This file contains all the parallel functions.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
******************************************************************************/

/*!***************************************************************************
 *  \defgroup parallel Parallel
 *  \brief Parallel implementation of GA
 *****************************************************************************/
/*!***************************************************************************
 *  \defgroup notimplemented Not yet implemented
 *  \brief Not yet implemented, mainly used for island/multiple demes.
 *****************************************************************************/

#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
#define DEBUG_EVAL 0
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief High-level routine to execute the genetic algorithm using the
           global model.
    \ingroup explicit
    \param  ctx       context variable
    \param  evaluate  a pointer to the user's evaluation function, which must
                      have the calling sequence shown in the example
    \param  comm      an MPI communicator
    \return None

    \rst

    Description
    -----------

    It is called after :c:func:`PGACreate` and :c:func:`PGASetUp` have
    been called. If a ``NULL`` communicator is given, a sequential
    execution method is used, otherwise, work is divided among the
    processors in the communicator.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      double f (PGAContext *ctx, int p, int pop, double *aux);

      ...
      PGARunGM (ctx, f, MPI_COMM_WORLD);

    \endrst

******************************************************************************/
void PGARunGM
    ( PGAContext *ctx
    , double (*evaluate)(PGAContext *, int, int, double *)
    , MPI_Comm comm
    )
{
    int       rank, Restarted;
    void    (*CreateNewGeneration) (PGAContext *, int, int) = NULL;

    /*  Let this be warned:
     *  The communicator is NOT duplicated.  There might be problems with
     *  PGAPack and the user program using the same communicator.
     */
    PGADebugEntered ("PGARunGM");

    rank = PGAGetRank (ctx, comm);

    if (rank == 0) {
        if (ctx->fops.PreEval) {
            int pop = PGA_OLDPOP;
            (*ctx->fops.PreEval)(&ctx, &pop);
        }
        if (ctx->cops.PreEval) {
            (*ctx->cops.PreEval)(ctx, PGA_OLDPOP);
        }
    }

    PGAEvaluate (ctx, PGA_OLDPOP, evaluate, comm);
    if (rank == 0) {
        int st = PGAGetSelectType (ctx);
        /* If epsilon constraints are used */
        if (ctx->ga.NumConstraint && ctx->ga.EpsilonGeneration) {
            int idx;
            PGAIndividual *ind;
            /* Sort population by auxiliary eval */
            /* No need to init the indeces, filled in by PGAEvalSort */
            PGAEvalSort (ctx, PGA_OLDPOP, ctx->scratch.intscratch);
            idx = ctx->scratch.intscratch [ctx->ga.EpsilonTheta];
            ind = PGAGetIndividual (ctx, idx, PGA_OLDPOP);
            assert (ind->auxtotalok);
            ctx->ga.Epsilon = ctx->ga.Epsilon_0 = ind->auxtotal;
            if (!ctx->ga.EpsilonExponent) {
                double l10 = log (10);
                assert (ctx->ga.EffEpsExponent == 0);
                ctx->ga.EffEpsExponent =
                    (-5 - log (ctx->ga.Epsilon) / l10) / (log (0.05) / l10);
                if (ctx->ga.EffEpsExponent < PGA_EPSILON_EXPONENT_MIN) {
                    ctx->ga.EffEpsExponent = PGA_EPSILON_EXPONENT_MIN;
                }
                if (ctx->ga.EffEpsExponent > PGA_EPSILON_EXPONENT_MAX) {
                    ctx->ga.EffEpsExponent = PGA_EPSILON_EXPONENT_MAX;
                }
            }
        }
        PGAUpdateBest (ctx, PGA_OLDPOP);
        if (st == PGA_SELECT_SUS || st == PGA_SELECT_PROPORTIONAL) {
            PGAFitness (ctx, PGA_OLDPOP);
        }
    }

    switch (PGAGetMixingType (ctx)) {
    case PGA_MIX_MUTATE_OR_CROSS:
        CreateNewGeneration = PGARunMutationOrCrossover;
        break;
    case PGA_MIX_MUTATE_AND_CROSS:
        CreateNewGeneration = PGARunMutationAndCrossover;
        break;
    case PGA_MIX_MUTATE_ONLY:
        CreateNewGeneration = PGARunMutationOnly;
        break;
    case PGA_MIX_TRADITIONAL:
        CreateNewGeneration = PGARunMutationAndCrossover;
        break;
    default:
        assert (0);
    }

    while (!PGADone(ctx, comm)) {
        if (rank == 0) {
            Restarted = PGA_FALSE;
            if ((ctx->ga.restart == PGA_TRUE) &&
                (ctx->ga.ItersOfSame % ctx->ga.restartFreq == 0)) {
                ctx->ga.ItersOfSame++;
                Restarted = PGA_TRUE;
                PGARestart(ctx, PGA_OLDPOP, PGA_NEWPOP);
            } else {
                PGASelect(ctx, PGA_OLDPOP);
                CreateNewGeneration(ctx, PGA_OLDPOP, PGA_NEWPOP);
                if (ctx->fops.PreEval) {
                    int pop = PGA_NEWPOP;
                    (*ctx->fops.PreEval)(&ctx, &pop);
                }
                if (ctx->cops.PreEval) {
                    (*ctx->cops.PreEval)(ctx, PGA_NEWPOP);
                }
            }
        }
        MPI_Bcast (&Restarted, 1, MPI_INT, 0, comm);

        PGAEvaluate (ctx, PGA_NEWPOP, evaluate, comm);
        if (rank == 0) {
            int st = PGAGetSelectType (ctx);
            if (st == PGA_SELECT_SUS || st == PGA_SELECT_PROPORTIONAL) {
                PGAFitness (ctx, PGA_NEWPOP);
            }
            /* If epsilon constraints are used */
            if (ctx->ga.NumConstraint && ctx->ga.EpsilonGeneration) {
                /* Smaller Exponent after generation EpsTLambda */
                if (ctx->ga.EpsTLambda && ctx->ga.iter == ctx->ga.EpsTLambda) {
                    assert (ctx->ga.EpsilonExponent == 0);
                    ctx->ga.EffEpsExponent = 0.3 * ctx->ga.EffEpsExponent
                                           + 0.7 * PGA_EPSILON_EXPONENT_MIN;
                }
                if (ctx->ga.iter >= ctx->ga.EpsilonGeneration) {
                    ctx->ga.Epsilon = 0;
                } else {
                    ctx->ga.Epsilon =
                        ( ctx->ga.Epsilon_0
                        * pow ( 1.0
                              - (double)ctx->ga.iter / ctx->ga.EpsilonGeneration
                              , ctx->ga.EffEpsExponent
                              )
                        );
                }
            }
        }

        /*  If the GA wasn't restarted, update the generation and print
         *  stuff.  We do this because a restart is NOT counted as a
         *  complete generation.
         */
        if (!Restarted) {
            PGAUpdateGeneration (ctx, comm);
            if (rank == 0) {
                PGAPrintReport (ctx, ctx->ga.OutputFile, PGA_OLDPOP);
            }
        }
    }

    if (rank == 0) {
        int pop = PGA_OLDPOP;
        int numaux = PGAGetNumAuxEval (ctx);
        int numcon = PGAGetNumConstraint (ctx);
        if (numaux == numcon) {
            int best_p = PGAGetBestIndex (ctx, pop);
            fprintf
                ( ctx->ga.OutputFile
                , "The Best Evaluation: %e"
                , _PGAGetEvaluation (ctx, best_p, pop, NULL)
                );
            if (numaux) {
                fprintf
                    ( ctx->ga.OutputFile
                    , " Constraints: %e", PGAGetAuxTotal (ctx, best_p, pop)
                    );
            }
            fprintf (ctx->ga.OutputFile, ".\n");
            #if 0
            /* Maybe make this an option? */
            fprintf
                ( ctx->ga.OutputFile
                , "Evaluations: %d\n"
                , PGAGetEvalCount (ctx)
                );
            #endif
            fprintf (ctx->ga.OutputFile, "The Best String:\n");
            PGAPrintString (ctx, ctx->ga.OutputFile, best_p, pop);
        } else {
            int i, k;
            char s [34];
            int p = ctx->rep.MOPrecision + 6;
            sprintf (s, "F %%5d %%%d.%de\n", p, ctx->rep.MOPrecision);
            PGAIndividual *ind = PGAGetIndividual (ctx, 0, pop);
            for (k=0; k<numaux+1; k++) {
                fprintf
                    ( ctx->ga.OutputFile
                    , "The Best (%d) evaluation: %e\n", k, ctx->rep.Best [k]
                    );
            }
            fprintf (ctx->ga.OutputFile, "The Nondominated Strings:\n");
            for (i=0; i<ctx->ga.PopSize; i++, ind++) {
                int j;
                if (ind->rank == 0) {
                    for (j=0; j<numcon; j++) {
                        if (ind->auxeval [numaux - numcon + j] > 0) {
                            break;
                        }
                    }
                    if (j < numcon) {
                        continue;
                    }
                    for (k=0; k<numaux+1; k++) {
                        double e = (k==0) ? ind->evalue : ind->auxeval [k-1];
                        fprintf (ctx->ga.OutputFile, s, k, e);
                    }
                    PGAPrintString (ctx, ctx->ga.OutputFile, i, pop);
                }
            }
        }
        fflush (ctx->ga.OutputFile);
    }
    PGADebugExited ("PGARunGM");
}


/*!****************************************************************************
    \brief Sequential internal evalution function.
    \ingroup internal
    \param   ctx       context variable
    \param   pop       symbolic constant of the population to be evaluated
    \param   evaluate  a pointer to a function to evaluate a string.

    \rst

    Description
    -----------

    Evaluates all strings that need to be evaluated using one processor.

    \endrst

******************************************************************************/
static void PGAEvaluateSeq
    ( PGAContext *ctx
    , int pop
    , double (*evaluate)(PGAContext *, int, int, double *)
    )
{
    int     p;
    double  e;

    PGADebugEntered ("PGAEvaluateSeq");

    /*  Standard sequential evaluation.  */
    for (p=0; p<ctx->ga.PopSize; p++) {
        if (!PGAGetEvaluationUpToDateFlag (ctx, p, pop)) {
            double *aux = PGAGetAuxEvaluation (ctx, p, pop);
            if (ctx->sys.UserFortran) {
                int fp = p + 1;
                e = (*((double(*)(void *, void *, void *, void *))evaluate))
                    (&ctx, &fp, &pop, aux);
            } else {
                e = (*evaluate)(ctx, p, pop, aux);
            }
            PGASetEvaluation (ctx, p, pop, e, aux);
            ctx->rep.nevals++;
        }
    }
    PGADebugExited ("PGAEvaluateSeq");
}
/*!****************************************************************************
    \brief Build an MPI datatype for eval and auxeval.
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
       dt = PGABuildEvaluation (ctx, p, pop);

    \endrst

******************************************************************************/
static
MPI_Datatype PGABuildEvaluation (PGAContext *ctx, int p, int pop)
{
    int            n = 1;
    int            counts[2];      /* Number of elements in each
                                      block (array of integer) */
    MPI_Aint       displs[2];      /* byte displacement of each
                                      block (array of integer) */
    MPI_Datatype   types[2];       /* type of elements in each block (array
                                      of handles to datatype objects) */
    MPI_Datatype   individualtype; /* new datatype (handle) */
    PGAIndividual *traveller;      /* address of individual in question */

    traveller = PGAGetIndividual(ctx, p, pop);
    MPI_Get_address(&traveller->evalue, &displs[0]);
    counts[0] = 1;
    types[0]  = MPI_DOUBLE;

    if (ctx->ga.NumAuxEval) {
        MPI_Get_address (traveller->auxeval, &displs[1]);
        counts[1] = ctx->ga.NumAuxEval;
        types[1]  = MPI_DOUBLE;
        n += 1;
    }

    MPI_Type_create_struct (n, counts, displs, types, &individualtype);
    MPI_Type_commit (&individualtype);

    return (individualtype);
}

/*!****************************************************************************
    \brief Transmit evaluation and aux eval to another process.
    \ingroup internal

    \param  ctx   context variable
    \param  p     index of an individual
    \param  pop   symbolic constant of the population
    \param  dest  ID of the process where this is going
    \param  tag   MPI tag to send with the individual
    \param  comm  MPI communicator

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int p, dest;

      ...
      dest = SelectAFreeProcessor ();
      PGASendEvaluation (ctx, p, PGA_NEWPOP, dest, PGA_COMM_EVALOFSTRING, comm);

    \endrst

******************************************************************************/
static
void PGASendEvaluation (PGAContext *ctx, int p, int pop, int dest, int tag,
                        MPI_Comm comm)
{
    MPI_Datatype individualtype;

    individualtype = PGABuildEvaluation (ctx, p, pop);
    MPI_Send (MPI_BOTTOM, 1, individualtype, dest, tag, comm);
    MPI_Type_free (&individualtype);
}

/*!****************************************************************************
    \brief Receive evaluation and aux eval from another process
    \ingroup internal
    \param  ctx     contex variable
    \param  p       index of an individual
    \param  pop     symbolic constant of the population
    \param  source  ID of the process from which to receive
    \param  tag     MPI tag to look for
    \param  comm    an MPI communicator
    \param  status  pointer to an MPI status structure
    \return string p in population pop is changed by side-effect

    \rst

    Example
    -------

    Receive evaluation from sub-process and place it into the first
    temporary location in :c:macro:`PGA_NEWPOP`.

    .. code-block:: c

      PGAContext *ctx;
      MPI_Comm    comm;
      MPI_Status  status;

      ...
      PGAReceiveEvaluation
        (ctx, PGA_TEMP1, PGA_NEWPOP, 0, PGA_COMM_EVALOFSTRING, comm, &status);

    \endrst

******************************************************************************/
static
void PGAReceiveEvaluation
    ( PGAContext *ctx, int p, int pop
    , int source, int tag, MPI_Comm comm, MPI_Status *status
    )
{
    MPI_Datatype individualtype;

    individualtype = PGABuildEvaluation (ctx, p, pop);
    MPI_Recv (MPI_BOTTOM, 1, individualtype, source, tag, comm, status);
    MPI_Type_free (&individualtype);
}

/*!****************************************************************************
    \brief Cooperative internal evaluation function.
    \ingroup internal
    \param   ctx      context variable
    \param   pop      symbolic constant of the population to be evaluated
    \param   evaluate a pointer to a function to evaluate a string.
    \param   comm     an MPI communicator
    \return  None

    \rst

    Description
    -----------

    Evaluates all strings that need to be evaluated using two processors
    cooperatively.  The first being the rank-0 process will send a
    string to the second for evaluation.  While the second is
    evaluating, the rank-0 process will *also* evaluate a string.

    \endrst

******************************************************************************/
static void PGAEvaluateCoop
    ( PGAContext *ctx
    , int pop
    , double (*evaluate)(PGAContext *, int, int, double *)
    , MPI_Comm comm
    )
{
    MPI_Status      stat;
    int             p, fp, q;
    double          e;
    PGAIndividual  *ind;

    PGADebugEntered ("PGAEvaluateCoop");

    q = -1;

    ind = PGAGetIndividual (ctx, 0, pop);

    for (p=0; p<ctx->ga.PopSize;) {
        while ((p<ctx->ga.PopSize) && (ind+p)->evaluptodate) {
            p++;
        }
        if (p<ctx->ga.PopSize) {
            PGASendIndividual (ctx, p, pop, 1, PGA_COMM_STRINGTOEVAL, comm);
            q = p;
        }
        p++;

        while ((p<ctx->ga.PopSize) && (ind+p)->evaluptodate) {
            p++;
        }
        if (p<ctx->ga.PopSize) {
            double *aux = PGAGetAuxEvaluation (ctx, p, pop);
            if (ctx->sys.UserFortran == PGA_TRUE) {
                fp = p+1;
                e = (*((double(*)(void *, void *, void *, void *))evaluate))
                    (&ctx, &fp, &pop, aux);
            } else {
                e = (*evaluate)(ctx, p, pop, aux);
            }
            PGASetEvaluation (ctx, p, pop, e, aux);
            ctx->rep.nevals++;
#if DEBUG_EVAL
            fprintf (stdout, "%4d: %10.8e Local\n", p, e);
            fflush  (stdout);
#endif
        }

        if (q >= 0) {
            PGAReceiveEvaluation
                (ctx, q, pop, 1, PGA_COMM_EVALOFSTRING, comm, &stat);
            PGASetEvaluationUpToDateFlag (ctx, q, pop, PGA_TRUE);
            ctx->rep.nevals++;
#if DEBUG_EVAL
            fprintf (stdout, "%4d: %10.8e Worker %d\n", p, e, 1);
            fflush  (stdout);
#endif
            q = -1;
        }
    }

    /*  Release the Worker  */
    MPI_Send (&q, 1, MPI_INT, 1, PGA_COMM_DONEWITHEVALS, comm);

    PGADebugExited ("PGAEvaluateCoop");
}



/*!****************************************************************************
    \brief Internal evaluation function, multiprocessing version.
    \ingroup internal
    \param   ctx      context variable
    \param   pop      symbolic constant of the population to be evaluated
    \param   comm     an MPI communicator
    \return  None

    \rst

    Description
    -----------

    Evaluates all strings that need evaluating using three or more
    processors.  The rank-0 process sends individuals to evaluate to all
    other processes and collects the results. The rank-0 process
    performs no evaluations itself.

    \endrst

******************************************************************************/
static void PGAEvaluateMP (PGAContext *ctx, int pop, MPI_Comm comm)
{
    int    *work;
    int     i, k, s, p, size, sentout;
    MPI_Status stat;
    PGAIndividual *ind;
    PGAIndividual *tmp1 = PGAGetIndividual (ctx, PGA_TEMP1, pop);

    PGADebugEntered ("PGAEvaluateMP");

    size = PGAGetNumProcs (ctx, comm);

    work = malloc (size *sizeof(int));
    if (work == NULL) {
        PGAError
            ( ctx, "PGAEvaluateMP:  Couldn't allocate work array"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }

    sentout = 0;
    s = 1;
    ind = PGAGetIndividual (ctx, 0, pop);

    /*  Send strings to all processes, since they are all unused.  */
    for (k=0; ((k<ctx->ga.PopSize) && (s<size)); k++) {
        if ((ind+k)->evaluptodate == PGA_FALSE) {
            work[s] = k;
            PGASendIndividual (ctx, k, pop, s, PGA_COMM_STRINGTOEVAL, comm);
#if DEBUG_EVAL
            fprintf (stdout, "%4d: Sent to worker %d.\n", k, s);
            fflush  (stdout);
#endif
            sentout++;
            s++;
        }
    }

    /*  Move to the next string to be evaluated.  Notice that all we need
     *  to do is skip any strings that are already evaluated, unlike
     *  below, where we need to _first_ go to the next string, then
     *  skip any that are up to date.
     */
    while ((k<ctx->ga.PopSize) && (ind+k)->evaluptodate) {
        k++;
    }

    /*  While there are still unevaluated individuals, receive whatever
     *  is waiting, then immediately send a new string to it.  This
     *  implicitly will balance the load across the machines, as we
     *  initially sent a string to _each_ process, so _each_ process
     *  will return an evaluation and get a new one immediately.
     */
    while(k<ctx->ga.PopSize) {
        /*  Receive the next evaluated string.  */
        PGAReceiveEvaluation
            ( ctx, PGA_TEMP1, pop
            , MPI_ANY_SOURCE, PGA_COMM_EVALOFSTRING, comm, &stat
            );
        p = work [stat.MPI_SOURCE];
        PGASetEvaluation (ctx, p, pop, tmp1->evalue, tmp1->auxeval);
        ctx->rep.nevals++;
#if DEBUG_EVAL
        fprintf
            ( stdout
            , "%4d: %10.8e Worker %d  Sent %d\n"
            , work[stat.MPI_SOURCE], e, stat.MPI_SOURCE, k
            );
        fflush (stdout);
#endif
        /*  Immediately send another string to be evaluated.  */
        work [stat.MPI_SOURCE] = k;
        PGASendIndividual
            (ctx, k, pop, stat.MPI_SOURCE, PGA_COMM_STRINGTOEVAL, comm);

        /*  Find the next unevaluated individual  */
        k++;
        while ((k<ctx->ga.PopSize) && (ind+k)->evaluptodate) {
            k++;
        }
    }

    /*  All strings have been sent out.  Wait for them to be done.  */
    while (sentout > 0) {
        PGAReceiveEvaluation
            ( ctx, PGA_TEMP1, pop
            , MPI_ANY_SOURCE, PGA_COMM_EVALOFSTRING, comm, &stat
            );
        p = work [stat.MPI_SOURCE];
        PGASetEvaluation (ctx, p, pop, tmp1->evalue, tmp1->auxeval);
        ctx->rep.nevals++;
        sentout--;
#if DEBUG_EVAL
        fprintf
            ( stdout
            , "%4d: %10.8e Worker %d\n"
            , work [stat.MPI_SOURCE], e, stat.MPI_SOURCE
            );
        fflush (stdout);
#endif
    }
    free (work);

    /* Release the workers */
    for (i=1; i<size; i++) {
        MPI_Send (&i, 1, MPI_INT, i, PGA_COMM_DONEWITHEVALS, comm);
    }

    PGADebugExited ("PGAEvaluateMP");
}


/*!****************************************************************************
    \brief Evaluation worker routine.
    \ingroup internal
    \param   ctx      context variable
    \param   pop      symbolic constant of the population to be evaluated
    \param   evaluate a pointer to a function to evaluate a string.
    \param   comm     an MPI communicator
    \return  None

    \rst

    Description
    -----------

    Sit around and wait for a string to eval to show up, then evaluate
    it and return the evaluation.  Terminates when it receives
    :c:macro:`PGA_COMM_DONEWITHEVALS`.

    \endrst

******************************************************************************/
static void PGAEvaluateWorker
    ( PGAContext *ctx
    , int pop
    , double (*evaluate)(PGAContext *, int, int, double *)
    , MPI_Comm comm
    )
{
    MPI_Status  stat;
    int         k;
    double      e;

    PGADebugEntered ("PGAEvaluateWorker");

    k = PGA_TEMP1;

    MPI_Probe (0, MPI_ANY_TAG, comm, &stat);
    while (  stat.MPI_TAG == PGA_COMM_STRINGTOEVAL
          || stat.MPI_TAG == PGA_COMM_SERIALIZE_SIZE
          )
    {
        double *aux;
        PGAReceiveIndividual
            (ctx, PGA_TEMP1, pop, 0, PGA_COMM_STRINGTOEVAL, comm, &stat);

        aux = PGAGetAuxEvaluation (ctx, PGA_TEMP1, pop);
        if (ctx->sys.UserFortran == PGA_TRUE) {
            e = (*((double(*)(void *, void *, void *, void *))evaluate))
                (&ctx, &k, &pop, aux);
        } else {
            e = (*evaluate)(ctx, PGA_TEMP1, pop, aux);
        }
        PGASetEvaluation (ctx, PGA_TEMP1, pop, e, aux);

        PGASendEvaluation (ctx, PGA_TEMP1, pop, 0, PGA_COMM_EVALOFSTRING, comm);
        MPI_Probe (0, MPI_ANY_TAG, comm, &stat);
    }
    MPI_Recv (&k, 1, MPI_INT, 0, PGA_COMM_DONEWITHEVALS, comm, &stat);

    PGADebugExited ("PGAEvaluateWorker");
}


/*!****************************************************************************
    \brief Call a user-specified function to return an evaluation of
           each string in the population.
    \ingroup explicit
    \param   ctx      context variable
    \param   pop      symbolic constant of the population to be evaluated
    \param   evaluate a pointer to a function to evaluate a string.
    \param   comm     an MPI communicator
    \return  Evaluates the population via side effect

    \rst

    Description
    -----------

    The user-specified function is only called if the string has been
    changed (e.g., by crossover or mutation) or the user has explicitly
    signaled the string's evaluation is out-of-date by a call to
    :c:func:`PGASetEvaluationUpToDateFlag`.

    The user-specified function will be called once for each string in
    population ``pop`` that requires evaluation.  This function must return
    a ``double`` (the evaluation function value) and must fit the prototype::

        double evaluate (PGAContext *c, int p, int pop, double *aux);

    Example
    -------

    Evaluate all strings in population :c:macro:`PGA_NEWPOP` using the
    user-defined evaluation function ``Energy``.

    .. code-block:: c

       double Energy (PGAContext *ctx, int p, int pop, double *aux) {
           ...
       };

       PGAContext *ctx;

       ...
       PGAEvaluate (ctx, PGA_NEWPOP, Energy, MPI_COMM_WORLD);

    \endrst

******************************************************************************/
void PGAEvaluate
    ( PGAContext *ctx
    , int pop
    , double (*evaluate)(PGAContext *, int, int, double *)
    , MPI_Comm comm
    )
{
    int  rank, size;

    PGADebugEntered ("PGAEvaluate");

    rank = PGAGetRank (ctx, comm);
    size = PGAGetNumProcs (ctx, comm);

    if (rank == 0) {
        if (size == 1) {
            PGAEvaluateSeq (ctx, pop, evaluate);
        } else if (size == 2) {
            PGAEvaluateCoop (ctx, pop, evaluate, comm);
        } else if (size > 2) {
            PGAEvaluateMP (ctx, pop, comm);
        }
    } else {
        PGAEvaluateWorker (ctx, pop, evaluate, comm);
    }

    PGADebugExited ("PGAEvaluate");
}


/*!****************************************************************************
    \brief Build an MPI datatype for string p in population pop.
    \ingroup explicit
    \param  ctx      context variable
    \param  p        index of an individual
    \param  pop      symbolic constant of the population
    \return An MPI datatype for member p of population pop

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int p;
      MPI_Datatype dt;

      ...
      dt = PGABuildDatatype (ctx, p, PGA_NEWPOP);

    \endrst

******************************************************************************/
MPI_Datatype PGABuildDatatype(PGAContext *ctx, int p, int pop)
{
    PGADebugEntered("PGABuildDatatype");

    PGADebugExited("PGABuildDatatype");

    return((*ctx->cops.BuildDatatype)(ctx, p, pop));
}


/*!****************************************************************************
    \brief Common part for building an MPI datatype for an individual
    \ingroup explicit
    \param  ctx       context variable
    \param  p         index of string
    \param  pop       symbolic constant of the population string p is in
    \param  counts    Number of elements in each block
    \param  displs    byte displacement array pointer
    \param  types     type elements
    \return The index of the next-to-be-filled counts, displs, types

    \rst

    Description
    -----------

    Returns by side effect a partially filled array of ``counts``,
    ``displs``, ``types``. The returned index will be max
    :c:macro:`PGA_MPI_HEADER_ELEMENTS`.
    This means callers may use a statically allocated buffer.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int idx = 0;
      int counts [PGA_MPI_HEADER_ELEMENTS + ...];
      int displs [PGA_MPI_HEADER_ELEMENTS + ...];
      MPI_Datatype types [PGA_MPI_HEADER_ELEMENTS + ...];

      ...
      idx = PGABuildDatatypeHeader (ctx, counts, displs, types);
      // Fill rest of counts, displs, types here and build datatype
      counts [idx] =
      displs [idx] =
      types  [idx] =
      idx++;
      ...

    \endrst

******************************************************************************/
int PGABuildDatatypeHeader
    ( PGAContext *ctx, int p, int pop
    , int *counts, MPI_Aint *displs, MPI_Datatype *types
    )
{
    int n = 6;
    PGAIndividual *traveller = PGAGetIndividual (ctx, p, pop);
    MPI_Get_address (&traveller->evalue, &displs [0]);
    counts [0] = 1;
    types  [0] = MPI_DOUBLE;
    MPI_Get_address (&traveller->fitness, &displs [1]);
    counts [1] = 1;
    types  [1] = MPI_DOUBLE;
    MPI_Get_address (&traveller->evaluptodate, &displs [2]);
    counts [2] = 1;
    types  [2] = MPI_INT;
    MPI_Get_address (&traveller->auxtotal, &displs [3]);
    counts [3] = 1;
    types  [3] = MPI_DOUBLE;
    MPI_Get_address (&traveller->auxtotalok, &displs [4]);
    counts [4] = 1;
    types  [4] = MPI_INT;
    MPI_Get_address (&traveller->index, &displs [5]);
    counts [5] = 1;
    types  [5] = MPI_INT;
    if (ctx->ga.NumAuxEval) {
        MPI_Get_address (traveller->auxeval, &displs [6]);
        counts [6] = ctx->ga.NumAuxEval;
        types  [6] = MPI_DOUBLE;
        n += 1;
    }
    assert (n <= PGA_MPI_HEADER_ELEMENTS);
    return n;
}


/*!****************************************************************************
    \brief Build datatype from serialized data.
    \ingroup internal
    \param  ctx  context variable
    \param  p    index of string
    \param  pop  symbolic constant of the population string p is in
    \return An MPI_Datatype

    \rst

    Description
    -----------

    This is the default function for the BuildDatatype user function if
    serialization is in use.

    \endrst

******************************************************************************/
MPI_Datatype PGASerializedBuildDatatype (PGAContext *ctx, int p, int pop)
{
    int            idx = 0;
    /* Number of elements in each block (array of integer) */
    int            counts [PGA_MPI_HEADER_ELEMENTS + 1];
    /* byte displacement of each block (array of integer) */
    MPI_Aint       displs [PGA_MPI_HEADER_ELEMENTS + 1];
    /* type of elements in each block (array of handles to datatype objects) */
    MPI_Datatype   types  [PGA_MPI_HEADER_ELEMENTS + 1];
    MPI_Datatype   individualtype; /* new datatype (handle) */

    idx = PGABuildDatatypeHeader (ctx, p, pop, counts, displs, types);

    MPI_Get_address (ctx->scratch.serialized, &displs [idx]);
    counts [idx] = ctx->scratch.serialization_size;
    types  [idx] = MPI_BYTE;
    idx++;

    MPI_Type_create_struct (idx, counts, displs, types, &individualtype);
    MPI_Type_commit (&individualtype);

    return individualtype;
}


/*!****************************************************************************
    \brief Transmit an individual to another process.
    \ingroup explicit
    \param  ctx   context variable
    \param  p     index of an individual
    \param  pop   symbolic constant of the population
    \param  dest  ID of the process where this is going
    \param  tag   MPI tag to send with the individual
    \param  comm  MPI communicator
    \return None

    \rst

    Description
    -----------

    This is usually called by one of the functions called via
    :c:func:`PGAEvaluate`.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      int p, dest;

      ...
      dest = SelectAFreeProcessor ();
      PGASendIndividual (ctx, p, PGA_NEWPOP, dest, PGA_SR_STRINGTOEVAL, comm);

    \endrst

******************************************************************************/
void PGASendIndividual
    (PGAContext *ctx, int p, int pop, int dest, int tag, MPI_Comm comm)
{
    MPI_Datatype individualtype;

    PGADebugEntered ("PGASendIndividual");

    /* If we do serialization of user defined type we must first send
     * the size
     */
    if (ctx->cops.Serialize) {
        unsigned long size = 0;
        assert (ctx->ga.datatype == PGA_DATATYPE_USER);
        assert (ctx->scratch.serialized == NULL);
        ctx->scratch.serialization_size = ctx->cops.Serialize
            (ctx, p, pop, (const void **)&ctx->scratch.serialized);
        assert (ctx->scratch.serialization_size != 0);
        size = (unsigned long)ctx->scratch.serialization_size;
        MPI_Send
            ( &size, 1, MPI_UNSIGNED_LONG
            , dest, PGA_COMM_SERIALIZE_SIZE, comm
            );
    }

    individualtype = PGABuildDatatype (ctx, p, pop);
    MPI_Send (MPI_BOTTOM, 1, individualtype, dest, tag, comm);
    MPI_Type_free (&individualtype);

    if (ctx->cops.Serialize) {
        ctx->cops.SerializeFree (ctx->scratch.serialized);
        ctx->scratch.serialized = NULL;
        ctx->scratch.serialization_size = 0;
    }

    PGADebugExited ("PGASendIndividual");
}

/*!****************************************************************************
    \brief Receive an individual from another process.
    \ingroup explicit
    \param  ctx     contex variable
    \param  p       index of an individual
    \param  pop     symbolic constant of the population
    \param  source  ID of the process from which to receive
    \param  tag     MPI tag to look for
    \param  comm    MPI communicator
    \param  status  pointer to an MPI status structure
    \return Status and string p in population pop are changed by side-effect

    \rst

    Description
    -----------

    This is usually called by one of the functions called via
    :c:func:`PGAEvaluate`.

    Example
    -------

    Receive a string from the rank-0 process with tag
    PGA_SR_STRINGTOEVAL, and place it into the first temporary location
    in PGA_NEWPOP.

    .. code-block:: c

      PGAContext *ctx;
      MPI_Comm    comm;
      MPI_Status  status;

      ...
      PGAReceiveIndividual
        (ctx, PGA_TEMP1, PGA_NEWPOP, 0, PGA_SR_STRINGTOEVAL, comm, &status);

    \endrst

******************************************************************************/
void PGAReceiveIndividual
    ( PGAContext *ctx, int p, int pop
    , int source, int tag, MPI_Comm comm, MPI_Status *status
    )
{
    MPI_Datatype individualtype;

    PGADebugEntered ("PGAReceiveIndividual");

    /* If we do serialization of user defined type we must first receive
     * the size
     */
    if (ctx->cops.Serialize) {
        unsigned long size = 0;
        assert (ctx->ga.datatype == PGA_DATATYPE_USER);
        assert (ctx->scratch.serialized == NULL);
        MPI_Recv
            ( &size, 1, MPI_UNSIGNED_LONG, source
            , PGA_COMM_SERIALIZE_SIZE, comm, status
            );
        ctx->scratch.serialization_size = size;
        assert (ctx->scratch.serialization_size != 0);
        ctx->scratch.serialized = malloc (ctx->scratch.serialization_size);
        if (ctx->scratch.serialized == NULL) {
            PGAErrorPrintf
                (ctx, PGA_FATAL, "Cannot allocate serialization buffer");
        }
    }

    individualtype = PGABuildDatatype (ctx, p, pop);
    MPI_Recv (MPI_BOTTOM, 1, individualtype, source, tag, comm, status);
    MPI_Type_free (&individualtype);

    if (ctx->cops.Serialize) {
        ctx->cops.Deserialize
            ( ctx, p, pop
            , ctx->scratch.serialized, ctx->scratch.serialization_size
            );
        free (ctx->scratch.serialized);
        ctx->scratch.serialized = NULL;
        ctx->scratch.serialization_size = 0;
    }

    PGADebugExited ("PGAReceiveIndividual");
}

/*!****************************************************************************
    \brief Send an individual to a process, while receiving a different
           individual from a different process.
    \ingroup explicit
    \param  ctx        context variable
    \param  send_p     index of string to send
    \param  send_pop   symbolic constant of population to send from
    \param  dest       destination process
    \param  send_tag   tag to send with
    \param  recv_p     index of string to receive
    \param  recv_pop   symbolic constant of population to receive from
    \param  source     process to receive from
    \param  recv_tag   tag to receive with
    \param  comm       an MPI communicator
    \param  status     pointer to the MPI status structure
    \return status and string recv_p in population recv_pop are modified by
            side-effect

    \rst

    Example
    -------

    A dedicated process is being used to perform an optimization algorithm
    on the strings.  Send a new string, ``s``, to the process, while
    receiving an optimized string, ``r``, from it.

    .. code-block:: c

      PGAContext *ctx;
      MPI_Comm    comm;
      MPI_Status  status;
      int  s, r;

      ...
      PGASendReceiveIndividual
        ( ctx
        , s, PGA_NEWPOP, 1, PGA_SR_STRINGTOMODIFY
        , r, PGA_NEWPOP, 1, PGA_SR_MODIFIEDSTRING
        , comm, &status
        );

    \endrst

******************************************************************************/
void PGASendReceiveIndividual
    ( PGAContext *ctx
    , int send_p, int send_pop, int dest, int send_tag
    , int recv_p, int recv_pop, int source, int recv_tag
    , MPI_Comm comm, MPI_Status *status
    )
{
    MPI_Datatype individualsendtype;
    MPI_Datatype individualrecvtype;

    PGADebugEntered ("PGASendReceiveIndividual");

    individualsendtype = PGABuildDatatype (ctx, send_p, send_pop);
    individualrecvtype = PGABuildDatatype (ctx, recv_p, recv_pop);

    MPI_Sendrecv
        ( MPI_BOTTOM, 1, individualsendtype, dest,   send_tag
        , MPI_BOTTOM, 1, individualrecvtype, source, recv_tag
        , comm, status
        );

    MPI_Type_free (&individualsendtype);
    MPI_Type_free (&individualrecvtype);

    PGADebugExited ("PGASendReceiveIndividual");
}


/*!****************************************************************************
    \brief Execute the island model genetic algorithm
    \ingroup notimplemented
    \param  ctx       context variable
    \param  evaluate  a pointer to the user's evaluation function, which must
                      have the calling sequence shown in the example.
    \param  comm      the MPI communicator to use
    \return None

    \rst

    Description
    -----------

    Not yet implemented.
    Based on ``ctx->par.topology`` this routine will need to create the
    appropriate communicator out of ``comm``.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      double f (PGAContext *ctx, int p, int pop, double *aux);
      MPI_Comm comm;

      ...
      PGARunIM (ctx, f, comm);

    \endrst

******************************************************************************/
void PGARunIM
    ( PGAContext *ctx
    , double (*evaluate)(PGAContext *c, int p, int pop, double *aux)
    , MPI_Comm comm
    )
{
     PGADebugEntered ("PGARunIM");
     PGAError
        ( ctx, "PGARunIM: Island model not implemented"
        , PGA_FATAL, PGA_VOID, NULL
        );
     PGADebugExited ("PGARunIM");
}


/*!****************************************************************************
    \brief Execute a neighborhood model genetic algorithm
    \ingroup notimplemented
    \param  ctx       context variable
    \param  evaluate  a pointer to the user's evaluation function, which must
                      have the calling sequence shown in the example.
    \param  comm      the MPI communicator to use
    \return None

    \rst

    Description
    -----------

    Not yet implemented.
    Based on ``ctx->par.topology`` this routine will need to create the
    appropriate communicator out of ``comm``.

    Example
    -------

    .. code-block:: c

      PGAContext *ctx,
      MPI_Comm comm;
      double f (PGAContext *ctx, int p, int pop);

      ...
      PGARunNM (ctx, f, comm);

    \endrst

******************************************************************************/
void PGARunNM
    ( PGAContext *ctx
    , double (*evaluate)(PGAContext *c, int p, int pop, double *aux)
    , MPI_Comm comm
    )
{
    PGADebugEntered ("PGARunNM");
    PGAError
        ( ctx, "PGARunNM: Island model not implemented"
        , PGA_FATAL, PGA_VOID, NULL
        );
    PGADebugExited ("PGARunNM");
}



/*!****************************************************************************
    \brief Return the rank of the processor in communicator comm.
    \ingroup parallel
    \param   ctx   context variable structure pointer
    \param   comm  an MPI communicator
    \return  The rank of this processor

    \rst

    Description
    -----------

    If ``comm`` is ``NULL`` or a sequential version of PGAPack is used,
    0 is returned.

    Example
    -------

    .. code-block:: c

        PGAContext  *ctx;
        int          rank;

        ...
        rank = PGAGetRank (ctx, MPI_COMM_WORLD);
        if (rank == 0) {
            LetRank0DoSomething ();
        }

    \endrst

******************************************************************************/
int PGAGetRank (PGAContext *ctx, MPI_Comm comm)
{
    int rank;

    PGADebugEntered ("PGAGetRank");

    if (comm == MPI_COMM_NULL) {
        rank = 0;
    } else {
        MPI_Comm_rank (comm, &rank);
    }

    PGADebugExited ("PGAGetRank");

    return rank;
}


/*!****************************************************************************
    \brief Return the size of communicator comm in processes.
    \ingroup parallel
    \param   ctx   context variable structure pointer
    \param   comm  an MPI communicator
    \return  The numbers of processors in communicator comm

    \rst

    Description
    -----------

    If ``comm`` is ``NULL`` or a sequential version of PGAPack is used,
    1 is returned.

    Example
    -------

    .. code-block:: c

        PGAContext  *ctx;

        ...
        if (PGAGetNumProcs (ctx, MPI_COMM_WORLD) < 4) {
            printf ("Too few processors for decent performance!\n");
            exit (-1);
        }

    \endrst

******************************************************************************/
int PGAGetNumProcs (PGAContext *ctx, MPI_Comm comm)
{
    int size;

    PGADebugEntered ("PGAGetNumProcs");

    if (comm == MPI_COMM_NULL) {
        size = 1;
    } else {
        MPI_Comm_size (comm, &size);
    }

    PGADebugExited ("PGAGetNumProcs");

    return size;
}


/*!****************************************************************************
    \brief Set the number of islands to use in an island model GA.
    \ingroup notimplemented
    \param   ctx  context variable
    \param   n    number of islands
    \return  None

    \rst

    Description
    -----------

    The default is one.  Currently must be the same as the number of
    processes in the default communicator.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx,
       double f (PGAContext *ctx, int p, int pop, double *aux);

       ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
       PGASetNumIslands (ctx, 10);
       PGASetUp (ctx);
       PGARun (ctx, f);
       PGADestroy (ctx);

    \endrst

******************************************************************************/
void PGASetNumIslands( PGAContext *ctx, int n)
{

    PGADebugEntered ("PGASetNumIslands");

    if (n < 1) {
        PGAError
            ( ctx, "PGASetNumIslands: Invalid value of n:"
            , PGA_FATAL, PGA_INT, (void *) &n
            );
    }

    ctx->par.NumIslands = n;

    PGADebugExited ("PGASetNumIslands");
}


/*!***************************************************************************
   \brief Returns the number of islands to use in an island model
   \ingroup notimplemented
   \param   ctx  context variable
   \return  The number of islands to use in an island model

   \rst

   Example
   -------

    .. code-block:: c

      PGAContext *ctx;
      int npop;

      ...
      npop = PGAGetNumIslands (ctx);

    \endrst

*****************************************************************************/
int PGAGetNumIslands (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetNumIslands");
    PGAFailIfNotSetUp ("PGAGetNumIslands");

    PGADebugExited ("PGAGetNumIslands");

    return ctx->par.NumIslands;
}

/*!****************************************************************************
    \brief Set the number of demes to use in a neighborhood model GA.
    \ingroup notimplemented

    \param   ctx           context variable
    \param   numdemes      number of demes
    \return  None

    \rst

    Description
    -----------

    Currently must be the same as the number of processes in the default
    communicator.  The default is one.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx,
       double f (PGAContext *ctx, int p, int pop, double *aux);

       ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
       PGASetNumDemes (ctx, 4);
       PGASetUp (ctx);
       PGARun (ctx, f);
       PGADestroy (ctx);

    \endrst

******************************************************************************/
void PGASetNumDemes( PGAContext *ctx, int numdemes)
{
    PGADebugEntered ("PGASetNumDemes");

    if (numdemes < 1) {
        PGAError
            ( ctx, "PGASetNumDemes: Invalid value of numdemes:"
            , PGA_FATAL, PGA_INT, (void *) &numdemes
            );
    }

    ctx->par.NumDemes = numdemes;

    PGADebugExited ("PGASetNumDemes");
}


/*!***************************************************************************
    \brief Returns the number of demes to use in a neighborhood model.
    \ingroup notimplemented
    \param   ctx  context variable
    \return  The number of demes to use in a neighborhood model

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int npop;

       ...
       npop = PGAGetNumDemes (ctx);

    \endrst

*****************************************************************************/
int PGAGetNumDemes (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetNumDemes");
    PGAFailIfNotSetUp ("PGAGetNumDemes");

    PGADebugExited ("PGAGetNumDemes");

    return ctx->par.NumDemes;
}


/*!****************************************************************************
    \brief Set the default communicator to use when PGARun is called.
    \ingroup init
    \param   ctx     context variable
    \param   comm    communicator to use
    \return  None

    \rst

    Description
    -----------

    Does not necessarily need to be ``MPI_COMM_WORLD`` (which is the
    default).

    Example
    -------

    .. code-block:: c

       MPI_Comm mycomm;
       PGAContext *ctx,
       double f (PGAContext *ctx, int p, int pop, double *aux);

       ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
       PGASetCommunicator (ctx, mycomm);
       PGASetUp (ctx);
       PGARun (ctx, f);
       PGADestroy (ctx);

    \endrst

******************************************************************************/
void PGASetCommunicator (PGAContext *ctx, MPI_Comm comm)
{

    PGADebugEntered ("PGASetCommunicator");

    ctx->par.DefaultComm = comm;

    PGADebugExited ("PGASetCommunicator");
}


/*!****************************************************************************
    \brief Returns the default communicator used when PGARun is called.
    \ingroup query
    \param   ctx     context variable
    \return  The default communicator

    \rst

    Example
    -------

    .. code-block:: c

       MPI_Comm comm;
       PGAContext *ctx,
       double f (PGAContext *ctx, int p, int pop, double *aux);

       ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
       PGASetUp (ctx);
       comm = PGAGetCommunicator (ctx);

    \endrst

******************************************************************************/
MPI_Comm PGAGetCommunicator (PGAContext *ctx)
{

    PGADebugEntered ("PGAGetCommunicator");

    PGADebugExited  ("PGAGetCommunicator");

    return ctx->par.DefaultComm;
}
