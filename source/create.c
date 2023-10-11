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
* This file contains functions to create and initialize data structures
* and populations.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include <stdint.h>
#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* Helper function to free builtin datatype */
static void ChromFree (PGAIndividual *ind)
{
    free (ind->chrom);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Create an uninitialized context variable.
    \ingroup standard-api
    \param   argc      Address of the count of the number of command
                       line arguments
    \param   argv      Array of command line arguments
    \param   datatype  The data type used for the strings
    \param   len       The string length (number of genes)
    \param   maxormin  The direction of optimization.
    \return  A pointer to the context variable

    \rst

    Description
    -----------

    The Fortran version of this function call contains only the last
    three arguments.

    The ``datatype`` must be one of :c:macro:`PGA_DATATYPE_BINARY`,
    :c:macro:`PGA_DATATYPE_CHARACTER`, :c:macro:`PGA_DATATYPE_INTEGER`,
    :c:macro:`PGA_DATATYPE_REAL`, or :c:macro:`PGA_DATATYPE_USER` to
    specify binary-valued, character-valued, integer-valued,
    real-valued, or a user-defined datatype, respectively. See
    :ref:`group:const-datatype` for the constants and types.

    The ``maxormin`` parameter must be one of :c:macro:`PGA_MAXIMIZE` or
    :c:macro:`PGA_MINIMIZE` for maximization or minimization, respectively.
    See :ref:`group:const-opt-dir` for the constants.

    Example
    -------

    In C:

    .. code-block:: c

      void main (int argc, char **argv) {
          PGAContext *ctx;

          ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
          //  Set options here
          PGASetUp (ctx);
          //  Run the GA here
          PGADestroy (ctx);
      }

    In Fortran:

    .. code-block:: fortran

               integer ctx

               ctx = PGACreate(PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE)
       c       Set options here
               call PGASetUp(ctx)
       c       Run the GA here
               call PGADestroy(ctx)
               stop
               end

    \endrst

******************************************************************************/
PGAContext *PGACreate
    (int *argc, char **argv, int datatype, int len, int maxormin)
{
    int i;
    PGAContext *ctx;

    ctx = (PGAContext *) malloc (sizeof (PGAContext));

    /*  We cannot make PGA calls until we sort the FuncNameIndex below,
     *  so we just manually print the (rather severe) error message.
     */
    if (ctx == NULL) {
        fprintf (stderr, "PGACreate: No room to allocate ctx\n");
        exit (-1);
    }
    memset (ctx, 0, sizeof (*ctx));


    /*  We use this (indirectly) in PGAReadCmdLine -- in processing
     *  -pgahelp and -pgahelp debug.
     */
    MPI_Initialized (&ctx->par.MPIAlreadyInit);

    /* Initialize MPI, only if it isn't already running (fortran)  */
    if (!ctx->par.MPIAlreadyInit) {
        MPI_Init (argc, &argv);
    }


    /*  Sort the FuncNameIndex.  This allows us to use a binary search
     *  for finding the function names.
     */
    PGASortFuncNameIndex (ctx);

    /* Initialize debug flags, then parse command line arguments.  */
    for (i=0; i<PGA_DEBUG_MAXFLAGS; i++) {
        ctx->debug.PGADebugFlags [i] = PGA_FALSE;
    }
    PGAReadCmdLine (ctx, argc, argv);


    /*  The context variable is now initialized enough to allow this
     *  call to complete successfully.
     */
    PGADebugEntered ("PGACreate");

    /* required parameter 1: abstract data type */
    switch (datatype)
    {
    case PGA_DATATYPE_BINARY:
    case PGA_DATATYPE_INTEGER:
    case PGA_DATATYPE_REAL:
    case PGA_DATATYPE_CHARACTER:
    case PGA_DATATYPE_USER:
        ctx->ga.datatype  = datatype;
        break;
    default:
        PGAError( ctx, "PGACreate: Invalid value of datatype:"
                , PGA_FATAL, PGA_INT, (void *) &datatype
                );
    };

    /* required parameter 2: string length */
    if (len <= 0) {
        PGAError ( ctx, "PGACreate: Invalid value of len:"
                 , PGA_FATAL, PGA_INT, (void *) &len
                 );
    } else {
        ctx->ga.StringLen = len;
    }


    /* required parameter 3: optimization direction */
    switch (maxormin) {
        case PGA_MAXIMIZE:
        case PGA_MINIMIZE:
            ctx->ga.optdir = maxormin;
            break;
        default:
            PGAError ( ctx, "PGACreate: Invalid value of optdir:"
                     , PGA_FATAL, PGA_INT, (void *) &maxormin
                     );
    };


    /*  For datatype == PGA_DATATYPE_BINARY, set how many full words
     *  are used in the packed representation, and how many extra bits
     *  this leaves us with.  Finally, set how many total words are used;
     *  if there are no extra bits, this is just the number of full words,
     *  else, there is one more word used than the number of full words.
     */
    switch (datatype) {
    case PGA_DATATYPE_BINARY:
        ctx->ga.fw = ctx->ga.StringLen/WL;
        ctx->ga.eb = ctx->ga.StringLen%WL;
        if (ctx->ga.eb == 0) {
            ctx->ga.tw = ctx->ga.fw;
        } else {
            ctx->ga.tw = ctx->ga.fw+1;
        }
        break;
    default:
        ctx->ga.fw = PGA_UNINITIALIZED_INT;
        ctx->ga.eb = PGA_UNINITIALIZED_INT;
        ctx->ga.tw = PGA_UNINITIALIZED_INT;
        break;
    }

    /*  Clear all the setting.  Later on, PGASetUp() will be called, and then
     *  it will notice which setting are uninitialized, and set them to the
     *  default value.
     */
    ctx->ga.PopSize            = PGA_UNINITIALIZED_INT;
    ctx->ga.NumAuxEval         = PGA_UNINITIALIZED_INT;
    ctx->ga.NumConstraint      = PGA_UNINITIALIZED_INT;
    ctx->ga.SumConstraints     = PGA_UNINITIALIZED_INT;
    ctx->ga.Epsilon            = 0.0;
    ctx->ga.Epsilon_0          = 0.0;
    ctx->ga.EpsilonGeneration  = PGA_UNINITIALIZED_INT;
    ctx->ga.EpsilonExponent    = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.EffEpsExponent     = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.EpsTLambda         = PGA_UNINITIALIZED_INT;
    ctx->ga.EpsilonTheta       = PGA_UNINITIALIZED_INT;
    ctx->ga.StoppingRule       = PGA_STOP_MAXITER;
    ctx->ga.MaxIter            = PGA_UNINITIALIZED_INT;
    ctx->ga.MaxNoChange        = PGA_UNINITIALIZED_INT;
    ctx->ga.MaxSimilarity      = PGA_UNINITIALIZED_INT;
    ctx->ga.NumReplace         = PGA_UNINITIALIZED_INT;
    ctx->ga.CrossoverType      = PGA_UNINITIALIZED_INT;
    ctx->ga.CrossBoundedFlag   = PGA_UNINITIALIZED_INT;
    ctx->ga.CrossBounceFlag    = PGA_UNINITIALIZED_INT;
    ctx->ga.CrossSBXEta        = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.CrossSBXOnce       = PGA_UNINITIALIZED_INT;
    ctx->ga.SelectType         = PGA_UNINITIALIZED_INT;
    ctx->ga.TournamentSize     = PGA_UNINITIALIZED_INT;
    ctx->ga.TournamentWithRepl = PGA_UNINITIALIZED_INT;
    ctx->ga.RandomizeSelect    = PGA_UNINITIALIZED_INT;
    ctx->ga.TruncProportion    = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.RTRWindowSize      = PGA_UNINITIALIZED_INT;
    ctx->ga.FitnessType        = PGA_UNINITIALIZED_INT;
    ctx->ga.FitnessMinType     = PGA_UNINITIALIZED_INT;
    ctx->ga.MutationType       = PGA_UNINITIALIZED_INT;
    ctx->ga.MixingType         = PGA_UNINITIALIZED_INT;
    ctx->ga.MutateRealValue    = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.MutateIntegerValue = PGA_UNINITIALIZED_INT;
    ctx->ga.MutateBoundedFlag  = PGA_UNINITIALIZED_INT;
    ctx->ga.MutateBounceFlag   = PGA_UNINITIALIZED_INT;
    ctx->ga.MutatePolyEta      = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.MutatePolyValue    = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.DEVariant          = PGA_UNINITIALIZED_INT;
    ctx->ga.DENumDiffs         = PGA_UNINITIALIZED_INT;
    ctx->ga.DECrossoverType    = PGA_UNINITIALIZED_INT;
    ctx->ga.DEDitherPerIndividual = PGA_UNINITIALIZED_INT;
    ctx->ga.DEScaleFactor      = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.DEAuxFactor        = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.DECrossoverProb    = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.DEJitter           = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.DEProbabilityEO    = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.DEDither           = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.NoDuplicates       = PGA_UNINITIALIZED_INT;
    ctx->ga.MutationProb       = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.CrossoverProb      = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.UniformCrossProb   = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.PTournamentProb    = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.FitnessRankMax     = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.FitnessCmaxValue   = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.PopReplace         = PGA_UNINITIALIZED_INT;
    ctx->ga.iter               = 0;
    ctx->ga.last_iter          = -1;
    ctx->ga.perm_idx           = 0;
    ctx->ga.ItersOfSame        = 0;
    ctx->ga.PercentSame        = 0;
    ctx->ga.selected           = NULL;
    ctx->ga.SelectIndex        = 0;
    ctx->ga.restart            = PGA_UNINITIALIZED_INT;
    ctx->ga.restartFreq        = PGA_UNINITIALIZED_INT;
    ctx->ga.restartAlleleProb  = PGA_UNINITIALIZED_DOUBLE;
    ctx->ga.OutputFile         = stdout;
    ctx->ga.OutFileName        = NULL;
    ctx->ga.CustomData         = NULL;
    ctx->ga.NAMWindow          = PGA_UNINITIALIZED_INT;

    /* Fixed edges for Edge Crossover */
    ctx->ga.n_edges            = 0;
    ctx->ga.edges              = NULL;
    ctx->ga.r_edge             = NULL;
    ctx->ga.symmetric          = PGA_FALSE;

    /* NSGA-III */
    ctx->ga.nrefdirs       = 0;
    ctx->ga.nrefpoints     = 0;
    ctx->ga.refdirs        = NULL;
    ctx->ga.refpoints      = NULL;
    ctx->ga.ndir_npart     = 0;
    ctx->ga.dirscale       = 0;
    ctx->ga.extreme        = NULL;
    ctx->ga.extreme_valid  = PGA_FALSE;
    ctx->ga.utopian        = NULL;
    ctx->ga.utopian_valid  = PGA_FALSE;
    ctx->ga.nadir          = NULL;
    ctx->ga.worst          = NULL;
    ctx->ga.worst_valid    = PGA_FALSE;
    ctx->ga.normdirs       = NULL;
    ctx->ga.ndpoints       = 0;

    /* Operations */
    ctx->cops.CreateString      = NULL;
    ctx->cops.Mutation          = NULL;
    ctx->cops.Crossover         = NULL;
    ctx->cops.PrintString       = NULL;
    ctx->cops.CopyString        = NULL;
    ctx->cops.Duplicate         = NULL;
    ctx->cops.InitString        = NULL;
    ctx->cops.BuildDatatype     = NULL;
    ctx->cops.StopCond          = NULL;
    ctx->cops.EndOfGen          = NULL;
    ctx->cops.GeneDistance      = NULL;
    ctx->cops.PreEval           = NULL;
    ctx->cops.Hash              = NULL;
    ctx->cops.Serialize         = NULL;
    ctx->cops.Deserialize       = NULL;
    ctx->cops.SerializeFree     = NULL;
    ctx->cops.ChromFree         = NULL;

    ctx->fops.Mutation          = NULL;
    ctx->fops.Crossover         = NULL;
    ctx->fops.PrintString       = NULL;
    ctx->fops.CopyString        = NULL;
    ctx->fops.Duplicate         = NULL;
    ctx->fops.InitString        = NULL;
    ctx->fops.StopCond          = NULL;
    ctx->fops.EndOfGen          = NULL;
    ctx->fops.GeneDistance      = NULL;
    ctx->fops.PreEval           = NULL;
    ctx->fops.Hash              = NULL;

    /* Parallel */
    ctx->par.NumIslands        = PGA_UNINITIALIZED_INT;
    ctx->par.NumDemes          = PGA_UNINITIALIZED_INT;
    ctx->par.DefaultComm       = MPI_COMM_NULL;
#ifdef FAKE_MPI
    ctx->par.MPIStubLibrary    = PGA_TRUE;
#else
    ctx->par.MPIStubLibrary    = PGA_FALSE;
#endif

    /* Reporting */
    ctx->rep.PrintFreq         = PGA_UNINITIALIZED_INT;
    ctx->rep.PrintOptions      = 0;
    ctx->rep.MOPrecision       = 14;
    ctx->rep.Offline           = NULL;
    ctx->rep.Online            = NULL;
    ctx->rep.Average           = NULL;
    ctx->rep.Best              = NULL;
    ctx->rep.BestIdx           = NULL;
    ctx->rep.starttime         = PGA_UNINITIALIZED_INT;
    ctx->rep.validcount        = 0;
    ctx->rep.validonline       = 0;
    ctx->rep.validoffline      = 0;
    ctx->rep.nevals            = 0;

    /* System
     *
     *  If ctx->sys.UserFortran is not set to PGA_TRUE in pgacreate_ (the
     *  fortran stub to PGACreate), the user program is in C.
     */
    if (ctx->sys.UserFortran != PGA_TRUE) {
        ctx->sys.UserFortran  = PGA_FALSE;
    }
    ctx->sys.SetUpCalled       = PGA_FALSE;
    ctx->sys.PGAMaxInt         = INT_MAX;
    ctx->sys.PGAMinInt         = INT_MIN;
    ctx->sys.PGAMaxDouble      = DBL_MAX;
    ctx->sys.PGAMinDouble      = DBL_MIN;

    /* Debug */
    /* Set above before parsing command line arguments */

    /* Initialization */
    ctx->init.RandomInit        = PGA_UNINITIALIZED_INT;
    ctx->init.BinaryProbability = PGA_UNINITIALIZED_DOUBLE;
    ctx->init.RealType          = PGA_UNINITIALIZED_INT;
    ctx->init.IntegerType       = PGA_UNINITIALIZED_INT;
    ctx->init.CharacterType     = PGA_UNINITIALIZED_INT;
    ctx->init.RandomSeed        = PGA_UNINITIALIZED_INT;

    /*  Allocate and clear arrays to define the minimum and maximum values
     *  allowed by integer and real datatypes.
     */
    switch (datatype)
    {
    case PGA_DATATYPE_INTEGER:
         ctx->init.IntegerMax = (int *) malloc (len * sizeof (PGAInteger));
         if (!ctx->init.IntegerMax) {
             PGAError ( ctx, "PGACreate: No room to allocate:", PGA_FATAL
                      , PGA_CHAR, (void *) "ctx->init.IntegerMax"
                      );
         }
         ctx->init.IntegerMin = (int *) malloc (len * sizeof (PGAInteger));
         if (!ctx->init.IntegerMin) {
             PGAError ( ctx, "PGACreate: No room to allocate:", PGA_FATAL
                      , PGA_CHAR, (void *) "ctx->init.IntegerMin"
                      );
         }
         ctx->init.RealMax = NULL;
         ctx->init.RealMin = NULL;
         for (i=0; i<len; i++) {
              ctx->init.IntegerMin [i] = PGA_UNINITIALIZED_INT;
              ctx->init.IntegerMax [i] = PGA_UNINITIALIZED_INT;
         }
         break;
    case PGA_DATATYPE_REAL:
         ctx->init.RealMax = (PGAReal *) malloc (len * sizeof (PGAReal));
         if (!ctx->init.RealMax) {
             PGAError ( ctx, "PGACreate: No room to allocate:", PGA_FATAL
                      , PGA_CHAR, (void *) "ctx->init.RealMax"
                      );
         }
         ctx->init.RealMin = (PGAReal *) malloc(len * sizeof(PGAReal));
         if (!ctx->init.RealMin) {
              PGAError ( ctx, "PGACreate: No room to allocate:", PGA_FATAL
                       , PGA_CHAR, (void *) "ctx->init.RealMin"
                       );
         }
         ctx->init.IntegerMax = NULL;
         ctx->init.IntegerMin = NULL;
         for (i=0; i<len; i++) {
              ctx->init.RealMin [i] = PGA_UNINITIALIZED_DOUBLE;
              ctx->init.RealMax [i] = PGA_UNINITIALIZED_DOUBLE;
         }
         break;
    default:
         ctx->init.RealMax = NULL;
         ctx->init.RealMin = NULL;
         ctx->init.IntegerMax = NULL;
         ctx->init.IntegerMin = NULL;
         break;
    }

    PGADebugExited ("PGACreate");

    return (ctx);
}

/*!****************************************************************************
    \brief Initialize to zero various data structures of an individual
           and call the appropriate function to create and initialize
           the string for the specific data type.
    \ingroup internal
    \param   ctx       Context variable
    \param   p         String index
    \param   pop       Symbolic constant of the population string p is in
    \param   initflag  If the value is PGA_TRUE, the string is randomly
                       initialized, otherwise it is set to zero
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p;

       ...
       PGACreateIndividual (ctx, p, PGA_NEWPOP, PGA_TRUE);
    \endrst

******************************************************************************/
static
void PGACreateIndividual (PGAContext *ctx, int p, int pop, int initflag)
{
    PGAIndividual *ind = PGAGetIndividual(ctx, p, pop);

    ind->ctx              = ctx;
    ind->index            = p;
    ind->pop              = pop == PGA_OLDPOP ? ctx->ga.oldpop : ctx->ga.newpop;
    ind->evalue           = 0.0;
    ind->fitness          = 0.0;
    ind->evaluptodate     = PGA_FALSE;
    ind->auxtotal         = 0.0;
    ind->auxtotalok       = PGA_FALSE;
    ind->rank             = UINT_MAX;
    ind->crowding         = 0;
    ind->funcidx          = 0;
    ind->distance         = 0;
    ind->point_idx        = 0;
    ind->next_hash        = NULL;
    if (ctx->ga.NumAuxEval) {
        ind->auxeval = malloc (sizeof (double) * ctx->ga.NumAuxEval);
        if (ind->auxeval == NULL) {
            PGAError(ctx, "PGACreateIndividual: Failed to allocate auxeval",
                     PGA_FATAL, PGA_VOID, NULL);
        }
    } else {
        ind->auxeval = NULL;
    }
    if (ctx->ga.PopReplace == PGA_POPREPL_NSGA_III) {
        int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
        ind->normalized = malloc (sizeof (double) * dim);
        if (ind->normalized == NULL) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate normalized point");
        }
        memset (ind->normalized, 0, sizeof (double) * dim);
    } else {
        ind->normalized = NULL;
    }

    (*ctx->cops.CreateString)(ctx, p, pop, initflag);
    /* For these we need to set up the hash and avoid duplicates in the
     * first generation:
     */
    if (  initflag && ctx->ga.NoDuplicates
       && (  ctx->ga.PopReplace == PGA_POPREPL_RTR
          || ctx->ga.PopReplace == PGA_POPREPL_PAIRWISE_BEST
          )
       )
    {
        while (PGADuplicate (ctx, p, pop, pop)) {
            PGAChange (ctx, p, pop);
        }
        PGAHashIndividual (ctx, p, pop);
    }
}

/*!****************************************************************************
    \brief Allocate a population of individuals and
           calls \ref PGACreateIndividual to set up each one.
    \ingroup internal

    \param   ctx  context variable
    \param   pop  symbolic constant of the population to create
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGACreatePop (ctx, PGA_NEWPOP);

    \endrst

******************************************************************************/
static
void PGACreatePop (PGAContext *ctx, int pop)
{
    int p, flag = 0;

    switch (pop) {
    case PGA_OLDPOP:
        ctx->ga.oldpop = (PGAIndividual *)malloc
            (sizeof(PGAIndividual) * (ctx->ga.PopSize + 2));
        if (ctx->ga.oldpop == NULL) {
            PGAError
                ( ctx, "PGACreatePop: No room to allocate ctx->ga.oldpop"
                , PGA_FATAL, PGA_VOID, NULL
                );
        }
        memset
            (ctx->ga.oldpop, 0, sizeof(PGAIndividual) * (ctx->ga.PopSize + 2));
        flag = ctx->init.RandomInit;
        break;
    case PGA_NEWPOP:
        ctx->ga.newpop = (PGAIndividual *)malloc
            (sizeof (PGAIndividual) * (ctx->ga.PopSize + 2));
        if (ctx->ga.newpop == NULL) {
            PGAError
                ( ctx, "PGACreatePop: No room to allocate ctx->ga.newpop"
                , PGA_FATAL, PGA_VOID, NULL
                );
        }
        memset
            (ctx->ga.newpop, 0, sizeof(PGAIndividual) * (ctx->ga.PopSize + 2));
        flag = PGA_FALSE;
        break;
    default:
        PGAError
            ( ctx, "PGACreatePop: Invalid value of pop:"
            , PGA_FATAL, PGA_INT, (void *) &pop
            );
        break;
    };
    for (p=0; p<ctx->ga.PopSize; p++) {
        PGACreateIndividual (ctx, p, pop, flag);
    }
    PGACreateIndividual (ctx, PGA_TEMP1, pop, PGA_FALSE);
    PGACreateIndividual (ctx, PGA_TEMP2, pop, PGA_FALSE);
}

/*!****************************************************************************
    \brief Set all uninitialized variables to default values and
           initialize some internal arrays.
    \ingroup standard-api
    \param   ctx  context variable
    \return  Uninitialized values in the context variable are set to
             defaults, and set values are checked for legality

    \rst

    Description
    -----------

    Must be called after :c:func:`PGACreate` and before the GA is started.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       PGACreate (ctx, ...);
       //  Set options here
       PGASetUp (ctx);

    \endrst

******************************************************************************/
void PGASetUp (PGAContext *ctx)
{
    /*  These are for temporary storage of datatype specific functions.
     *  They allow some (understatement of the year!!) cleaning of the
     *  code below.
     */
    void    (*CreateString)(PGAContext *, int, int, int) = NULL;
    int     (*Mutation)(PGAContext *, int, int, double) = NULL;
    void    (*Crossover)(PGAContext *, int, int, int, int, int, int) = NULL;
    void    (*PrintString)(PGAContext *, FILE *, int, int) = NULL;
    void    (*CopyString)(PGAContext *, int, int, int, int) = NULL;
    int     (*Duplicate)(PGAContext *, int, int, int, int) = NULL;
    PGAHash (*Hash)(PGAContext *, int, int) = NULL;
    void    (*InitString)(PGAContext *, int, int) = NULL;
    double  (*GeneDist)(PGAContext *, int, int, int, int) = NULL;
    MPI_Datatype (*BuildDatatype)(PGAContext *, int, int) = NULL;
    int err=0, i;

    PGADebugEntered("PGASetUp");
    PGAFailIfSetUp("PGASetUp");

    if (ctx->ga.datatype != PGA_DATATYPE_USER
       && (  ctx->cops.Serialize
          || ctx->cops.Deserialize
          || ctx->cops.SerializeFree
          )
       )
    {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetUp: Serialize/Deserialize only with user defined datatype"
            );
    }
    if (ctx->ga.datatype != PGA_DATATYPE_USER && ctx->cops.ChromFree) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetUp: ChromFree only with user defined datatype"
            );
    }

    if (  ctx->ga.datatype == PGA_DATATYPE_BINARY
       && ctx->ga.tw       == PGA_UNINITIALIZED_INT
       )
    {
        PGAError( ctx,
                  "PGASetUp: Binary: Total Words (ctx->ga.tw) == UNINITIALIZED?"
                , PGA_FATAL, PGA_INT, (void *) &ctx->ga.tw
                );
    }

    if (  ctx->ga.datatype == PGA_DATATYPE_BINARY
       && ctx->ga.fw       == PGA_UNINITIALIZED_INT
       )
    {
        PGAError ( ctx,
                   "PGASetUp: Binary: Full Words (ctx->ga.fw) == UNINITIALIZED?"
                 , PGA_FATAL, PGA_INT,  (void *) &ctx->ga.fw
                 );
    }

    if (  ctx->ga.datatype == PGA_DATATYPE_BINARY
       && ctx->ga.eb       == PGA_UNINITIALIZED_INT
       )
    {
        PGAError ( ctx,
                   "PGASetUp: Binary: Empty Bits (ctx->ga.eb) == UNINITIALIZED?"
                 , PGA_FATAL, PGA_INT, (void *) &ctx->ga.eb
                 );
    }

    if (ctx->ga.NumAuxEval == PGA_UNINITIALIZED_INT) {
        ctx->ga.NumAuxEval = 0;
    }

    if (ctx->ga.NumConstraint == PGA_UNINITIALIZED_INT) {
        ctx->ga.NumConstraint = ctx->ga.NumAuxEval;
    }

    if (  ctx->ga.NumConstraint > ctx->ga.NumAuxEval
       || ctx->ga.NumConstraint < 0
       )
    {
         PGAError ( ctx, "PGASetUp: We need 0 <= NumConstraint <= NumAuxEval"
                  , PGA_FATAL, PGA_VOID, NULL
                  );
    }

    if (ctx->ga.PopReplace        == PGA_UNINITIALIZED_INT) {
        if (ctx->ga.NumAuxEval - ctx->ga.NumConstraint > 0) {
            ctx->ga.PopReplace    = PGA_POPREPL_NSGA_II;
        } else {
            ctx->ga.PopReplace    = PGA_POPREPL_BEST;
        }
    }

    if (ctx->ga.nrefpoints > 0 && ctx->ga.PopReplace != PGA_POPREPL_NSGA_III) {
        PGAErrorPrintf
            (ctx, PGA_FATAL, "PGASetUp: Reference points only for NSGA-III");
    }
    if (ctx->ga.nrefdirs > 0 && ctx->ga.PopReplace != PGA_POPREPL_NSGA_III) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetUp: Reference directions only for NSGA-III"
            );
    }

    if (  ctx->ga.NumAuxEval - ctx->ga.NumConstraint > 0
       && ctx->ga.PopReplace != PGA_POPREPL_NSGA_II
       && ctx->ga.PopReplace != PGA_POPREPL_NSGA_III
       )
    {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetUp: NumAuxEval=%d, NumConstraint=%d: Population "
              "replacement with multi-objective optimization must be NSGA-II"
              " or NSGA-III"
            , ctx->ga.NumAuxEval, ctx->ga.NumConstraint
            );
    }

    ctx->ga.extreme_valid = PGA_FALSE;
    ctx->ga.utopian_valid = PGA_FALSE;
    ctx->ga.worst_valid   = PGA_FALSE;
    if (ctx->ga.PopReplace == PGA_POPREPL_NSGA_III) {
        int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
        ctx->ga.extreme = malloc (sizeof (double) * dim * dim);
        if (ctx->ga.extreme == NULL) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate extreme point");
        }
        memset (ctx->ga.extreme, 0, sizeof (double) * dim * dim);
        ctx->ga.utopian = malloc (sizeof (double) * dim);
        if (ctx->ga.utopian == NULL) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate utopian point");
        }
        memset (ctx->ga.utopian, 0, sizeof (double) * dim);
        ctx->ga.nadir = malloc (sizeof (double) * dim);
        if (ctx->ga.nadir == NULL) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate nadir point");
        }
        memset (ctx->ga.nadir, 0, sizeof (double) * dim);
        ctx->ga.worst = malloc (sizeof (double) * dim);
        if (ctx->ga.worst == NULL) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate worst point");
        }
        memset (ctx->ga.worst, 0, sizeof (double) * dim);
        if (ctx->ga.nrefdirs == 0) {
            assert (ctx->ga.refdirs == NULL);
            if (ctx->ga.nrefpoints == 0) {
                assert (ctx->ga.refpoints == NULL);
                (void)LIN_dasdennis (dim, 2, &ctx->ga.refpoints, 0, 1, NULL);
                ctx->ga.nrefpoints = LIN_binom (dim + 2 - 1, 2);
                if (ctx->ga.refpoints == NULL) {
                    PGAErrorPrintf
                        (ctx, PGA_FATAL, "Cannot allocate ref points");
                }
            }
        } else {
            /* Allocate space for normalized reference directions */
            int npart = ctx->ga.ndir_npart;
            size_t lb = LIN_binom (dim + npart - 1, npart);
            size_t n = lb * ctx->ga.nrefdirs;
            /* Check that there is/was no overflow */
            assert (lb <= n);
            assert (n < SIZE_MAX / (sizeof (double) * dim));
            ctx->ga.normdirs = malloc (sizeof (double) * dim * n);
            if (ctx->ga.normdirs == NULL) {
                PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate normdirs");
            }
            ctx->ga.ndpoints = lb;
            if (ctx->ga.nrefpoints == 0) {
                assert (ctx->ga.refpoints == NULL);
                ctx->ga.nrefpoints = LIN_dasdennis
                    (dim, 1, &ctx->ga.refpoints, 0, 1, NULL);
                if (ctx->ga.refpoints == NULL) {
                    PGAErrorPrintf
                        (ctx, PGA_FATAL, "Cannot allocate ref points");
                }
            }
        }
    }

    /* Init PopSize from refdirs/refpoints if applicable, otherwise use
     * default of 100
     */
    if (ctx->ga.PopSize == PGA_UNINITIALIZED_INT) {
        int dim = ctx->ga.NumAuxEval - ctx->ga.NumConstraint + 1;
        int mod;
        size_t ps = 0;
        if (ctx->ga.nrefpoints > 0) {
            ps += ctx->ga.nrefpoints;
        }
        if (ctx->ga.nrefdirs > 0) {
            int npart = ctx->ga.ndir_npart;
            size_t lb = LIN_binom (dim + npart - 1, npart);
            size_t n = lb * ctx->ga.nrefdirs;
            ps += n;
        }
        /* Next multiple of 4, see NSGA-III papers */
        mod = ps % 4;
        if (mod) {
            ps += 4 - mod;
        }

        assert (ps < INT_MAX);
        if (ps == 0) {
            ps = 100;
        }
        ctx->ga.PopSize = ps;
    }

    if (ctx->ga.SumConstraints == PGA_UNINITIALIZED_INT) {
        ctx->ga.SumConstraints = PGA_TRUE;
    }

    if (ctx->ga.MaxIter == PGA_UNINITIALIZED_INT) {
        ctx->ga.MaxIter = 1000;
    }

    if (ctx->ga.MaxNoChange == PGA_UNINITIALIZED_INT) {
        ctx->ga.MaxNoChange = 100;
    }

    if (ctx->ga.MaxSimilarity == PGA_UNINITIALIZED_INT) {
        ctx->ga.MaxSimilarity  = 95;
    }

    if (ctx->ga.NumReplace == PGA_UNINITIALIZED_INT) {
        ctx->ga.NumReplace = (int) ceil (ctx->ga.PopSize * 0.1);
    }

    if (ctx->ga.NumReplace > ctx->ga.PopSize) {
        PGAError ( ctx, "PGASetUp: NumReplace > PopSize"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }

    if (ctx->ga.EpsilonGeneration == PGA_UNINITIALIZED_INT) {
        ctx->ga.EpsilonGeneration = 0;
    }
    if (ctx->ga.EpsilonGeneration > ctx->ga.MaxIter) {
        PGAError ( ctx, "PGASetUp: EpsilonGeneration > MaxIter"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }

    if (ctx->ga.EpsilonExponent == PGA_UNINITIALIZED_DOUBLE) {
        /* This is later used to initialize it to a dynamic value */
        ctx->ga.EpsilonExponent = 0;
        ctx->ga.EffEpsExponent  = 0;
        ctx->ga.EpsTLambda      = round (0.95 * ctx->ga.EpsilonGeneration);
        if (ctx->ga.EpsTLambda < 1) {
            ctx->ga.EpsTLambda = 1;
        }
    } else {
        ctx->ga.EffEpsExponent  = ctx->ga.EpsilonExponent;
        ctx->ga.EpsTLambda      = 0;
    }

    if (ctx->ga.EpsilonTheta == PGA_UNINITIALIZED_INT) {
        ctx->ga.EpsilonTheta = round (0.2 * ctx->ga.PopSize);
        if (!ctx->ga.EpsilonGeneration) {
            ctx->ga.EpsilonTheta = 0;
        }
    }

    if (ctx->ga.EpsilonTheta >= ctx->ga.PopSize - 1) {
        PGAError ( ctx, "PGASetUp: EpsilonTheta >= PopSize - 1"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }

    if (ctx->ga.CrossoverType == PGA_UNINITIALIZED_INT) {
        ctx->ga.CrossoverType = PGA_CROSSOVER_TWOPT;
    }

    if (  ctx->ga.CrossoverType == PGA_CROSSOVER_SBX
       && ctx->ga.datatype != PGA_DATATYPE_INTEGER
       && ctx->ga.datatype != PGA_DATATYPE_REAL
       )
    {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetUp: SBX crossover only for Integer and Real datatypes"
            );
    }

    if (ctx->ga.CrossBoundedFlag == PGA_UNINITIALIZED_INT) {
        ctx->ga.CrossBoundedFlag = PGA_FALSE;
    }

    if (ctx->ga.CrossBounceFlag == PGA_UNINITIALIZED_INT) {
        ctx->ga.CrossBounceFlag = PGA_FALSE;
    }

    if (ctx->ga.CrossSBXEta == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.CrossSBXEta = 2;
    }

    if (ctx->ga.CrossSBXOnce == PGA_UNINITIALIZED_INT) {
        ctx->ga.CrossSBXOnce = PGA_FALSE;
    }

    if (ctx->ga.SelectType == PGA_UNINITIALIZED_INT) {
        ctx->ga.SelectType = PGA_SELECT_TOURNAMENT;
    }

    if (ctx->ga.TournamentSize    == PGA_UNINITIALIZED_INT) {
        ctx->ga.TournamentSize     = 2;
    }

    if (ctx->ga.TournamentWithRepl == PGA_UNINITIALIZED_INT) {
        ctx->ga.TournamentWithRepl  = PGA_TRUE;
    }

    if (ctx->ga.RandomizeSelect == PGA_UNINITIALIZED_INT) {
        ctx->ga.RandomizeSelect  = PGA_FALSE;
    }

    if (ctx->ga.TruncProportion == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.TruncProportion  = 0.5;
    }

    if ( ctx->ga.RTRWindowSize     == PGA_UNINITIALIZED_INT) {
        int v = ctx->ga.PopSize / 20;
        if (ctx->ga.StringLen < ctx->ga.PopSize / 20) {
            v = ctx->ga.StringLen;
        }
        if (v < 1) {
            v = 1;
        }
        ctx->ga.RTRWindowSize = v;
    } else if (  ctx->ga.RTRWindowSize > ctx->ga.PopSize
              || ctx->ga.RTRWindowSize <= 0
              )
    {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetUp: required: 0 < RTR window size <= popsize"
            );
    }

    if ( ctx->ga.NumAuxEval > 0
       && (  ctx->ga.SelectType == PGA_SELECT_PROPORTIONAL
          || ctx->ga.SelectType == PGA_SELECT_SUS
          )
       )
    {
        if (ctx->ga.SelectType == PGA_SELECT_SUS) {
            PGAError ( ctx
                     , "PGASetUp: Auxiliary evaluation with default"
                       " string compare is incompatible with "
                       "Stochastic universal selection"
                     , PGA_FATAL, PGA_VOID, NULL
                     );
        } else {
            PGAError ( ctx
                     , "PGASetUp: Auxiliary evaluation with default"
                       " string compare is incompatible with "
                       "Proportional selection"
                     , PGA_FATAL, PGA_VOID, NULL
                     );
        }
    }

    if (ctx->ga.FitnessType       == PGA_UNINITIALIZED_INT) {
        ctx->ga.FitnessType        = PGA_FITNESS_RAW;
    }

    if (ctx->ga.FitnessMinType    == PGA_UNINITIALIZED_INT) {
        ctx->ga.FitnessMinType     = PGA_FITNESSMIN_CMAX;
    }

    if (ctx->ga.MixingType        == PGA_UNINITIALIZED_INT) {
        ctx->ga.MixingType         = PGA_MIX_MUTATE_OR_CROSS;
    }

    /* Require a gene length of at least two if we're doing crossover */
    if (  (ctx->ga.MixingType != PGA_MIX_MUTATE_ONLY && ctx->ga.StringLen <= 1)
       || ctx->ga.StringLen <= 0
       )
    {
        PGAError ( ctx, "PGACreate: Invalid value of StringLen:"
                 , PGA_FATAL, PGA_INT, (void *) &ctx->ga.StringLen
                 );
    }

    if (  ctx->ga.CrossoverType == PGA_CROSSOVER_TWOPT
       && ctx->ga.StringLen == 2 && ctx->ga.MixingType != PGA_MIX_MUTATE_ONLY
       )
    {
        PGAError ( ctx
                 , "PGASetUp: Invalid Crossover type for string of length 2"
                 , PGA_FATAL, PGA_INT, (void *) &ctx->ga.CrossoverType
                 );
    }

    if (ctx->ga.MutationProb      == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.MutationProb       = 1. / ctx->ga.StringLen;
    }

    if ( ctx->ga.MutationType      == PGA_UNINITIALIZED_INT) {
        switch (ctx->ga.datatype) {
        case PGA_DATATYPE_BINARY:
        case PGA_DATATYPE_CHARACTER:
        case PGA_DATATYPE_USER:
             /* leave PGA_UNINITIALIZED_INT for these data types */
             break;
        case PGA_DATATYPE_REAL:
             ctx->ga.MutationType   = PGA_MUTATION_GAUSSIAN;
             break;
        case PGA_DATATYPE_INTEGER:
             switch (ctx->init.IntegerType) {
                 case PGA_UNINITIALIZED_INT:
                 case PGA_IINIT_PERMUTE:
                     ctx->ga.MutationType   = PGA_MUTATION_PERMUTE;
                     break;
                 case PGA_IINIT_RANGE:
                     ctx->ga.MutationType   = PGA_MUTATION_RANGE;
                     break;
             }
             break;
        default:
            PGAError ( ctx, "PGASetup: Invalid value of ctx->ga.datatype:"
                     , PGA_FATAL, PGA_INT, (void *) &(ctx->ga.datatype)
                     );
        }
    }

    if (ctx->ga.MutateRealValue   == PGA_UNINITIALIZED_DOUBLE) {
        switch (ctx->ga.MutationType) {
        case PGA_MUTATION_GAUSSIAN:
            ctx->ga.MutateRealValue   = 0.1;
            break;
        case PGA_MUTATION_UNIFORM:
            ctx->ga.MutateRealValue   = 0.1;
            break;
        case PGA_MUTATION_CONSTANT:
            ctx->ga.MutateRealValue   = 0.01;
            break;
        case PGA_MUTATION_RANGE:
        default:
            ctx->ga.MutateRealValue   = 0.0;
        }
    }

    if (ctx->ga.MutateIntegerValue == PGA_UNINITIALIZED_INT) {
        ctx->ga.MutateIntegerValue  = 1;
    }

    if (ctx->ga.MutateBoundedFlag == PGA_UNINITIALIZED_INT) {
        ctx->ga.MutateBoundedFlag  = PGA_FALSE;
    }

    if (ctx->ga.MutateBounceFlag  == PGA_UNINITIALIZED_INT) {
        ctx->ga.MutateBounceFlag   = PGA_FALSE;
    }

    if (ctx->ga.MutatePolyEta  == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.MutatePolyEta   = 100;
    }

    if (ctx->ga.MutatePolyValue  == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.MutatePolyValue   = -1.0;
    }

    if (ctx->ga.DEVariant         == PGA_UNINITIALIZED_INT) {
        ctx->ga.DEVariant          = PGA_DE_VARIANT_RAND;
    }

    if (ctx->ga.DENumDiffs        == PGA_UNINITIALIZED_INT) {
        ctx->ga.DENumDiffs         = 1;
    }

    if (ctx->ga.DECrossoverType   == PGA_UNINITIALIZED_INT) {
        ctx->ga.DECrossoverType    = PGA_DE_CROSSOVER_BIN;
    }

    if (ctx->ga.DEDitherPerIndividual == PGA_UNINITIALIZED_INT) {
        ctx->ga.DEDitherPerIndividual  = PGA_FALSE;
    }

    if (ctx->ga.DEScaleFactor     == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.DEScaleFactor      = 0.9;
    }

    if (ctx->ga.DEAuxFactor       == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.DEAuxFactor        = 0.5 * (ctx->ga.DEScaleFactor + 1.0);
    }

    if (ctx->ga.DECrossoverProb   == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.DECrossoverProb    = 0.9;
    }

    if (ctx->ga.DEJitter          == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.DEJitter           = 0.0;
    }

    if (ctx->ga.DEDither          == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.DEDither           = 0.0;
    }

    if (ctx->ga.DEProbabilityEO   == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.DEProbabilityEO    = 0.5;
    }

    if (ctx->ga.NoDuplicates      == PGA_UNINITIALIZED_INT) {
        ctx->ga.NoDuplicates       = PGA_FALSE;
    }

    if ( ctx->ga.NoDuplicates
       && ((ctx->ga.StoppingRule & PGA_STOP_TOOSIMILAR) == PGA_STOP_TOOSIMILAR)
       )
    {
        PGAError ( ctx
                 , "PGASetUp: No Duplicates inconsistent with Stopping Rule:"
                 , PGA_FATAL, PGA_INT, (void *) &ctx->ga.StoppingRule
                 );
    }

    if (ctx->ga.CrossoverProb     == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.CrossoverProb      = 0.85;
    }

    if (ctx->ga.UniformCrossProb  == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.UniformCrossProb   = 0.6;
    }

    if (ctx->ga.PTournamentProb   == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.PTournamentProb    = 0.6;
    }

    if (ctx->ga.FitnessRankMax    == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.FitnessRankMax     = 1.2;
    }

    if (ctx->ga.FitnessCmaxValue  == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.FitnessCmaxValue   = 1.01;
    }

    if (ctx->ga.restart           == PGA_UNINITIALIZED_INT) {
        ctx->ga.restart            = PGA_FALSE;
    }

    if (ctx->ga.restartFreq       == PGA_UNINITIALIZED_INT) {
        ctx->ga.restartFreq        = 50;
    }

    if (ctx->ga.restartAlleleProb == PGA_UNINITIALIZED_DOUBLE) {
        ctx->ga.restartAlleleProb = 0.5;
    }

    /* Negative assortative mating */
    if (ctx->ga.NAMWindow == PGA_UNINITIALIZED_INT) {
        ctx->ga.NAMWindow = 1;
    }

    if (ctx->ga.NAMWindow > ctx->ga.PopSize - 2) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetUp: NAM window size must be <= PopSize - 2"
            );
    }


/* ops */
    /*  If no user supplied "done" function, use the built in one.
     *  No need to check EndOfGen; they only get called if they
     *  are defined.
     */
    if (  ((void *)ctx->cops.StopCond == (void *)PGADone)
       || ((void *)ctx->fops.StopCond == (void *)PGADone)
       )
    {
        PGAError ( ctx
                 , "PGASetUp: Using PGADone as the user stopping condition will"
                   " result in an infinite loop!"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }

    switch (ctx->ga.datatype) {
    case PGA_DATATYPE_BINARY:
        CreateString  = PGABinaryCreateString;
        BuildDatatype = PGABinaryBuildDatatype;
        Mutation      = PGABinaryMutation;

        switch (ctx->ga.CrossoverType) {
          case PGA_CROSSOVER_ONEPT:
            Crossover  = PGABinaryOneptCrossover;
            break;
          case PGA_CROSSOVER_TWOPT:
            Crossover  = PGABinaryTwoptCrossover;
            break;
          case PGA_CROSSOVER_UNIFORM:
            Crossover  = PGABinaryUniformCrossover;
            break;
        }
        PrintString    = PGABinaryPrintString;
        CopyString     = PGABinaryCopyString;
        Duplicate      = PGABinaryDuplicate;
        InitString     = PGABinaryInitString;
        GeneDist       = PGABinaryGeneDistance;
        Hash           = PGABinaryHash;
        break;
      case PGA_DATATYPE_INTEGER:
        CreateString   = PGAIntegerCreateString;
        BuildDatatype  = PGAIntegerBuildDatatype;
        Mutation       = PGAIntegerMutation;
        switch (ctx->ga.CrossoverType) {
          case PGA_CROSSOVER_ONEPT:
            Crossover  = PGAIntegerOneptCrossover;
            break;
          case PGA_CROSSOVER_TWOPT:
            Crossover  = PGAIntegerTwoptCrossover;
            break;
          case PGA_CROSSOVER_UNIFORM:
            Crossover  = PGAIntegerUniformCrossover;
            break;
          case PGA_CROSSOVER_SBX:
            Crossover  = PGAIntegerSBXCrossover;
            break;
          case PGA_CROSSOVER_EDGE:
            Crossover  = PGAIntegerEdgeCrossover;
            break;
        }
        PrintString    = PGAIntegerPrintString;
        CopyString     = PGAIntegerCopyString;
        Duplicate      = PGAIntegerDuplicate;
        InitString     = PGAIntegerInitString;
        GeneDist       = PGAIntegerGeneDistance;
        Hash           = PGAIntegerHash;
        break;
      case PGA_DATATYPE_REAL:
        CreateString   = PGARealCreateString;
        BuildDatatype  = PGARealBuildDatatype;
        Mutation       = PGARealMutation;
        switch (ctx->ga.CrossoverType) {
          case PGA_CROSSOVER_ONEPT:
            Crossover  = PGARealOneptCrossover;
            break;
          case PGA_CROSSOVER_TWOPT:
            Crossover  = PGARealTwoptCrossover;
            break;
          case PGA_CROSSOVER_UNIFORM:
            Crossover  = PGARealUniformCrossover;
            break;
          case PGA_CROSSOVER_SBX:
            Crossover  = PGARealSBXCrossover;
            break;
        }
        PrintString   = PGARealPrintString;
        CopyString    = PGARealCopyString;
        Duplicate     = PGARealDuplicate;
        InitString    = PGARealInitString;
        GeneDist      = PGARealGeneDistance;
        Hash          = PGARealHash;
        break;
      case PGA_DATATYPE_CHARACTER:
        CreateString  = PGACharacterCreateString;
        BuildDatatype = PGACharacterBuildDatatype;
        Mutation      = PGACharacterMutation;
        switch (ctx->ga.CrossoverType) {
          case PGA_CROSSOVER_ONEPT:
            Crossover  = PGACharacterOneptCrossover;
            break;
          case PGA_CROSSOVER_TWOPT:
            Crossover  = PGACharacterTwoptCrossover;
            break;
          case PGA_CROSSOVER_UNIFORM:
            Crossover  = PGACharacterUniformCrossover;
            break;
        }
        PrintString = PGACharacterPrintString;
        CopyString  = PGACharacterCopyString;
        Duplicate   = PGACharacterDuplicate;
        InitString  = PGACharacterInitString;
        GeneDist    = PGACharacterGeneDistance;
        Hash        = PGACharacterHash;
        break;
      case PGA_DATATYPE_USER:
        if (ctx->cops.CreateString == NULL) {
            PGAError
                ( ctx, "PGASetUp: User datatype needs CreateString function:"
                , PGA_FATAL, PGA_INT, (void *) &err
                );
        }
        if (ctx->cops.Mutation == NULL) {
            PGAError
                ( ctx, "PGASetUp: User datatype needs Mutation function:"
                , PGA_WARNING, PGA_INT, (void *) &err
                );
        }
        if (ctx->cops.Crossover == NULL) {
            PGAError
                ( ctx, "PGASetUp: User datatype needs Crossover function:"
                , PGA_WARNING, PGA_INT, (void *) &err
                );
        }
        if (ctx->cops.PrintString == NULL) {
            PGAError
                ( ctx, "PGASetUp: User datatype needs PrintString function:"
                , PGA_WARNING, PGA_INT, (void *) &err
                );
        }
        if (ctx->cops.Duplicate == NULL && ctx->ga.NoDuplicates) {
            PGAError
                ( ctx, "PGASetUp: User datatype needs Duplicate function:"
                , PGA_FATAL, PGA_INT, (void *) &err
                );
        }
        if (ctx->cops.CopyString == NULL) {
            PGAError
                ( ctx, "PGASetUp: User datatype needs CopyString function:"
                , PGA_FATAL, PGA_INT, (void *) &err
                );
        }
        if (  ( ctx->cops.Serialize && !ctx->cops.Deserialize)
           || (!ctx->cops.Serialize &&  ctx->cops.Deserialize)
           )
        {
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "PGASetUp: Serialize/Deserialize must be specified together"
                );
        }
        if (ctx->cops.Serialize && ctx->cops.BuildDatatype) {
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "PGASetUp: Serialize/Deserialize must not "
                  "specify BuildDatatype"
                );
        }
        if (ctx->cops.Serialize) {
            ctx->cops.BuildDatatype = PGASerializedBuildDatatype;
            /* If not set we asume it's allocated by malloc and friends */
            if (ctx->cops.SerializeFree == NULL) {
                ctx->cops.SerializeFree = &free;
            }
        }
        if (ctx->cops.BuildDatatype == NULL) {
             PGAError
                ( ctx, "PGASetUp: User datatype needs BuildDatatype function:"
                , PGA_FATAL, PGA_INT, (void *) &err
                );
        }
        if (  ctx->cops.GeneDistance == NULL
           && ctx->ga.PopReplace == PGA_POPREPL_RTR
           )
        {
             PGAError
                ( ctx, "PGASetUp: User datatype needs GeneDistance function:"
                , PGA_FATAL, PGA_INT, (void *) &err
                );
        }
        if (ctx->cops.Hash == NULL && ctx->ga.NoDuplicates) {
            PGAError
                ( ctx, "PGASetUp: User datatype needs Hash function:"
                , PGA_FATAL, PGA_INT, (void *) &err
                );
        }
        break;
    }
    if (ctx->cops.ChromFree == NULL) {
        ctx->cops.ChromFree = &ChromFree;
    }
    if ((ctx->cops.Mutation     == NULL) && (ctx->fops.Mutation    == NULL)) {
        ctx->cops.Mutation      = Mutation;
    }
    if ((ctx->cops.Crossover    == NULL) && (ctx->fops.Crossover   == NULL)) {
        ctx->cops.Crossover     = Crossover;
        if (Crossover == NULL) {
            PGAErrorPrintf (ctx, PGA_FATAL, "PGASetUp: No crossover specified");
        }
    }
    if ((ctx->cops.PrintString  == NULL) && (ctx->fops.PrintString == NULL)) {
        ctx->cops.PrintString   = PrintString;
    }
    if ((ctx->cops.Duplicate    == NULL) && (ctx->fops.Duplicate   == NULL)) {
        ctx->cops.Duplicate     = Duplicate;
    }
    if ((ctx->cops.InitString   == NULL) && (ctx->fops.InitString  == NULL)) {
        ctx->cops.InitString    = InitString;
    }
    if ((ctx->cops.GeneDistance == NULL) && (ctx->fops.GeneDistance == NULL)) {
        ctx->cops.GeneDistance  = GeneDist;
    }
    if (ctx->cops.CreateString  == NULL) {
        ctx->cops.CreateString  = CreateString;
    }
    if (ctx->cops.CopyString    == NULL) {
        ctx->cops.CopyString    = CopyString;
    }
    if (ctx->cops.BuildDatatype == NULL) {
        ctx->cops.BuildDatatype = BuildDatatype;
    }
    if ((ctx->cops.Hash == NULL) && (ctx->fops.Hash == NULL)) {
        ctx->cops.Hash = Hash;
    }

/* par */
    if (ctx->par.NumIslands == PGA_UNINITIALIZED_INT) {
        ctx->par.NumIslands = 1;
    }
    if (ctx->par.NumDemes == PGA_UNINITIALIZED_INT) {
        ctx->par.NumDemes = 1;
    }
    if (ctx->par.DefaultComm == MPI_COMM_NULL) {
        ctx->par.DefaultComm = MPI_COMM_WORLD;
    }



/* rep */
    if (ctx->rep.PrintFreq == PGA_UNINITIALIZED_INT) {
        ctx->rep.PrintFreq  = 10;
    }

/* sys */
    /* no more sets necessary here. */

/* debug */

/* init */
    if (ctx->init.RandomInit == PGA_UNINITIALIZED_INT) {
        ctx->init.RandomInit  = PGA_TRUE;
    }

    if (ctx->init.BinaryProbability == PGA_UNINITIALIZED_DOUBLE) {
        ctx->init.BinaryProbability  = 0.5;
    }

    if (ctx->init.RealType == PGA_UNINITIALIZED_INT) {
        ctx->init.RealType  = PGA_RINIT_RANGE;
    }
    if (ctx->init.IntegerType == PGA_UNINITIALIZED_INT) {
        ctx->init.IntegerType  = PGA_IINIT_PERMUTE;
    }
    if (ctx->init.CharacterType == PGA_UNINITIALIZED_INT) {
        ctx->init.CharacterType = PGA_CINIT_LOWER;
    }

    switch (ctx->ga.datatype)
    {
    case PGA_DATATYPE_INTEGER:
         for (i=0; i<ctx->ga.StringLen; i++) {
              if (ctx->init.IntegerMin[i] == PGA_UNINITIALIZED_INT) {
                  ctx->init.IntegerMin[i] = 0;
              }
              if (ctx->init.IntegerMax[i] == PGA_UNINITIALIZED_INT) {
                  ctx->init.IntegerMax[i] = ctx->ga.StringLen - 1;
              }
         }
         break;
    case PGA_DATATYPE_REAL:
         for (i=0; i<ctx->ga.StringLen; i++) {
              if (ctx->init.RealMin[i] == PGA_UNINITIALIZED_DOUBLE) {
                  ctx->init.RealMin[i] = 0.;
              }
              if (ctx->init.RealMax[i] == PGA_UNINITIALIZED_DOUBLE) {
                  ctx->init.RealMax[i] = 1.;
              }
         }
         break;
    }

    /* If a seed was not specified, get one from a time of day call */
    if (ctx->init.RandomSeed == PGA_UNINITIALIZED_INT) {
        ctx->init.RandomSeed = (int)time(NULL);
    }

    /* seed random number generator with this process' unique seed */
    ctx->init.RandomSeed += PGAGetRank (ctx, MPI_COMM_WORLD);
    PGARandom01 (ctx, ctx->init.RandomSeed);

    ctx->ga.selected = (int *)malloc( sizeof(int) * ctx->ga.PopSize );
    if (ctx->ga.selected == NULL) {
        PGAError ( ctx, "PGASetUp: No room to allocate ctx->ga.selected"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }

    ctx->ga.sorted = (int *)malloc( sizeof(int) * ctx->ga.PopSize );
    if (ctx->ga.sorted == NULL) {
        PGAError ( ctx, "PGASetUp: No room to allocate ctx->ga.sorted"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }

    ctx->scratch.intscratch = malloc( sizeof(int) * ctx->ga.PopSize );
    if (ctx->scratch.intscratch == NULL) {
        PGAError
            ( ctx, "PGASetUp: No room to allocate ctx->scratch.intscratch"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    ctx->scratch.permute = NULL;
    ctx->ga.perm_idx = 0;
    if (  ctx->ga.SelectType == PGA_SELECT_TOURNAMENT
       || ctx->ga.SelectType == PGA_SELECT_TRUNCATION
       )
    {
        ctx->scratch.permute = malloc (sizeof (int) * ctx->ga.PopSize);
        if (ctx->scratch.permute == NULL) {
            PGAError
                ( ctx, "PGASetUp: No room to allocate ctx->scratch.permute"
                , PGA_FATAL, PGA_VOID, NULL
                );
        }
        /* This forces a first shuffle */
        ctx->ga.perm_idx = ctx->ga.PopSize;
    }

    ctx->scratch.dblscratch = malloc (sizeof (double) * ctx->ga.PopSize);

    if (ctx->scratch.dblscratch == NULL) {
        PGAError
            ( ctx, "PGASetUp: No room to allocate ctx->scratch.dblscratch"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }

    /* If we're doing non-dominated sorting */
    if (ctx->ga.NumAuxEval - ctx->ga.NumConstraint >= 1) {
        int intsfor2pop = (ctx->ga.PopSize * 2 + WL - 1) / WL;
        ctx->scratch.dominance = malloc
            (sizeof (PGABinary) * intsfor2pop * 2 * ctx->ga.PopSize);
        if (ctx->scratch.dominance == NULL) {
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "PGASetUp: No room to allocate ctx->scratch.dominance"
                );
        }
    } else {
        ctx->scratch.dominance = NULL;
    }

    if (ctx->ga.NoDuplicates) {
        size_t hashsize = sizeof (PGAIndividual *) * ctx->ga.PopSize;
        ctx->scratch.hashed = malloc (hashsize);
        if (ctx->scratch.hashed == NULL) {
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "PGASetUp: No room to allocate ctx->scratch.hashed"
                );
        }
        memset (ctx->scratch.hashed, 0, hashsize);
    } else {
        ctx->scratch.hashed = NULL;
    }

    /* If the crossover type is Edge crossover */
    if (ctx->ga.CrossoverType == PGA_CROSSOVER_EDGE) {
        ctx->scratch.edgemap = malloc
            (sizeof (PGAInteger) * 4 * ctx->ga.StringLen);
        if (ctx->scratch.edgemap == NULL) {
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "PGASetUp: No room to allocate ctx->scratch.edgemap"
                );
        }
    } else {
        ctx->scratch.edgemap = NULL;
        if (ctx->ga.n_edges) {
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "PGASetUp: Fixed edges only for edge crossover"
                );
        }
    }
    if (ctx->ga.n_edges && ctx->init.IntegerMin [0] != 0) {
        PGAErrorPrintf
            (ctx, PGA_FATAL , "PGASetUp: Fixed edges only with IntegerMin=0");
    }
    ctx->scratch.serialization_size = 0;
    ctx->scratch.serialized = NULL;

    PGACreatePop (ctx, PGA_OLDPOP);
    PGACreatePop (ctx, PGA_NEWPOP);

    /* Now that we know how many eval functions we have, we can set up
     * the remaining fields of the report structure
     */
    ctx->rep.Offline = malloc (sizeof (double) * (1 + ctx->ga.NumAuxEval));
    if (ctx->rep.Offline == NULL) {
        PGAError ( ctx, "PGASetUp: No room to allocate rep.Offline"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }
    memset (ctx->rep.Offline, 0, sizeof (double) * (1 + ctx->ga.NumAuxEval));
    ctx->rep.Online = malloc (sizeof (double) * (1 + ctx->ga.NumAuxEval));
    if (ctx->rep.Online == NULL) {
        PGAError ( ctx, "PGASetUp: No room to allocate rep.Online"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }
    memset (ctx->rep.Online, 0, sizeof (double) * (1 + ctx->ga.NumAuxEval));
    ctx->rep.Average = malloc (sizeof (double) * (1 + ctx->ga.NumAuxEval));
    if (ctx->rep.Average == NULL) {
        PGAError ( ctx, "PGASetUp: No room to allocate rep.Average"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }
    for (i=0; i<=ctx->ga.NumAuxEval; i++) {
        ctx->rep.Average [i] = PGA_UNINITIALIZED_DOUBLE;
    }
    ctx->rep.Best = malloc (sizeof (double) * (1 + ctx->ga.NumAuxEval));
    if (ctx->rep.Best == NULL) {
        PGAError ( ctx, "PGASetUp: No room to allocate rep.Best"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }
    for (i=0; i<=ctx->ga.NumAuxEval; i++) {
        ctx->rep.Best [i] = PGA_UNINITIALIZED_DOUBLE;
    }
    ctx->rep.BestIdx = malloc (sizeof (int) * (1 + ctx->ga.NumAuxEval));
    if (ctx->rep.BestIdx == NULL) {
        PGAError ( ctx, "PGASetUp: No room to allocate rep.BestIdx"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }
    memset (ctx->rep.BestIdx, 0, sizeof (int) * (1 + ctx->ga.NumAuxEval));

    if (ctx->ga.OutFileName != NULL && PGAGetRank (ctx, MPI_COMM_WORLD) == 0) {
        ctx->ga.OutputFile = fopen (ctx->ga.OutFileName, "w");
        if (ctx->ga.OutputFile == NULL) {
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "Cannot open output file: %s", strerror (errno)
                );
        }
    }

    ctx->rep.starttime = time (NULL);

    /* This is done at the end to avoid trying to free unallocated memory */
    ctx->sys.SetUpCalled = PGA_TRUE;

    PGADebugExited ("PGASetUp");
}

/*!****************************************************************************
    \brief A boolean flag to indicate whether to randomly initialize
           alleles.
    \ingroup init
    \param   ctx   context variable
    \param   flag  indicates whether random initialization should be
                   performed
    \return  None

    \rst

    Description
    -----------

    Legal values are :c:macro:`PGA_TRUE` and :c:macro:`PGA_FALSE`.  Default
    is :c:macro:`PGA_TRUE`: randomly initialize alleles.

    Example
    -------

    Set the initialization routine to initialize all alleles to zero:

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetRandomInitFlag (ctx,PGA_FALSE);
    \endrst

******************************************************************************/
void PGASetRandomInitFlag (PGAContext *ctx, int flag)
{
    PGADebugEntered ("PGASetRandomInitFlag");
    PGAFailIfSetUp  ("PGASetRandomInitFlag");

    switch (flag) {
    case PGA_TRUE:
    case PGA_FALSE:
      ctx->init.RandomInit = flag;
      break;
    default:
      PGAError
        ( ctx, "PGASetRandomInitFlag: Invalid value of flag:"
        , PGA_FATAL, PGA_INT, (void *) &flag
        );
      break;
    }
    PGADebugExited ("PGASetRandomInitFlag");
}

/*!***************************************************************************
    \brief Return true/false to indicate whether or not alleles are
           randomly initialized.
    \ingroup query
    \param   ctx  context variable
    \return  Return true if alleles are randomly initialized

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       if (PGAGetRandomInitFlag (ctx)) {
           printf ("Population is randomly initialized\n");
       } else {
           printf ("Population initialized to zero\n");
       }

    \endrst

*****************************************************************************/
int PGAGetRandomInitFlag (PGAContext *ctx)
{
    PGADebugEntered("PGAGetRandomInitFlag");

    PGAFailIfNotSetUp("PGAGetRandomInitFlag");

    PGADebugExited("PGAGetRandomInitFlag");

    return(ctx->init.RandomInit);
}

/*!****************************************************************************
    \brief Initialize the number of auxiliary evaluations
    \ingroup init
    \param   ctx       context variable
    \param   n         Number of auxiliary evaluations
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetNumAuxEval (ctx, 5);
    \endrst

******************************************************************************/
void PGASetNumAuxEval (PGAContext *ctx, int n)
{
    PGADebugEntered("PGASetNumAuxEval");

    if (n <= 0) {
        PGAError(ctx, "PGASetNumAuxEval: Parameter needs to be positive",
                 PGA_FATAL, PGA_VOID, NULL);
    } else {
        ctx->ga.NumAuxEval = n;
    }

    PGADebugExited("PGASetNumAuxEval");
}

/*!****************************************************************************
    \brief Get the number of auxiliary evaluations
    \ingroup query

    \param   ctx       context variable
    \return  Number of auxiliary evaluations

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int num;

       ...
       num = PGAGetNumAuxEval (ctx);

    \endrst

******************************************************************************/
int PGAGetNumAuxEval (PGAContext *ctx)
{
    return ctx->ga.NumAuxEval;
}

/*!****************************************************************************
    \brief Initialize the number of constraints.
    \ingroup init
    \param   ctx       context variable
    \param   n         Number of constraints
    \return  None

    \rst

    Description
    -----------

    The maximum number of constraints (and the default) is the number of
    Auxiliary evaluations, see :c:func:`PGASetNumAuxEval`.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetNumConstraint (ctx, 5);

    \endrst

******************************************************************************/
void PGASetNumConstraint (PGAContext *ctx, int n)
{
    PGADebugEntered("PGASetNumConstraint");

    if (n < 0) {
        PGAError(ctx, "PGASetNumConstraint: Parameter needs to be positive",
                 PGA_FATAL, PGA_VOID, NULL);
    } else {
        ctx->ga.NumConstraint = n;
    }

    PGADebugExited("PGASetNumConstraint");
}

/*!****************************************************************************
    \brief Get the number of constraints.
    \ingroup query
    \param   ctx       context variable
    \return  Number of constraints

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int num;

       ...
       num = PGAGetNumConstraint (ctx);

    \endrst

******************************************************************************/
int PGAGetNumConstraint (PGAContext *ctx)
{
    return ctx->ga.NumConstraint;
}

/*!****************************************************************************
    \brief Configure if constraints are summed for minimization or use
           nondominated sorting for optimization.
    \ingroup init

    \param   ctx       context variable
    \param   n         PGA_TRUE or PGA_FALSE
    \return  None

    \rst

    Description
    -----------

    This only has an effect if the NSGA-II or NSGA-III replacement
    scheme is configured. In that case, using nondominated sorting for
    constraint optimization can be turned on by this option. By default
    constraints are summed and the sum in minimized. This is also done
    if no NSGA replacement scheme is in use. If summing is disabled by
    setting this configuration to :c:macro:`PGA_FALSE`, constraints are
    minimized using nondominated sorting. Note that nondominated sorting
    for constrains may not work very well on many problems.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetSumConstraintsFlag (ctx, PGA_FALSE);

    \endrst

******************************************************************************/
void PGASetSumConstraintsFlag (PGAContext *ctx, int n)
{
    PGADebugEntered("PGASetSumConstraintsFlag");

    if (n != PGA_FALSE && n != PGA_TRUE) {
        PGAError
            ( ctx, "PGASetSumConstraints: PGA_TRUE or PGA_FALSE required"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    ctx->ga.SumConstraints = n;

    PGADebugExited("PGASetSumConstraintsFlag");
}

/*!****************************************************************************
    \brief Query if constraints are summed or optimized by NSGA-II
           nondominated sorting
    \ingroup query
    \param   ctx       context variable
    \return  sum constraints flag

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int n;

       ...
       n = PGAGetSumConstraintsFlag (ctx);

    \endrst

******************************************************************************/
int PGAGetSumConstraintsFlag (PGAContext *ctx)
{
    return ctx->ga.SumConstraints;
}

/*!****************************************************************************
    \brief Configure the generation until which constraints are relaxed
           via the Epsilon Contraint method.
    \ingroup init


    \param   ctx       context variable
    \param   gen       Epsilon contraint generation
    \return  None

    \rst

    Description
    -----------

    The default is 0 (no Epsilon Contraint method is used).
    The parameter ``gen`` must be below the value set with
    :c:func:`PGASetMaxGAIterValue`.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetEpsilonGeneration (ctx, 50);

    \endrst

******************************************************************************/
void PGASetEpsilonGeneration (PGAContext *ctx, int gen)
{
    ctx->ga.EpsilonGeneration = gen;
}

/*!****************************************************************************
    \brief Get value of the generation until which constraints are
           relaxed via the Epsilon Contraint method.
    \ingroup query
    \param   ctx       context variable
    \return  The epsilon generation

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int n;

       ...
       n = PGAGetEpsilonGeneration (ctx);

    \endrst

******************************************************************************/
int PGAGetEpsilonGeneration (PGAContext *ctx)
{
    return ctx->ga.EpsilonGeneration;
}

/*!****************************************************************************
    \brief Configure the exponent of the term computing the epsilon
           value in each generation.
    \ingroup init

    \param   ctx       context variable
    \param   e         Exponent
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetEpsilonExponent (ctx, 5);
    \endrst

******************************************************************************/
void PGASetEpsilonExponent (PGAContext *ctx, double e)
{
    /* Note that in some papers the recommended minimum is 3 while in
     * others it is 2. We allow 2.
     */
    if (e < 2 || e > PGA_EPSILON_EXPONENT_MAX) {
        PGAError
            ( ctx, "PGASetEpsilonExponent: 2 <= e <= 10 required"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    ctx->ga.EpsilonExponent = e;
}

/*!****************************************************************************
    \brief Get the exponent used for epsilon constraints
    \ingroup query

    \param   ctx       context variable
    \return  The epsilon exponent

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double e;

       ...
       e = PGAGetEpsilonExponent (ctx);
    \endrst

******************************************************************************/
double PGAGetEpsilonExponent (PGAContext *ctx)
{
    return ctx->ga.EpsilonExponent;
}

/*!****************************************************************************
    \brief Set the theta generation value that is used for initializing
           the epsilon constraint.
    \ingroup init

    \param   ctx       context variable
    \param   theta     population index
    \return  None

    \rst

    Description
    -----------

    The initial population is sorted by constraint violation and the
    individual with index ``theta`` is used for initializing the initial
    epsilon for the epsilon constraint method.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetEpsilonTheta (ctx, n);
    \endrst

******************************************************************************/
void PGASetEpsilonTheta (PGAContext *ctx, int theta)
{
    if (theta < 1) {
        PGAError ( ctx, "PGASetUp: EpsilonTheta must be >= 1"
                 , PGA_FATAL, PGA_VOID, NULL
                 );
    }
    ctx->ga.EpsilonTheta = theta;
}

/*!****************************************************************************
    \brief Query the population index for initializing epsilon for the
           epsilon constraint method.
    \ingroup query

    \param   ctx       context variable
    \return  The epsilon theta

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int n;

       ...
       n = PGAGetEpsilonTheta (ctx);
    \endrst

******************************************************************************/
int PGAGetEpsilonTheta (PGAContext *ctx)
{
    return ctx->ga.EpsilonTheta;
}

/*!****************************************************************************
    \brief Set output file name for printing statistics etc.
    \ingroup init
    \param   ctx       context variable
    \param   name      output filename
    \return  None

    \rst

    Description
    -----------

    Note that the file is not immediately opened, instead it is later
    opened in the rank 0 individual.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       char *name = "output.file";

       ...
       PGASetOutputFile (ctx, name);
    \endrst

******************************************************************************/
void PGASetOutputFile (PGAContext *ctx, const char *name)
{
    char *n = malloc (strlen (name) + 1);
    if (n == NULL) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGASetOutputFile: Cannot allocate name"
            );
    }
    strcpy (n, name);
    ctx->ga.OutFileName = n;
}
