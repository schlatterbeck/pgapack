/******************************************************************************
*     FILE: pgapack.h: This file contains all constant and structure
*                      definitions definitions for PGAPack as well as all
*                      function declarations.
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
******************************************************************************/
/* Microsoft choses to arbitrarily deprecate some standard C-Library functions
 * (CRT stands for C runtime library) Disable deprecation warning
 * And Microsoft ist stuck in the 1980s with their C Compiler (well
 * technically it isn't a C-Compiler because it doesn't support the
 * latest version of the standard) because they do not support
 * dynamically allocated variable arrays (in the C standard since 1999).
 * So this hack uses _alloca (which doesn't return a NULL return value
 * but throws an exception if allocation fails). We do not bother to
 * catch the exception in the error-case.
 */
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define USE_ALLOCA
#define ALLOCA _alloca
#endif

#ifdef USE_ALLOCA
#ifndef ALLOCA
#define ALLOCA alloca
#endif
#define DECLARE_DYNARRAY(type, name, size) \
        type *name = ALLOCA (sizeof (type) * (size))
#define DECLARE_DYNARRAY2(type, name, size1, size2) \
        type *name = ALLOCA (sizeof (type) * (size1) * (size2))
#define DECLARE_DYNPTR(type, name, size) type *name
#define DEREF1_DYNPTR(name, size, idx) &(name[(idx) * (size)])
#define DEREF2_DYNPTR(name, size, idx1, idx2) name[(idx1) * (size) + (idx2)]
#else /* !USE_ALLOCA */
#define DECLARE_DYNARRAY(type, name, size) type name [size]
#define DECLARE_DYNARRAY2(type, name, size1, size2) type name [size1][size2]
#define DECLARE_DYNPTR(type, name, size) type (*name)[size]
#define DEREF1_DYNPTR(name, size, idx) name[idx]
#define DEREF2_DYNPTR(name, size, idx1, idx2) name[idx1][idx2]
#endif /* !USE_ALLOCA */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <mpi.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Hack for regression testing */
#ifndef STATIC
#define STATIC static
#endif


/*  If OPTIMIZED, remove various sanity checks, and debug output.
 *
 *  PGADebugEntered(a)      - Print a debug message about entering "a"
 *  PGADebugExited(a)       - Print a debug message about exiting "a"
 *  PGAFailIfNotSetup(a)    - Fail fatally if PGASetUp has not been called
 *  PGAFailIfSetup(a)       - Fail fatally if PGASetUp has been called
 *  PGACheckDataType(a, D)  - Fail fatally if the datatype is not D
 */
#ifndef OPTIMIZE
#define OPTIMIZE 0
#define PGADebugEntered(a) \
  PGADebugPrint(ctx, PGA_DEBUG_ENTERED, a, "Entered", PGA_VOID, NULL)
#define PGADebugExited(a) \
  PGADebugPrint(ctx, PGA_DEBUG_EXIT, a, "Exited", PGA_VOID, NULL)
#define PGAFailIfNotSetUp(Name)  \
  if (ctx->sys.SetUpCalled == PGA_FALSE) \
     PGAError(ctx, "PGASetUp must be called before " Name, \
              PGA_FATAL, PGA_VOID, NULL)
#define PGAFailIfSetUp(Name)  \
  if (ctx->sys.SetUpCalled == PGA_TRUE) \
     PGAError(ctx, Name " must be called before PGASetUp", PGA_FATAL, \
	      PGA_VOID, NULL)
#define PGACheckDataType(Name, DataType) \
  if (ctx->ga.datatype != DataType) \
     PGAError(ctx, "DataType is incorrect for " Name,PGA_FATAL,PGA_VOID,NULL)
#else
#undef OPTIMIZE
#define OPTIMIZE 1
#define PGADebugPrint(a,b,c,x,y,z) ;
#define PGADebugEntered(a) ;
#define PGADebugExited(a) ;
#define PGAFailIfNotSetUp(Name)  ;
#define PGAFailIfSetUp(Name)  ;
#define PGACheckDataType(Name, DataType) ;
#endif

/*****************************************
*           BINARY   MACROS              *
*****************************************/
/* Note: WL used to be a macro defined on the command-line.
 *       Since it is used only to represent the size of the PGABinary
 *       data type (an unsigned long) we can safely use sizeof here.
 */
#ifndef WL
#define WL (sizeof(PGABinary) * 8)
#endif

#define ONEL        ((PGABinary)1)
#define BIT(x,y)    (y&(ONEL<<((WL-1)-(x))))    /* true if bit is 1,         */
#define SET(x,y)    (y|=(ONEL<<((WL-1)-(x))))   /* set a bit to 1            */
#define UNSET(x,y)  (y&=(~(ONEL<<((WL-1)-(x)))))/* set a bit to 0, clear     */
#define TOGGLE(x,y) (y^=(ONEL<<((WL-1)-(x))))   /* complement a bits value   */
#define INDEX(ix,bx,bit,WL) ix=bit/WL;bx=bit%WL /* map global column (bit)   */
                                                /* to word (ix) and bit (bx) */

/* Used in PGAEvalCompare and others */
static inline int CMP (const double a, const double b)
{
    return (a < b ? -1 : (a > b ? 1 : 0));
}

/*****************************************
*       ABSTRACT DATA TYPES              *
*****************************************/
#define PGA_DATATYPE_BINARY      1    /* Array of unsigned ints            */
                                      /* parsed into bits    : binary.c    */
#define PGA_DATATYPE_INTEGER     2    /* Array of ints       : integer.c   */
#define PGA_DATATYPE_REAL        3    /* Array of doubles    : real.c      */
#define PGA_DATATYPE_CHARACTER   4    /* Array of characters : character.c */
#define PGA_DATATYPE_USER        5    /*  --user defined--                 */

#define PGABinary                unsigned long
#define PGAInteger               signed long int
#define PGAReal                  double
#define PGACharacter             signed char

#define PGA_INT                   1
#define PGA_DOUBLE                2
#define PGA_CHAR                  3
#define PGA_VOID                  4


/*****************************************
*              BOOLEANS                  *
*****************************************/
#define PGA_TRUE                   1
#define PGA_FALSE                  0

/* And some functions for manipulating bit arrays */
static inline void SET_BIT (PGABinary *bitptr, int idx)
{
    int iidx   = idx / WL;
    PGABinary ishift = 1lu << (idx % WL);
    bitptr [iidx] |= ishift;
}

static inline int GET_BIT (PGABinary *bitptr, int idx)
{
    int iidx   = idx / WL;
    PGABinary ishift = 1lu << (idx % WL);
    return bitptr [iidx] & ishift;
}

static inline void CLEAR_BIT (PGABinary *bitptr, int idx)
{
    int iidx   = idx / WL;
    PGABinary ishift = 1lu << (idx % WL);
    bitptr [iidx] &= ~ishift;
}


/*****************************************
*                FLAGS                   *
*****************************************/
#define PGA_FATAL                 1
#define PGA_WARNING               2

/*****************************************
*             MISC CONSTANT              *
*****************************************/
#define PGA_TEMP1                -1138
#define PGA_TEMP2                -4239

#define PGA_OLDPOP               -6728
#define PGA_NEWPOP               -8376

#define PGA_UNINITIALIZED_INT    -3827
#define PGA_UNINITIALIZED_DOUBLE -968.3827

/*****************************************
*        DEBUG LEVELS                    *
*****************************************/
#define PGA_DEBUG_ENTERED        12
#define PGA_DEBUG_EXIT           13
#define PGA_DEBUG_MALLOC         80
#define PGA_DEBUG_PRINTVAR       82
#define PGA_DEBUG_SEND           22
#define PGA_DEBUG_RECV           23
#define PGA_DEBUG_MAXFLAGS       1000

/*****************************************
*           DIRECTION                    *
*****************************************/
#define PGA_MAXIMIZE            1    /* specify direction for fitness calc  */
#define PGA_MINIMIZE            2    /* specify direction for fitness calc  */

/*****************************************
*         STOPPING CRITERIA              *
*****************************************/
#define PGA_STOP_MAXITER        1    /* Stop: for maximum iterations      */
#define PGA_STOP_NOCHANGE       2    /* Stop: no change in best string    */
#define PGA_STOP_TOOSIMILAR     4    /* Stop: homogeneous population      */

/*****************************************
*            CROSSOVER                   *
*****************************************/
#define PGA_CROSSOVER_ONEPT     1    /* One point crossover                */
#define PGA_CROSSOVER_TWOPT     2    /* Two point crossover                */
#define PGA_CROSSOVER_UNIFORM   3    /* Uniform   crossover                */
#define PGA_CROSSOVER_SBX       4    /* Simulated binary crossover (SBX)   */
#define PGA_CROSSOVER_EDGE      5    /* Edge Recombination                 */

/*****************************************
*            SELECTION                   *
*****************************************/
#define PGA_SELECT_PROPORTIONAL 1    /* proportional selection              */
#define PGA_SELECT_SUS          2    /* stochastic universal selection      */
#define PGA_SELECT_TOURNAMENT   3    /* tournament selection                */
#define PGA_SELECT_PTOURNAMENT  4    /* probabilistic tournament selection  */
#define PGA_SELECT_TRUNCATION   5    /* truncation selection                */
#define PGA_SELECT_LINEAR       6    /* linear selection                    */

/*****************************************
*            FITNESS                     *
*****************************************/
#define PGA_FITNESS_RAW         1    /* use raw fitness (evaluation)        */
#define PGA_FITNESS_NORMAL      2    /* linear normalization fitness        */
#define PGA_FITNESS_RANKING     3    /* linear ranking fitness              */

/*****************************************
*            FITNESS (MINIMIZATION)      *
*****************************************/
#define PGA_FITNESSMIN_RECIPROCAL  1 /* reciprocal fitness                  */
#define PGA_FITNESSMIN_CMAX        2 /* cmax fitness                        */

/*****************************************
*               MUTATION                 *
*****************************************/
#define PGA_MUTATION_CONSTANT   1    /* Real/Integer: Fixed value           */
#define PGA_MUTATION_RANGE      2    /* Real/Integer: Uniform range         */
#define PGA_MUTATION_UNIFORM    3    /* Real: +- Uniform random no.         */
#define PGA_MUTATION_GAUSSIAN   4    /* Real: +- Gaussian random no.        */
#define PGA_MUTATION_PERMUTE    5    /* Integer: Permutation (swap)         */
#define PGA_MUTATION_DE         6    /* Differential Evolution (only real)  */
#define PGA_MUTATION_POLY       7    /* Polynomial mutation                 */

/*****************************************
*               MIXING                   *
*****************************************/
/* This defines how mutation/crossover are combined (or not)
 * The MUTATE_AND_CROSS variant performs mutation only if crossover was
 * also performed. The TRADITIONAL variant performs mutation with the
 * configured probability and then mutates with the given probability
 * regardless if crossover was performed or not (this is the way all
 * traditional implementations of GA are handling it).
 * Note: This replaces the previous flags (PGASetMutationOrCrossoverFlag
 * and friends) which are still supported for legacy reasons.
 * The default is PGA_MIX_MUTATE_OR_CROSS also for legacy reasons.
 */
#define PGA_MIX_MUTATE_OR_CROSS   1    /* Either mutation or crossover      */
#define PGA_MIX_MUTATE_AND_CROSS  2    /* Mutation only if crossover        */
#define PGA_MIX_MUTATE_ONLY       3    /* Only mutation                     */
#define PGA_MIX_TRADITIONAL       4    /* Mutation after crossover          */

/*****************************************
* Differential Evolution Variant         *
*****************************************/
#define PGA_DE_VARIANT_RAND      1   /* Standard DE */
#define PGA_DE_VARIANT_BEST      2   /* Derive from best string */
#define PGA_DE_VARIANT_EITHER_OR 3   /* Either-or variant */

/*******************************************
* Differential Evolution Crossover Variant *
*******************************************/
#define PGA_DE_CROSSOVER_BIN      1  /* Standard DE binomial crossover */
#define PGA_DE_CROSSOVER_EXP      2  /* Exponential crossover */

/*****************************************
*        POPULATION REPLACEMENT          *
*****************************************/
#define PGA_POPREPL_BEST          1  /* Select best   string                */
#define PGA_POPREPL_RANDOM_NOREP  2  /* Select random string w/o replacement*/
#define PGA_POPREPL_RANDOM_REP    3  /* Select random string w/  replacement*/
#define PGA_POPREPL_RTR           4  /* Restricted tournament replacement   */
#define PGA_POPREPL_PAIRWISE_BEST 5  /* Pairwise compare old/newpop         */
#define PGA_POPREPL_NSGA_II       6  /* NSGA-II non-dominated sorting       */
#define PGA_POPREPL_NSGA_III      7  /* NSGA-III non-dominated sorting      */

/****************************************
 *       REPORT OPTIONS                 *
 ****************************************/
#define PGA_REPORT_ONLINE        1    /* Print the online analysis           */
#define PGA_REPORT_OFFLINE       2    /* Print the offline analysis          */
#define PGA_REPORT_HAMMING       4    /* Print the Hamming distance          */
#define PGA_REPORT_STRING        8    /* Print the string                    */
#define PGA_REPORT_WORST         16   /* Print the worst individual          */
#define PGA_REPORT_AVERAGE       32   /* Print the average of the population */

/*****************************************
*            RANDOMIZER                  *
*****************************************/
#define PGA_RINIT_PERCENT        1  /* real percent offset                   */
#define PGA_RINIT_RANGE          2  /* real range                            */
#define PGA_IINIT_PERMUTE        1  /* integer permutation                   */
#define PGA_IINIT_RANGE          2  /* integer range (nonunique)             */
#define PGA_CINIT_LOWER          1  /* all lowercase letters                 */
#define PGA_CINIT_UPPER          2  /* all uppercase letters                 */
#define PGA_CINIT_MIXED          3  /* both upper and lower case letters     */

/*****************************************
*         SET USER FUNCTION              *
*****************************************/
#define PGA_USERFUNCTION_CREATESTRING            1
#define PGA_USERFUNCTION_MUTATION                2
#define PGA_USERFUNCTION_CROSSOVER               3
#define PGA_USERFUNCTION_PRINTSTRING             4
#define PGA_USERFUNCTION_COPYSTRING              5
#define PGA_USERFUNCTION_DUPLICATE               6
#define PGA_USERFUNCTION_INITSTRING              7
#define PGA_USERFUNCTION_BUILDDATATYPE           8
#define PGA_USERFUNCTION_STOPCOND                9
#define PGA_USERFUNCTION_ENDOFGEN                10
#define PGA_USERFUNCTION_GEN_DIFFERENCE          11
#define PGA_USERFUNCTION_PRE_EVAL                12
#define PGA_NUM_USERFUNCTIONS                    12

/*****************************************
*           MPI SEND/RECV TAGS           *
*****************************************/
#define PGA_COMM_STRINGTOEVAL        1  /* MPI tag for sending string       */
#define PGA_COMM_EVALOFSTRING        2  /* MPI tag for returning evaluation */
#define PGA_COMM_DONEWITHEVALS       3  /* MPI tag for ending parallel eval */

/*****************************************
*           EPSILON CONSTRAINTS          *
*****************************************/
#define PGA_EPSILON_EXPONENT_MIN   3.0  /* minimum exponent cp from paper   */
#define PGA_EPSILON_EXPONENT_MAX  10.0  /* maximum exponent cp from paper   */


/*****************************************
*       INDIVIDUAL STRUCTURE             *
*****************************************/

typedef struct PGAIndividual {         /* primary population data structure */
  int                   index;         /* index of this indiv of the pop    */
  double                evalue;        /* evaluation function value         */
  double                fitness;       /* fitness    function value         */
  int                   evaluptodate;  /* flag whether evalue is current    */
  void                 *chrom;         /* pointer to the GA string          */
  double               *auxeval;       /* Auxiliary evaluations             */
  double                auxtotal;      /* Total aux evaluation              */
  int                   auxtotalok;    /* flag wether auxtotal is current*/
  unsigned int          rank;          /* Rank for dominance-sorting        */
  /* The following are not transmitted via MPI */
  struct PGAContext    *ctx;           /* Pointer to our PGAContext         */
  struct PGAIndividual *pop;           /* The population of this indiv.     */
  double                crowding;      /* Crowding metric for NSGA-II,-III  */
  int                   funcidx;       /* Temporary function index          */
  /* The following are for NSGA-III only */
  double               *normalized;    /* Normalized point for NSGA-III     */
  double                distance;      /* Distance to associated point      */
  int                   point_idx;     /* Index of associated point         */
} PGAIndividual;

/*****************************************
*      Fixed edges data structure        *
*****************************************/
typedef struct PGAFixedEdge_s {
    PGAInteger             lhs;
    PGAInteger             rhs;
    struct PGAFixedEdge_s *next;
    struct PGAFixedEdge_s *prev;
} PGAFixedEdge;


/*****************************************
*          GA ALGORITHM STRUCTURE        *
*****************************************/
typedef struct {
    int datatype;            /* data type: binary, integer, or real       */
    int optdir;              /* direction of optimization                 */
    int tw;                  /* total number of words, full + partial     */
    int fw;                  /* number of full (WL length) words          */
    int eb;                  /* number of extra bits in last NOT full word*/
    int PopSize;             /* Number of strings to use                  */
    int NumAuxEval;          /* Number of auxiliary evaluation values     */
    int NumConstraint;       /* Number of constraints                     */
    int SumConstraints;      /* PGA_TRUE if no dominance-sorting for
                                constraints                               */
    double Epsilon;          /* Current epsilon for eps constraints       */
    double Epsilon_0;        /* Initial Epsilon for eps constraints       */
    int EpsilonGeneration;   /* Max Generation for epsilon constraints    */
    double EpsilonExponent;  /* Exponent for tightening epsilon           */
    double EffEpsExponent;   /* Effective Exponent for tightening epsilon */
    int EpsTLambda;          /* Generation lambda for dynamic exponent    */
    int EpsilonTheta;        /* Theta best individual of epsilon init     */
    int StringLen;           /* string lengths                            */
    int StoppingRule;        /* Termination Criteria                      */
    int MaxIter;             /* Maximum number of iterations to run       */
    int MaxNoChange;         /* # of iters with no change before stopping */
    int MaxSimilarity;       /* % of pop the same before stopping         */
    int NumReplace;          /* Number of string to replace each gen      */
    int PopReplace;          /* Method of choosing ind.s to copy to newpop*/
    int iter;                /* iteration (generation) counter            */
    int ItersOfSame;         /* # iters with no change in best            */
    int PercentSame;         /* % of pop that is homogeneous              */
    int NoDuplicates;        /* Don't allow duplicate strings             */
    int CrossoverType;       /* Type of crossover for genetic algorithm   */
    int CrossBoundedFlag;    /* Confine alleles to given range (bound)    */
    int CrossBounceFlag;     /* Confine alleles to given range (bounce)   */
    double CrossSBXEta;      /* eta value for SBX                         */
    int CrossSBXOnce;        /* SBX probability once for whole string     */
    int SelectType;          /* Type of selection for genetic algorithm   */
    int SelectIndex;         /* index of Select for next two individuals  */
    int FitnessType;         /* Type of fitness transformation used       */
    int FitnessMinType;      /* Transformation for minimization problems  */
    int MixingType;          /* Combination of crossover/mutation         */
    int MutationType;        /* Type of mutation used                     */
    int MutateIntegerValue;  /* Multiplier to mutate Integer strings with */
    int MutateBoundedFlag;   /* Confine alleles to given range (bound)    */
    int MutateBounceFlag;    /* Confine alleles to given range (random)   */
    double MutatePolyEta;    /* Eta for polynomial mutation               */
    double MutatePolyValue;  /* Value for polynomial mutation             */
    double TournamentSize;   /* Number of participants in tournament      */
    int RTRWindowSize;       /* Window for restricted tournament select   */
    int TournamentWithRepl;  /* Tournament with / without replacement     */
    int RandomizeSelect;     /* Additional randomisation during select    */
    int DEVariant;           /* Differential evolution (DE) variant       */
    int DENumDiffs;          /* Number of differences for DE (1 or 2)     */
    int DECrossoverType;     /* Crossover type DE                         */
    int DEDitherPerIndividual; /* Per indidivual or per generation        */
    double DEDither;         /* Dither value centered around F            */
    double DEScaleFactor;    /* Scale Factor F for DE                     */
    double DEAuxFactor;      /* Auxiliary Factor K for DE                 */
    double DECrossoverProb;  /* Crossover probability Cr for DE           */
    double DEJitter;         /* Jitter interval DE (uniform dist)         */
    double DEProbabilityEO;  /* Either-Or-Probability for DE              */
    double MutateRealValue;  /* Multiplier to mutate Real strings with    */
    double MutationProb;     /* Starting mutation probability             */
    double CrossoverProb;    /* Crossover probability                     */
    double UniformCrossProb; /* Prob of bit select in uniform crossover   */
    double PTournamentProb;  /* Prob of selection in Prob. Tournament     */
    double FitnessRankMax;   /* MAX value for use in ranking              */
    double FitnessCmaxValue; /* Cmax value used to convert minimizations  */
    double restartAlleleProb;/* prob of changing an allele in a restart   */
    double TruncProportion;  /* proportion for truncation selection       */
    int restart;             /* whether to use the restart operator       */
    int restartFreq;         /* frequency with which to restart           */
    int *selected;           /* array of indices for selection            */
    int *sorted;             /* array of sorted individual indices        */
    size_t nrefdirs;         /* Number of reference directions            */
    void *refdirs;           /* Reference directions for NSGA-III         */
    void *normdirs;          /* normalized reference directions           */
    size_t ndpoints;         /* Number of points in refdir point cloud    */
    double dirscale;         /* Scale factor for reference directions     */
    int ndir_npart;          /* Number Das Dennis partitions for refdir   */
    size_t nrefpoints;       /* Number of reference points                */
    void *refpoints;         /* Ref points on normalized hyperplane       */
    void *extreme;           /* Extreme vector for NSGA-III               */
    int extreme_valid;       /* PGA_TRUE of above is valid                */
    double *utopian;         /* Utopian vector for NSGA-III               */
    int utopian_valid;       /* PGA_TRUE of above is valid                */
    double *nadir;           /* nadir point for NSGA-III                  */
    double *worst;           /* Worst point discovered so far             */
    int worst_valid;         /* PGA_TRUE of above is valid                */
    size_t n_edges;          /* Number of fixed edges                     */
    int symmetric;           /* Fixed edges are symmetric?                */
    PGAFixedEdge *edges;     /* Fixed edges for edge crossover            */
    PGAInteger (*r_edge)[2]; /* Right node + index into edges             */
    PGAIndividual *oldpop;   /* pointer to population (old)               */
    PGAIndividual *newpop;   /* pointer to population (new)               */
} PGAAlgorithm;



/*****************************************
*        OPERATIONS STRUCTURES           *
*****************************************/
typedef struct PGAContext PGAContext;
typedef struct {
    void         (*CreateString)(PGAContext *, int, int, int);
    int          (*Mutation)(PGAContext *, int, int, double);
    void         (*Crossover)(PGAContext *, int, int, int, int, int, int);
    void         (*PrintString)(PGAContext *, FILE *, int, int);
    void         (*CopyString)(PGAContext *, int, int, int, int);
    int          (*Duplicate)(PGAContext *, int, int, int, int);
    void         (*InitString)(PGAContext *, int, int);
    MPI_Datatype (*BuildDatatype)(PGAContext *, int, int);
    int          (*StopCond)(PGAContext *);
    void         (*EndOfGen)(PGAContext *);
    double       (*GeneDistance)(PGAContext *, int, int, int, int);
    void         (*PreEval)(PGAContext *, int);
    int          (*EvalCompare)(PGAContext *, int, int, int, int);
} PGACOperations;

typedef struct {
    int          (*Mutation)(void *, void *, void *, void *);
    void         (*Crossover)(void *, void *, void *, void *, void *, void *, void *);
    void         (*PrintString)(void *, void *, void *, void *);
    void         (*CopyString)(void *, void *, void *, void *, void *);
    int          (*Duplicate)(void *, void *, void *, void *, void *);
    void         (*InitString)(void *, void *, void *);
    int          (*StopCond)(void *);
    void         (*EndOfGen)(void *);
    double       (*GeneDistance)(void *, void *, void *, void *, void *);
    void         (*PreEval)(void *, void *);
    int          (*EvalCompare)(void *, void *, void *, void *, void *);
} PGAFortranOperations;

typedef struct sample_state_s {
    int         n;
    int         k;
    int         idx;
    PGAContext *ctx;
} PGASampleState;


/*****************************************
*          PARALLEL STRUCTURE            *
*****************************************/
typedef struct {
    int      MPIAlreadyInit;   /* Flag whether MPI was previously initialized*/
    int      NumIslands;       /* Number of islands in island model          */
    int      NumDemes;         /* Number of demes in neighborhood model      */
    MPI_Comm DefaultComm;      /* Default communicator for PGARun            */
    int      MPIStubLibrary;   /* Boolean: real or stub version of MPI       */
} PGAParallel;

/*****************************************
*          REPORT STRUCTURE              *
*****************************************/
typedef struct {
     int             PrintFreq;    /* How often to print statistics reports */
     int             PrintOptions;
     int             MOPrecision;  /* Precision of multi objective eval     */
     double         *Offline;      /* One value for each function           */
     double         *Online;       /* One value for each function           */
     double         *Average;      /* One value for each function           */
     double         *Best;         /* One value for each function           */
     double          MinSumConstr; /* Min Sum of violated constraints       */
     int            *BestIdx;      /* Indices of best individuals           */
     int             validcount;   /* # indivs w/o constraint-violations    */
     int             validonline;  /* cumulated validcount                  */
     int             validoffline; /* cumulated best indiv valid count      */
     int             nevals;       /* Number of evaluations                 */
     time_t          starttime;
} PGAReport;


/*****************************************
*          SYSTEM STRUCTURE              *
*****************************************/
typedef struct {
    int    UserFortran;             /* user routines in Fortran or C?        */
    int    SetUpCalled;             /* has PGASetUp been called?             */
    int    PGAMaxInt;               /* largest  int     of machine           */
    int    PGAMinInt;               /* smallest int     of machine           */
    double PGAMaxDouble;            /* largest  double  of machine           */
    double PGAMinDouble;            /* smallest double  of machine           */
} PGASystem;


/*****************************************
*          DEBUG STRUCTURE               *
*****************************************/
typedef struct {
    int PGADebugFlags[PGA_DEBUG_MAXFLAGS];
} PGADebug;

/*****************************************
*      INITIALIZATION STRUCTURE          *
*****************************************/
typedef struct {
    int    RandomInit;             /* flag whether to randomize strings    */
    double BinaryProbability;      /* probability that a Bit will be 1     */
    int    RealType;               /* type of real      initialization     */
    int    IntegerType;            /* type of integer   initialization     */
    int    CharacterType;          /* type of character initialization     */
    int    *IntegerMin;            /* minimum of range of integers         */
    int    *IntegerMax;            /* maximum of range of integers         */
    double *RealMin;               /* minimum of range of reals            */
    double *RealMax;               /* maximum of range of reals            */
    int    RandomSeed;             /* integer to seed random numbers with  */
} PGAInitialize;

/*****************************************
*      SCRATCH DATA STRUCTURES           *
*****************************************/
typedef struct {
    int          *intscratch;            /* integer-scratch space          */
    double       *dblscratch;            /* double- scratch space          */
    PGABinary    *dominance;             /* for dominance sorting          */
    PGAInteger   (*edgemap)[4];          /* For Edge Crossover             */
} PGAScratch;

/*****************************************
*          CONTEXT STRUCTURE             *
*****************************************/
struct PGAContext {
    PGAAlgorithm           ga;
    PGACOperations         cops;
    PGAFortranOperations   fops;
    PGAParallel            par;
    PGAReport              rep;
    PGASystem              sys;
    PGADebug               debug;
    PGAInitialize          init;
    PGAScratch             scratch;
};

/*****************************************
*          binary.c
*****************************************/

void PGASetBinaryAllele ( PGAContext *ctx, int p, int pop, int i, int val );
int PGAGetBinaryAllele ( PGAContext *ctx, int p, int pop, int i );
void PGASetBinaryInitProb ( PGAContext *ctx, double probability );
double PGAGetBinaryInitProb (PGAContext *ctx);
void PGABinaryCreateString(PGAContext *ctx, int p, int pop, int initflag);
int PGABinaryMutation( PGAContext *ctx, int p, int pop, double mr );
void PGABinaryOneptCrossover(PGAContext *ctx, int p1, int p2, int pop1, int c1,
                             int c2, int pop2);
void PGABinaryTwoptCrossover(PGAContext *ctx, int p1, int p2, int pop1, int c1,
                             int c2, int pop2);
void PGABinaryUniformCrossover(PGAContext *ctx, int p1, int p2, int pop1,
                               int c1, int c2, int pop2);
void PGABinaryPrintString( PGAContext *ctx, FILE *fp, int p, int pop );
void PGABinaryCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
int PGABinaryDuplicate( PGAContext *ctx, int p1, int pop1, int p2, int pop2);
void PGABinaryInitString(PGAContext *ctx, int p, int pop);
MPI_Datatype PGABinaryBuildDatatype(PGAContext *ctx, int p, int pop);
int PGABinaryHammingDistance ( PGAContext *ctx, PGABinary *s1, PGABinary *s2 );
void PGABinaryPrint( PGAContext *ctx, FILE *fp, PGABinary *chrom, int nb );
double PGABinaryGeneDistance (PGAContext *ctx, int p1, int pop1, int p2,
                              int pop2);

/*****************************************
*          char.c
*****************************************/

void PGASetCharacterAllele (PGAContext *ctx, int p, int pop, int i, char value);
char PGAGetCharacterAllele (PGAContext *ctx, int p, int pop, int i);
void PGASetCharacterInitType(PGAContext *ctx, int value);
void PGACharacterCreateString (PGAContext *ctx, int p, int pop, int InitFlag);
int PGACharacterMutation( PGAContext *ctx, int p, int pop, double mr );
void PGACharacterOneptCrossover(PGAContext *ctx, int p1, int p2, int pop1,
                                int c1, int c2, int pop2);
void PGACharacterTwoptCrossover( PGAContext *ctx, int p1, int p2, int pop1,
                              int c1, int c2, int pop2);
void PGACharacterUniformCrossover(PGAContext *ctx, int p1, int p2, int pop1,
                                int c1, int c2, int pop2);
void PGACharacterPrintString ( PGAContext *ctx, FILE *fp, int p, int pop);
void PGACharacterCopyString (PGAContext *ctx, int p1, int pop1, int p2,
                             int pop2);
int PGACharacterDuplicate( PGAContext *ctx, int p1, int pop1, int p2, int pop2);
void PGACharacterInitString(PGAContext *ctx, int p, int pop);
MPI_Datatype PGACharacterBuildDatatype(PGAContext *ctx, int p, int pop);
double PGACharacterGeneDistance (PGAContext *ctx, int p1, int pop1, int p2,
                                 int pop2);

/*****************************************
*          cmdline.c
*****************************************/

void PGAReadCmdLine( PGAContext *ctx, int *argc, char **argv );
void PGAParseDebugArg(PGAContext *ctx, char *st);
void PGAStripArgs(char **argv, int *argc, int *c, int num);

/*****************************************
*          create.c
*****************************************/

PGAContext *PGACreate ( int *argc, char **argv,
                        int datatype, int len, int maxormin );
void PGASetUp ( PGAContext *ctx );
void PGASetRandomInitFlag(PGAContext *ctx, int RandomBoolean);
int PGAGetRandomInitFlag (PGAContext *ctx);
void PGACreatePop (PGAContext *ctx, int pop);
void PGACreateIndividual (PGAContext *ctx, int p, int pop, int initflag);
void PGASetNumAuxEval (PGAContext *ctx, int n);
int PGAGetNumAuxEval (PGAContext *ctx);
void PGASetNumConstraint (PGAContext *ctx, int n);
int PGAGetNumConstraint (PGAContext *ctx);
void PGASetSumConstraintsFlag (PGAContext *ctx, int n);
int PGAGetSumConstraintsFlag (PGAContext *ctx);
void PGASetEpsilonGeneration (PGAContext *ctx, int generation);
int PGAGetEpsilonGeneration (PGAContext *ctx);
void PGASetEpsilonExponent (PGAContext *ctx, double exponent);
double PGAGetEpsilonExponent (PGAContext *ctx);
void PGASetEpsilonTheta (PGAContext *ctx, int theta);
int PGAGetEpsilonTheta (PGAContext *ctx);

/*****************************************
*          cross.c
*****************************************/

void PGACrossover ( PGAContext *ctx, int p1, int p2, int pop1,
                    int c1, int c2, int pop2 );
int PGAGetCrossoverType (PGAContext *ctx);
double PGAGetCrossoverProb (PGAContext *ctx);
double PGAGetUniformCrossoverProb (PGAContext *ctx);
void PGASetCrossoverType (PGAContext *ctx, int crossover_type);
void PGASetCrossoverProb( PGAContext *ctx, double crossover_prob);
void PGASetUniformCrossoverProb( PGAContext *ctx, double uniform_cross_prob);
void PGASetCrossoverBoundedFlag(PGAContext *ctx, int val);
int PGAGetCrossoverBoundedFlag(PGAContext *ctx);
void PGASetCrossoverBounceBackFlag(PGAContext *ctx, int val);
int PGAGetCrossoverBounceBackFlag(PGAContext *ctx);
void PGASetCrossoverSBXEta(PGAContext *ctx, double eta);
double PGAGetCrossoverSBXEta(PGAContext *ctx);
void PGASetCrossoverSBXOncePerString (PGAContext *ctx, int val);
int PGAGetCrossoverSBXOncePerString (PGAContext *ctx);
void PGACrossoverSBX
    (PGAContext *ctx, double p1, double p2, double u, double *c1, double *c2);

/*****************************************
*          debug.c
*****************************************/

void PGASortFuncNameIndex(PGAContext *ctx);
#if OPTIMIZE==0
void PGADebugPrint( PGAContext *ctx, int level, char *funcname,
                   char *msg, int datatype, void *data );
#endif
void PGASetDebugLevel(PGAContext *ctx, int level);
void PGAClearDebugLevel(PGAContext *ctx, int level);
void PGASetDebugLevelByName(PGAContext *ctx, char *funcname);
void PGAClearDebugLevelByName(PGAContext *ctx, char *funcname);
int PGAGetDebugLevelOfName(PGAContext *ctx, char *funcname);
int PGAGetDebugFlag(PGAContext *ctx, char *funcname);
void PGASetDebugFlag11(PGAContext *ctx, int Flag);
void PGASetDebugFlag20(PGAContext *ctx, int Flag);
void PGASetDebugFlag21(PGAContext *ctx, int Flag);
void PGASetDebugFlag30(PGAContext *ctx, int Flag);
void PGASetDebugFlag32(PGAContext *ctx, int Flag);
void PGASetDebugFlag34(PGAContext *ctx, int Flag);
void PGASetDebugFlag36(PGAContext *ctx, int Flag);
void PGASetDebugFlag40(PGAContext *ctx, int Flag);
void PGASetDebugFlag42(PGAContext *ctx, int Flag);
void PGASetDebugFlag44(PGAContext *ctx, int Flag);
void PGASetDebugFlag46(PGAContext *ctx, int Flag);
void PGASetDebugFlag48(PGAContext *ctx, int Flag);
void PGASetDebugFlag50(PGAContext *ctx, int Flag);
void PGASetDebugFlag52(PGAContext *ctx, int Flag);
void PGASetDebugFlag54(PGAContext *ctx, int Flag);
void PGASetDebugFlag56(PGAContext *ctx, int Flag);
void PGASetDebugFlag58(PGAContext *ctx, int Flag);
void PGASetDebugFlag60(PGAContext *ctx, int Flag);
void PGASetDebugFlag62(PGAContext *ctx, int Flag);
void PGASetDebugFlag64(PGAContext *ctx, int Flag);
void PGASetDebugFlag66(PGAContext *ctx, int Flag);
void PGAPrintDebugOptions(PGAContext *ctx);

/*****************************************
*          duplcate.c
*****************************************/

int PGADuplicate(PGAContext *ctx, int p, int pop1, int pop2, int n);
void PGAChange( PGAContext *ctx, int p, int pop );
void PGASetNoDuplicatesFlag( PGAContext *ctx, int no_dup);
int PGAGetNoDuplicatesFlag (PGAContext *ctx);

/*****************************************
*          evaluate.c
*****************************************/

/* Macro trickery for optional aux argument */
/* See https://stackoverflow.com/questions/1472138/c-default-arguments */
struct opt_ptr_ptr { const double **d; };
struct opt_ptr     { const double *d;  };
void _PGASetEvaluation
    (PGAContext *ctx, int p, int pop, double v, const double *aux);
static inline void _PGASetEvaluation_s
    (PGAContext *ctx, int p, int pop, double v, const struct opt_ptr opt_ptr)
{
    _PGASetEvaluation (ctx, p, pop, v, opt_ptr.d);
}
#define PGASetEvaluation(a, b, c, d, ...) \
    _PGASetEvaluation_s(a, b, c, d, (struct opt_ptr){__VA_ARGS__})
double _PGAGetEvaluation (PGAContext *ctx, int p, int pop, const double **aux);
static inline double _PGAGetEvaluation_s
    (PGAContext *ctx, int p, int pop, const struct opt_ptr_ptr opt_ptr_ptr)
{
    return _PGAGetEvaluation (ctx, p, pop, opt_ptr_ptr.d);
}
#define PGAGetEvaluation(a, b, c, ...) \
    _PGAGetEvaluation_s(a, b, c, (struct opt_ptr_ptr){__VA_ARGS__})
double *PGAGetAuxEvaluation (PGAContext *ctx, int p, int pop);
void PGASetEvaluationUpToDateFlag (PGAContext *ctx, int p, int pop, int status);
int PGAGetEvaluationUpToDateFlag (PGAContext *ctx, int p, int pop);
double PGAGetRealFromBinary
    ( PGAContext *ctx, int p, int pop
    , int start, int end, double lower, double upper
    );
double PGAGetRealFromGrayCode
    ( PGAContext *ctx, int p, int pop
    , int start, int end, double lower, double upper
    );
void PGAEncodeRealAsBinary
    ( PGAContext *ctx, int p, int pop
    , int start, int end, double low, double high, double val
    );
void PGAEncodeRealAsGrayCode
    ( PGAContext *ctx, int p, int pop
    , int start, int end, double low, double high, double val
    );
unsigned int PGAGetIntegerFromBinary
    (PGAContext *ctx, int p, int pop, int start, int end);
unsigned int PGAGetIntegerFromGrayCode
    (PGAContext *ctx, int p, int pop, int start, int end);
void PGAEncodeIntegerAsBinary
    (PGAContext *ctx, int p, int pop, int start, int end, unsigned int val);
void PGAEncodeIntegerAsGrayCode
    (PGAContext *ctx, int p, int pop, int start, int end, unsigned int val);
double PGAMapIntegerToReal
    (PGAContext *ctx, int v, int a, int b, double l, double u);
int PGAMapRealToInteger
    (PGAContext *ctx, double r, double l, double u, int a, int b);

/*****************************************
*          fitness.c
*****************************************/

void PGAFitness (PGAContext *ctx, int popindex);
int PGARank (PGAContext *ctx, int p, int *order, int n);
double PGAGetFitness (PGAContext *ctx, int p, int pop);
int PGAGetFitnessType (PGAContext *ctx);
int PGAGetFitnessMinType (PGAContext *ctx);
double PGAGetMaxFitnessRank (PGAContext *ctx);
void PGASetFitnessType (PGAContext *ctx, int fitness_type);
void PGASetFitnessMinType (PGAContext *ctx, int fitness_type);
void PGASetMaxFitnessRank (PGAContext *ctx, double fitness_rank_max);
void PGASetFitnessCmaxValue (PGAContext *ctx, double val);
double PGAGetFitnessCmaxValue (PGAContext *ctx);

/*****************************************
*          hamming.c
*****************************************/

double PGAHammingDistance( PGAContext *ctx, int popindex);

/*****************************************
*          integer.c
*****************************************/

void PGASetIntegerAllele (PGAContext *ctx, int p, int pop, int i, int value);
int PGAGetIntegerAllele (PGAContext *ctx, int p, int pop, int i);
void PGASetIntegerInitPermute (PGAContext *ctx, int min, int max);
void PGASetIntegerInitRange (PGAContext *ctx, const int *min, const int *max);
int PGAGetIntegerInitType (PGAContext *ctx);
int PGAGetMinIntegerInitValue (PGAContext *ctx, int i);
int PGAGetMaxIntegerInitValue (PGAContext *ctx, int i);
void PGAIntegerCreateString (PGAContext *ctx, int p, int pop, int InitFlag);
int PGAIntegerMutation (PGAContext *ctx, int p, int pop, double mr);
void PGAIntegerOneptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGAIntegerTwoptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGAIntegerUniformCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGAIntegerSBXCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGAIntegerEdgeCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGAIntegerSetFixedEdges
    (PGAContext *ctx, size_t n, PGAInteger (*edge)[2], int symmetric);
void PGAIntegerPrintString (PGAContext *ctx, FILE *fp, int p, int pop);
void PGAIntegerCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
int PGAIntegerDuplicate (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
void PGAIntegerInitString (PGAContext *ctx, int p, int pop);
MPI_Datatype PGAIntegerBuildDatatype (PGAContext *ctx, int p, int pop);
double PGAIntegerGeneDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
double PGAIntegerEuclidianDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2);

/*****************************************
*          linalg.c
*****************************************/

int LIN_solve (int n, void *a, double *b);
void LIN_print_matrix (int n, void *a);
void LIN_print_vector (int n, double *v);
int LIN_gcd (int a, int b);
size_t LIN_binom (int a, int b);
void LIN_normalize_to_refplane (int dim, double *v);
void LIN_dasdennis_allocated
    (int dim, int npart, double scale, double *dir, int npoints, void *mem);
int LIN_dasdennis
    (int dim, int npart, void *result, int nexist, double scale, double *dir);
double LIN_euclidian_distance (int dim, double *v1, double *v2);
double LIN_2norm (int dim, double *v);

/*****************************************
*          mpi_stub.c
*****************************************/


/*****************************************
*          mutation.c
*****************************************/

int PGAMutate(PGAContext *ctx, int p, int pop);
void PGASetMutationType( PGAContext *ctx, int mutation_type);
int PGAGetMutationType (PGAContext *ctx);
void PGASetMutationRealValue( PGAContext *ctx, double val);
double PGAGetMutationRealValue (PGAContext *ctx);
void PGASetMutationIntegerValue( PGAContext *ctx, int val);
int PGAGetMutationIntegerValue (PGAContext *ctx);
void PGASetMutationBoundedFlag(PGAContext *ctx, int val);
int PGAGetMutationBoundedFlag(PGAContext *ctx);
void PGASetMutationBounceBackFlag(PGAContext *ctx, int val);
int PGAGetMutationBounceBackFlag(PGAContext *ctx);
void PGASetMutationProb(PGAContext *ctx, double mutation_prob);
double PGAGetMutationProb (PGAContext *ctx);
void PGASetDEVariant (PGAContext *ctx, int val);
int PGAGetDEVariant (PGAContext *ctx);
void PGASetDENumDiffs (PGAContext *ctx, int val);
int PGAGetDENumDiffs (PGAContext *ctx);
void PGASetDEScaleFactor (PGAContext *ctx, double val);
double PGAGetDEScaleFactor (PGAContext *ctx);
void PGASetDEAuxFactor (PGAContext *ctx, double val);
double PGAGetDEAuxFactor (PGAContext *ctx);
void PGASetDECrossoverProb (PGAContext *ctx, double val);
double PGAGetDECrossoverProb (PGAContext *ctx);
void PGASetDEJitter (PGAContext *ctx, double val);
double PGAGetDEJitter (PGAContext *ctx);
void PGASetDEProbabilityEO (PGAContext *ctx, double val);
double PGAGetDEProbabilityEO (PGAContext *ctx);
void PGASetDECrossoverType (PGAContext *ctx, int val);
int PGAGetDECrossoverType (PGAContext *ctx);
void PGASetDEDither (PGAContext *ctx, double val);
double PGAGetDEDither (PGAContext *ctx);
void PGASetDEDitherPerIndividual (PGAContext *ctx, int val);
int PGAGetDEDitherPerIndividual (PGAContext *ctx);
void PGASetMutationPolyEta (PGAContext *ctx, double eta);
double PGAGetMutationPolyEta (PGAContext *ctx);
void PGASetMutationPolyValue (PGAContext *ctx, double c);
double PGAGetMutationPolyValue (PGAContext *ctx);

/*****************************************
*          parallel.c
*****************************************/

void PGARunGM(PGAContext *ctx, double (*f)(PGAContext *, int, int, double *),
	      MPI_Comm comm);
void PGAEvaluateSeq(PGAContext *ctx, int pop,
		    double (*f)(PGAContext *, int, int, double *));
void PGAEvaluateCoop(PGAContext *ctx, int pop,
		     double (*f)(PGAContext *, int, int, double *),
                     MPI_Comm comm);
void PGAEvaluateMS(PGAContext *ctx, int pop,
		   double (*f)(PGAContext *c, int p, int pop, double *),
                   MPI_Comm comm);
void PGAEvaluateSlave(PGAContext *ctx, int pop,
		      double (*f)(PGAContext *, int, int, double *),
                      MPI_Comm comm);
void PGAEvaluate(PGAContext *ctx, int pop,
		 double (*f)(PGAContext *, int, int, double *), MPI_Comm comm);
MPI_Datatype PGABuildDatatype(PGAContext *ctx, int p, int pop);
void PGASendIndividual(PGAContext *ctx, int p, int pop, int dest, int tag,
                       MPI_Comm comm);
void PGAReceiveIndividual(PGAContext *ctx, int p, int pop, int source, int tag,
                          MPI_Comm comm, MPI_Status *status);
void PGASendReceiveIndividual(PGAContext *ctx, int send_p, int send_pop,
                              int dest, int send_tag, int recv_p, int recv_pop,
                              int source, int recv_tag, MPI_Comm comm,
                              MPI_Status *status);
void PGARunIM(PGAContext *ctx, double (*f)
              (PGAContext *c, int p, int pop, double *), MPI_Comm tcomm);
void PGARunNM(PGAContext *ctx,
              double (*f)(PGAContext *c, int p, int pop, double *),
              MPI_Comm tcomm);
int PGAGetRank (PGAContext *ctx, MPI_Comm comm);
int PGAGetNumProcs (PGAContext *ctx, MPI_Comm comm);
void PGASetNumIslands( PGAContext *ctx, int n);
int PGAGetNumIslands (PGAContext *ctx);
void PGASetNumDemes( PGAContext *ctx, int numdemes);
int PGAGetNumDemes (PGAContext *ctx);
void PGASetCommunicator( PGAContext *ctx, MPI_Comm comm);
MPI_Comm PGAGetCommunicator( PGAContext *ctx);

/*****************************************
*          pga.c
*****************************************/

void PGARun (PGAContext *ctx,
    double (*evaluate)(PGAContext *c, int p, int pop, double *auxeval));
void PGARunMutationAndCrossover (PGAContext *ctx, int oldpop, int newpop);
void PGARunMutationOrCrossover (PGAContext *ctx, int oldpop, int newpop );
void PGARunMutationOnly (PGAContext *ctx, int oldpop, int newpop );
void PGAUpdateGeneration (PGAContext *ctx, MPI_Comm comm);
int PGAGetDataType (PGAContext *ctx);
int PGAGetOptDirFlag (PGAContext *ctx);
int PGAGetStringLength (PGAContext *ctx);
int PGAGetVariableStringLength (PGAContext *ctx, int p, int pop);
int PGAGetGAIterValue (PGAContext *ctx);
int PGAGetEvalCount (PGAContext *ctx);
void PGASetMutationOrCrossoverFlag (PGAContext *ctx, int flag);
void PGASetMutationAndCrossoverFlag (PGAContext *ctx, int flag);
void PGASetMutationOnlyFlag (PGAContext *ctx, int flag);
void PGASetMixingType (PGAContext *ctx, int type);
int PGAGetMutationOrCrossoverFlag (PGAContext *ctx);
int PGAGetMutationAndCrossoverFlag (PGAContext *ctx);
int PGAGetMutationOnlyFlag (PGAContext *ctx);
int PGAGetMixingType (PGAContext *ctx);

/*****************************************
*          pop.c
*****************************************/

void PGASortPop ( PGAContext *ctx, int pop );
int PGAGetPopSize (PGAContext *ctx);
int PGAGetNumReplaceValue (PGAContext *ctx);
int PGAGetPopReplaceType (PGAContext *ctx);
int PGAGetRTRWindowSize (PGAContext *ctx);
int PGAGetSortedPopIndex ( PGAContext *ctx, int n );
void PGASetPopSize (PGAContext *ctx, int popsize);
void PGASetNumReplaceValue( PGAContext *ctx, int pop_replace);
void PGASetPopReplaceType( PGAContext *ctx, int pop_replace);
void PGASetRTRWindowSize (PGAContext *ctx, int window);
void PGARestrictedTournamentReplacement (PGAContext *ctx);
void PGAPairwiseBestReplacement (PGAContext *ctx);
void PGA_NSGA_II_Replacement (PGAContext *ctx);
void PGA_NSGA_III_Replacement (PGAContext *ctx);
void PGASetReferencePoints (PGAContext *ctx, size_t npoints, void *points);
void PGASetReferenceDirections
    (PGAContext *ctx, size_t ndirs, void *dirs, int npart, double scale);


/*****************************************
*          random.c
*****************************************/

int PGARandomFlip ( PGAContext *ctx, double p );
int PGARandomInterval( PGAContext *ctx, int start, int end);
double PGARandom01( PGAContext *ctx, int newseed );
double PGARandomUniform( PGAContext *ctx, double start, double end);
double PGARandomGaussian( PGAContext *ctx, double mean, double sigma);
int PGAGetRandomSeed(PGAContext *ctx);
void PGASetRandomSeed(PGAContext *ctx, int seed);
void PGARandomSampleInit(PGAContext *ctx, PGASampleState *state, int k, int n);
int PGARandomNextSample(PGASampleState *state);

/*****************************************
*          real.c
*****************************************/

void PGASetRealAllele (PGAContext *ctx, int p, int pop, int i, double value);
double PGAGetRealAllele (PGAContext *ctx, int p, int pop, int i);
void PGASetRealInitPercent (PGAContext *ctx, double *median, double *percent);
void PGASetRealInitRange
    (PGAContext *ctx, const double *min, const double *max);
double PGAGetMinRealInitValue (PGAContext *ctx, int i);
double PGAGetMaxRealInitValue (PGAContext *ctx, int i);
int PGAGetRealInitType (PGAContext *ctx);
void PGARealCreateString (PGAContext *ctx, int p, int pop, int initflag);
int PGARealMutation (PGAContext *ctx, int p, int pop, double mr);
void PGARealOneptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGARealTwoptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGARealUniformCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGARealSBXCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2);
void PGARealPrintString (PGAContext *ctx, FILE *fp, int p, int pop);
void PGARealCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
int PGARealDuplicate (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
void PGARealInitString (PGAContext *ctx, int p, int pop);
MPI_Datatype PGARealBuildDatatype (PGAContext *ctx, int p, int pop);
double PGARealGeneDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
double PGARealEuclidianDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2);

/*****************************************
*          report.c
*****************************************/

void PGAPrintReport(PGAContext *ctx, FILE *fp, int pop);
void PGASetPrintOptions (PGAContext *ctx, int option);
void PGASetPrintFrequencyValue( PGAContext *ctx, int print_freq);
int PGAGetPrintFrequencyValue (PGAContext *ctx);
void PGASetMultiObjPrecision (PGAContext *ctx, int prec);
int PGAGetMultiObjPrecision (PGAContext *ctx);
void PGAPrintPopulation ( PGAContext *ctx, FILE *fp, int pop );
void PGAPrintIndividual ( PGAContext *ctx, FILE *fp, int p, int pop );
void PGAPrintString ( PGAContext *ctx, FILE *file, int p, int pop );
void PGAPrintContextVariable ( PGAContext *ctx, FILE *fp );

/*****************************************
*          restart.c
*****************************************/

void PGARestart(PGAContext *ctx, int source_pop, int dest_pop);
void PGASetRestartFlag(PGAContext *ctx, int val);
int PGAGetRestartFlag(PGAContext *ctx);
void PGASetRestartFrequencyValue(PGAContext *ctx, int numiter);
int PGAGetRestartFrequencyValue(PGAContext *ctx);
void PGASetRestartAlleleChangeProb(PGAContext *ctx, double prob);
double PGAGetRestartAlleleChangeProb(PGAContext *ctx);

/*****************************************
*          select.c
*****************************************/

void PGASelect (PGAContext *ctx, int popix);
int PGASelectNextIndex (PGAContext *ctx, int popix);
void PGASetSelectType (PGAContext *ctx, int select_type);
int PGAGetSelectType (PGAContext *ctx);
void PGASetPTournamentProb (PGAContext *ctx, double ptournament_prob);
double PGAGetPTournamentProb (PGAContext *ctx);
int PGASelectProportional (PGAContext *ctx, PGAIndividual *pop);
void PGASelectSUS (PGAContext *ctx, PGAIndividual *pop);
int PGASelectTournament (PGAContext *ctx, int pop);
int PGASelectPTournament (PGAContext *ctx, int pop);
int PGASelectTruncation (PGAContext *ctx, int pop);
int PGASelectLinear (PGAContext *ctx, PGAIndividual *pop);
void PGASetTournamentSize (PGAContext *ctx, double tournament_size);
double PGAGetTournamentSize (PGAContext *ctx);
void PGASetTournamentWithReplacement (PGAContext *ctx, int value);
int PGAGetTournamentWithReplacement (PGAContext *ctx);
void PGASetTruncationProportion (PGAContext *ctx, double proportion);
double PGAGetTruncationProportion (PGAContext *ctx);
void PGASetRandomizeSelect (PGAContext *ctx, int value);
int PGAGetRandomizeSelect (PGAContext *ctx);
double INDGetAuxTotal (PGAIndividual *ind);
double PGAGetAuxTotal (PGAContext *ctx, int p, int pop);

/*****************************************
*          stop.c
*****************************************/

int PGADone(PGAContext *ctx, MPI_Comm comm);
int PGACheckStoppingConditions( PGAContext *ctx);
void PGASetStoppingRuleType (PGAContext *ctx, int stoprule);
int PGAGetStoppingRuleType (PGAContext *ctx);
void PGASetMaxGAIterValue(PGAContext *ctx, int maxiter);
int PGAGetMaxGAIterValue (PGAContext *ctx);
void PGASetMaxNoChangeValue(PGAContext *ctx, int max_no_change);
void PGASetMaxSimilarityValue(PGAContext *ctx, int max_similarity);

/*****************************************
*          system.c
*****************************************/

void PGAError (PGAContext *ctx, char *msg,
               int level, int datatype, void *data);
void PGAErrorPrintf (PGAContext *ctx, int level, char *fmt, ...);
void PGADestroy (PGAContext *ctx);
int PGAGetMaxMachineIntValue (PGAContext *ctx);
int PGAGetMinMachineIntValue (PGAContext *ctx);
double PGAGetMaxMachineDoubleValue (PGAContext *ctx);
double PGAGetMinMachineDoubleValue (PGAContext *ctx);
void PGAUsage (PGAContext *ctx);
void PGAPrintVersionNumber (PGAContext *ctx);

/*****************************************
*          user.c
*****************************************/

void PGASetUserFunction(PGAContext *ctx, int constant, void *f);

/*****************************************
*          utility.c
*****************************************/

double PGAMean (PGAContext *ctx, double *a, int n);
double PGAStddev (PGAContext *ctx, double *a, int n, double mean);
int PGARound (PGAContext *ctx, double x);
void INDCopyIndividual (PGAIndividual *src, PGAIndividual *dst);
void PGACopyIndividual (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
int PGACheckSum (PGAContext *ctx, int p, int pop);
int PGAGetWorstIndex (PGAContext *ctx, int pop);
int PGAGetBestIndex (PGAContext *ctx, int pop);
double PGAGetBestReport (PGAContext *ctx, int pop, int idx);
int PGAGetBestReportIndex (PGAContext *ctx, int pop, int idx);
PGAIndividual *PGAGetIndividual (PGAContext *ctx, int p, int pop);
void PGAUpdateAverage (PGAContext *ctx, int pop);
void PGAUpdateBest (PGAContext *ctx, int pop);
void PGAUpdateOnline (PGAContext *ctx, int pop);
void PGAUpdateOffline (PGAContext *ctx, int pop);
int PGAComputeSimilarity (PGAContext *ctx, int popix);
int INDEvalCompare (PGAIndividual *ind1, PGAIndividual *ind2);
int PGAEvalCompare (PGAContext *ctx, int p1, int pop1, int p2, int pop2);
void PGAEvalSort (PGAContext *ctx, int pop, int *idx);
int PGAEvalSortHelper (const void *i1, const void *i2);
void PGAShuffle (PGAContext *ctx, int *list, int n);

#ifdef __cplusplus
}
#endif

