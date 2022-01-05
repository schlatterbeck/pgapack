c******************************************************************************
c     FILE: pgapackf.h
c
c     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
c              Brian P. Walenz
c*****************************************************************************/
c *** I/O FLAGS
CBARF !! is NULL ok?
      integer STDIN, STDOUT, STDERR, NULL
      parameter(STDIN=5, STDOUT=6, STDERR=6, NULL=0)


c *** ABSTRACT DATA TYPES
      integer PGA_DATATYPE_BINARY
      parameter( PGA_DATATYPE_BINARY =      1)
      integer PGA_DATATYPE_INTEGER
      parameter( PGA_DATATYPE_INTEGER =     2)
      integer PGA_DATATYPE_REAL
      parameter( PGA_DATATYPE_REAL =        3)
      integer PGA_DATATYPE_CHARACTER
      parameter( PGA_DATATYPE_CHARACTER =   4)
      integer PGA_DATATYPE_USER
      parameter( PGA_DATATYPE_USER =        5)

      integer PGA_INT
      parameter( PGA_INT =                   1)
      integer PGA_DOUBLE
      parameter( PGA_DOUBLE =                2)
      integer PGA_CHAR
      parameter( PGA_CHAR =                  3)
      integer PGA_VOID
      parameter( PGA_VOID =                  4)


c *** BOOLEANS &  FLAGS
      integer PGA_TRUE
      parameter ( PGA_TRUE =                   1)
      integer PGA_FALSE
      parameter ( PGA_FALSE =                  0)

      integer PGA_FATAL
      parameter ( PGA_FATAL =                 1)
      integer PGA_WARNING
      parameter ( PGA_WARNING =               2)


      integer PGA_UNINITIALIZED_INT
      parameter ( PGA_UNINITIALIZED_INT =    -3827)
      double precision PGA_UNINITIALIZED_DOUBLE
      parameter ( PGA_UNINITIALIZED_DOUBLE = -968.3827)

c *** TEMP & POP REFERENT CONSTANTS
      integer PGA_TEMP1
      parameter ( PGA_TEMP1 =                -1138)
      integer PGA_TEMP2
      parameter ( PGA_TEMP2 =                -4239)

      integer PGA_OLDPOP
      parameter ( PGA_OLDPOP =               -6728)
      integer PGA_NEWPOP
      parameter ( PGA_NEWPOP =               -8376)

c *** DEBUG LEVELS
      integer PGA_DEBUG_ENTERED
      parameter ( PGA_DEBUG_ENTERED =          12)
      integer PGA_DEBUG_EXIT
      parameter ( PGA_DEBUG_EXIT =             13)
      integer PGA_DEBUG_MALLOC
      parameter ( PGA_DEBUG_MALLOC =           80)
      integer PGA_DEBUG_PRINTVAR
      parameter ( PGA_DEBUG_PRINTVAR =         82)
      integer PGA_DEBUG_SEND
      parameter ( PGA_DEBUG_SEND =             22)
      integer PGA_DEBUG_RECV
      parameter ( PGA_DEBUG_RECV =             23)
      integer PGA_DEBUG_MAXFLAGS
      parameter ( PGA_DEBUG_MAXFLAGS =       1000)

c *** DIRECTION
      integer PGA_MAXIMIZE
      parameter ( PGA_MAXIMIZE =            1)
      integer PGA_MINIMIZE
      parameter ( PGA_MINIMIZE =            2)

c *** STOPPING CRITERIA
      integer PGA_STOP_MAXITER
      parameter ( PGA_STOP_MAXITER =        1)
      integer PGA_STOP_NOCHANGE
      parameter ( PGA_STOP_NOCHANGE =       2)
      integer PGA_STOP_TOOSIMILAR
      parameter ( PGA_STOP_TOOSIMILAR =     4)

c *** CROSSOVER
      integer PGA_CROSSOVER_ONEPT
      parameter ( PGA_CROSSOVER_ONEPT =     1)
      integer PGA_CROSSOVER_TWOPT
      parameter ( PGA_CROSSOVER_TWOPT =     2)
      integer PGA_CROSSOVER_UNIFORM
      parameter ( PGA_CROSSOVER_UNIFORM =   3)
      integer PGA_CROSSOVER_SBX
      parameter ( PGA_CROSSOVER_SBX =       4)

c *** SELECTION
      integer PGA_SELECT_PROPORTIONAL
      parameter ( PGA_SELECT_PROPORTIONAL = 1)
      integer PGA_SELECT_SUS
      parameter ( PGA_SELECT_SUS =          2)
      integer PGA_SELECT_TOURNAMENT
      parameter ( PGA_SELECT_TOURNAMENT =   3)
      integer PGA_SELECT_PTOURNAMENT
      parameter ( PGA_SELECT_PTOURNAMENT =  4)
      integer PGA_SELECT_TRUNCATION
      parameter ( PGA_SELECT_TRUNCATION =   5)
      integer PGA_SELECT_LINEAR
      parameter ( PGA_SELECT_LINEAR =       6)

c *** FITNESS
      integer PGA_FITNESS_RAW
      parameter ( PGA_FITNESS_RAW =         1)
      integer PGA_FITNESS_NORMAL
      parameter ( PGA_FITNESS_NORMAL =      2)
      integer PGA_FITNESS_RANKING
      parameter ( PGA_FITNESS_RANKING =     3)

c *** FITNESS (MINIMIZATION)
      integer PGA_FITNESSMIN_RECIPROCAL
      parameter ( PGA_FITNESSMIN_RECIPROCAL =  1)
      integer PGA_FITNESSMIN_CMAX
      parameter ( PGA_FITNESSMIN_CMAX =        2)

c *** MUTATION
      integer PGA_MUTATION_CONSTANT
      parameter ( PGA_MUTATION_CONSTANT =  1)
      integer PGA_MUTATION_RANGE
      parameter ( PGA_MUTATION_RANGE    =  2)
      integer PGA_MUTATION_UNIFORM
      parameter ( PGA_MUTATION_UNIFORM  =  3)
      integer PGA_MUTATION_GAUSSIAN
      parameter ( PGA_MUTATION_GAUSSIAN =  4)
      integer PGA_MUTATION_PERMUTE
      parameter ( PGA_MUTATION_PERMUTE  =  5)
      integer PGA_MUTATION_DE
      parameter ( PGA_MUTATION_DE       =  6)

c *** Differential Evolution Variant
      integer PGA_DE_VARIANT_RAND
      parameter ( PGA_DE_VARIANT_RAND      =  1)
      integer PGA_DE_VARIANT_BEST
      parameter ( PGA_DE_VARIANT_BEST      =  2)
      integer PGA_DE_VARIANT_EITHER_OR
      parameter ( PGA_DE_VARIANT_EITHER_OR =  3)

c *** Differential Evolution Crossover Variant
      integer PGA_DE_CROSSOVER_BIN
      parameter ( PGA_DE_CROSSOVER_BIN =  1)
      integer PGA_DE_CROSSOVER_EXP
      parameter ( PGA_DE_CROSSOVER_EXP =  2)

c *** POPULATION REPLACEMENT
      integer PGA_POPREPL_BEST
      parameter ( PGA_POPREPL_BEST =            1)
      integer PGA_POPREPL_RANDOM_NOREP
      parameter ( PGA_POPREPL_RANDOM_NOREP =    2)
      integer PGA_POPREPL_RANDOM_REP
      parameter ( PGA_POPREPL_RANDOM_REP =      3)
      integer PGA_POPREPL_RTR
      parameter ( PGA_POPREPL_RTR =             4)
      integer PGA_POPREPL_PAIRWISE_BEST
      parameter ( PGA_POPREPL_PAIRWISE_BEST =   5)
      integer PGA_POPREPL_NSGA_II
      parameter ( PGA_POPREPL_NSGA_II =         6)

c *** REPORT OPTIONS
      integer PGA_REPORT_ONLINE
      parameter ( PGA_REPORT_ONLINE =   1 )
      integer PGA_REPORT_OFFLINE
      parameter ( PGA_REPORT_OFFLINE =  2 )
      integer PGA_REPORT_HAMMING
      parameter ( PGA_REPORT_HAMMING =  4 )
      integer PGA_REPORT_STRING
      parameter ( PGA_REPORT_STRING =   8 )
      integer PGA_REPORT_WORST
      parameter ( PGA_REPORT_WORST =   16 )
      integer PGA_REPORT_AVERAGE
      parameter ( PGA_REPORT_AVERAGE = 32 )

c *** RANDOMIZER
      integer PGA_IINIT_PERMUTE
      parameter ( PGA_IINIT_PERMUTE =             1)
      integer PGA_IINIT_RANGE
      parameter ( PGA_IINIT_RANGE =               2)
      integer PGA_CINIT_LOWER
      parameter ( PGA_CINIT_LOWER =               1)
      integer PGA_CINIT_UPPER
      parameter ( PGA_CINIT_UPPER =               2)
      integer PGA_CINIT_MIXED
      parameter ( PGA_CINIT_MIXED =               3)

c *** SET USER FUNCTION
      integer PGA_USERFUNCTION_CREATESTRING
      parameter ( PGA_USERFUNCTION_CREATESTRING =       1)
      integer PGA_USERFUNCTION_MUTATION
      parameter ( PGA_USERFUNCTION_MUTATION =           2)
      integer PGA_USERFUNCTION_CROSSOVER
      parameter ( PGA_USERFUNCTION_CROSSOVER =          3)
      integer PGA_USERFUNCTION_PRINTSTRING
      parameter ( PGA_USERFUNCTION_PRINTSTRING  =       4)
      integer PGA_USERFUNCTION_COPYSTRING
      parameter ( PGA_USERFUNCTION_COPYSTRING =         5)
      integer PGA_USERFUNCTION_DUPLICATE
      parameter ( PGA_USERFUNCTION_DUPLICATE =          6)
      integer PGA_USERFUNCTION_INITSTRING
      parameter ( PGA_USERFUNCTION_INITSTRING =         7)
      integer PGA_USERFUNCTION_BUILDDATATYPE
      parameter ( PGA_USERFUNCTION_BUILDDATATYPE =      8)
      integer PGA_USERFUNCTION_STOPCOND
      parameter ( PGA_USERFUNCTION_STOPCOND =           9)
      integer PGA_USERFUNCTION_ENDOFGEN
      parameter ( PGA_USERFUNCTION_ENDOFGEN =          10)
      integer PGA_USERFUNCTION_GEN_DIFFERENCE
      parameter ( PGA_USERFUNCTION_GEN_DIFFERENCE =    11)
      integer PGA_USERFUNCTION_PRE_EVAL
      parameter ( PGA_USERFUNCTION_PRE_EVAL =          12)

c *** TAGS
      integer PGA_COMM_STRINGTOEVAL
      parameter ( PGA_COMM_STRINGTOEVAL =              1)
      integer PGA_COMM_EVALOFSTRING
      parameter ( PGA_COMM_EVALOFSTRING =              2)
      integer PGA_COMM_DONEWITHEVALS
      parameter ( PGA_COMM_DONEWITHEVALS =             3)
c *** binary
      integer PGAGetBinaryAllele
      external PGAGetBinaryAllele
      double precision PGAGetBinaryInitProb
      external PGAGetBinaryInitProb
c *** char
      character PGAGetCharacterAllele
      external PGAGetCharacterAllele
c *** create
      integer(8) PGACreate
      external PGACreate
      integer PGAGetRandomInitFlag
      external PGAGetRandomInitFlag
      integer PGAGetNumAuxEval
      external PGAGetNumAuxEval
      integer PGAGetNumConstraint
      external PGAGetNumConstraint
      integer PGAGetSumConstraintsFlag
      external PGAGetSumConstraintsFlag
c *** cross
      integer PGAGetCrossoverType
      external PGAGetCrossoverType
      double precision PGAGetCrossoverProb
      external PGAGetCrossoverProb
      double precision PGAGetUniformCrossoverProb
      external PGAGetUniformCrossoverProb
      integer PGAGetCrossoverBoundedFlag
      external PGAGetCrossoverBoundedFlag
      integer PGAGetCrossoverBounceBackFlag
      external PGAGetCrossoverBounceBackFlag
      integer PGASetCrossoverSBXOncePerString
      external PGASetCrossoverSBXOncePerString
      double precision PGAGetCrossoverSBXNu
      external PGAGetCrossoverSBXNu
c *** duplcate
      integer PGADuplicate
      external PGADuplicate
      integer PGAGetNoDuplicatesFlag
      external PGAGetNoDuplicatesFlag
c *** evaluate
      double precision PGAGetEvaluation
      external PGAGetEvaluation
      double precision PGAGetEvaluationAux
      external PGAGetEvaluationAux
      integer PGAGetEvaluationUpToDateFlag
      external PGAGetEvaluationUpToDateFlag
      double precision PGAGetRealFromBinary
      external PGAGetRealFromBinary
      double precision PGAGetRealFromGrayCode
      external PGAGetRealFromGrayCode
      integer PGAGetIntegerFromBinary
      external PGAGetIntegerFromBinary
      integer PGAGetIntegerFromGrayCode
      external PGAGetIntegerFromGrayCode
c *** fitness
      integer PGARank
      external PGARank
      double precision PGAGetFitness
      external PGAGetFitness
      integer PGAGetFitnessType
      external PGAGetFitnessType
      integer PGAGetFitnessMinType
      external PGAGetFitnessMinType
      double precision PGAGetMaxFitnessRank
      external PGAGetMaxFitnessRank
      double precision PGAGetFitnessCmaxValue
      external PGAGetFitnessCmaxValue
c *** hamming
      double precision PGAHammingDistance
      external PGAHammingDistance
c *** integer
      integer PGAGetIntegerAllele
      external PGAGetIntegerAllele
      integer PGAGetIntegerInitType
      external PGAGetIntegerInitType
      integer PGAGetMinIntegerInitValue
      external PGAGetMinIntegerInitValue
      integer PGAGetMaxIntegerInitValue
      external PGAGetMaxIntegerInitValue
c *** mutation
      integer PGAMutate
      external PGAMutate
      integer PGAGetMutationType
      external PGAGetMutationType
      double precision PGAGetMutationRealValue
      external PGAGetMutationRealValue
      integer PGAGetMutationIntegerValue
      external PGAGetMutationIntegerValue
      integer PGAGetMutationBoundedFlag
      external PGAGetMutationBoundedFlag
      integer PGAGetMutationBounceBackFlag
      external PGAGetMutationBounceBackFlag
      double precision PGAGetMutationProb
      external PGAGetMutationProb
      integer PGAGetDEVariant
      external PGAGetDEVariant
      integer PGAGetDENumDiffs
      external PGAGetDENumDiffs
      double precision PGAGetDEScaleFactor
      external PGAGetDEScaleFactor
      double precision PGAGetDEAuxFactor
      external PGAGetDEAuxFactor
      double precision PGAGetDECrossoverProb
      external PGAGetDECrossoverProb
      double precision PGAGetDEJitter
      external PGAGetDEJitter
      double precision PGAGetDEProbabilityEO
      external PGAGetDEProbabilityEO
      integer PGAGetDECrossoverType
      external PGAGetDECrossoverType
      double precision PGAGetDEDither
      external PGAGetDEDither
      integer PGAGetDEDitherPerIndividual
      external PGAGetDEDitherPerIndividual
c *** parallel
      integer(8) PGABuildDatatype
      external PGABuildDatatype
      integer PGAGetRank
      external PGAGetRank
      integer PGAGetNumProcs
      external PGAGetNumProcs
      integer(8) PGAGetCommunicator
      external PGAGetCommunicator
c *** pga
      integer PGAGetDataType
      external PGAGetDataType
      integer PGAGetOptDirFlag
      external PGAGetOptDirFlag
      integer PGAGetStringLength
      external PGAGetStringLength
      integer PGAGetGAIterValue
      external PGAGetGAIterValue
      integer PGAGetEvalCount
      external PGAGetEvalCount
      integer PGAGetMutationOrCrossoverFlag
      external PGAGetMutationOrCrossoverFlag
      integer PGAGetMutationAndCrossoverFlag
      external PGAGetMutationAndCrossoverFlag
      integer PGAGetMutationOnlyFlag
      external PGAGetMutationOnlyFlag
c *** pop
      integer PGAGetPopSize
      external PGAGetPopSize
      integer PGAGetNumReplaceValue
      external PGAGetNumReplaceValue
      integer PGAGetPopReplaceType
      external PGAGetPopReplaceType
      integer PGAGetRTRWindowSize
      external PGAGetRTRWindowSize
      integer PGAGetSortedPopIndex
      external PGAGetSortedPopIndex
c *** random
      integer PGARandomFlip
      external PGARandomFlip
      integer PGARandomInterval
      external PGARandomInterval
      double precision PGARandom01
      external PGARandom01
      double precision PGARandomUniform
      external PGARandomUniform
      double precision PGARandomGaussian
      external PGARandomGaussian
      integer PGAGetRandomSeed
      external PGAGetRandomSeed
      integer PGARandomNextSample
      external PGARandomNextSample
c *** real
      double precision PGAGetRealAllele
      external PGAGetRealAllele
      double precision PGAGetMinRealInitValue
      external PGAGetMinRealInitValue
      double precision PGAGetMaxRealInitValue
      external PGAGetMaxRealInitValue
      integer PGAGetRealInitType
      external PGAGetRealInitType
c *** report
      integer PGAGetPrintFrequencyValue
      external PGAGetPrintFrequencyValue
c *** restart
      integer PGAGetRestartFlag
      external PGAGetRestartFlag
      integer PGAGetRestartFrequencyValue
      external PGAGetRestartFrequencyValue
      double precision PGAGetRestartAlleleChangeProb
      external PGAGetRestartAlleleChangeProb
c *** select
      integer PGASelectNextIndex
      external PGASelectNextIndex
      integer PGAGetSelectType
      external PGAGetSelectType
      double precision PGAGetPTournamentProb
      external PGAGetPTournamentProb
      double precision PGAGetTournamentSize
      external PGAGetTournamentSize
      integer PGAGetTournamentWithReplacement
      external PGAGetTournamentWithReplacement
      double precision PGAGetTruncationProportion
      external PGAGetTruncationProportion
      integer PGAGetRandomizeSelect
      external PGAGetRandomizeSelect
      double precision PGAGetAuxTotal
      external PGAGetAuxTotal
c *** stop
      integer PGADone
      external PGADone
      integer PGACheckStoppingConditions
      external PGACheckStoppingConditions
      integer PGAGetStoppingRuleType
      external PGAGetStoppingRuleType
      integer PGAGetMaxGAIterValue
      external PGAGetMaxGAIterValue
c *** system
      integer PGAGetMaxMachineIntValue
      external PGAGetMaxMachineIntValue
      integer PGAGetMinMachineIntValue
      external PGAGetMinMachineIntValue
      double precision PGAGetMaxMachineDoubleValue
      external PGAGetMaxMachineDoubleValue
      double precision PGAGetMinMachineDoubleValue
      external PGAGetMinMachineDoubleValue
c *** utility
      double precision PGAMean
      external PGAMean
      double precision PGAStddev
      external PGAStddev
      integer PGARound
      external PGARound
      integer PGACheckSum
      external PGACheckSum
      integer PGAGetWorstIndex
      external PGAGetWorstIndex
      integer PGAGetBestIndex
      external PGAGetBestIndex
      double precision PGAGetBestReport
      external PGAGetBestReport
      integer PGAGetBestReportIndex
      external PGAGetBestReportIndex
      integer PGAEvalCompare
      external PGAEvalCompare
