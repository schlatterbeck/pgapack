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

/*****************************************************************************
*     File: mutation.c: This file contains the data structure neutral mutation
*                       routines
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
*****************************************************************************/

#include "pgapack.h"

/*U****************************************************************************
  PGAMutate - This routine performs mutation on a string.  The type of mutation
  depends on the data type.  Refer to the user guide for data-specific
  examples.

  Category: Operators

  Inputs:
      ctx  - context variable
      p   - index of string to mutate
      pop - symbolic constant of the population containing p

  Output:
      The number of mutations performed.  Member p in population pop is
      mutated by side-effect.

  Example:
      Mutate the best string in the population, until 10 or more mutations
      have occured.

      PGAContext *ctx;
      int p, count = 0;
      :
      p = PGAGetBestIndex(ctx, PGA_NEWPOP);
      while (count < 10) {
          count += PGAMutate(ctx, p, PGA_NEWPOP);
      }

****************************************************************************U*/
int PGAMutate(PGAContext *ctx, int p, int pop)
{
    double mr;
    int count;
    int fp;
    PGADebugEntered("PGAMutate");
    
    mr    = ctx->ga.MutationProb;
    if (ctx->fops.Mutation) {
	fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
        count = (*ctx->fops.Mutation)(&ctx, &fp, &pop, &mr);
    } else {
	count = (*ctx->cops.Mutation)( ctx, p, pop, mr );
    }
    
    if ( count > 0 )
	PGASetEvaluationUpToDateFlag(ctx, p, pop, PGA_FALSE);
    
    PGADebugExited("PGAMutate");
    
    return(count);
}

/*U****************************************************************************
   PGASetMutationType - set type of mutation to use. Only effects integer-
   and real-valued strings.  Binary-valued strings are always complemented.
   In character-valued strings, one alphabetic character is replaced with
   another chosen uniformly randomly.  The alphabetic characters will be lower,
   upper, or mixed case depending on how the strings were initialized.

   Valid choices are PGA_MUTATION_CONSTANT (Real/Integer), PGA_MUTATION_RANGE
   (Real/Integer), PGA_MUTATION_UNIFORM (Real), PGA_MUTATION_GAUSSIAN (Real),
   and PGA_MUTATION_PERMUTE (Integer).  The default for integer-valued strings
   conforms to how the strings were initialized.  The default for real-valued
   strings is PGA_MUTATION_GAUSSIAN.  See the user guide for more details.

   Category: Operators

   Inputs:
      ctx           - context variable
      mutation_type - symbolic constant to specify the mutation type

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetMutationType(ctx, PGA_MUTATION_UNIFORM);

****************************************************************************U*/
void PGASetMutationType( PGAContext *ctx, int mutation_type)
{
    PGADebugEntered("PGASetMutationType");

     switch (mutation_type)
     {
     case PGA_MUTATION_CONSTANT:
     case PGA_MUTATION_RANGE:
     case PGA_MUTATION_UNIFORM:
     case PGA_MUTATION_GAUSSIAN:
     case PGA_MUTATION_PERMUTE:
     case PGA_MUTATION_DE:
          ctx->ga.MutationType = mutation_type;
          break;
     default:
          PGAError ( ctx,
                    "PGASetMutationType: Invalid value of mutation_type:",
                    PGA_FATAL, PGA_INT, (void *) &mutation_type);
          break;
     }

    PGADebugExited("PGASetMutationType");
}

/*U***************************************************************************
   PGAGetMutationType - Returns the type of mutation used

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      Returns the integer corresponding to the symbolic constant
      used to specify the type of mutation specified

   Example:
      PGAContext *ctx;
      int mutatetype;
      :
      mutatetype = PGAGetMutationType(ctx);
      switch (mutatetype) {
      case PGA_MUTATION_CONSTANT:
          printf("Mutation Type = PGA_MUTATION_CONSTANT\n");
          break;
      case PGA_MUTATION_RANGE:
          printf("Mutation Type = PGA_MUTATION_RANGE\n");
          break;
      case PGA_MUTATION_UNIFORM:
          printf("Mutation Type = PGA_MUTATION_UNIFORM\n");
          break;
      case PGA_MUTATION_GAUSSIAN:
          printf("Mutation Type = PGA_MUTATION_GAUSSIAN\n");
          break;
      case PGA_MUTATION_PERMUTE:
          printf("Mutation Type = PGA_MUTATION_PERMUTE\n");
          break;
      }

***************************************************************************U*/
int PGAGetMutationType (PGAContext *ctx)
{
    PGADebugEntered("PGAGetMutationType");
    PGAFailIfNotSetUp("PGAGetMutationType");
    PGADebugExited("PGAGetMutationType");
    return(ctx->ga.MutationType);
}

/*U****************************************************************************
   PGASetMutationRealValue - Set multiplier to mutate PGA_DATATYPE_REAL
   strings with.  The use of this value depends on the type of mutation
   being used.  The default value is 0.1.  See the user guide for more details.

   Category: Operators

   Inputs:
      ctx - context variable
      val - the mutation value to use for Real mutation

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetMutationRealValue(ctx,50.0);

****************************************************************************U*/
void PGASetMutationRealValue( PGAContext *ctx, double val)
{
    PGADebugEntered("PGASetMutationRealValue");

    if (val < 0.0)
        PGAError ( ctx,
                  "PGASetMutationRealValue: Invalid value of val:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &val);
    else
        ctx->ga.MutateRealValue = val;

    PGADebugExited("PGASetMutationRealValue");
}

/*U***************************************************************************
   PGAGetMutationRealValue - Returns the value of the multiplier used to
   mutate PGA_DATATYPE_REAL strings with.

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the multiplier used to mutate PGA_DATATYPE_REAL strings with

   Example:
      PGAContext *ctx;
      double val;
      :
      val = PGAGetMutationRealValue(ctx);

***************************************************************************U*/
double PGAGetMutationRealValue (PGAContext *ctx)
{
    PGADebugEntered("PGAGetMutationRealValue");
    PGAFailIfNotSetUp("PGAGetMutationRealValue");

    PGADebugExited("PGAGetMutationRealValue");

    return(ctx->ga.MutateRealValue);
}

/*U****************************************************************************
   PGASetMutationIntegerValue - Set multiplier to mutate PGA_DATATYPE_INTEGER
   strings with.  The use of this value depends on the type of mutation
   being used.  The default value is 1.  See the user guide for more details.

   Category: Operators

   Inputs:
      ctx - context variable
      val - the mutation value to use for Integer mutation

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetMutationIntegerValue(ctx, 5);

****************************************************************************U*/
void PGASetMutationIntegerValue( PGAContext *ctx, int val)
{
    PGADebugEntered("PGASetMutationIntegerValue");

    if (val < 0.0)
        PGAError ( ctx,
                  "PGASetMutationIntegerValue: Invalid value of val:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &val);
    else
        ctx->ga.MutateIntegerValue = val;

    PGADebugExited("PGASetMutationIntegerValue");
}


/*U***************************************************************************
  PGAGetMutationIntegerValue - Returns the value of the multiplier
  used to mutate PGA_DATATYPE_INTEGER strings with.

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the multiplier used to mutate PGA_DATATYPE_INTEGER
      strings with

   Example:
      PGAContext *ctx;
      int ival;
      :
      ival = PGAGetMutationIntegerValue(ctx);

***************************************************************************U*/
int PGAGetMutationIntegerValue (PGAContext *ctx)
{
    PGADebugEntered("PGAGetMutationIntegerValue");
    PGAFailIfNotSetUp("PGAGetMutationIntegerValue");

    PGADebugExited("PGAGetMutationIntegerValue");

    return(ctx->ga.MutateIntegerValue);
}

/*U****************************************************************************
   PGASetMutationBoundedFlag - If this flag is set to PGA_TRUE, then for
   Integer and Real strings whenever a gene is mutated, if it underflows
   (overflows) the lower (upper)bound it is reset to the lower (upper) bound.
   In this way all allele values remain within the range the integer strings
   were initialized on.  If this flag is PGA_FALSE (the default), the alleles
   may take any values.

   Category: Operators

   Inputs:
      ctx  - context variable
      flag - either PGA_TRUE or PGA_FALSE

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetMutationBoundedFlag(ctx, PGA_TRUE);

****************************************************************************U*/
void PGASetMutationBoundedFlag(PGAContext *ctx, int val)
{
    PGADebugEntered("PGASetMutationBoundedFlag");

    switch (val)
    {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.MutateBoundedFlag = val;
         break;
    default:
         PGAError(ctx, "PGASetMutationBoundedFlag: Invalid value:",
                  PGA_FATAL, PGA_INT, (void *) &val);
         break;
    }

    PGADebugExited("PGASetMutationBoundedFlag");
}


/*U****************************************************************************
   PGAGetMutationBoundedFlag - returns PGA_TRUE or PGA_FALSE to indicate
   whether mutated integer strings remain in the range specified when
   initialized with PGASetIntegerInitRange.

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      PGA_TRUE if restricted to the given range, otherwise PGA_FALSE.

   Example:
      PGAContext *ctx;
      int val;
      :
      val = PGAGetMutationBoundedFlag(ctx);

****************************************************************************U*/
int PGAGetMutationBoundedFlag(PGAContext *ctx)
{
    PGADebugEntered  ("PGAGetMutationBoundedFlag");
    PGAFailIfNotSetUp("PGAGetMutationBoundedFlag");
    PGADebugExited   ("PGAGetMutationBoundedFlag");
    return (ctx->ga.MutateBoundedFlag);
}

/*U****************************************************************************
   PGASetMutationBounceBackFlag - If this flag is set to PGA_TRUE, then for
   Integer and Real strings whenever a gene is mutated, if it underflows
   (overflows) the lower (upper)bound it is reset to a random value between
   the old value and the violated bound.
   In this way all allele values remain within the range the strings
   were initialized on.  If this flag is PGA_FALSE (the default), the alleles
   may take any values. See also PGASetMutationBoundedFlag.

   Category: Operators

   Inputs:
      ctx  - context variable
      flag - either PGA_TRUE or PGA_FALSE

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetMutationBounceBackFlag(ctx, PGA_TRUE);

****************************************************************************U*/
void PGASetMutationBounceBackFlag(PGAContext *ctx, int val)
{
    PGADebugEntered("PGASetMutationBounceBackFlag");

    switch (val)
    {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.MutateBounceFlag = val;
         break;
    default:
         PGAError(ctx, "PGASetMutationBounceBackFlag: Invalid value:",
                  PGA_FATAL, PGA_INT, (void *) &val);
         break;
    }

    PGADebugExited("PGASetMutationBounceBackFlag");
}


/*U****************************************************************************
   PGAGetMutationBounceBackFlag - returns PGA_TRUE or PGA_FALSE to indicate
   whether mutated strings remain within the initialization range

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      PGA_TRUE if restricted to the given range, otherwise PGA_FALSE.

   Example:
      PGAContext *ctx;
      int val;
      :
      val = PGAGetMutationBounceBackFlag(ctx);

****************************************************************************U*/
int PGAGetMutationBounceBackFlag(PGAContext *ctx)
{
    PGADebugEntered  ("PGAGetMutationBounceBackFlag");
    PGAFailIfNotSetUp("PGAGetMutationBounceBackFlag");
    PGADebugExited   ("PGAGetMutationBounceBackFlag");
    return (ctx->ga.MutateBounceFlag);
}



/*U****************************************************************************
   PGASetMutationProb - Specifies the probability that a given allele will
   be mutated.  If this is called without calling PGASetMutationType(), the
   default mutation type is PGA_MUTATION_FIXED.  The default probability is
   the reciprocal of the string length.

   Category: Operators

   Inputs:
      ctx - context variable
      p   - the mutation probability

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetMutationProb(ctx,0.001);

****************************************************************************U*/
void PGASetMutationProb(PGAContext *ctx, double mutation_prob)
{
    PGADebugEntered("PGASetMutationProb");

    if ((mutation_prob < 0.0) || (mutation_prob > 1.0))
        PGAError ( ctx,
                  "PGASetMutationProb: Invalid value of mutation_prob:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &mutation_prob);
    else
        ctx->ga.MutationProb = mutation_prob;

    PGADebugExited("PGASetMutationProb");
}

/*U***************************************************************************
   PGAGetMutationProb - Returns the probability of mutation.

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The mutation probability

   Example:
      PGAContext *ctx;
      double pm;
      :
      pm = PGAGetMutateProb(ctx);

***************************************************************************U*/
double PGAGetMutationProb (PGAContext *ctx)
{
    PGADebugEntered("PGAGetMutationProb");
    PGAFailIfNotSetUp("PGAGetMutationProb");
    PGADebugExited("PGAGetMutationProb");
    return(ctx->ga.MutationProb);
}


/*U****************************************************************************
  PGASetDEVariant - sets the variant used for Differential Evolution.
  Only used if the mutation type is PGA_MUTATION_DE.

  Category: Initialization

  Inputs:
     ctx - context variable
     v   - symbolic constant, currently one of PGA_DE_VARIANT_RAND,
           PGA_DE_VARIANT_BEST, PGA_DE_VARIANT_EITHER_OR

  Outputs:

  Example:

     PGAContext *ctx;
     :
     PGASetDEVariant(ctx, PGA_DE_VARIANT_BEST);

****************************************************************************U*/
void PGASetDEVariant (PGAContext *ctx, int variant)
{
    PGAFailIfSetUp("PGASetDEVariant");
    switch (variant)
    {
    case PGA_DE_VARIANT_BEST:
    case PGA_DE_VARIANT_RAND:
    case PGA_DE_VARIANT_EITHER_OR:
        ctx->ga.DEVariant = variant;
        break;
    default:
        PGAError
            ( ctx, "PGASetDEVariant: Invalid value of DE variant:"
            , PGA_FATAL, PGA_INT, (void *) &variant
            );
        break;
    }
}

/*U***************************************************************************
   PGAGetDEVariant - Returns the variant of Differential Evolution

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      Returns the integer corresponding to the symbolic constant
      used to specify the variant of differential evolution

   Example:
      PGAContext *ctx;
      int variant;
      :
      variant = PGAGetDEVariant(ctx);
      switch (variant) {
      case PGA_DE_VARIANT_RAND:
          printf("DE Variant = PGA_DE_VARIANT_RAND\n");
          break;
      case PGA_DE_VARIANT_BEST:
          printf("DE Variant = PGA_DE_VARIANT_BEST\n");
          break;
      }

***************************************************************************U*/
int PGAGetDEVariant (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDEVariant");
    return(ctx->ga.DEVariant);
}

/*U****************************************************************************
   PGASetDEScaleFactor - Set the scale factor F for DE

   Category: Operators

   Inputs:
      ctx - context variable
      val - the scale factor

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDEScaleFactor (ctx, 0.75);

****************************************************************************U*/
void PGASetDEScaleFactor (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 2.0)
        PGAError ( ctx,
                  "PGASetDEScaleFactor: Invalid value of F:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &val);
    else
        ctx->ga.DEScaleFactor = val;
}

/*U***************************************************************************
   PGAGetDEScaleFactor - Returns the value of the scale factor F for DE

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the scale factor

   Example:
      PGAContext *ctx;
      double val;
      :
      val = PGAGetDEScaleFactor(ctx);

***************************************************************************U*/
double PGAGetDEScaleFactor (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDEScaleFactor");
    return(ctx->ga.DEScaleFactor);
}

/*U****************************************************************************
   PGASetDEAuxFactor - Set the auxiliary factor K for DE

   Category: Operators

   Inputs:
      ctx - context variable
      val - the auxiliary factor

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDEAuxFactor (ctx, 0.75);

****************************************************************************U*/
void PGASetDEAuxFactor (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 2.0)
        PGAError ( ctx,
                  "PGASetDEAuxFactor: Invalid value of K:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &val);
    else
        ctx->ga.DEAuxFactor = val;
}

/*U***************************************************************************
   PGAGetDEAuxFactor - Returns the value of the auxiliary factor K for DE

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the auxiliary factor

   Example:
      PGAContext *ctx;
      double val;
      :
      val = PGAGetDEAuxFactor(ctx);

***************************************************************************U*/
double PGAGetDEAuxFactor (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDEAuxFactor");
    return(ctx->ga.DEAuxFactor);
}

/*U****************************************************************************
   PGASetDECrossoverProb - Set the crossover probability for DE

   Category: Operators

   Inputs:
      ctx - context variable
      val - the crossover probability

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDECrossoverProb (ctx, 0.75);

****************************************************************************U*/
void PGASetDECrossoverProb (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 1.0)
        PGAError ( ctx,
                  "PGASetDECrossoverProb: Invalid value of crossover:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &val);
    else
        ctx->ga.DECrossoverProb = val;
}

/*U***************************************************************************
   PGAGetDECrossoverProb - Returns the crossover probability for DE

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the DE crossover probability

   Example:
      PGAContext *ctx;
      double val;
      :
      val = PGAGetDECrossoverProb(ctx);

***************************************************************************U*/
double PGAGetDECrossoverProb (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDECrossoverProb");
    return(ctx->ga.DECrossoverProb);
}

/*U****************************************************************************
   PGASetDEJitter - Set the jitter for DE

   Category: Operators

   Inputs:
      ctx - context variable
      val - the jitter for DE

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDEJitter (ctx, 0.75);

****************************************************************************U*/
void PGASetDEJitter (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 2.0)
        PGAError ( ctx,
                  "PGASetDEJitter: Invalid value of jitter:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &val);
    else
        ctx->ga.DEJitter = val;
}

/*U***************************************************************************
   PGAGetDEJitter - Returns the jitter for DE

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the DE jitter

   Example:
      PGAContext *ctx;
      double val;
      :
      val = PGAGetDEJitter(ctx);

***************************************************************************U*/
double PGAGetDEJitter (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDEJitter");
    return(ctx->ga.DEJitter);
}

/*U****************************************************************************
   PGASetDEProbabilityEO - Set the either-or probability for
   PGA_DE_VARIANT_EITHER_OR of Differential Evolution

   Category: Operators

   Inputs:
      ctx - context variable
      val - the either-or probability

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDEProbabilityEO (ctx, 0.75);

****************************************************************************U*/
void PGASetDEProbabilityEO (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 1.0)
        PGAError ( ctx,
                  "PGASetDEProbabilityEO: Invalid value of EO probabilty:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &val);
    else
        ctx->ga.DEProbabilityEO = val;
}

/*U***************************************************************************
   PGAGetDEProbabilityEO - Returns the probability of the either-or
   variant of Differential Evolution

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the DE EO-Probability

   Example:
      PGAContext *ctx;
      double val;
      :
      val = PGAGetDEProbabilityEO(ctx);

***************************************************************************U*/
double PGAGetDEProbabilityEO (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDEProbabilityEO");
    return(ctx->ga.DEProbabilityEO);
}

/*U****************************************************************************
   PGASetDENumDiffs - Set the number of differences for DE

   Category: Operators

   Inputs:
      ctx - context variable
      val - the number of differences

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDENumDiffs (ctx, 2);

****************************************************************************U*/
void PGASetDENumDiffs (PGAContext *ctx, int val)
{
    if (val < 1 || val > 2)
        PGAError ( ctx,
                  "PGASetDENumDiffs: Invalid value of num diffs:",
                   PGA_FATAL, PGA_INT, (void *) &val);
    else
        ctx->ga.DENumDiffs = val;
}

/*U***************************************************************************
   PGAGetDENumDiffs - Returns the number of differences for DE

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the DE number of differences

   Example:
      PGAContext *ctx;
      int val;
      :
      val = PGAGetDENumDiffs(ctx);

***************************************************************************U*/
int PGAGetDENumDiffs (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDENumDiffs");
    return(ctx->ga.DENumDiffs);
}

/*U****************************************************************************
   PGASetDECrossoverType - Set the crossover type for DE

   Category: Operators

   Inputs:
      ctx - context variable
      val - the crossover type

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDECrossoverType (ctx, PGA_DE_CROSSOVER_EXP);

****************************************************************************U*/
void PGASetDECrossoverType (PGAContext *ctx, int val)
{
    switch (val)
    {
    case PGA_DE_CROSSOVER_BIN:
    case PGA_DE_CROSSOVER_EXP:
        ctx->ga.DECrossoverType = val;
        break;
    default:
        PGAError
            ( ctx, "PGASetDECrossoverType: Invalid crossover type:"
            , PGA_FATAL, PGA_INT, (void *) &val
            );
        break;
    }
}

/*U***************************************************************************
   PGAGetDECrossoverType - Returns the DE crossover type

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the DE crossover type

   Example:
      PGAContext *ctx;
      int val;
      :
      val = PGAGetDECrossoverType(ctx);

***************************************************************************U*/
int PGAGetDECrossoverType (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDECrossoverType");
    return(ctx->ga.DECrossoverType);
}

/*U****************************************************************************
   PGASetDEDither - Set the DE dither range (+/-)

   Category: Operators

   Inputs:
      ctx - context variable
      val - the dither range

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDEDither (ctx, 0.25);

****************************************************************************U*/
void PGASetDEDither (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 1.0)
        PGAError ( ctx,
                  "PGASetDEProbabilityEO: Invalid value of Dither:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &val);
    else
        ctx->ga.DEDither = val;
}

/*U***************************************************************************
   PGAGetDEDither - Returns the DE dither value

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The value of the DE dither

   Example:
      PGAContext *ctx;
      double val;
      :
      val = PGAGetDEDither(ctx);

***************************************************************************U*/
double PGAGetDEDither (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetDEDither");
    return(ctx->ga.DEDither);
}

/*U****************************************************************************
   PGASetDEDitherPerIndividual - If this is set to PGA_TRUE, then for
   Differential Evolution if the Dither value is non-zero we produce a
   new random value to add to the scale factor F *for each individual*.
   Otherwise if the flag is not set (PGA_FALSE), the we produce a new 
   value in each generation, the same value for *all* individuals.

   Category: Operators

   Inputs:
      ctx  - context variable
      flag - either PGA_TRUE or PGA_FALSE

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetDEDitherPerIndividual(ctx, PGA_TRUE);

****************************************************************************U*/
void PGASetDEDitherPerIndividual(PGAContext *ctx, int val)
{
    switch (val)
    {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.DEDitherPerIndividual = val;
         break;
    default:
         PGAError(ctx, "PGASetDEDitherPerIndividual: Invalid value:",
                  PGA_FATAL, PGA_INT, (void *) &val);
         break;
    }
}


/*U****************************************************************************
   PGAGetDEDitherPerIndividual - returns PGA_TRUE or PGA_FALSE to indicate
   whether the dither is applied anew for each individual or if the
   value is re-used for all individuals in one generation.

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      PGA_TRUE if dither is applied for each individual.

   Example:
      PGAContext *ctx;
      int val;
      :
      val = PGAGetDEDitherPerIndividual(ctx);

****************************************************************************U*/
int PGAGetDEDitherPerIndividual(PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetMutationBounceBackFlag");
    return (ctx->ga.DEDitherPerIndividual);
}

