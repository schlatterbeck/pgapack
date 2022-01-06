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
*     FILE: crossover.c: This file contains the data structure neutral
*                        crossover routines.
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
*****************************************************************************/

#include "pgapack.h"

/*U****************************************************************************
   PGACrossover - performs crossover on two parent strings to create two
   child strings (via side-effect).  The type of crossover performed is
   either the default or that specified by PGASetCrossoverType

   Category: Operators

   Inputs:
      ctx  - context variable
      p1   - the first parent string
      p2   - the second parent string
      pop1 - symbolic constant of the population containing string p1 and p2
      c1   - the first child string
      c2   - the second child string
      pop2 - symbolic constant of the population to contain string c1 and c2

   Outputs:
      c1 and c2 in pop2 are children of p1 and p2 in pop1.  p1 and p2 are not
      modified.

   Example:
      Perform crossover on the two parent strings mom and dad in population
      PGA_OLDPOP, and insert the child strings, child1 and child1, in
      population PGA_NEWPOP.

      PGAContext *ctx;
      int mom, dad, child1, child2;
      :
      PGACrossover(ctx, mom, dad, PGA_OLDPOP, child1, child2, PGA_NEWPOP);

****************************************************************************U*/
void PGACrossover ( PGAContext *ctx, int p1, int p2, int pop1,
                    int c1, int c2, int pop2 )
{
    int fp1, fp2, fc1, fc2;

    PGADebugEntered("PGACrossover");
    PGADebugPrint( ctx, PGA_DEBUG_PRINTVAR, "PGACrossover", " p1 = ",
		  PGA_INT, (void *) &p1 );
    PGADebugPrint( ctx, PGA_DEBUG_PRINTVAR, "PGACrossover", " p2 = ",
		  PGA_INT, (void *) &p2 );
    PGADebugPrint( ctx, PGA_DEBUG_PRINTVAR, "PGACrossover", " c1 = ",
		  PGA_INT, (void *) &c1 );
    PGADebugPrint( ctx, PGA_DEBUG_PRINTVAR, "PGACrossover", " c2 = ",
		  PGA_INT, (void *) &c2 );

    if (ctx->fops.Crossover) {
	fp1 = ((p1 == PGA_TEMP1) || (p1 == PGA_TEMP2)) ? p1 : p1+1;
	fp2 = ((p2 == PGA_TEMP1) || (p2 == PGA_TEMP2)) ? p2 : p2+1;
	fc1 = ((c1 == PGA_TEMP1) || (c1 == PGA_TEMP2)) ? c1 : c1+1;
	fc2 = ((c2 == PGA_TEMP1) || (c2 == PGA_TEMP2)) ? c2 : c2+1;
	(*ctx->fops.Crossover)(&ctx, &fp1, &fp2, &pop1, &fc1, &fc2, &pop2);
    } else {
	(*ctx->cops.Crossover)(ctx, p1, p2, pop1, c1, c2, pop2);
    }

    PGASetEvaluationUpToDateFlag(ctx, c1, pop2, PGA_FALSE);
    PGASetEvaluationUpToDateFlag(ctx, c2, pop2, PGA_FALSE);
    
    PGADebugExited("PGACrossover");
}

/*U***************************************************************************
   PGAGetCrossoverType - Returns the type of crossover selected

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      Returns the integer corresponding to the symbolic constant
      used to specify the crossover type

   Example:
      PGAContext *ctx;
      int crosstype;
      :
      crosstype = PGAGetCrossoverType(ctx);
      switch (crosstype) {
      case PGA_CROSSOVER_ONEPT:
          printf("Crossover Type = PGA_CROSSOVER_ONEPT\n");
          break;
      case PGA_CROSSOVER_TWOPT:
          printf("Crossover Type = PGA_CROSSOVER_TWOPT\n");
          break;
      case PGA_CROSSOVER_UNIFORM:
          printf("Crossover Type = PGA_CROSSOVER_UNIFORM\n");
          break;
      case PGA_CROSSOVER_SBX:
          printf("Crossover Type = PGA_CROSSOVER_SBX\n");
          break;
      }

***************************************************************************U*/
int PGAGetCrossoverType (PGAContext *ctx)
{
    PGADebugEntered("PGAGetCrossoverType");

    PGAFailIfNotSetUp("PGAGetCrossoverType");

    PGADebugExited("PGAGetCrossoverType");

    return(ctx->ga.CrossoverType);
}

/*U***************************************************************************
   PGAGetCrossoverProb - Returns the crossover probability

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The crossover probability

   Example:
      PGAContext *ctx;
      double pc;
      :
      pc = PGAGetCrossoverProb(ctx);

***************************************************************************U*/
double PGAGetCrossoverProb (PGAContext *ctx)
{
    PGADebugEntered("PGAGetCrossoverProb");

    PGAFailIfNotSetUp("PGAGetCrossoverProb");

    PGADebugExited("PGAGetCrossoverProb");

    return(ctx->ga.CrossoverProb);
}

/*U***************************************************************************
   PGAGetUniformCrossoverProb - returns the probability of a bit being
   selected from the first child string in uniform crossover

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The uniform crossover probability

   Example:
      PGAContext *ctx;
      double pu;
      :
      pu = PGAGetUniformCrossoverProb(ctx);

***************************************************************************U*/
double PGAGetUniformCrossoverProb (PGAContext *ctx)
{
    PGADebugEntered("PGAGetUniformCrossoverProb");

    PGAFailIfNotSetUp("PGAGetUniformCrossoverProb");

    PGADebugExited("PGAGetUniformCrossoverProb");

    return(ctx->ga.UniformCrossProb);
}

/*U****************************************************************************
   PGASetCrossoverType - specify the type of crossover to use.
   The default is PGA_CROSSOVER_TWOPT.

   Category: Operators

   Inputs:
      ctx            - context variable
      crossover_type - symbolic constant to specify crossover type

   Outputs:
      None

   Example:
      Use uniform crossover when crossingover strings.

      PGAContext *ctx;
      :
      PGASetCrossoverType(ctx, PGA_CROSSOVER_UNIFORM);

****************************************************************************U*/
void PGASetCrossoverType (PGAContext *ctx, int crossover_type)
{

    PGADebugEntered("PGASetCrossoverType");

    switch (crossover_type) {
        case PGA_CROSSOVER_ONEPT:
        case PGA_CROSSOVER_TWOPT:
        case PGA_CROSSOVER_UNIFORM:
        case PGA_CROSSOVER_SBX:
            ctx->ga.CrossoverType = crossover_type;
            break;
        default:
            PGAError( ctx,
                     "PGASetCrossoverType: Invalid value of crossover_type:",
                      PGA_FATAL, PGA_INT, (void *) &crossover_type );
    };

    PGADebugExited("PGASetCrossoverType");
}


/*U****************************************************************************
   PGASetCrossoverProb - Probability that a selected string will undergo
   crossover.  The default is 0.85.

   Category: Operators

   Inputs:
      ctx - context variable
      p   - the crossover probability

   Outputs:
      None

   Example:
      Make crossover happen infrequently.

      PGAContext *ctx;
      :
      PGASetCrossoverProb(ctx,0.001);

****************************************************************************U*/
void PGASetCrossoverProb( PGAContext *ctx, double crossover_prob)
{
    PGADebugEntered("PGASetCrossoverProb");

    if ((crossover_prob < 0.0) || (crossover_prob > 1.0))
        PGAError ( ctx,
                  "PGASetCrossoverProb: Invalid value of crossover_prob:",
                   PGA_FATAL, PGA_DOUBLE, (void *) &crossover_prob);
    else
        ctx->ga.CrossoverProb = crossover_prob;

    PGADebugExited("PGASetCrossoverProb");
}

/*U****************************************************************************
   PGASetUniformCrossoverProb - Probability used in uniform crossover
   to specify that an allele value value be selected from a particular
   parent. The default is 0.6.  The crossover type must have been set
   to PGA_CROSSOVER_UNIFORM with PGASetCrossoverType for this function
   call to have any effect.

   Category: Operators

   Inputs:
      ctx - context variable
      p   - the crossover probability

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetUniformCrossoverProb(ctx,0.9);

****************************************************************************U*/
void PGASetUniformCrossoverProb( PGAContext *ctx, double uniform_cross_prob)
{
    PGADebugEntered("PGASetUniformCrossoverProb");

    if ((uniform_cross_prob < 0.0) || (uniform_cross_prob > 1.0))
        PGAError ( ctx,
                  "PGASetUniformCrossoverProb: Invalid value of "
                  "uniform_cross_prob:", PGA_FATAL, PGA_DOUBLE,
                  (void *) &uniform_cross_prob);
    else
        ctx->ga.UniformCrossProb = uniform_cross_prob;

    PGADebugExited("PGASetUniformCrossoverProb");
}

/*U****************************************************************************
   PGASetCrossoverBoundedFlag - If this flag is set to PGA_TRUE, then for
   Integer and Real strings with simulated binary crossover (SBX)
   crossed over values that exceed the bounds are confined to the
   bounds by setting them to the boundary.

   Category: Operators

   Inputs:
      ctx  - context variable
      flag - either PGA_TRUE or PGA_FALSE

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetCrossoverBoundedFlag (ctx, PGA_TRUE);

****************************************************************************U*/
void PGASetCrossoverBoundedFlag (PGAContext *ctx, int val)
{
    switch (val)
    {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.CrossBoundedFlag = val;
         break;
    default:
         PGAError(ctx, "PGASetCrossoverBoundedFlag: Invalid value:",
                  PGA_FATAL, PGA_INT, (void *) &val);
         break;
    }
}

/*U****************************************************************************
   PGAGetCrossoverBoundedFlag - returns PGA_TRUE or PGA_FALSE to indicate
   whether crossed over strings remain in the range specified.

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      PGA_TRUE if restricted to the given range, otherwise PGA_FALSE.

   Example:
      PGAContext *ctx;
      int val;
      :
      val = PGAGetCrossoverBoundedFlag (ctx);

****************************************************************************U*/
int PGAGetCrossoverBoundedFlag (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetCrossoverBoundedFlag");
    return (ctx->ga.CrossBoundedFlag);
}

/*U****************************************************************************
   PGASetCrossoverBounceBackFlag - If this flag is set to PGA_TRUE, then
   for Integer and Real strings with simulated binary crossover (SBX)
   crossed over values that exceed the bounds are confined to the
   bounds by bouncing them back to a random value between the boundary
   and the neares parent.

   Category: Operators

   Inputs:
      ctx  - context variable
      flag - either PGA_TRUE or PGA_FALSE

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetCrossoverBounceBackFlag (ctx, PGA_TRUE);

****************************************************************************U*/
void PGASetCrossoverBounceBackFlag (PGAContext *ctx, int val)
{
    switch (val)
    {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.CrossBounceFlag = val;
         break;
    default:
         PGAError(ctx, "PGASetCrossoverBounceBackFlag: Invalid value:",
                  PGA_FATAL, PGA_INT, (void *) &val);
         break;
    }
}

/*U****************************************************************************
   PGAGetCrossoverBounceBackFlag - returns PGA_TRUE or PGA_FALSE to indicate
   whether crossed over strings are bounced back when exceeding the
   bounds.

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      PGA_TRUE if restricted to the given range, otherwise PGA_FALSE.

   Example:
      PGAContext *ctx;
      int val;
      :
      val = PGAGetCrossoverBounceBackFlag (ctx);

****************************************************************************U*/
int PGAGetCrossoverBounceBackFlag (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetCrossoverBounceBackFlag");
    return (ctx->ga.CrossBounceFlag);
}

/*U****************************************************************************
   PGASetCrossoverSBXEta - Set the eta parameter for simulated binary
   crossover (SBX)

   Category: Operators

   Inputs:
      ctx  - context variable
      eta  - eta >= 0

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetCrossoverSBXEta (ctx, 10);

****************************************************************************U*/
void PGASetCrossoverSBXEta (PGAContext *ctx, double eta)
{
    if (eta < 0) {
        PGAError(ctx, "PGASetCrossoverSBXEta: Invalid value:",
                 PGA_FATAL, PGA_DOUBLE, (void *) &eta);
    }
    ctx->ga.CrossSBXEta = eta;
}

/*U****************************************************************************
   PGAGetCrossoverSBXEta - returns SBX eta value

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The SBX eta value

   Example:
      PGAContext *ctx;
      double eta;
      :
      eta = PGAGetCrossoverSBXEta (ctx);

****************************************************************************U*/
double PGAGetCrossoverSBXEta (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetCrossoverSBXEta");
    return (ctx->ga.CrossSBXEta);
}

/*U****************************************************************************
   PGASetCrossoverSBXOncePerString - Compute random number for SBX
   polynomial distribution only once per string/individual. If set to
   PGA_TRUE all alleles will use the same value which means that the
   resulting string will point into the same direction as the vector
   between both parents.

   Category: Operators

   Inputs:
      ctx  - context variable
      val  - PGA_TRUE or PGA_FALSE, default is PGA_FALSE

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetCrossoverSBXOncePerString (ctx, PGA_TRUE);

****************************************************************************U*/
void PGASetCrossoverSBXOncePerString (PGAContext *ctx, int val)
{
    if (val != PGA_TRUE && val != PGA_FALSE) {
        PGAError(ctx, "PGASetCrossoverSBXOncePerString: Invalid value:",
                 PGA_FATAL, PGA_INT, (void *) &val);
    }
    ctx->ga.CrossSBXOnce = val;
}

/*U****************************************************************************
   PGAGetCrossoverSBXOncePerString - returns SBX setting if random
   number for SBX polynomial distribution is computed once per string.

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      The SBX once-per-string value

   Example:
      PGAContext *ctx;
      int r;
      :
      r = PGAGetCrossoverSBXOncePerString (ctx);

****************************************************************************U*/
int PGAGetCrossoverSBXOncePerString (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetCrossoverSBXOncePerString");
    return (ctx->ga.CrossSBXOnce);
}


/*U****************************************************************************
   PGACrossoverSBX - Cross over two parent alleles with simulated binary
   crossover. This uses double for both parent alleles but is re-used in
   both, integer and real SBX crossover. The probability is used to
   compute the new alleles from the polynomial distribution.

   Category: Operators

   Inputs:
      ctx - context variable
      p1  - (double) Allele of first string
      p2  - (double) Allele of second string
      u   - Random value between 0 and 1
      c1  - pointer to new first child allele
      c2  - pointer to new second child allele


   Outputs:
      None

   Example:
      PGAContext *ctx;
      double p1, p2, u;
      double c1, c2;
      double result;
      :
      u = PGARandom01 (ctx, 0);
      result = PGACrossoverSBX (ctx, p1, p2, u, &c1, &c2);

****************************************************************************U*/
void PGACrossoverSBX
    (PGAContext *ctx, double p1, double p2, double u, double *c1, double *c2)
{
    double beta;
    double eta = PGAGetCrossoverSBXEta (ctx) + 1.0;
    if (u < 0.5) {
        beta = pow (2.0 * u, 1.0 / eta);
    } else {
        beta = pow (2.0 - 2.0 * u, -1.0 / eta);
    }
    *c1 = 0.5 * ((p1 + p2) - beta * fabs (p1 - p2));
    *c2 = 0.5 * ((p1 + p2) + beta * fabs (p1 - p2));
}
