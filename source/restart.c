/*
 *  
 *  ********************************************************************* 
 *  (C) COPYRIGHT 1995 UNIVERSITY OF CHICAGO 
 *  *********************************************************************
 *  
 *  This software was authored by
 *  
 *  D. Levine
 *  Mathematics and Computer Science Division Argonne National Laboratory
 *  Argonne IL 60439
 *  levine@mcs.anl.gov
 *  (708) 252-6735
 *  (708) 252-5986 (FAX)
 *  
 *  with programming assistance of participants in Argonne National 
 *  Laboratory's SERS program.
 *  
 *  This program contains material protectable under copyright laws of the 
 *  United States.  Permission is hereby granted to use it, reproduce it, 
 *  to translate it into another language, and to redistribute it to 
 *  others at no charge except a fee for transferring a copy, provided 
 *  that you conspicuously and appropriately publish on each copy the 
 *  University of Chicago's copyright notice, and the disclaimer of 
 *  warranty and Government license included below.  Further, permission 
 *  is hereby granted, subject to the same provisions, to modify a copy or 
 *  copies or any portion of it, and to distribute to others at no charge 
 *  materials containing or derived from the material.
 *  
 *  The developers of the software ask that you acknowledge its use in any 
 *  document referencing work based on the  program, such as published 
 *  research.  Also, they ask that you supply to Argonne National 
 *  Laboratory a copy of any published research referencing work based on 
 *  the software.
 *  
 *  Any entity desiring permission for further use must contact:
 *  
 *  J. Gleeson
 *  Industrial Technology Development Center Argonne National Laboratory
 *  Argonne IL 60439
 *  gleesonj@smtplink.eid.anl.gov
 *  (708) 252-6055
 *  
 *  ******************************************************************** 
 *  DISCLAIMER
 *  
 *  THIS PROGRAM WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY AN AGENCY 
 *  OF THE UNITED STATES GOVERNMENT.  NEITHER THE UNIVERSITY OF CHICAGO, 
 *  THE UNITED STATES GOVERNMENT NOR ANY OF THEIR EMPLOYEES MAKE ANY 
 *  WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR 
 *  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY 
 *  INFORMATION OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT 
 *  INFRINGE PRIVATELY OWNED RIGHTS.
 *  
 *  ********************************************************************** 
 *  GOVERNMENT LICENSE
 *  
 *  The Government is granted for itself and others acting on its behalf a 
 *  paid-up, non-exclusive, irrevocable worldwide license in this computer 
 *  software to reproduce, prepare derivative works, and perform publicly 
 *  and display publicly.
 */

/*****************************************************************************
*     FILE: restart.c: This file contains the routines needed to handle
*                      the restart operator, and restarting the GA.
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
*****************************************************************************/

#include "pgapack.h"

/*U****************************************************************************
   PGARestart - reseeds a population from the best string

   Category: Operators

   Inputs:
      val         - the probability of changing an allele when copying the
                    best string to the new population
      source_pop  - the source population
      dest_pop    - symbolic constant of the destination population

   Outputs:
      dest_pop is modified by side-effect.

   Example:
      Perform an unspecified test to determine if the current evolution is
      not evolving fast enough, and if so, restart the evolution.

      PGAContext *ctx;	    PGAEvaluateMS(ctx, PGA_OLDPOP, f, comm);
	    PGAFitness   (ctx, PGA_OLDPOP);
	    }

      :
      if (StagnantEvolution()) {
          PGARestart(ctx, PGA_OLDPOP, PGA_NEWPOP);
          PGAEvaluate(ctx, PGA_NEWPOP, EvalFunc);
          PGAUpdateGeneration(ctx);
      }

****************************************************************************U*/
void PGARestart(PGAContext *ctx, int source_pop, int dest_pop)
{
    /* For integers and reals, the amount by which to change is set with
       PGASetMutationIntegerValue and PGASetMutationRealValue, respectively.
       For binary strings, the bits are complemented. */

    int dest_p, old_mut_type, source_p;
    double val;
    
    PGADebugEntered("PGARestart");
    
    printf("Restarting the algorithm . . . \n");
    fflush(stdout);
    source_p = PGAGetBestIndex(ctx, source_pop);
    if (source_p != 0 || source_pop != dest_pop)
	PGACopyIndividual(ctx, source_p, source_pop, 0, dest_pop);
    PGASetEvaluationUpToDateFlag(ctx, 0, dest_pop, PGA_FALSE);
    old_mut_type = PGAGetMutationType(ctx);
    ctx->ga.MutationType = PGA_MUTATION_UNIFORM;
    val = ctx->ga.restartAlleleProb;
    
    if (ctx->fops.Mutation) {
	for (dest_p = 2; dest_p <= ctx->ga.PopSize; dest_p++) {
	    PGACopyIndividual(ctx, 0, dest_pop, dest_p-1, dest_pop);
	    (*ctx->fops.Mutation)(&ctx, &dest_p, &dest_pop, &val);
	    PGASetEvaluationUpToDateFlag(ctx, dest_p-1, dest_pop, PGA_FALSE);
	}
    } else {
	for (dest_p = 1; dest_p < ctx->ga.PopSize; dest_p++) {
	    PGACopyIndividual(ctx, 0, dest_pop, dest_p, dest_pop);
	    (*ctx->cops.Mutation)(ctx, dest_p, dest_pop, val);
	    PGASetEvaluationUpToDateFlag(ctx, dest_p, dest_pop, PGA_FALSE);
	}
    }
    ctx->ga.MutationType = old_mut_type;
    
    PGADebugExited("PGARestart");
}

/*U****************************************************************************
  PGASetRestartFlag - specifies whether the algorithm should employ
  the restart operator

   Category: Operators

   Inputs:
      ctx - context variable
      val - boolean variable

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetRestartFlag(ctx, PGA_TRUE);

****************************************************************************U*/
void PGASetRestartFlag(PGAContext *ctx, int val)
{
    PGADebugEntered("PGASetRestartFlag");

    switch (val)
    {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.restart = val;
         break;
    default:
         PGAError(ctx, "PGASetRestartFlag: Invalid value for restart:",
                  PGA_FATAL, PGA_INT, (void *) &val);
         break;
    }

    PGADebugExited("PGASetRestartFlag");
}

/*U****************************************************************************
   PGAGetRestartFlag - returns whether the algorithm should employ the
   restart operator

   Category: Operators

   Inputs:
      ctx - context variable

   Outputs:
      PGA_TRUE if restarting is enabled, otherwise PGA_FALSE.

   Example:
      PGAContext *ctx;
      int val;
      :
      val = PGAGetRestartFlag(ctx);

****************************************************************************U*/
int PGAGetRestartFlag(PGAContext *ctx)
{
    PGADebugEntered("PGAGetRestartFlag");
    PGAFailIfNotSetUp("PGAGetRestartFlag");

    PGADebugExited("PGAGetRestartFlag");

    return (ctx->ga.restart);
}

/*U****************************************************************************
  PGASetRestartFrequencyValue - specifies the number of iterations of no
  change in the best string after which the algorithm should restart

  Category: Operators

  Inputs:
      ctx - context variable
      numiter - number of changeless iterations

  Outputs:
      None

  Example:
      PGAContext *ctx;
      :
      PGASetRestartFrequencyValue(ctx, 100);

****************************************************************************U*/
void PGASetRestartFrequencyValue(PGAContext *ctx, int numiter)
{
    PGADebugEntered("PGASetRestartFrequencyValue");

    if (numiter > 0)
         ctx->ga.restartFreq = numiter;
    else
         PGAError(ctx, "PGASetRestartFrequencyValue: Invalid value for "
                  "restart freqency:", PGA_FATAL, PGA_INT, (void *) &numiter);

    PGADebugExited("PGASetRestartFrequencyValue");
}

/*U****************************************************************************
  PGAGetRestartFrequencyValue - returns the number of iterations of no
  change in the best string after which the algorithm should restart

  Category: Operators

  Inputs:
      ctx     - context variable
      numiter - number of changeless iterations

  Outputs:
      The number of iteration of no change required for a restart.

  Example:
      PGAContext *ctx;
      :
      numiter = PGAGetRestartFrequencyValue(ctx);

****************************************************************************U*/
int PGAGetRestartFrequencyValue(PGAContext *ctx)
{
    PGADebugEntered("PGAGetRestartFrequencyValue");
    PGAFailIfNotSetUp("PGAGetRestartFrequencyValue");

    PGADebugExited("PGAGetRestartFrequencyValue");

    return (ctx->ga.restartFreq);
}

/*U****************************************************************************
  PGASetRestartAlleleChangeProb - specifies the probability with which
  an allele will be mutated during a restart

  Category: Operators

  Inputs:
      ctx - context variable
      prob - probability of mutation

  Outputs:
      None

  Example:
      PGAContext *ctx;
      :
      PGASetRestartAlleleChangeProb(ctx, 0.5);

****************************************************************************U*/
void PGASetRestartAlleleChangeProb(PGAContext *ctx, double prob)
{
    PGADebugEntered("PGASetRestartAlleleChangeProb");

    if (prob >= 0.0 && prob <= 1.0)
         ctx->ga.restartAlleleProb = prob;
    else
         PGAError(ctx, "PGASetRestartAlleleChangeProb: Invalid probability:",
                  PGA_FATAL, PGA_DOUBLE, (void *) &prob);

    PGADebugExited("PGASetRestartAlleleChangeProb");
}

/*U****************************************************************************
  PGAGetRestartAlleleChangeProb - returns the probability with which
  an allele will be mutated during a restart

  Category: Operators

  Inputs:
      ctx - context variable

  Outputs:
      The probability of mutating an allele during a restart.

  Example:
      PGAContext *ctx;
      :
      prob = PGASetRestartAlleleChangeProb(ctx);

****************************************************************************U*/
double PGAGetRestartAlleleChangeProb(PGAContext *ctx)
{
    PGADebugEntered("PGAGetRestartAlleleChangeProb");
    PGAFailIfNotSetUp("PGAGetRestartAlleleChangeProb");

    PGADebugExited("PGAGetRestartAlleleChangeProb");

    return (ctx->ga.restartAlleleProb);
}

