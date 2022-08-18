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
*     FILE: integer.c: This file contains the routines specific to the integer
*                      data structure
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
*****************************************************************************/

#include <stdint.h>
#include "pgapack.h"

/* Helper for bounds/bounce check */
static void bouncheck
    ( PGAContext *ctx, int idx, int boundflag, int bounceflag
    , PGAInteger *child, PGAInteger minp, PGAInteger maxp
    )
{
    if (boundflag || bounceflag) {
        if (child [idx] < ctx->init.IntegerMin [idx]) {
            if (bounceflag) {
                child [idx] = PGARandomInterval
                    (ctx, ctx->init.IntegerMin [idx], minp);
            } else {
                child [idx] = ctx->init.IntegerMin [idx];
            }
        }
        if (child [idx] > ctx->init.IntegerMax [idx]) {
            if (bounceflag) {
                child [idx] = PGARandomInterval
                    (ctx, maxp, ctx->init.IntegerMax [idx]);
            } else {
                child [idx] = ctx->init.IntegerMax [idx];
            }
        }
    }
}
/* Helper for sorting / searching */
int intcmp (const void *v1, const void *v2)
{
    const PGAInteger *i1 = v1;
    const PGAInteger *i2 = v2;
    if (*i1 < *i2) {
        return -1;
    }
    if (*i1 > *i2) {
        return 1;
    }
    return 0;
}
/* Helpers for checking an allele against fixed edges */

/* Search for node in rev edges and return pointer or NULL if not found */
static PGAInteger (*rev_edge (PGAContext *ctx, PGAInteger node))[2]
{
    if (ctx->ga.n_edges == 0) {
        return NULL;
    }
    return bsearch
        ( &node
        , ctx->ga.r_edge
        , ctx->ga.n_edges
        , 2 * sizeof (PGAInteger)
        , intcmp
        );
}

/* Search for node in forward edges and return pointer or NULL if not found */
static PGAFixedEdge *get_edge (PGAContext *ctx, PGAInteger node)
{
    if (ctx->ga.n_edges == 0) {
        return NULL;
    }
    return bsearch
        ( &node
        , ctx->ga.edges
        , ctx->ga.n_edges
        , sizeof (PGAFixedEdge)
        , intcmp
        );
}

#ifdef DEBUG
static void assert_has_edges (PGAContext *ctx, PGAInteger *a)
{
    int i;
    int l = ctx->ga.StringLen;
    if (!ctx->ga.n_edges) {
        return;
    }
    for (i=0; i<l; i++) {
        PGAFixedEdge *p = get_edge (ctx, a [i]);
        if (p && !p->prev) {
            int j = 0, inc = 0;
            while (p) {
                /* Construct to easily set breakpoint on failing assert */
                if (a [(i + j + l) % l] != p->lhs) {
                    assert (a [(i + j + l) % l] == p->lhs);
                }
                /* First item */
                if (inc == 0) {
                    if (a [(i + 1) % l] == p->rhs) {
                        inc = 1;
                    } else if (a [(i - 1 + l) % l] == p->rhs) {
                        inc = -1;
                    } else {
                        assert (0);
                    }
                } else { /* Subsequent items */
                    assert (a [(i + j + inc + l) % l] == p->rhs);
                }
                p = p->next;
                j += inc;
            }
            if (j > 0) {
                i += j;
            }
        }
    }
}
#endif /* DEBUG */


/*U****************************************************************************
   PGASetIntegerAllele - sets the value of a (integer) allele.

   Category: Fitness & Evaluation

   Inputs:
      ctx - context variable
      p   - string index
      pop - symbolic constant of the population the string is in
      i   - allele index
      val - integer value to set the allele to

   Outputs:

   Example:
      Set the value of the ith allele of string p in population PGA_NEWPOP
      to 64.

      PGAContext *ctx;
      int p, i;
      :
      PGASetIntegerAllele (ctx, p, PGA_NEWPOP, i, 64)

****************************************************************************U*/
void PGASetIntegerAllele (PGAContext *ctx, int p, int pop, int i, int value)
{
    PGAIndividual *ind;
    PGAInteger     *chrom;

    PGADebugEntered("PGASetIntegerAllele");
    PGACheckDataType("PGASetIntegerAllele", PGA_DATATYPE_INTEGER);

    ind = PGAGetIndividual ( ctx, p, pop );
    chrom = (PGAInteger *)ind->chrom;
    chrom[i] = value;

    PGADebugExited("PGASetIntegerAllele");
}

/*U****************************************************************************
   PGAGetIntegerAllele - Returns the value of allele i of member p in
   population pop.  Assumes the data type is PGA_DATATYPE_INTEGER.

   Category: Fitness & Evaluation

   Inputs:
      ctx - context variable
      p   - string index
      pop - symbolic constant of the population the string is in
      i   - allele index

   Outputs:

   Example:
      Returns the value of the ith integer allele of string p
      in population PGA_NEWPOP.

      PGAContext *ctx;
      int p, i, k;
      :
      k =  PGAGetIntegerAllele ( ctx, p, PGA_NEWPOP, i )

****************************************************************************U*/
int PGAGetIntegerAllele (PGAContext *ctx, int p, int pop, int i)
{
    PGAIndividual *ind;
    PGAInteger     *chrom;

    PGADebugEntered("PGAGetIntegerAllele");

    PGACheckDataType("PGAGetIntegerAllele", PGA_DATATYPE_INTEGER);

    ind = PGAGetIndividual ( ctx, p, pop );
    chrom = (PGAInteger *)ind->chrom;

    PGADebugExited("PGAGetIntegerAllele");

    return( (int) chrom[i] );
}

/*U****************************************************************************
  PGASetIntegerInitPermute - sets a flag to tell the initialization routines
  to set each integer-valued gene to a random permutation of the values given
  by an upper and lower bound.  The length of the interval must be the same
  as the string length.  This is the default strategy for initializing
  integer-valued strings. The default interval is [0,L-1] where L is the
  string length.  No string initialization is done by this call.

  Category: Initialization

  Inputs:
     ctx - context variable
     min - the lower bound of numbers used in the permutation
     max - the upper bound of numbers used in the permutation

  Outputs:

  Example:
      Set the initialization routines to set each gene to a random and
      unique value from the interval $[500,599]$.

      PGAContext *ctx;
      :
      PGASetIntegerInitPermute(ctx, 500, 599)}

****************************************************************************U*/
void PGASetIntegerInitPermute ( PGAContext *ctx, int min, int max)
{
     int i, range;

    PGADebugEntered("PGASetIntegerInitPermute");
    PGAFailIfSetUp("PGASetIntegerInitPermute");
    PGACheckDataType("PGASetIntegerInitPermute", PGA_DATATYPE_INTEGER);

     range = max - min + 1;
     if (max <= min)
          PGAError(ctx, "PGASetIntegerInitPermute: max does not exceed min:",
                   PGA_FATAL, PGA_INT, (void *) &max);
     else if (range != ctx->ga.StringLen) {
          PGAError(ctx, "PGASetIntegerInitPermute: range of:",
                   PGA_FATAL, PGA_INT, (void *) &range);
          PGAError(ctx, "PGASetIntegerInitPermute: does not equal "
                   "string length:", PGA_FATAL, PGA_INT,
                    (void *) &(ctx->ga.StringLen));
     }
     else
     {
          ctx->init.IntegerType = PGA_IINIT_PERMUTE;
          for (i = 0; i < ctx->ga.StringLen; i++)
          {
               ctx->init.IntegerMin[i] = min;
               ctx->init.IntegerMax[i] = max;
          }
     }

    PGADebugExited("PGASetIntegerInitPermute");
}

/*U****************************************************************************
  PGASetIntegerInitRange - sets a flag to tell the initialization routines to
  set each integer-valued gene to a value chosen randomly from the interval
  given by an upper and lower bound.  No string initialization is done by
  this call.

  Category: Initialization

  Inputs:
     ctx - context variable
     min - array of lower bounds that define the interval the gene is
           initialized from
     max - array of upper bounds that define the interval the gene is
           initialized from

  Outputs:

  Example:
      Set the initialization routines to select a value for gene i
      uniformly randomly from the interval [0,i].  Assumes all strings
      are of the same length.

      PGAContext *ctx;
      int *low, *high, stringlen, i;
      :
      stringlen = PGAGetStringLength(ctx);
      low  = (int *) malloc(stringlen*sizeof(int));
      high = (int *) malloc(stringlen*sizeof(int));
      for(i=0;i<stringlen;i++) {
          low[i]  = 0;
          high[i] = i
      }
      PGASetIntegerInitRange(ctx, low, high);

****************************************************************************U*/
void PGASetIntegerInitRange (PGAContext *ctx, const int *min, const int *max)
{
     int i;

     PGADebugEntered("PGASetIntegerInitRange");
     PGAFailIfSetUp("PGASetIntegerInitRange");
     PGACheckDataType("PGASetIntegerInitRange", PGA_DATATYPE_INTEGER);

     for (i = 0; i < ctx->ga.StringLen; i++)
     {
        if (max[i] < min[i])
            PGAError(ctx, "PGASetIntegerInitRange: Lower bound exceeds upper "
                    "bound for allele #", PGA_FATAL, PGA_INT, (void *) &i);
        else {
            ctx->init.IntegerMin[i] = min[i];
            ctx->init.IntegerMax[i] = max[i];
        }
     }
     ctx->init.IntegerType = PGA_IINIT_RANGE;

     PGADebugExited("PGASetIntegerInitRange");
}

/*U***************************************************************************
  PGAGetIntegerInitType - returns the type of scheme used to randomly
  initialize strings of data type PGA_DATATYPE_INTEGER.

   Category: Initialization

   Inputs:
      ctx - context variable

   Outputs:
      Returns the integer corresponding to the symbolic constant
      used to specify the scheme used to initialize integer strings

   Example:
      PGAContext *ctx;
      int inittype;
      :
      inittype = PGAGetIntegerInitType(ctx);
      switch (inittype) {
      case PGA_IINIT_PERMUTE:
          printf("Data Type = PGA_IINIT_PERMUTE\n");
          break;
      case PGA_IINIT_RANGE:
          printf("Data Type = PGA_IINIT_RANGE\n");
          break;
      }

***************************************************************************U*/
int PGAGetIntegerInitType (PGAContext *ctx)
{
    PGADebugEntered("PGAGetIntegerInitType");
    PGAFailIfNotSetUp("PGAGetIntegerInitType");
    PGACheckDataType("PGAGetIntegerInitType", PGA_DATATYPE_INTEGER);

    PGADebugExited("PGAGetIntegerInitType");

    return(ctx->init.IntegerType);
}

/*U***************************************************************************
   PGAGetMinIntegerInitValue - returns the minimum of the range of integers
   used to randomly initialize integer strings.

   Category: Initialization

   Inputs:
      ctx - context variable

   Outputs:
      The minimum of the range of integers used to randomly initialize
      integer strings

   Example:
      PGAContext *ctx;
      int min;
      :
      min = PGAGetMinIntegerInitValue(ctx);

***************************************************************************U*/
int PGAGetMinIntegerInitValue (PGAContext *ctx, int i)
{
    PGADebugEntered("PGAGetMinIntegerInitValue");
    PGAFailIfNotSetUp("PGAGetMinIntegerInitValue");
    PGACheckDataType("PGASetIntegerAllele", PGA_DATATYPE_INTEGER);

    if (i < 0 || i >= ctx->ga.StringLen)
         PGAError(ctx, "PGAGetMinIntegerInitValue: Index out of range:",
                  PGA_FATAL, PGA_INT, (int *) &i);

    PGADebugExited("PGAGetMinIntegerInitValue");

    return(ctx->init.IntegerMin[i]);
}

/*U***************************************************************************
   PGAGetMaxIntegerInitValue - returns the maximum of the range of integers
   used to randomly initialize integer strings.

   Category: Initialization

   Inputs:
      ctx - context variable

   Outputs:
      The maximum of the range of integers used to randomly initialize
      integer strings.

   Example:
      PGAContext *ctx;
      int max;
      :
      max = PGAGetMaxIntegerInitValue(ctx);

***************************************************************************U*/
int PGAGetMaxIntegerInitValue (PGAContext *ctx, int i)
{
    PGADebugEntered("PGAGetMaxIntegerInitValue");
    PGAFailIfNotSetUp("PGAGetMaxIntegerInitValue");
    PGACheckDataType("PGAGetMaxIntegerInitValue", PGA_DATATYPE_INTEGER);

    if (i < 0 || i >= ctx->ga.StringLen)
         PGAError(ctx, "PGAGetMaxIntegerInitValue: Index out of range:",
                  PGA_FATAL, PGA_INT, (int *) &i);

    PGADebugExited("PGAGetMaxIntegerInitValue");

    return(ctx->init.IntegerMax[i]);
}


/*I****************************************************************************
   PGAIntegerCreateString - Allocate memory for a string of type PGAInteger,
   and initializes or clears the string according to initflag.

   Inputs:
      ctx      - context variable
      p        - string index
      pop      - symbolic constant of the population string p is in
      initflag - A true/false flag used in conjunction with ctx->ga.RandomInit
                 to initialize the string either randomly or set to zero

   Outputs:
      new      - a pointer set to the address of the allocated memory

   Example:
      Allocates and clears memory and assigns the address of the allocated
      memory to the string field (ind->chrom) of the individual.

      PGAContext *ctx;
      PGAIndividual *ind;
      :
      PGAIntegerCreateString( ctx, ind, PGA_FALSE );

****************************************************************************I*/
void PGAIntegerCreateString (PGAContext *ctx, int p, int pop, int InitFlag)
{
    int i, fp;
    PGAInteger *c;
    PGAIndividual *new = PGAGetIndividual(ctx, p, pop);

    PGADebugEntered("PGAIntegerCreateString");

    new->chrom = (void *)malloc(ctx->ga.StringLen * sizeof(PGAInteger));
    if (new->chrom == NULL)
	PGAError(ctx, "PGAIntegerCreateString: No room to allocate "
		 "new->chrom", PGA_FATAL, PGA_VOID, NULL);
    c = (PGAInteger *)new->chrom;
    if (InitFlag)
	if (ctx->fops.InitString) {
	    fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
	    (*ctx->fops.InitString)(&ctx, &fp, &pop);
	} else {
	    (*ctx->cops.InitString)(ctx, p, pop);
	}
    else
	for (i=0; i<ctx->ga.StringLen; i++)
	    c[i] = 0;
    
    PGADebugExited("PGAIntegerCreateString");
}

/*I****************************************************************************
   PGAIntegerMutation - randomly mutates an integer-valued gene with a
   specified probability. This routine is called from PGAMutation and must
   cast the void string pointer it is passed as the second argument.

   Inputs:
      ctx      - context variable
      p        - string index
      pop      - symbolic constant of the population string p is in
      mr       - probability of mutating an integer-valued gene

   Outputs:
      Returns the number of mutations

   Example:

****************************************************************************I*/
int PGAIntegerMutation( PGAContext *ctx, int p, int pop, double mr )
{
    PGAInteger *c;
    int i, j, temp;
    int count = 0;

    PGADebugEntered("PGAIntegerMutation");

    c = (PGAInteger *)PGAGetIndividual(ctx, p, pop)->chrom;
    for(i=0; i<ctx->ga.StringLen; i++) {
        int old_value = c [i];
        /* Do not permute fixed edges */
        if (  ctx->ga.n_edges
           && (get_edge (ctx, c [i]) || rev_edge (ctx, c [i]))
           )
        {
            continue;
        }

        /* randomly choose an allele   */
        if ( PGARandomFlip(ctx, mr) ) {
            /* apply appropriate mutation operator */
            switch (ctx->ga.MutationType) {
            case PGA_MUTATION_CONSTANT:
                /* add or subtract from allele */
                if ( PGARandomFlip(ctx, .5) )
                    c[i] += ctx->ga.MutateIntegerValue;
                else
                    c[i] -= ctx->ga.MutateIntegerValue;
                break;
            case PGA_MUTATION_PERMUTE:
            {
                /* could check for j == i if we were noble */
	        /* edd: 16 Jun 2007  applying patch from Debian bug
                 * report #333381 correcting an 'off-by-one' here
		 * bu reducing StringLen by 1
                 */
                j = PGARandomInterval(ctx, 0, ctx->ga.StringLen - 1);
                if (ctx->ga.n_edges) {
                    while (get_edge (ctx, c [j]) || rev_edge (ctx, c [j])) {
                        j = PGARandomInterval(ctx, 0, ctx->ga.StringLen - 1);
                    }
                }
                temp = c[i];
                c[i] = c[j];
                c[j] = temp;
                break;
            }
            case PGA_MUTATION_RANGE:
                c[i] = PGARandomInterval(ctx, ctx->init.IntegerMin[i],
                                              ctx->init.IntegerMax[i]);
                break;
            case PGA_MUTATION_POLY:
              {
		double u = PGARandom01 (ctx, 0);
		double eta = PGAGetMutationPolyEta (ctx) + 1;
                double delta, val;
		if (u < 0.5) {
		    delta = pow (2 * u, 1.0 / eta) - 1.0;
		} else {
		    delta = 1.0 - pow (2 * (1 - u), 1.0 / eta);
		}
		if (ctx->ga.MutatePolyValue >= 0) {
		    c [i] += (int)round (delta * ctx->ga.MutatePolyValue);
		} else {
		    if (delta < 0) {
			val = fabs (c [i] - ctx->init.IntegerMin [i] + 0.4999);
		    } else {
			val = fabs (ctx->init.IntegerMax [i] - c [i] + 0.4999);
		    }
		    c [i] += (int)round (delta * val);
		}
		break;
              }
            default:
                PGAError(ctx, "PGAIntegerMutation: Invalid value of "
                         "ga.MutationType:", PGA_FATAL, PGA_INT,
                         (void *) &(ctx->ga.MutationType));
                break;
            }

            /* reset to min/max or bounce if outside range */
	    bouncheck
                ( ctx, i, ctx->ga.MutateBoundedFlag, ctx->ga.MutateBounceFlag
                , c, old_value, old_value
                );
            count++;
        }
    }
    PGADebugExited("PGAIntegerMutation");
    return(count);
}

/*I****************************************************************************
   PGAIntegerOneptCrossover - performs one-point crossover on two parent
   strings producing two children via side-effect

   Inputs:
      ctx  - context variable
      p1   - the first parent string
      p2   - the second parent string
      pop1 - symbolic constant of the population containing string p1 and p2
      c1   - the first child string
      c2   - the second child string
      pop2 - symbolic constant of the population to contain string c1 and c2

   Outputs:

   Example:
      Performs crossover on the two parent strings m and d, producing
      children s and b.

      PGAContext *ctx;
      int m, d, s, b;
      :
      PGAIntegerOneptCrossover(ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

****************************************************************************I*/
void PGAIntegerOneptCrossover(PGAContext *ctx, int p1, int p2, int pop1,
                              int c1, int c2, int pop2)
{
     PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual(ctx, p1,
                                                          pop1)->chrom;
     PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual(ctx, p2,
                                                          pop1)->chrom;
     PGAInteger *child1  = (PGAInteger *)PGAGetIndividual(ctx, c1,
                                                          pop2)->chrom;
     PGAInteger *child2  = (PGAInteger *)PGAGetIndividual(ctx, c2,
                                                          pop2)->chrom;
     int i, xsite;

    PGADebugEntered("PGAIntegerOneptCrossover");

    xsite = PGARandomInterval(ctx, 1,ctx->ga.StringLen-1);

    for(i=0;i<xsite;i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }

    for(i=xsite;i<ctx->ga.StringLen;i++) {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }

    PGADebugExited("PGAIntegerOneptCrossover");
}


/*I****************************************************************************
   PGAIntegerTwoptCrossover - performs two-point crossover on two parent
   strings producing two children via side-effect

   Inputs:
      ctx  - context variable
      p1   - the first parent string
      p2   - the second parent string
      pop1 - symbolic constant of the population containing string p1 and p2
      c1   - the first child string
      c2   - the second child string
      pop2 - symbolic constant of the population to contain string c1 and c2

   Outputs:

   Example:
      Performs crossover on the two parent strings m and d, producing
      children s and b.

      PGAContext *ctx;
      int m, d, s, b;
      :
      PGAIntegerTwoptCrossover(ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

****************************************************************************I*/
void PGAIntegerTwoptCrossover( PGAContext *ctx, int p1, int p2, int pop1,
                              int c1, int c2, int pop2)
{
     PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual(ctx, p1,
                                                          pop1)->chrom;
     PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual(ctx, p2,
                                                          pop1)->chrom;
     PGAInteger *child1  = (PGAInteger *)PGAGetIndividual(ctx, c1,
                                                          pop2)->chrom;
     PGAInteger *child2  = (PGAInteger *)PGAGetIndividual(ctx, c2,
                                                          pop2)->chrom;
     int i, temp, xsite1, xsite2;

    PGADebugEntered("PGAIntegerTwoptCrossover");

    /* pick two cross sites such that xsite2 > xsite1 */
    xsite1 = PGARandomInterval(ctx, 1,ctx->ga.StringLen-1);
    xsite2 = xsite1;
    while ( xsite2 == xsite1 )
        xsite2 = PGARandomInterval(ctx, 1,ctx->ga.StringLen-1);
    if ( xsite1 > xsite2 ) {
        temp   = xsite1;
        xsite1 = xsite2;
        xsite2 = temp;
    }

    for(i=0;i<xsite1;i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }

    for(i=xsite1;i<xsite2;i++) {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }

    for(i=xsite2;i<ctx->ga.StringLen;i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }

    PGADebugExited("PGAIntegerTwoptCrossover");
}


/*I****************************************************************************
   PGAIntegerUniformCrossover - performs uniform crossover on two parent
   strings producing two children via side-effect

   Inputs:
      ctx  - context variable
      p1   - the first parent string
      p2   - the second parent string
      pop1 - symbolic constant of the population containing string p1 and p2
      c1   - the first child string
      c2   - the second child string
      pop2 - symbolic constant of the population to contain string c1 and c2

   Outputs:

   Example:
      Performs crossover on the two parent strings m and d, producing
      children s and b.

      PGAContext *ctx;
      int m, d, s, b;
      :
      PGAIntegerUniformCrossover( ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

****************************************************************************I*/
void PGAIntegerUniformCrossover(PGAContext *ctx, int p1, int p2, int pop1,
                                int c1, int c2, int pop2)
{
     PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual(ctx, p1,
                                                          pop1)->chrom;
     PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual(ctx, p2,
                                                          pop1)->chrom;
     PGAInteger *child1  = (PGAInteger *)PGAGetIndividual(ctx, c1,
                                                          pop2)->chrom;
     PGAInteger *child2  = (PGAInteger *)PGAGetIndividual(ctx, c2,
                                                          pop2)->chrom;
     int i;

    PGADebugEntered("PGAIntegerUniformCrossover");

    for(i=0;i<ctx->ga.StringLen;i++) {
        if ( parent1[i] == parent2[i] ) {
            child1[i] = parent1[i];
            child2[i] = parent2[i];
        }
        else {
            if(PGARandomFlip(ctx, ctx->ga.UniformCrossProb)) {
                child1[i] = parent1[i];
                child2[i] = parent2[i];
            }
            else {
                child1[i] = parent2[i];
                child2[i] = parent1[i];
            }
        }
    }

    PGADebugExited("PGAIntegerUniformCrossover");
}

/*I****************************************************************************
   PGAIntegerSBXCrossover - performs simulated binary crossover (SBX)
   on two parent strings producing two children via side-effect

   Inputs:
      ctx  - context variable
      p1   - the first parent string
      p2   - the second parent string
      pop1 - symbolic constant of the population containing string p1 and p2
      c1   - the first child string
      c2   - the second child string
      pop2 - symbolic constant of the population to contain string c1 and c2

   Outputs:
      c1 and c2 in population pop2 are modified by side-effect.

   Example:
      Performs crossover on the two parent strings m and d, producing
      children s and b.

      PGAContext *ctx;
      int m, d, s, b;
      :
      PGAIntegerSBXCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

****************************************************************************I*/
void PGAIntegerSBXCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent1 = (PGAInteger *)PGAGetIndividual
        (ctx, p1, pop1)->chrom;
    PGAInteger *parent2 = (PGAInteger *)PGAGetIndividual
        (ctx, p2, pop1)->chrom;
    PGAInteger *child1  = (PGAInteger *)PGAGetIndividual
        (ctx, c1, pop2)->chrom;
    PGAInteger *child2  = (PGAInteger *)PGAGetIndividual
        (ctx, c2, pop2)->chrom;
    int i;
    double u = 0;

    if (ctx->ga.CrossSBXOnce) {
        u = PGARandom01 (ctx, 0);
    }

    for (i=0; i<ctx->ga.StringLen; i++) {
        if (  parent1[i] == parent2[i]
           || !PGARandomFlip (ctx, ctx->ga.UniformCrossProb)
           )
        {
            child1 [i] = parent1 [i];
            child2 [i] = parent2 [i];
        } else {
            int j;
            double c1, c2;
            PGAInteger minp =
                parent1 [i] < parent2 [i] ? parent1 [i] : parent2 [i];
            PGAInteger maxp =
                parent1 [i] > parent2 [i] ? parent1 [i] : parent2 [i];
            if (!ctx->ga.CrossSBXOnce) {
                u = PGARandom01 (ctx, 0);
            }
            PGACrossoverSBX (ctx, parent1 [i], parent2 [i], u, &c1, &c2);
            child1 [i] = (int)round (c1);
            child2 [i] = (int)round (c2);
            for (j=0; j<2; j++) {
                bouncheck
                    ( ctx, i, ctx->ga.CrossBoundedFlag, ctx->ga.CrossBounceFlag
                    , j ? child2 : child1, minp, maxp
                    );
            }
        }
    }
}

/*I****************************************************************************
   PGAIntegerEdgeCrossover - performs Edge Recombination
   on two parent strings producing two children via side-effect

   Inputs:
      ctx  - context variable
      p1   - the first parent string
      p2   - the second parent string
      pop1 - symbolic constant of the population containing string p1 and p2
      c1   - the first child string
      c2   - the second child string
      pop2 - symbolic constant of the population to contain string c1 and c2

   Outputs:
      c1 and c2 in population pop2 are modified by side-effect.

   Example:
      Performs crossover on the two parent strings m and d, producing
      children s and b.

      PGAContext *ctx;
      int m, d, s, b;
      :
      PGAIntegerEdgeCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

****************************************************************************I*/

static void append_edge (PGAContext *ctx, PGAInteger n1, PGAInteger n2)
{
    int i;
    PGAInteger *em = ctx->scratch.edgemap [n1];
    for (i=0; i<4; i++) {
        if (em [i] == 0) {
            em [i] = n2 + 1;
            return;
        } else if (em [i] == n2 + 1) {
            em [i] = -em [i];
            return;
        } else {
            assert (-em [i] != n2 + 1);
        }
    }
}

static void build_edge_map (PGAContext *ctx, PGAInteger **parent)
{
    int j;
    PGAInteger i, l = ctx->ga.StringLen;
    memset (ctx->scratch.edgemap, 0, sizeof (PGAInteger) * 4 * l);
    unsigned long long s [2] = {0, 0};
    for (i=0; i<l; i++) {
        for (j=0; j<2; j++) {
            PGAInteger n1 = parent [j][i];
            PGAInteger n2 = parent [j][(i + 1) % l];
            if (n1 < 0 || n1 >= l || n2 < 0 || n2 >= l) {
                PGAErrorPrintf
                    (ctx, PGA_FATAL, "Parent gene is no permutation");
            }
            append_edge (ctx, n1, n2);
            append_edge (ctx, n2, n1);
            s [j] += parent [j][i];
        }
    }
    for (j=0; j<2; j++) {
        if (s [j] != (unsigned long long)(l) * (unsigned long long)(l-1) / 2) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Parent gene is no permutation");
        }
    }
}

/* Remove new edge from all right sides, note that only the edges in
 * the edge table for the new edge have this edge on the right side
 */
void remove_edge_from_right (PGAContext *ctx, PGAInteger cidx)
{
    int i, j;
    for (j=0; j<4; j++) {
        PGAInteger v = labs (ctx->scratch.edgemap [cidx][j]) - 1;
        for (i=0; i<4; i++) {
            if (labs (ctx->scratch.edgemap [v][i]) - 1 == cidx) {
                ctx->scratch.edgemap [v][i] = 0;
                break;
            }
        }
    }
}

static void fix_edge_map (PGAContext *ctx)
{
    size_t i, j, k;
    PGAFixedEdge *e;
    /* We only need to do something if we have fixed edges */
    if (ctx->ga.n_edges == 0) {
        return;
    }
    for (i=0; i<ctx->ga.n_edges; i++) {
        e = ctx->ga.edges + i;
        if (ctx->ga.symmetric) {
            /* If intermediate edge we don't need to do anything */
            if (e->prev && e->next) {
                continue;
            }
            /* Ensure that fixed edge is always preferred */
            for (j=0; j<2; j++) {
                PGAInteger node  = j == 0 ? e->lhs : e->rhs;
                PGAInteger other = j == 0 ? e->rhs : e->lhs;
                PGAInteger *em = ctx->scratch.edgemap [node];
                /* Ignore warning if compiling without assertions */
                #pragma GCC diagnostic push
                #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
                int found = 0;
                #pragma GCC diagnostic pop
                if (j == 0 && e->prev) {
                    continue;
                }
                if (j == 1 && e->next) {
                    continue;
                }
                for (k=0; k<4; k++) {
                    if (abs (em [k]) - 1 == other) {
                        /* This assertion fails if not both of the
                         * parents have the fixed edge
                         */
                        assert (em [k] < 0);
                        found = 1;
                    } else {
                        em [k] = abs (em [k]);
                    }
                }
                assert (found);
            }
        } else {
            /* Asymmetric case, edge must have given orientation */
            /* right side may never occur in the right side of edge table */
            remove_edge_from_right (ctx, e->rhs);
            /* Entry in edgemap can only have fixed right side */
            ctx->scratch.edgemap [e->lhs][0] = e->rhs;
            ctx->scratch.edgemap [e->lhs][1]
                = ctx->scratch.edgemap [e->lhs][2]
                = ctx->scratch.edgemap [e->lhs][3] = 0;
        }
    }
}

int count_edges (PGAContext *ctx, PGAInteger idx)
{
    int j;
    int c = 0;
    PGAInteger *em = ctx->scratch.edgemap [idx];
    for (j=0; j<4; j++) {
        if (em [j]) {
            c++;
        }
    }
    return c;
}

void next_edge (PGAContext *ctx, PGAInteger *child, PGAInteger idx)
{
    PGAInteger i, j;
    PGAInteger *em = ctx->scratch.edgemap [child [idx]];
    assert (idx < ctx->ga.StringLen - 1);
    /* Prefer common edges */
    if (em [0] < 0 && em [1] >= 0) {
        child [idx + 1] = -em [0] - 1;
    } else if (em [1] < 0 && em [0] >= 0) {
        child [idx + 1] = -em [1] - 1;
    } else if (em [2] < 0) {
        child [idx + 1] = -em [2] - 1;
    } else {
        PGAInteger idxm = 0;
        PGAInteger emin [4];
        int minv = -1;
        for (j=0; j<4; j++) {
            int v;
            if (em [j] == 0) {
                continue;
            }
            v = count_edges (ctx, labs (em [j]) - 1);
            if (idxm == 0 || v < minv) {
                minv = v;
                emin [0] = labs (em [j]) - 1;
                idxm = 1;
            } else if (v == minv) {
                emin [idxm++] = labs (em [j]) - 1;
            }
        }
        if (idxm == 1) {
            child [idx + 1] = emin [0];
        } else if (idxm > 1) {
            child [idx + 1] = emin [PGARandomInterval (ctx, 0, idxm - 1)];
        } else {
            PGAInteger used [idx + 1];
            PGAInteger mini = -1;
            PGAInteger lastu = 0;
            memcpy (used, child, sizeof (PGAInteger) * (idx + 1));
            qsort (used, idx + 1, sizeof (PGAInteger), intcmp);
            for (j=0; j<idx+1; j++) {
                for (i=lastu; i<used [j]; i++) {
                    int ec = count_edges (ctx, i);
                    PGAFixedEdge *p = get_edge (ctx, i);
                    /* May not start in the middle of a run of fixed edges */
                    if (p && p->prev) {
                        continue;
                    }
                    /* If we have asymmetric (directed) fixed edges and
                     * this is the last of a run (intermediate cases
                     * above) we may not start with that edge
                     */
                    if (!ctx->ga.symmetric && NULL != rev_edge (ctx, i)) {
                        continue;
                    }
                    if (minv < 0 || ec < minv) {
                        minv = ec;
                        mini = i;
                    }
                }
                lastu = used [j] + 1;
            }
            /* don't randomize different indexes with same count */
            if (mini < 0) {
                /* The current (sorted) list is tight
                 * Or did contain only invalid starting points
                 * (right sides of fixed edges)
                 */
                for (j=idx+1; j<ctx->ga.StringLen; j++) {
                    PGAFixedEdge *p = get_edge (ctx, j);
                    /* Don't use intermediate fixed edge */
                    if (p && p->prev) {
                        continue;
                    }
                    if (ctx->ga.symmetric || NULL == rev_edge (ctx, j)) {
                        child [idx + 1] = j;
                        break;
                    }
                }
            } else {
                child [idx + 1] = mini;
            }
        }
    }
    remove_edge_from_right (ctx, child [idx + 1]);
}

void PGAIntegerEdgeCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAInteger *parent [2];
    PGAInteger *child  [2];
    int j, oc;
    PGAInteger i, ci;
    PGAInteger l = ctx->ga.StringLen;
    int p_ok [] = {1, 1};
    PGAInteger ok [2];

    parent [0] = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent [1] = (PGAInteger *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child  [0] = (PGAInteger *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child  [1] = (PGAInteger *)PGAGetIndividual (ctx, c2, pop2)->chrom;

    /* Build Edge-Map, 0 is unallocated, a node is represented by index+1 */
    build_edge_map (ctx, parent);
    fix_edge_map   (ctx);
    /* Find a node with edge-count < 4 for each child
     * If not possible (e.g. for the 2nd child) we use first node of one
     * of the parents.
     */
    ci = 0;
    oc = 0;
    for (i=0; i<l; i++) {
        PGAFixedEdge *p = get_edge (ctx, i);
        /* May not start with intermediate edge of fixed-edge run */
        if (p && p->prev) {
            continue;
        }
        /* May not start with end node if asymmetric fixed edges */
        if (!ctx->ga.symmetric && rev_edge (ctx, i)) {
            continue;
        }
        for (j=0; j<4; j++) {
            if (!ctx->scratch.edgemap [i][j]) {
                child [ci++][0] = i;
                break;
            }
        }
        if (oc < 2 && j == 4) {
            ok [oc++] = i;
        }
        if (ci >= 2) {
            break;
        }
    }
    if (ci < 2) {
        for (i=0; i<2; i++) {
            PGAFixedEdge *p = get_edge (ctx, parent [i][0]);
            /* May not start with intermediate edge of fixed-edge run */
            if (p && p->prev) {
                p_ok [i] = 0;
                continue;
            }
            /* May not start with end node if asymmetric fixed edges */
            if (!ctx->ga.symmetric && rev_edge (ctx, i)) {
                p_ok [i] = 0;
                continue;
            }
        }
        assert (oc == 2);
        if (ci == 0) {
            if (p_ok [0]) {
                child [0][0] = parent [0][0];
            } else {
                child [0][0] = ok [0];
            }
            if (p_ok [1]) {
                child [1][0] = parent [1][0];
            } else {
                child [1][0] = ok [1];
            }
        } else {
            if (child [0][0] == parent [0][0] && p_ok [1]) {
                child [1][0] = parent [1][0];
            } else if (child [0][0] == parent [1][0] && p_ok [0]) {
                child [1][0] = parent [0][0];
            } else if (!(p_ok [0] && p_ok [1])) {
                child [1][0] = ok [0];
            } else {
                child [1][0] = parent [PGARandomFlip (ctx, 0.5)][0];
            }
        }
    }
    remove_edge_from_right (ctx, child [0][0]);
    /* Now do the actual crossover */
    for (ci=0; ci<2; ci++) {
        for (i=0; i<l-1; i++) {
            next_edge (ctx, child [ci], i);
        }
        #ifdef DEBUG
        assert_has_edges (ctx, child [ci]);
        #endif /* DEBUG */

        /* Need to rebuild, edge-map is consumed above */
        build_edge_map (ctx, parent);
        fix_edge_map   (ctx);
        remove_edge_from_right (ctx, child [1][0]);
    }
}

/*I****************************************************************************
   PGAIntegerSetFixedEdges - Set edges that have to be present
   This is used only in Edge Crossover
   Note: The edges data structure is copied and must be freed by the
   caller. It is admissible that the edges data is in automatic
   variables allocated on the stack.

   Inputs:
      ctx       - context variable
      n         - Number of edges
      edge      - Pointer to edges, each edge consists of two indices
      symmetric - Flag that indicates if edges are allowed in reverse
                  direction

   Outputs:

   Example:
      Set edges (0, 213) and (7, 11) as a fixed edges

      PGAContext *ctx;
      PGAInteger edge[][2] = {(PGAInteger []){0, 213}, (PGAInteger []){7, 11}};
      int  n = sizeof (edge) / (2 * sizeof (PGAInteger));
      :
      PGAIntegerSetFixedEdges (ctx, n, edge, PGA_TRUE);

****************************************************************************I*/
void PGAIntegerSetFixedEdges
    (PGAContext *ctx, size_t n, PGAInteger (*edge)[2], int symmetric)
{
    size_t i, j;
    PGAInteger prev;
    /* Do not allocate twice */
    assert (ctx->ga.edges == NULL);
    ctx->ga.edges = malloc (sizeof (PGAFixedEdge) * n);
    if (ctx->ga.edges == NULL) {
        PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate edges");
    }
    ctx->ga.r_edge = malloc (2 * sizeof (PGAInteger) * n);
    if (ctx->ga.r_edge == NULL) {
        PGAErrorPrintf (ctx, PGA_FATAL, "Cannot allocate rev-edges");
    }
    for (i=0; i<n; i++) {
        ctx->ga.edges  [i].lhs  = edge [i][0];
        ctx->ga.edges  [i].rhs  = edge [i][1];
        ctx->ga.edges  [i].next = ctx->ga.edges [i].prev = NULL;
    }
    qsort (ctx->ga.edges,  n, sizeof (PGAFixedEdge),   intcmp);
    for (i=0; i<n; i++) {
        ctx->ga.r_edge [i][0] = ctx->ga.edges [i].rhs;
        ctx->ga.r_edge [i][1] = i;
    }
    qsort (ctx->ga.r_edge, n, 2 * sizeof (PGAInteger), intcmp);
    ctx->ga.n_edges = n;
    ctx->ga.symmetric = symmetric ? PGA_TRUE : PGA_FALSE;
    /* predecessors, cycle check */
    prev = -1;
    for (i=0; i<n; i++) {
        PGAFixedEdge *p;
        if (ctx->ga.edges [i].lhs == prev) {
            PGAErrorPrintf (ctx, PGA_FATAL, "Fixed edges may not share nodes");
        }
        prev = ctx->ga.edges [i].lhs;
        /* Try finding predecessor */
        p = get_edge (ctx, ctx->ga.edges [i].rhs);
        if (p != NULL) {
            if (p->prev) {
                PGAErrorPrintf
                    (ctx, PGA_FATAL, "Fixed edges may not share nodes");
            } else {
                p->prev = ctx->ga.edges + i;
            }
        }
    }
    /* Now loop again over all nodes and perform cycle check:
     * We follow the predecessor indexes until we encounter a stop (-1)
     * or we have exhausted n.
     */
    for (i=0; i<n; i++) {
        if (ctx->ga.edges [i].prev) {
            PGAFixedEdge *p = ctx->ga.edges + i;
            for (j=0; j<n; j++) {
                if (p->prev) {
                    PGAFixedEdge *np = p->prev;
                    if (np->next) {
                        assert (np->next == p);
                    } else {
                        np->next = p;
                    }
                    p = np;
                } else {
                    break;
                }
            }
            if (j == n) {
                PGAErrorPrintf (ctx, PGA_FATAL, "Fixed edges contain a cycle");
            }
        }
    }
}

/*I****************************************************************************
   PGAIntegerPrintString - writes an integer-valued string to a file.

   Inputs:
      ctx - context variable
      fp  - file pointer to file to write the string to
      p   - index of the string to write out
      pop - symbolic constant of the population string p is in

   Outputs:

   Example:
      Write member p in population PGA_NEWPOP to stdout.

      PGAContext *ctx;
      int  p;
      :
      PGAIntegerPrintString(ctx, stdout, p, PGA_NEWPOP);

****************************************************************************I*/
void PGAIntegerPrintString ( PGAContext *ctx, FILE *fp, int p, int pop)
{
    PGAInteger *c = (PGAInteger *)PGAGetIndividual(ctx, p, pop)->chrom;
    int i;

    PGADebugEntered("PGAIntegerPrintString");

    for(i = 0; i < ctx->ga.StringLen; i++)
    {
        switch ( i % 6 )
        {
        case 0:
            fprintf ( fp, "#%5d: [%8ld]",i,c[i]);
            break;
        case 1:
        case 2:
        case 3:
        case 4:
            fprintf ( fp, ", [%8ld]",c[i]);
            break;
        case 5:
            fprintf ( fp, ", [%8ld]",c[i]);
            if (i+1 < ctx->ga.StringLen)
                fprintf ( fp, "\n");
            break;
        }
    }
    fprintf ( fp, "\n" );

    PGADebugExited("PGAIntegerPrintString");
}

/*I****************************************************************************
   PGAIntegerCopyString - Copy one integer-valued string to another.

   Inputs:
      ctx - context variable
      p1   - string to copy
      pop1 - symbolic constant of population containing string p1
      p2   - string to copy p1 to
      pop2 - symbolic constant of population containing string p2

   Outputs:

   Example:

****************************************************************************I*/
void PGAIntegerCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAInteger *source = (PGAInteger *)PGAGetIndividual(ctx, p1, pop1)->chrom;
    PGAInteger *dest   = (PGAInteger *)PGAGetIndividual(ctx, p2, pop2)->chrom;
    int i;

    PGADebugEntered("PGAIntegerCopyString");

    for (i = 0; i < ctx->ga.StringLen; i++)
        dest[i] = source[i];

    PGADebugExited("PGAIntegerCopyString");
}

/*I****************************************************************************
   PGAIntegerDuplicate - Returns true if string a is a duplicate of
   string b, else returns false.

   Inputs:
      ctx - context variable
      p1   - string index of the first string to compare
      pop1 - symbolic constant of the population string p1 is in
      p2   - string index of the second string to compare
      pop2 - symbolic constant of the population string p2 is in

   Outputs:
      Returns true/false if strings are duplicates

   Example:

****************************************************************************I*/
int PGAIntegerDuplicate( PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAInteger *a = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAInteger *b = (PGAInteger *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int count = 0;

    PGADebugEntered ("PGAIntegerDuplicate");

    for (count=0; count<ctx->ga.StringLen; count++) {
        if (a [count] != b [count]) {
            break;
        }
    }

    PGADebugExited("PGAIntegerDuplicate");

    return count == ctx->ga.StringLen ? PGA_TRUE : PGA_FALSE;
}

/*I****************************************************************************
   PGAIntegerInitString - randomly initialize a string of type PGAInteger

   Inputs:
      ctx - context variable
      p   - index of string to randomly initialize
      pop - symbolic constant of the population string p is in

   Outputs:

   Example:

****************************************************************************I*/

/* Helper function for computing index into gene for given value if it
 * is part of a fixed edge
 */
static void compute_idx
    (PGAContext *ctx, size_t (*x)[2], size_t k, PGAInteger a)
{
    PGAInteger (*r)[2];
    PGAFixedEdge *p;

    p = get_edge (ctx, a);
    if (p != NULL) {
        size_t off = p - ctx->ga.edges;
        x [off][0] = k;
        if (p->prev) {
            assert (p->prev->rhs == a);
            off = p->prev - ctx->ga.edges;
            x [off][1] = k;
        }
    } else if (NULL != (r = rev_edge (ctx, a))) {
        size_t off = (*r) [1];
        assert ((ctx->ga.edges + off)->rhs == a);
        x [off][1] = k;
    }
}

void PGAIntegerInitString (PGAContext *ctx, int p, int pop)
{
    int *list;
    int len, i, j;
    PGAInteger *c = (PGAInteger *)PGAGetIndividual (ctx, p, pop)->chrom;

    PGADebugEntered ("PGAIntegerInitString");

    len = ctx->ga.StringLen;

    switch (ctx->init.IntegerType)
    {
    case PGA_IINIT_PERMUTE:
        list = (int *)malloc (sizeof(int) * len);
        if (list == NULL)
            PGAErrorPrintf
                ( ctx, PGA_FATAL
                , "PGAIntegerInitString: No room to allocate list"
                );
        j = ctx->init.IntegerMin [0];
        for (i=0; i<len; i++) {
             list [i] = j++;
        }
        for (i=0; i<len; i++) {
             j = PGARandomInterval (ctx, 0, len - i - 1);
             c [i] = list [j];
             list [j] = list [len - i - 1];
        }
        free (list);
        /* Ensure fixed edges if configured, all fixed edges are used in
         * forward direction regardless of the symmetric flag
         */
        if (ctx->ga.n_edges) {
            size_t (*x)[2] = malloc (2 * ctx->ga.n_edges * sizeof (size_t));
            size_t k, j = 0;
            if (x == NULL) {
                PGAErrorPrintf
                    (ctx, PGA_FATAL, "PGAIntegerInitString: Cannot allocate");
            }
            for (k=0; k<ctx->ga.n_edges; k++) {
                x [k][0] = x [k][1] = SIZE_MAX;
            }
            for (k=0; k<(size_t)len; k++) {
                compute_idx (ctx, x, k, c [k]);
            }
            for (k=0; k<ctx->ga.n_edges; k++) {
                assert (x [k][0] != SIZE_MAX && x [k][1] != SIZE_MAX);
            }
            for (k=0; k<ctx->ga.n_edges; k++) {
                PGAInteger v;
                PGAFixedEdge *p = ctx->ga.edges + k;
                size_t ix;
                if (p->prev) {
                    continue;
                }
                /* We start at 0, otherwise sequences may overlap and
                 * destroy each other
                 */
                ix = k;
                /* Relocate first item of sequence */
                v = c [j];
                c [j] = c [x [ix][0]];
                c [x [ix][0]] = v;
                compute_idx (ctx, x, x [ix][0], v);
                j++;
                do {
                    size_t tmp;
                    v = c [j];
                    c [j] = c [x [ix][1]];
                    c [x [ix][1]] = v;
                    compute_idx (ctx, x, x [ix][1], v);
                    j++;
                    p = p->next;
                    if (p) {
                        tmp = p - ctx->ga.edges;
                        assert  (x [tmp][0] == x [ix][1]);
                        x [tmp][0] = j;
                        ix = tmp;
                    }
                } while (p);
            }
            free (x);
            #ifdef DEBUG
            assert_has_edges (ctx, c);
            # endif /* DEBUG */
        }
        break;
    case PGA_IINIT_RANGE:
        for (i = 0; i < len; i++)
            c[i] = PGARandomInterval
                (ctx, ctx->init.IntegerMin[i], ctx->init.IntegerMax[i]);
        break;
    }

    PGADebugExited ("PGAIntegerInitString");
}

/*I****************************************************************************
  PGAIntegerBuildDatatype - Build an MPI datatype for a string of type
  PGA_DATATYPE_INTEGER.

  Inputs:
      ctx - context variable
      p   - index of string to randomly initialize
      pop - symbolic constant of the population string p is in

  Outputs:

  Example:

****************************************************************************I*/
MPI_Datatype PGAIntegerBuildDatatype(PGAContext *ctx, int p, int pop)
{
    int            n = 6;
    int            counts[7];      /* Number of elements in each
                                      block (array of integer) */
    MPI_Aint       displs[7];      /* byte displacement of each
                                      block (array of integer) */
    MPI_Datatype   types[7];       /* type of elements in each block (array
                                      of handles to datatype objects) */
    MPI_Datatype   individualtype; /* new datatype (handle) */
    PGAIndividual *traveller;      /* address of individual in question */

    PGADebugEntered("PGAIntegerBuildDatatype");

    traveller = PGAGetIndividual(ctx, p, pop);
    MPI_Get_address(&traveller->evalue, &displs[0]);
    counts[0] = 1;
    types[0]  = MPI_DOUBLE;

    MPI_Get_address(&traveller->fitness, &displs[1]);
    counts[1] = 1;
    types[1]  = MPI_DOUBLE;

    MPI_Get_address(&traveller->evaluptodate, &displs[2]);
    counts[2] = 1;
    types[2]  = MPI_INT;

    MPI_Get_address(traveller->chrom, &displs[3]);
    counts[3] = ctx->ga.StringLen;
    types[3]  = MPI_LONG;

    MPI_Get_address(&traveller->auxtotal, &displs[4]);
    counts[4] = 1;
    types[4]  = MPI_DOUBLE;

    MPI_Get_address(&traveller->auxtotalok, &displs[5]);
    counts[5] = 1;
    types[5]  = MPI_INT;

    if (ctx->ga.NumAuxEval) {
        MPI_Get_address(traveller->auxeval, &displs[6]);
        counts[6] = ctx->ga.NumAuxEval;
        types[6]  = MPI_DOUBLE;
        n += 1;
    }

    MPI_Type_create_struct(n, counts, displs, types, &individualtype);
    MPI_Type_commit(&individualtype);

    PGADebugExited("PGAIntegerBuildDatatype");

    return (individualtype);
}

/*I****************************************************************************
   PGAIntegerGeneDistance - Compute genetic difference of two strings.
   Sum of the absolute values of the differences of each allele.

   Inputs:
      ctx   - context variable
      p1    - first string index
      pop1  - symbolic constant of the population the first string is in
      p2    - second string index
      pop2  - symbolic constant of the population the second string is in

   Outputs:
      genetic distance of the two strings

   Example:
      Internal function.  Use PGAGeneDistance.

****************************************************************************I*/
double PGAIntegerGeneDistance (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAInteger *c1 = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAInteger *c2 = (PGAInteger *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int ret = 0;
    int i;

    PGADebugEntered("PGAIntegerGeneDistance");
    for (i=0; i<ctx->ga.StringLen; i++) {
        ret += labs (c1 [i] - c2 [i]);
    }
    PGADebugExited("PGAIntegerGeneDistance");
    return ret;
}

/*I****************************************************************************
   PGAIntegerEuclidianDistance - Compute genetic difference of two strings.
   This uses the Euclidian distance metric, the square-root of the sum
   of all squared differences of each allele.

   Inputs:
      ctx   - context variable
      p1    - first string index
      pop1  - symbolic constant of the population the first string is in
      p2    - second string index
      pop2  - symbolic constant of the population the second string is in

   Outputs:
      genetic euclidian distance of the two strings

   Example:
      Use in PGASetUserFunction.

****************************************************************************I*/
double PGAIntegerEuclidianDistance (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAInteger *c1 = (PGAInteger *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAInteger *c2 = (PGAInteger *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    double ret = 0.0;
    int i;

    PGADebugEntered("PGAIntegerGeneDistance");
    for (i=0; i<ctx->ga.StringLen; i++) {
        ret += pow (c1 [i] - c2 [i], 2);
    }
    PGADebugExited("PGAIntegerGeneDistance");
    return sqrt (ret);
}

