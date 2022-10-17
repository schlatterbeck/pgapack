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
*     FILE: real.c: This file contains the routines specific to the floating
*                   point data structure
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
*****************************************************************************/

#include <pgapack.h>
/* Helper for bounds/bounce check */
static void bouncheck
    ( PGAContext *ctx, int idx, int boundflag, int bounceflag
    , PGAReal *child, PGAReal minp, PGAReal maxp
    )
{
    if (boundflag || bounceflag) {
        if (child [idx] < ctx->init.RealMin [idx]) {
            if (bounceflag) {
                child [idx] = PGARandomUniform
                    (ctx, ctx->init.RealMin [idx], minp);
            } else {
                child [idx] = ctx->init.RealMin [idx];
            }
        }
        if (child [idx] > ctx->init.RealMax [idx]) {
            if (bounceflag) {
                child [idx] = PGARandomUniform
                    (ctx, maxp, ctx->init.RealMax [idx]);
            } else {
                child [idx] = ctx->init.RealMax [idx];
            }
        }
    }
}

/*U****************************************************************************
   PGASetRealAllele - sets the value of real-valued allele i in string p
   in population pop

   Category: Fitness & Evaluation

   Inputs:
      ctx - context variable
      p   - string index
      pop - symbolic constant of the population the string is in
      i   - allele index
      val - real value to set the allele to

   Outputs:
      The specified allele in p is modified by side-effect.

   Example:
      Sets the value of the ith allele of string p in population PGA_NEWPOP
      to 1.57

      PGAContext *ctx;
      int i, p;
      :
      PGASetRealAllele ( ctx, p, PGA_NEWPOP, i, 1.57)

****************************************************************************U*/
void PGASetRealAllele (PGAContext *ctx, int p, int pop, int i, double value)
{
    PGAIndividual *ind;
    PGAReal      *chrom;

    PGADebugEntered("PGASetRealAllele");
    PGACheckDataType("PGASetRealAllele", PGA_DATATYPE_REAL);

    ind = PGAGetIndividual ( ctx, p, pop );
    chrom = (PGAReal *)ind->chrom;
    chrom[i] = value;

    PGADebugExited("PGASetRealAllele");
}

/*U****************************************************************************
   PGAGetRealAllele - returns the value of real-valued allele i in string p
   in population pop

   Category: Fitness & Evaluation

   Inputs:
      ctx - context variable
      p   - string index
      pop - symbolic constant of the population the string is in
      i   - allele index

   Outputs:
      The value of allele i

   Example:
      Returns the value of the ith real-valued allele of string p
      in population PGA_NEWPOP

      PGAContext *ctx;
      int p, i, r;
      r =  PGAGetRealAllele (ctx, p, PGA_NEWPOP, i)

****************************************************************************U*/
double PGAGetRealAllele (PGAContext *ctx, int p, int pop, int i)
{
    PGAIndividual *ind;
    PGAReal      *chrom;

    PGADebugEntered("PGAGetRealAllele");
    PGACheckDataType("PGAGetRealAllele", PGA_DATATYPE_REAL);

    ind = PGAGetIndividual ( ctx, p, pop );
    chrom = (PGAReal *)ind->chrom;

    PGADebugExited("PGAGetRealAllele");

    return( (double) chrom[i] );
}

/*U****************************************************************************
  PGASetRealInitPercent - sets the upper and lower bounds for randomly
  initializing real-valued genes.  For each gene these bounds define an
  interval from which the initial allele value is selected uniformly randomly.
  With this routine the user specifies a median value and a percent offset
  for each allele.

  Category: Initialization

  Inputs:
     ctx     - context variable
     median  - an array containing the mean value of the interval
     percent - an array containing the percent offset to add and subtract to
               the median to define the interval

  Outputs:

  Example:
     Set the initialization routines to select a value for each real-valued
     gene i uniformly randomly from the interval [i-v,i+v], where $v = i/2$.
     Assumes all strings are the same length.

     PGAContext *ctx;
     double *median, *percent;
     int i, stringlen;
     :
     stringlen = PGAGetStringLength(ctx);
     median  = (double *) malloc(stringlen*sizeof(double));
     percent = (double *) malloc(stringlen*sizeof(double));
     for(i=0;i<stringlen;i++) {
        median[i]  = (double) i;
        percent[i] = 0.5;
     }
     PGASetRealInitPercent(ctx, median, percent);

****************************************************************************U*/
void PGASetRealInitPercent ( PGAContext *ctx, double *median, double *percent)
{
     int i;
     int stringlen;
     double offset;

    PGADebugEntered("PGASetRealInitPercent");
    PGAFailIfSetUp("PGASetRealInitPercent");
    PGACheckDataType("PGASetRealInitPercent", PGA_DATATYPE_REAL);

    stringlen = PGAGetStringLength(ctx);
    for (i=0; i<stringlen; i++) {
    }
    for (i=0; i<stringlen; i++) {
         offset = fabs(median[i] * percent[i]);
         ctx->init.RealMin[i] = median[i] - offset;
         ctx->init.RealMax[i] = median[i] + offset;
    }
    ctx->init.RealType = PGA_RINIT_PERCENT;

    PGADebugExited("PGASetRealInitPercent");
}

/*U****************************************************************************
  PGASetRealInitRange - sets the upper and lower bounds for randomly
  initializing real-valued genes.  For each gene these bounds define an
  interval from which the initial allele value is selected uniformly randomly.
  The user specifies two arrays containing lower and bound for each gene to
  define the interval.  This is the default strategy for initializing
  real-valued strings.  The default interval is $[0,1.0]$ for each gene.

  Category: Initialization

  Inputs:
     ctx - context variable
     min - array containing the lower bound of the interval for each gene
     mac - array containing the upper bound of the interval for each gene

  Outputs:

  Example:
     Set the initialization routines to select a value for each real-valued
     gene i uniformly randomly from the interval [-10.,i]
     Assumes all strings are of the same length.

     PGAContext *ctx;
     double *low, *high;
     int i, stringlen;
     :
     stringlen = PGAGetStringLength(ctx);
     low  = (double *) malloc(stringlen*sizeof(double));
     high = (double *) malloc(stringlen*sizeof(double));
     for(i=0;i<stringlen;i++) {
        low[i]  = -10.0;
        high[i] = i;
     }
     PGASetRealInitRange(ctx, low, high);

****************************************************************************U*/
void PGASetRealInitRange (PGAContext *ctx, const double *min, const double *max)
{
     int i;
    PGADebugEntered("PGASetRealInitRange");
    PGAFailIfSetUp("PGASetRealInitRange");
    PGACheckDataType("PGASetRealInitRange", PGA_DATATYPE_REAL);

    for (i=ctx->ga.StringLen-1; i>=0; i--) {
         if (max[i] < min[i])
              PGAError(ctx, "PGASetRealInitRange: Lower bound exceeds upper "
                       "bound for allele #", PGA_FATAL, PGA_INT, (void *) &i);
         else
         {
              ctx->init.RealMin[i] = min[i];
              ctx->init.RealMax[i] = max[i];
         }
    }
    ctx->init.RealType = PGA_RINIT_RANGE;

    PGADebugExited("PGASetRealInitRange");
}


/*U***************************************************************************
  PGAGetMinRealInitValue - returns the minimum value used to randomly
  initialize allele i in a real string

   Category: Initialization

   Inputs:
      ctx - context variable
      i   - an allele position

   Outputs:
      The minimum value used to randomly initialize allele i

   Example:
      PGAContext *ctx;
      int min;
      :
      min = PGAGetMinRealInitValue(ctx, 0);

***************************************************************************U*/
double PGAGetMinRealInitValue (PGAContext *ctx, int i)
{
    PGADebugEntered("PGAGetMinRealInitValue");
    PGAFailIfNotSetUp("PGAGetMinRealInitValue");
    PGACheckDataType("PGAGetMinRealInitValue", PGA_DATATYPE_REAL);

    if (i < 0 || i >= ctx->ga.StringLen)
         PGAError(ctx, "PGAGetMinRealInitValue: Index out of range:",
                  PGA_FATAL, PGA_INT, (int *) &i);

    PGADebugExited("PGAGetMinRealInitValue");

    return(ctx->init.RealMin[i]);
}

/*U***************************************************************************
  PGAGetMaxRealInitValue - returns the maximum value used to randomly
  initialize allele i in a real string

   Category: Initialization

   Inputs:
      ctx - context variable
      i   - an allele position

   Outputs:
      The maximum value used to randomly initialize allele i

   Example:
      PGAContext *ctx;
      int max;
      :
      max = PGAGetMaxRealInitValue(ctx, 0);

***************************************************************************U*/
double PGAGetMaxRealInitValue (PGAContext *ctx, int i)
{
    PGADebugEntered("PGAGetMaxRealInitValue");
    PGAFailIfNotSetUp("PGAGetMaxRealInitValue");
    PGACheckDataType("PGAGetMaxRealInitValue", PGA_DATATYPE_REAL);

    if (i < 0 || i >= ctx->ga.StringLen)
         PGAError(ctx, "PGAGetMaxRealInitValue: Index out of range:",
                  PGA_FATAL, PGA_INT, (int *) &i);

    PGADebugExited("PGAGetMaxRealInitValue");

    return(ctx->init.RealMax[i]);
}


/*U***************************************************************************
  PGAGetRealInitType - returns the type of scheme used to randomly
  initialize strings of data type PGA_DATATYPE_REAL.

   Category: Initialization

   Inputs:
      ctx - context variable

   Outputs:
      Returns the integer corresponding to the symbolic constant
      used to specify the scheme used to initialize real strings

   Example:
      PGAContext *ctx;
      int inittype;
      :
      inittype = PGAGetRealInitType(ctx);
      switch (inittype) {
      case PGA_RINIT_PERCENT:
          printf("Data Type = PGA_RINIT_PERCENT\n");
          break;
      case PGA_RINIT_RANGE:
          printf("Data Type = PGA_RINIT_RANGE\n");
          break;
      }

***************************************************************************U*/
int PGAGetRealInitType (PGAContext *ctx)
{
    PGADebugEntered("PGAGetRealInitType");
    PGAFailIfNotSetUp("PGAGetRealInitType");
    PGACheckDataType("PGAGetRealInitType", PGA_DATATYPE_REAL);

    PGADebugExited("PGAGetRealInitType");

    return(ctx->init.RealType);
}


/*I****************************************************************************
   PGARealCreateString - Allocate memory for a string of type PGAReal

   Inputs:
      ctx      - context variable
      p        - string index
      pop      - symbolic constant of the population string p is in
      initflag - A true/false flag used in conjunction with ctx->ga.RandomInit
                 to initialize the string either randomly or set to zero

   Outputs:

   Example:
      Allocates memory and assigns the address of the allocated memory to
      the real string field (ind->chrom) of the individual.  Also, clears
      the string.

      PGAContext *ctx;
      int p;
      :
      PGARealCreateString( ctx, p, PGA_NEWPOP, PGA_FALSE );

****************************************************************************I*/
void PGARealCreateString (PGAContext *ctx, int p, int pop, int initflag)
{
    PGAIndividual *new = PGAGetIndividual(ctx, p, pop);
    int i, fp;
    PGAReal *c;

    PGADebugEntered("PGARealCreateString");

    new->chrom = (void *) malloc (ctx->ga.StringLen * sizeof(PGAReal));
    if (new->chrom == NULL)
	PGAError(ctx, "PGARealCreateString: No room to allocate new->chrom",
		 PGA_FATAL, PGA_VOID, NULL);
    c = (PGAReal *)new->chrom;
    if (initflag)
	if (ctx->fops.InitString) {
	    fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
	    (*ctx->fops.InitString)(&ctx, &fp, &pop);
	} else {
	    (*ctx->cops.InitString)(ctx, p, pop);
	}
    else
	for (i=ctx->ga.StringLen-1; i>=0; i--)
	    c[i] = 0.0;

    PGADebugExited("PGARealCreateString");
}

/*I****************************************************************************
   PGARealMutation - randomly mutates a floating point string with probability
   mr.  Three of the four mutation operators are of the form v = v +- p*v.
   That is, the new value of v (allele i) is the old value + or - a percentage,
   p, of the old value. There are three possibilities for choosing p: (1)
   constant value (0.01 by default), (2) selected uniformly on (0,UB) (UB is
   .1 by default), and (3) selected from a Gaussian distribution (with mean 0
   and standard deviation .1 be default).  The change to an allele, p*v, is
   added or subtracted to the old value with a probability of .5. The fourth
   option is to replace v with a value selected uniformly random from the
   initialization range of that gene. Alleles to mutate are randomly selected.
   The value set by the routine PGASetMutationRealValue is used as p, UB, and
   sigma in cases 1,2, and 3, respectively.

   Inputs:
      ctx - context variable
      p        - string index
      pop      - symbolic constant of the population string p is in
      mr  - probability of mutating a real-valued gene

   Outputs: The number of mutations performed.

   Example:
      Sets the value of the ith gene of string p
      in population PGA_NEWPOP to one

      PGAContext *ctx;
      int NumMutations, p;
      :
      NumMutations = PGARealMutation( ctx, p, PGA_NEWPOP, .001 );

****************************************************************************I*/
int PGARealMutation (PGAContext *ctx, int p, int pop, double mr)
{
    PGAReal *c;
    int i, j;
    int count = 0;
    double val = 0.0;
    /* The following are used for DE variants */
    int midx = 0;
    int do_best = (ctx->ga.DEVariant == PGA_DE_VARIANT_BEST);
    int nrand  = 2 * ctx->ga.DENumDiffs + (!do_best);
    int maxidx = 2 * ctx->ga.DENumDiffs + 1;
    DECLARE_DYNARRAY (PGAReal *, indivs, maxidx);
    int do_crossover = 1;
    static double de_dither = 0.0;
    static int last_iter = -1;


    PGADebugEntered("PGARealMutation");
    if (ctx->ga.MutationType == PGA_MUTATION_DE) {
        DECLARE_DYNARRAY (int, idx, maxidx);
        PGASampleState sstate;
        int best = 0;
        int avoid = 1;
        if (do_best) {
            best = PGAGetBestIndex (ctx, PGA_OLDPOP);
        }

        if (ctx->ga.DEDither > 0) {
            if (ctx->ga.DEDitherPerIndividual || last_iter != ctx->ga.iter) {
                de_dither = ctx->ga.DEDither * (PGARandom01 (ctx, 0) - 0.5);
                last_iter = ctx->ga.iter;
            }
        }

        /* We rely on the fact that we operate on an individual of
         * PGA_NEWPOP and take data from PGA_OLDPOP. Assert this is the
         * case.
         */
        if (pop != PGA_NEWPOP) {
            PGAError
                ( ctx, "PGARealMutation: Invalid value of pop:"
                , PGA_FATAL, PGA_INT, (void *) &(pop)
                );
        }

        /* for BEST strategy we have to avoid collision with running
         * index p *and* the best individual
         */
        if (do_best && best != p) {
            avoid = 2;
        }
        /* Use (PopSize - avoid) to avoid collision with running index
         * and with the best index (unless this is the same)
         * We may use this form of selection of the indeces needed for
         * DE for linear selection, for other selection schemes we use
         * the selected individuals and re-sample if we get collisions.
         */
        if (ctx->ga.SelectType == PGA_SELECT_LINEAR) {
            PGARandomSampleInit (ctx, &sstate, nrand, ctx->ga.PopSize - avoid);
            for (i=0; i<nrand; i++) {
                int rawidx = PGARandomNextSample (&sstate);
                /* Avoid collision with p and optionally best, samples
                 * are drawn with reduced upper bound, see
                 * PGARandomSampleInit above
                 */
                if (rawidx >= p) {
                    rawidx += 1;
                }
                /* The best individual can be in the part of the old
                 * population that is copied verbatim to the new
                 * population. In that case never increment the index
                 */
                if (do_best && best != p && rawidx >= best) {
                    rawidx += 1;
                }
                idx [i] = rawidx;
            }
        } else {
            for (i=0; i<nrand; i++) {
                do {
                    idx [i] = PGASelectNextIndex (ctx, PGA_OLDPOP);
                } while
                    (  idx [i] == p
                    || (do_best && idx [i] == best)
                    || (i > 0 && idx [i] == idx [i-1])
                    || (i > 1 && idx [i] == idx [i-2])
                    || (i > 2 && idx [i] == idx [i-3])
                    || (i > 3 && idx [i] == idx [i-4])
                    );
            }
        }
        /* Since indices from PGARandomNextSample are
         * returned in order we need to shuffle
         */
        PGAShuffle (ctx, idx, nrand);
        /* Now we have a list of shuffled indexes that do not
         * collide with one-another or with p or best
         * Add best index as last to the list if do_best
         */
        if (do_best) {
            idx [maxidx - 1] = best;
        }
        for (i=0; i<maxidx; i++) {
            indivs [i] = (PGAReal *)
                PGAGetIndividual (ctx, idx [i], PGA_OLDPOP)->chrom;
        }
        /* Index of allele that is mutated in any case */
        midx = PGARandomInterval (ctx, 0, ctx->ga.StringLen - 1);
    }

    c = (PGAReal *)PGAGetIndividual(ctx, p, pop)->chrom;
    for(i=0; i<ctx->ga.StringLen; i++) {
        double old_value = c [i];
        int idx = i;

        switch (ctx->ga.MutationType) {
            case PGA_MUTATION_RANGE:
            case PGA_MUTATION_CONSTANT:
            case PGA_MUTATION_UNIFORM:
            case PGA_MUTATION_GAUSSIAN:
            case PGA_MUTATION_POLY:
                /* randomly choose an allele   */
                if ( PGARandomFlip(ctx, mr) ) {
                    /* generate on range, or calculate multiplier */
                    switch (ctx->ga.MutationType) {
                      case PGA_MUTATION_RANGE:
                        c [i] = PGARandomUniform
                            (ctx, ctx->init.RealMin [i], ctx->init.RealMax [i]);
                        break;
                      case PGA_MUTATION_CONSTANT:
                        val = ctx->ga.MutateRealValue;
                        break;
                      case PGA_MUTATION_UNIFORM:
                        val = PGARandomUniform
                            (ctx, 0.0, ctx->ga.MutateRealValue);
                        break;
                      case PGA_MUTATION_GAUSSIAN:
                        val = PGARandomGaussian
                            (ctx, 0.0, ctx->ga.MutateRealValue);
                        break;
                      case PGA_MUTATION_POLY:
                      {
                        double u = PGARandom01 (ctx, 0);
                        double eta = PGAGetMutationPolyEta (ctx) + 1;
                        double delta;
                        if (u < 0.5) {
                            delta = pow (2 * u, 1.0 / eta) - 1.0;
                        } else {
                            delta = 1.0 - pow (2 * (1 - u), 1.0 / eta);
                        }
                        if (ctx->ga.MutatePolyValue >= 0) {
                            c [i] += delta * ctx->ga.MutatePolyValue;
                        } else {
                            if (delta < 0) {
                                val = fabs (c [i] - ctx->init.RealMin [i]);
                            } else {
                                val = fabs (ctx->init.RealMax [i] - c [i]);
                            }
                            c [i] += delta * val;
                        }
                        break;
                      }
                    }
                    /* apply multiplier calculated in switch above */
                    if ( (ctx->ga.MutationType == PGA_MUTATION_CONSTANT) ||
                         (ctx->ga.MutationType == PGA_MUTATION_UNIFORM)  ||
                         (ctx->ga.MutationType == PGA_MUTATION_GAUSSIAN)
                       )
                    {
                         /* add/subtract from allele */
                        if ( PGARandomFlip(ctx, .5) )
                            c [i] += val * c [i];
                        else
                            c [i] -= val * c [i];
                    }
                    /* increment mutation count */
                    count++;
                }
                break;
            case PGA_MUTATION_DE:
                switch (ctx->ga.DECrossoverType) {
                case PGA_DE_CROSSOVER_BIN:
                    do_crossover =
                        (  idx == midx
                        || PGARandomFlip (ctx, ctx->ga.DECrossoverProb)
                        );
                    break;
                case PGA_DE_CROSSOVER_EXP:
                    /* The first index copied is midx, then all indices
                     * are copied while the coin flip is valid
                     */
                    if (do_crossover) {
                        idx = (midx + i) % ctx->ga.StringLen;
                        if (i > 0) {
                            do_crossover =
                                (PGARandomFlip (ctx, ctx->ga.DECrossoverProb));
                        }
                    }
                    break;
                default:
                    PGAError
                        ( ctx, "PGARealMutation: Invalid DE crossover type:"
                        , PGA_FATAL, PGA_INT
                        , (void *) &(ctx->ga.DECrossoverType)
                        );
                    break;
                }
                if (do_crossover){
                    double f = ctx->ga.DEScaleFactor + de_dither;
                    if (ctx->ga.DEJitter > 0) {
                        f += ctx->ga.DEJitter * (PGARandom01 (ctx, 0) - 0.5);
                    }
                    switch (ctx->ga.DEVariant) {
                    case PGA_DE_VARIANT_RAND:
                    case PGA_DE_VARIANT_BEST:
                        /* the last element is either a random individual
                         * or the best depending on variant
                         */
                        c [idx] = indivs [maxidx - 1][idx];
                        /* Add difference vectors */
                        for (j=0; j < (maxidx - 1); j+=2) {
                            c [idx] +=
                                f * (indivs [j][idx] - indivs [j+1][idx]);
                        }
                        break;
                    case PGA_DE_VARIANT_EITHER_OR:
                        /* We use only 1 difference and ignore DENumDiffs */
                        if (PGARandom01 (ctx, 0) < ctx->ga.DEProbabilityEO) {
                            c [idx] = indivs [0][idx]
                                  + f * (indivs [1][idx] - indivs [2][idx]);
                        } else {
                            double k = ctx->ga.DEAuxFactor;
                            c [idx] = indivs [0][idx]
                                  + k * ( indivs [1][idx] + indivs [2][idx]
                                        - 2 * indivs [0][idx]
                                        );
                        }
                        break;
                    default:
                        PGAError(ctx, "PGARealMutation: Invalid value of "
                                 "ga.DEVariant:", PGA_FATAL, PGA_INT,
                                 (void *) &(ctx->ga.DEVariant));
                        break;
                    }
                    count++;
                }
                break;
            default:
                PGAError(ctx, "PGARealMutation: Invalid value of "
                         "ga.MutationType:", PGA_FATAL, PGA_INT,
                         (void *) &(ctx->ga.MutationType));
                break;
        }

        /* reset to min/max or bounce if outside range */
        bouncheck
            ( ctx, idx, ctx->ga.MutateBoundedFlag, ctx->ga.MutateBounceFlag
            , c, old_value, old_value
            );
    }

    PGADebugExited("PGARealMutation");

    return(count);
}

/*I****************************************************************************
   PGARealOneptCrossover - this routine performs one point crossover on two
   parent strings, producing (via side effect) the crossed children child1 and
   child2

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
      PGARealOneptCrossover( ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP );

****************************************************************************I*/
void PGARealOneptCrossover( PGAContext *ctx, int p1, int p2, int pop1,
                           int c1, int c2, int pop2)
{
     PGAReal *parent1 = (PGAReal *)PGAGetIndividual(ctx, p1,
                                                    pop1)->chrom;
     PGAReal *parent2 = (PGAReal *)PGAGetIndividual(ctx, p2,
                                                    pop1)->chrom;
     PGAReal *child1  = (PGAReal *)PGAGetIndividual(ctx, c1,
                                                    pop2)->chrom;
     PGAReal *child2  = (PGAReal *)PGAGetIndividual(ctx, c2,
                                                    pop2)->chrom;
     int i, xsite;

    PGADebugEntered("PGARealOneptCrossover");

    xsite = PGARandomInterval(ctx, 1,ctx->ga.StringLen-1);

    for(i=0;i<xsite;i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }

    for(i=xsite;i<ctx->ga.StringLen;i++) {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }

    PGADebugExited("PGARealOneptCrossover");
}


/*I****************************************************************************
   PGARealTwoptCrossover - performs two-point crossover on two parent strings
   producing two children via side-effect

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
      PGARealTwoptCrossover( ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP );

****************************************************************************I*/
void PGARealTwoptCrossover( PGAContext *ctx, int p1, int p2, int pop1,
                           int c1, int c2, int pop2)
{
     PGAReal *parent1 = (PGAReal *)PGAGetIndividual(ctx, p1,
                                                    pop1)->chrom;
     PGAReal *parent2 = (PGAReal *)PGAGetIndividual(ctx, p2,
                                                    pop1)->chrom;
     PGAReal *child1  = (PGAReal *)PGAGetIndividual(ctx, c1,
                                                    pop2)->chrom;
     PGAReal *child2  = (PGAReal *)PGAGetIndividual(ctx, c2,
                                                    pop2)->chrom;
     int i, temp, xsite1, xsite2;

    PGADebugEntered("PGARealTwoptCrossover");

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

    PGADebugExited("PGARealTwoptCrossover");
}


/*I****************************************************************************
   PGARealUniformCrossover - performs uniform crossover on two parent strings
   producing two children via side-effect

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
      PGARealUniformCrossover( ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP );

****************************************************************************I*/
void PGARealUniformCrossover( PGAContext *ctx, int p1, int p2, int pop1,
                             int c1, int c2, int pop2)
{
     PGAReal *parent1 = (PGAReal *)PGAGetIndividual(ctx, p1,
                                                    pop1)->chrom;
     PGAReal *parent2 = (PGAReal *)PGAGetIndividual(ctx, p2,
                                                    pop1)->chrom;
     PGAReal *child1  = (PGAReal *)PGAGetIndividual(ctx, c1,
                                                    pop2)->chrom;
     PGAReal *child2  = (PGAReal *)PGAGetIndividual(ctx, c2,
                                                    pop2)->chrom;
     int i;

    PGADebugEntered("PGARealUniformCrossover");

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

    PGADebugExited("PGARealUniformCrossover");
}

/*I****************************************************************************
   PGARealSBXCrossover - performs simulated binary crossover (SBX)
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
      PGARealSBXCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

****************************************************************************I*/
void PGARealSBXCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGAReal *parent1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *parent2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    PGAReal *child1  = (PGAReal *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    PGAReal *child2  = (PGAReal *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    int i;
    double u = 0;

    if (ctx->ga.CrossSBXOnce) {
        u = PGARandom01 (ctx, 0);
    }

    for (i=0; i<ctx->ga.StringLen; i++) {
        if (  parent1 [i] == parent2 [i]
           || !PGARandomFlip (ctx, ctx->ga.UniformCrossProb)
           )
        {
            child1 [i] = parent1 [i];
            child2 [i] = parent2 [i];
        } else {
            int j;
            PGAReal minp =
                parent1 [i] < parent2 [i] ? parent1 [i] : parent2 [i];
            PGAReal maxp =
                parent1 [i] > parent2 [i] ? parent1 [i] : parent2 [i];
            if (!ctx->ga.CrossSBXOnce) {
                u = PGARandom01 (ctx, 0);
            }
            PGACrossoverSBX
                (ctx, parent1 [i], parent2 [i], u, child1 + i, child2 + i);
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
   PGARealPrintString - writes a real-valued string to a file.  This routine
   casts the void string pointer it is passed as the second argument.

   Inputs:
      ctx - context variable
      fp  - file pointer to file to write the string to
      p   - index of the string to write out
      pop - symbolic constant of the population string p is in

   Outputs:

   Example:
      Write string s to stdout.

      PGAContext *ctx;
      int s;
      :
      PGARealPrintString( ctx, stdout, s, PGA_NEWPOP );

****************************************************************************I*/
void PGARealPrintString (PGAContext *ctx, FILE *fp, int p, int pop)
{
    PGAReal *c = (PGAReal *)PGAGetIndividual(ctx, p, pop)->chrom;
    int i;

    PGADebugEntered("PGARealPrintString");

    for(i = 0; i < ctx->ga.StringLen; i++)
    {
        switch ( i % 5 )
        {
        case 0:
            fprintf ( fp, "#%4d: [%11.7g]",i,c[i]);
            break;
        case 1:
        case 2:
        case 3:
            fprintf ( fp, ", [%11.7g]",c[i]);
            break;
        case 4:
            fprintf ( fp, ", [%11.7g]",c[i]);
            if (i+1 < ctx->ga.StringLen)
                fprintf ( fp, "\n");
            break;
        }
    }
    fprintf ( fp, "\n" );

    PGADebugExited("PGARealPrintString");
}


/*I****************************************************************************
   PGARealCopyString - Copy one real-valued string string to another

   Inputs:
      ctx - context variable
      p1   - string to copy
      pop1 - symbolic constant of population containing string p1
      p2   - string to copy p1 to
      pop2 - symbolic constant of population containing string p2

   Outputs:
      String p2 in population pop2 is modified to be a copy of string
      p1 in population pop1.

   Example:
      Copy string x to y.

      PGAContext *ctx;
      int x, y;
      :
      PGARealCopyString (ctx, x, PGA_OLDPOP, y, PGA_NEWPOP);

****************************************************************************I*/
void PGARealCopyString ( PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAReal *source = (PGAReal *)PGAGetIndividual(ctx, p1, pop1)->chrom;
    PGAReal *dest   = (PGAReal *)PGAGetIndividual(ctx, p2, pop2)->chrom;
    int i;

    PGADebugEntered("PGARealCopyString");

    for (i=ctx->ga.StringLen-1; i>=0; i--)
        *(dest++) = *(source++);

    PGADebugExited("PGARealCopyString");
}


/*I****************************************************************************
   PGARealDuplicate - Returns true if real-valued string a is a duplicate of
   real-valued string b, else returns false.

   Inputs:
      ctx - context variable
      p1   - string index of the first string to compare
      pop1 - symbolic constant of the population string p1 is in
      p2   - string index of the second string to compare
      pop2 - symbolic constant of the population string p2 is in

   Outputs:
      Returns true/false if strings are duplicates

   Example:
      Compare strings x with y to see if they are duplicates

      PGAContext *ctx;
      int x, y;
      :
      if ( PGARealDuplicate( ctx, x, PGA_OLDPOP, y, PGA_OLDPOP ) )
          printf("strings are duplicates\n");

****************************************************************************I*/
int PGARealDuplicate( PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAReal *a = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *b = (PGAReal *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int count = 0;

    PGADebugEntered ("PGARealDuplicate");

    for (count=0; count<ctx->ga.StringLen; count++) {
        if (a [count] != b [count]) {
            break;
        }
    }

    PGADebugExited ("PGARealDuplicate");

    return count == ctx->ga.StringLen ? PGA_TRUE : PGA_FALSE;
}

/*I****************************************************************************
   PGARealInitString - randomly initialize a string of type PGAReal

   Inputs:
      ctx - context variable
      p   - index of string to randomly initialize
      pop - symbolic constant of the population string p is in

   Outputs:
      String p in population pop is randomly initialized by side-effect.

   Example:
      PGAContext *ctx;
      int p;
      :
      PGARealInitString (ctx, p, PGA_NEWPOP);

****************************************************************************I*/
void PGARealInitString ( PGAContext *ctx, int p, int pop)
{
     int i;
     PGAReal *c = (PGAReal *)PGAGetIndividual(ctx, p, pop)->chrom;

     PGADebugEntered("PGARealInitString");

     for (i = 0; i < ctx->ga.StringLen; i++)
          c[i] = PGARandomUniform(ctx, ctx->init.RealMin[i],
                                  ctx->init.RealMax[i]);

     PGADebugExited("PGARealInitString");
}

/*I****************************************************************************
  PGARealBuildDatatype - Build an MPI datatype for a string.

  Inputs:
     ctx   - context variable
     p     - index of string
     pop   - symbolic constant of population string p is in

  Outputs:
     An MPI_Datatype.

  Example:
     PGAContext   *ctx;
     int           p;
     MPI_Datatype  dt;
     :
     dt = PGARealBuildDatatype(ctx, p, pop);

****************************************************************************I*/
MPI_Datatype PGARealBuildDatatype(PGAContext *ctx, int p, int pop)
{
    int            n = 7;
    int            counts[8];      /* Number of elements in each
                                      block (array of integer) */
    MPI_Aint       displs[8];      /* byte displacement of each
                                      block (array of integer) */
    MPI_Datatype   types[8];       /* type of elements in each block (array
                                      of handles to datatype objects) */
    MPI_Datatype   individualtype; /* new datatype (handle) */
    PGAIndividual *traveller;      /* address of individual in question */

    PGADebugEntered("PGARealBuildDatatype");

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
    types[3]  = MPI_DOUBLE;

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
 
    if (ctx->ga.NumAuxEval){
        MPI_Get_address(&traveller->index, &displs[7]);
        counts[7] = 1;
        types[7]  = MPI_INT;
    }
    else
    {
        MPI_Get_address(&traveller->index, &displs[6]);
        counts[6] = 1;
        types[6]  = MPI_INT;
    }


    MPI_Type_create_struct(n, counts, displs, types, &individualtype);
    MPI_Type_commit(&individualtype);

    PGADebugExited("PGARealBuildDatatype");

    return (individualtype);
}

/*I****************************************************************************
   PGARealGeneDistance - Compute genetic difference of two strings.
   Sum of the absolute values of the differences of each allele.
   So this is a Manhattan distance (mainly for performance reasons).

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
double PGARealGeneDistance (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAReal *c1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *c2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    double ret = 0.0;
    int i;

    PGADebugEntered("PGARealGeneDistance");
    for (i=0; i<ctx->ga.StringLen; i++) {
        ret += fabs (c1 [i] - c2 [i]);
    }
    PGADebugExited("PGARealGeneDistance");
    return ret;
}

/*I****************************************************************************
   PGARealEuclidianDistance - Compute genetic difference of two strings.
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
double PGARealEuclidianDistance (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGAReal *c1 = (PGAReal *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGAReal *c2 = (PGAReal *)PGAGetIndividual (ctx, p2, pop2)->chrom;

    PGADebugEntered("PGARealGeneDistance");
    PGADebugExited("PGARealGeneDistance");
    return LIN_euclidian_distance (ctx->ga.StringLen, c1, c2);
}
