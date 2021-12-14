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
*     FILE: pop.c: This file contains systme routines that act on entire
*                  populations.
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
*****************************************************************************/

#include <assert.h>
#include "pgapack.h"

/*U****************************************************************************
   PGASortPop - Creates an (internal) array of indices according to one of
   three criteria.  If PGA_POPREPL_BEST is used (the default) the array is
   sorted from most fit to least fit.  If PGA_POPREPL_RANDOM_REP is
   used the indices in the array are selected randomly with replacement.
   If PGA_POPREPL_RANDOM_NOREP is used the indices in the array are selected
   randomly without replacement.  The function PGASetPopReplaceType() is used
   to specify which strategy is used.  The indices of the sorted population
   members may then be accessed from the internal array via
   PGAGetSortedPopIndex().  This routine is typically used during population
   replacement.

   Category: Generation

   Inputs:
       ctx      - context variable
       popindex - symbolic constant of the population from which to create
                  the srted array.

   Output:
      An inteneral array of indices sorted according to one of three
      criteria is created.

   Example:
      Copy the five best strings from the old population into the new
      population.  The rest of the new population will be created by
      recombination, and is not shown.

      PGAContext *ctx;
      int i,j;
      :
      PGASetPopReplaceType(ctx,PGA_POPREPL_BEST)
      :
      PGASortPop(ctx, PGA_OLDPOP);
      for ( i=0; i < 5; i++) {
          j = PGAGetSortedPopIndex(ctx, i);
      PGACopyIndividual (ctx, j, PGA_OLDPOP, i, PGA_NEWPOP);
      :

****************************************************************************U*/
void PGASortPop (PGAContext *ctx, int pop)
{
    int i,j;

    PGADebugEntered ("PGASortPop");
    if (pop != PGA_OLDPOP && pop != PGA_NEWPOP) {
        PGAError(ctx,
                 "PGASort: Invalid value of pop:",
                 PGA_FATAL,
                 PGA_INT,
                 (void *) &pop);
    }
    switch (ctx->ga.PopReplace) {
    case PGA_POPREPL_BEST:
        /* No need to init ga.sorted, done by PGAEvalSort */
        PGAEvalSort (ctx, pop, ctx->ga.sorted);
        break;
    case PGA_POPREPL_RANDOM_REP:
        for (i = 0; i < ctx->ga.PopSize; i++) {
            ctx->scratch.intscratch [i] = i;
        };
        for (i = 0; i < ctx->ga.PopSize; i++) {
            j = PGARandomInterval (ctx, 0, ctx->ga.PopSize-1);
            ctx->ga.sorted [i] = ctx->scratch.intscratch [j];
        };
        break;
    case PGA_POPREPL_RANDOM_NOREP:
        for (i = 0; i < ctx->ga.PopSize; i++) {
            ctx->scratch.intscratch [i] = i;
        };
        for (i = 0; i < ctx->ga.PopSize; i++) {
            j = PGARandomInterval (ctx, 0, ctx->ga.PopSize-i-1);
            ctx->ga.sorted[i] = ctx->scratch.intscratch[j];
            ctx->scratch.intscratch [j] =
                ctx->scratch.intscratch [ctx->ga.PopSize-i-1];
        };
        break;
    }
    PGADebugExited ("PGASortPop");
}


/*U***************************************************************************
   PGAGetPopSize - Returns the population size

   Category: Generation

   Inputs:
      ctx - context variable

   Outputs:
      The population size

   Example:
      PGAContext *ctx;
      int popsize;
      :
      popsize = PGAGetPopSize(ctx);

***************************************************************************U*/
int PGAGetPopSize (PGAContext *ctx)
{
    PGADebugEntered("PGAGetPopSize");
    PGAFailIfNotSetUp("PGAGetPopSize");

    PGADebugExited("PGAGetPopSize");

    return(ctx->ga.PopSize);
}

/*U***************************************************************************
   PGAGetNumReplaceValue - Returns the maximum number of strings to replace
   each generation.

   Category: Generation

   Inputs:
      ctx - context variable

   Outputs:
      The maximum number number of strings to replace each generation

   Example:
      PGAContext *ctx;
      int numreplace;
      :
      numreplace = PGAGetNumReplaceValue(ctx);

***************************************************************************U*/
int PGAGetNumReplaceValue (PGAContext *ctx)
{
    PGADebugEntered("PGAGetNumReplaceValue");
    PGAFailIfNotSetUp("PGAGetNumReplaceValue");

    PGADebugExited("PGAGetNumReplaceValue");

    return(ctx->ga.NumReplace);
}

/*U***************************************************************************
   PGAGetPopReplaceType - returns the symbolic constant used to determine
   which strings to copy from the old population to the new population.

   Category: Generation

   Inputs:
      ctx - context variable

   Outputs:
      The symbolic constant of the replacement strategy.

   Example:
      PGAContext *ctx;
      int popreplace;
      :
      popreplace = PGAGetPopReplaceType(ctx);
      switch (popreplace) {
      case PGA_POPREPL_BEST:
          printf("Replacement Strategy = PGA_POPREPL_BEST\n");
          break;
      case PGA_POPREPL_RANDOM_REP:
          printf("Replacement Strategy = PGA_POPREPL_RANDOM_REP\n");
          break;
      case PGA_POPREPL_RANDOM_NOREP:
          printf("Replacement Strategy = PGA_POPREPL_RANDOM_NOREP\n");
          break;
      }

****************************************************************************U*/
int PGAGetPopReplaceType (PGAContext *ctx)
{
    PGADebugEntered("PGAGetPopReplaceType");
    PGAFailIfNotSetUp("PGAGetPopRelaceType");

    PGADebugExited("PGAGetPopReplaceType");

    return(ctx->ga.PopReplace);
}

/*U****************************************************************************
   PGASetRTRWindowSize - Set window size used for restricted tournament
   selection. The window size must be smaller than the population size.
   The default is min (n, N/20) where n is the string length and N is
   the population size.

   Category: Generation

   Inputs:
      ctx         - context variable
      windowsize  - size of the window for restricted tournament
                    selection

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetRTRWindowSize(ctx, windowsize);

****************************************************************************U*/
void PGASetRTRWindowSize( PGAContext *ctx, int windowsize)
{
    PGADebugEntered("PGASetRTRWindowSize");

    if (windowsize > ctx->ga.PopSize) {
        PGAError ( ctx,
                  "PGASetRTRWindowSize: Invalid value of windowsize:",
                   PGA_FATAL, PGA_INT, (void *) &windowsize);
    }

    ctx->ga.RTRWindowSize = windowsize;

    PGADebugExited("PGASetRTRWindowSize");
}

/*U***************************************************************************
   PGAGetRTRWindowSize - Returns the window size for restricted
   tournamen replacement.

   Category: Generation

   Inputs:
      ctx - context variable

   Outputs:
      The size of the window for restricted tournament selection

   Example:
      PGAContext *ctx;
      int windowsize;
      :
      windowsize = PGAGetRTRWindowSize(ctx);

***************************************************************************U*/
int PGAGetRTRWindowSize (PGAContext *ctx)
{
    PGADebugEntered("PGAGetRTRWindowSize");

    PGADebugExited("PGAGetRTRWindowSize");

    return(ctx->ga.RTRWindowSize);
}

/*U****************************************************************************
   PGAGetSortedPopIndex - returns a population string index from the array
   created by PGASortPop().

   Category: Generation

   Inputs:
       ctx      - context variable
       n        - specified which index element is to be returned.

   Output:
       A population string index from the array created by PGASortPop

   Example:
      Copy the five best strings from the old population into the new
      population.  The rest of the new population will be created by
      recombination, and is not shown.

      PGAContext *ctx;
      int i,j;
      :
      PGASetPopReplaceType(ctx,PGA_POPREPL_BEST)
      PGASortPop(ctx, PGA_OLDPOP);
      for ( i=0; i < 5; i++) {
          j = PGAGetSortedPopIndex(ctx, i);
      PGACopyIndividual (ctx, j, PGA_OLDPOP, i, PGA_NEWPOP);
      :

****************************************************************************U*/
int PGAGetSortedPopIndex ( PGAContext *ctx, int n )
{
     int temp = 0;

    PGADebugEntered("PGAGetSortedPopIndex");
     if (n >= 0 && n < ctx->ga.PopSize )
          temp = ctx->ga.sorted[n];
     else
          PGAError( ctx, "PGAGetSorted: Invalid value of n:",
                   PGA_FATAL, PGA_INT, (void *) &n );

    PGADebugExited("PGAGetSortedPopIndex");

     return (temp);
}

/*U****************************************************************************
   PGASetPopSize - Specifies the size of the genetic algorithm population.
   The default population size is 100.

   Category: Generation

   Inputs:
      ctx     - context variable
      popsize - the genetic algorithm population size to use

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetPopSize(ctx, 200);

****************************************************************************U*/
void PGASetPopSize (PGAContext *ctx, int popsize)
{

    PGADebugEntered("PGASetPopSize");
    PGAFailIfSetUp("PGASetPopSize");

    if (popsize < 1 || popsize % 2)
        PGAError( ctx, "PGASetPopSize: Invalid value of popsize:",
                  PGA_FATAL, PGA_INT, (void *) &popsize );
    else
        ctx->ga.PopSize = popsize;

    PGADebugExited("PGASetPopSize");
}



/*U****************************************************************************
   PGASetNumReplaceValue - specifies the number of new strings to create each
   generation.  The default is ten percent of the population size

   Category: Generation

   Inputs:
      ctx         - context variable
      pop_replace - the genetic algorithm population size to use

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetNumReplaceValue(ctx, 35);

****************************************************************************U*/
void PGASetNumReplaceValue( PGAContext *ctx, int pop_replace)
{
    PGADebugEntered("PGASetNumReplaceValue");

    if (pop_replace < 0)
      PGAError( ctx,
               "PGASetNumReplaceValue: Invalid value of pop_replace:",
                PGA_FATAL, PGA_INT, (void *) &pop_replace );
    else
        ctx->ga.NumReplace = pop_replace;

    PGADebugExited("PGASetNumReplaceValue");
}




/*U****************************************************************************
   PGASetPopReplaceType - Choose method of sorting strings to copy from old
   population to new population.  Valid choices are PGA_POPREPL_BEST,
   PGA_POPREPL_RANDOM_NOREP, or PGA_POPREPL_RANDOM_REP for copying the best
   strings, or  random string, with or without replacement, respectively,
   from the old population into the new population. Additional
   replacement types are PGA_POPREPL_RTR for restricted tournament
   replacement, PGA_POPREPL_PAIRWISE_BEST for pairwise comparison of
   each individual in the old/new population, and PGA_POPREPL_NSGA_II
   for multiobjective optimization using the Nondominated Sorting
   Genetic Algorithm (NSGA-II). The default is PGA_POPREPL_BEST.

   Category: Generation

   Inputs:
      ctx         - context variable
      pop_replace - symbolic constant to specify the population replacement
                    strategy

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGASetPopReplaceType(ctx, PGA_POPREPL_RANDOM_NOREP);

****************************************************************************U*/
void PGASetPopReplaceType( PGAContext *ctx, int pop_replace)
{
    PGADebugEntered("PGASetPopReplaceType");

    switch (pop_replace) {
    case PGA_POPREPL_BEST:
    case PGA_POPREPL_RANDOM_NOREP:
    case PGA_POPREPL_RANDOM_REP:
    case PGA_POPREPL_RTR:
    case PGA_POPREPL_PAIRWISE_BEST:
    case PGA_POPREPL_NSGA_II:
        ctx->ga.PopReplace = pop_replace;
        break;
    default:
        PGAError ( ctx,
                  "PGASetPopReplaceType: Invalid value of pop_replace:",
                   PGA_FATAL, PGA_INT, (void *) &pop_replace);
        break;
    }

    PGADebugExited("PGASetPopReplaceType");
}

/*U****************************************************************************
   PGARestrictedTournamentReplacement - Perform restricted tournament
   replacement: for each individual in PGA_NEWPOP we select a window of
   individuals from PGA_OLDPOP, find the one genetically most like the
   new candidate and replace the individual if the new candidate has
   better evalutation. Note that we may not use the fitness here:
   Fitness from two different populations are uncompareable!
   After this populations are swapped (echange of PGA_NEWPOP and
   PGA_OLDPOP) for further processing.

   Category: Generation

   Inputs:
      ctx         - context variable

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGARestrictedTournamentReplacement(ctx);

****************************************************************************U*/
void PGARestrictedTournamentReplacement (PGAContext *ctx)
{
    int i, j;
    int popsize = PGAGetPopSize(ctx);
    int numreplace = PGAGetNumReplaceValue(ctx);
    PGASampleState state;
    PGAIndividual *temp;
    int oldpop = PGA_OLDPOP;
    int newpop = PGA_NEWPOP;

    PGADebugEntered("PGARestrictedTournamentReplacement");
    for (i=popsize - numreplace; i<popsize; i++) {
        double dist = -1.0;
        int closest = 0;
        PGARandomSampleInit
            (ctx, &state, ctx->ga.RTRWindowSize, ctx->ga.PopSize);
        for (j=0; j<ctx->ga.RTRWindowSize; j++) {
            double d;
            int idx = PGARandomNextSample (&state);
            if (ctx->fops.GeneDistance) {
                d = (*ctx->fops.GeneDistance)
                    (&ctx, &idx, &oldpop, &i, &newpop);
            } else {
                d = (*ctx->cops.GeneDistance)
                    (ctx, idx, PGA_OLDPOP, i, PGA_NEWPOP);
            }
            if (dist < 0 || d < dist) {
                dist = d;
                closest = idx;
            }
        }

        /* If new population individual is better */
        if (PGAEvalCompare (ctx, i, PGA_NEWPOP, closest, PGA_OLDPOP) <= 0) {
            /* Copy i in PGA_NEWPOP to closest in PGA_OLDPOP */
            PGACopyIndividual (ctx, i, PGA_NEWPOP, closest, PGA_OLDPOP);
        }
    }
    /* Exchange old/newpop, will be done again by PGAUpdateGeneration,
     * we just make sure that the result of above replacements is now
     * newpop.
     */
    temp           = ctx->ga.oldpop;
    ctx->ga.oldpop = ctx->ga.newpop;
    ctx->ga.newpop = temp;
    PGADebugExited("PGARestrictedTournamentReplacement");
}
/*U****************************************************************************

   PGAPairwiseBestReplacement - Perform pairwise best replacement:
   Compare individuals with same index in PGA_OLDPOP and PGA_NEWPOP and
   select the one with better evalutation. Note that we may not use the
   fitness here: Fitness from two different populations are
   uncompareable!
   This replacement strategy is used in evolutionary algorithms that
   modify a single individual and replace the parent if the offspring is
   better. A popular example is Differential Evolution (DE).
   After this populations are swapped (echange of PGA_NEWPOP and
   PGA_OLDPOP) for further processing.

   Category: Generation

   Inputs:
      ctx         - context variable

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGAPairwiseBestReplacement(ctx);

****************************************************************************U*/
void PGAPairwiseBestReplacement (PGAContext *ctx)
{
    int i;
    int popsize = PGAGetPopSize(ctx);
    int numreplace = PGAGetNumReplaceValue(ctx);
    PGAIndividual *temp;

    PGADebugEntered("PGAPairwiseBestReplacement");
    for (i=popsize - numreplace; i<popsize; i++) {
        /* Note the '<=' comparison, differential evolution can walk across
         * areas with equal evaluation this way
         */
        if (PGAEvalCompare (ctx, i, PGA_NEWPOP, i, PGA_OLDPOP) <= 0) {
            /* Copy i in PGA_NEWPOP to i in PGA_OLDPOP */
            PGACopyIndividual (ctx, i, PGA_NEWPOP, i, PGA_OLDPOP);
        }
    }
    /* Exchange old/newpop, will be done again by PGAUpdateGeneration,
     * we just make sure that the result of above replacements is now
     * newpop.
     */
    temp           = ctx->ga.oldpop;
    ctx->ga.oldpop = ctx->ga.newpop;
    ctx->ga.newpop = temp;
    PGADebugExited("PGAPairwiseBestReplacement");
}

/*U****************************************************************************

   PGA_NSGA_II_Replacement - Perform NSGA-II Replacement
   First compute a dominance matrix of N x N bits. The rows are the
   dominated-by relation. We loop over all n^1 pairs of individuals and
   fill the matrix. Initit all ranks with -1.
   Then starting with rank0:
   - Get all rows of the matrix which are 0 and where the individual has
     no rank yet: These are the currently non-dominated individuals,
     assign the current rank
   - Loop over all individuals with the current rank and remove their
     bits from the dominance matrix
   - Increment the rank counter

   Category: Generation

   Inputs:
      ctx         - context variable

   Outputs:
      None

   Example:
      PGAContext *ctx;
      :
      PGA_NSGA_II_Replacement(ctx);

****************************************************************************U*/
/* Helper functions for PGA_NSGA_II_Replacement */

#define NONNEGEVAL(v, is_ev) ((is_ev) ? (v) : (((v) < 0) ? 0 : (v)))
#define GETEVAL(ind, fidx, is_ev)                       \
    ( ((fidx) == 0)                                     \
    ? (ind)->evalue                                     \
    : NONNEGEVAL ((ind)->auxeval [(fidx) - 1], (is_ev)) \
    )

static int crowdsort_cmp (const void *a1, const void *a2)
{
    PGAIndividual **i1 = (void *)a1;
    PGAIndividual **i2 = (void *)a2;
    int is_ev = INDGetAuxTotal (*i1) ? 0 : 1;
    double e1 = GETEVAL(*i1, (*i1)->funcidx, is_ev);
    double e2 = GETEVAL(*i2, (*i2)->funcidx, is_ev);
    return CMP(e1, e2);
}

static int nondom_cmp (const void *a1, const void *a2)
{
    const PGAIndividual * const *i1 = a1;
    const PGAIndividual * const *i2 = a2;
    if ((*i1)->rank < (*i2)->rank) {
        return -1;
    }
    if ((*i1)->rank > (*i2)->rank) {
        return 1;
    }
    return CMP((*i2)->crowding, (*i1)->crowding);
}

/* Compute crowding distance over the given individuals */
static void crowding (PGAContext *ctx, PGAIndividual **start, int n, int rank)
{
    int i, k;
    int nrank = 0;
    int is_ev = INDGetAuxTotal (*start) ? 0 : 1;
    PGAIndividual *crowd [n];
    int nc = ctx->ga.NumConstraint;
    int na = ctx->ga.NumAuxEval;
    int base = is_ev ? 0 : (na - nc);
    int nfun = is_ev ? (na - nc + 1) : nc;
    double f_min [nfun], f_max [nfun];

    for (i=0; i<n; i++) {
        PGAIndividual *ind = start [i];
        if (ind->rank == rank) {
            for (k=0; k<nfun; k++) {
                double e = GETEVAL (ind, k + base, is_ev);
                if (nrank == 0 || f_min [k] > e) {
                    f_min [k] = e;
                }
                if (nrank == 0 || f_max [k] < e) {
                    f_max [k] = e;
                }
            }
            ind->crowding = 0;
            crowd [nrank++] = ind;
        }
    }
    assert (nrank > 0);
    for (k=0; k<nfun; k++) {
        double norm = f_max [k] - f_min [k];
        for (i=0; i<nrank; i++) {
            (crowd [i])->funcidx = k + base;
        }
        qsort (crowd, nrank, sizeof (crowd [0]), crowdsort_cmp);
        (crowd [0])->crowding = DBL_MAX;
        (crowd [nrank-1])->crowding = DBL_MAX;
        for (i=1; i<nrank-1; i++) {
            if ((crowd [i])->crowding != DBL_MAX) {
                (crowd [i])->crowding +=
                    ( GETEVAL(crowd [i+1], k + base, is_ev)
                    - GETEVAL(crowd [i-1], k + base, is_ev)
                    ) / norm;
            }
        }
    }
}

/* Dominance computation, return the maximum rank given or UINT_MAX if
 * goal was reached exactly (in which case no crowding is necessary)
 */
static int ranking (PGAContext *ctx, PGAIndividual **start, int n, int goal)
{
    int i, j, k;
    int is_ev = INDGetAuxTotal (*start) ? 0 : 1;
    int rank;
    int nranked = 0;
    int nc = ctx->ga.NumConstraint;
    int na = ctx->ga.NumAuxEval;
    int base = is_ev ? 0 : (na - nc);
    int nfun = is_ev ? (na - nc + 1) : nc;
    int nintbits = sizeof (int) * 8;
    int intsforn = (n + nintbits - 1) / nintbits;
    int (*dominance)[n][intsforn] =
        (int (*)[n][intsforn])(ctx->scratch.dominance);

    for (i=0; i<n; i++) {
        (start [i])->rank = UINT_MAX;
        for (j=0; j<intsforn; j++) {
            (*dominance) [i][j] = 0;
        }
    }
    for (i=0; i<n; i++) {
        int iidx   = i / nintbits;
        int ishift = 1 << (i % nintbits);
        for (j=i+1; j<n; j++) {
            int jidx   = j / nintbits;
            int jshift = 1 << (j % nintbits);
            int cmp = 0;
            for (k=0; k<nfun; k++) {
                double e1, e2;
                int ncmp;
                e1 = GETEVAL (start [i], k + base, is_ev);
                e2 = GETEVAL (start [j], k + base, is_ev);
                if (is_ev && ctx->ga.optdir == PGA_MAXIMIZE) {
                    ncmp = CMP (e2, e1);
                } else {
                    ncmp = CMP (e1, e2);
                }
                if (cmp && ncmp && ncmp != cmp) {
                    break;
                }
                cmp = ncmp;
            }
            /* Non-dominated? */
            if (!cmp || k<nfun) {
                continue;
            }
            /* j dominated by i */
            if (cmp < 0) {
                (*dominance) [j][iidx] |= ishift;
            /* i dominated by j */
            } else {
                (*dominance) [i][jidx] |= jshift;
            }
        }
    }
    /* Now repeatedly loop over individuals and establish rank */
    nranked = 0;
    for (rank=0; rank<n; rank++) {
        for (i=0; i<n; i++) {
            if ((*(start+i))->rank != UINT_MAX) {
                continue;
            }
            for (j=0; j<intsforn; j++) {
                if ((*dominance) [i][j]) {
                    break;
                }
            }
            /* Non-dominated in this rank */
            if (j==intsforn) {
                (*(start+i))->rank = rank;
                nranked ++;
            }
        }
        /* Need to rank only goal individuals */
        if (nranked >= goal) {
            break;
        }
        /* Remove dominance bits for this rank */
        for (i=0; i<n; i++) {
            int iidx   = i / nintbits;
            int ishift = 1 << (i % nintbits);
            if ((*(start+i))->rank != rank) {
                continue;
            }
            for (j=0; j<n; j++) {
                (*dominance) [j][iidx] &= ~ishift;
            }
        }
    }
    /* No need for crowding computation if we hit goal exactly */
    if (nranked == goal) {
        return UINT_MAX;
    }
    return rank;
}

void PGA_NSGA_II_Replacement (PGAContext *ctx)
{
    int i;
    int n_unc_ind, n_con_ind;
    int popsize = PGAGetPopSize (ctx);
    int numreplace = PGAGetNumReplaceValue (ctx);
    PGAIndividual *all_individuals [popsize + numreplace];
    PGAIndividual **constrained = all_individuals + popsize + numreplace;
    PGAIndividual **unconstrained = all_individuals;
    PGAIndividual *oldpop = ctx->ga.oldpop;
    PGAIndividual *newpop = ctx->ga.newpop;
    PGAIndividual *temp;

    PGADebugEntered("PGA_NSGA_II_Replacement");

    /* We keep two pointers into the all_individuals array. One with
     * constrained individuals starts from the end. The other with
     * unconstrained individuals starts from the start.
     */

    /* First loop over all old individuals and put them into the
     * all_individuals array.
     */
    for (i=0; i<popsize; i++) {
        if (INDGetAuxTotal (oldpop + i) > 0) {
            constrained--;
            *constrained = oldpop + i;
        } else {
            *unconstrained = oldpop + i;
            unconstrained++;
        }
    }
    /* Now put all the new individuals into the same array */
    for (i=popsize - numreplace; i<popsize; i++) {
        if (INDGetAuxTotal (newpop + i) > 0) {
            constrained--;
            *constrained = newpop + i;
        } else {
            *unconstrained = newpop + i;
            unconstrained++;
        }
    }
    /* When array has been filled both pointers point to the same
     * individual somewhere in the middle
     */
    assert (constrained == unconstrained);
    n_unc_ind = unconstrained - all_individuals;
    n_con_ind = popsize + numreplace - n_unc_ind;

    /* First perform non-dominated sorting on unconstrained individuals
     * Normal sorting if only one eval function */
    if (ctx->ga.NumConstraint == ctx->ga.NumAuxEval) {
        qsort
            ( all_individuals
            , n_unc_ind
            , sizeof (all_individuals [0])
            , PGAEvalSortHelper
            );
    } else {
        int rank;
        rank = ranking (ctx, all_individuals, n_unc_ind, popsize);
        if (n_unc_ind >= popsize && rank != UINT_MAX) {
            crowding (ctx, all_individuals, n_unc_ind, rank);
        }
        qsort \
            ( all_individuals
            , n_unc_ind
            , sizeof (all_individuals [0])
            , nondom_cmp
            );
    }

    /* Only sort constrained individuals if some of them go into next
     * generation, i.e. if number of unconstrained individuals is below
     * popsize.
     */
    if (n_unc_ind < popsize) {
        assert (ctx->ga.NumConstraint);
        /* Normal sorting if only one constraint function or
         * SumConstraints is PGA_TRUE
         */
        if (ctx->ga.SumConstraints || ctx->ga.NumConstraint == 1) {
            qsort
                ( all_individuals + n_unc_ind
                , n_con_ind
                , sizeof (all_individuals [0])
                , PGAEvalSortHelper
                );
        } else {
            int rank;
            rank = ranking
                ( ctx
                , all_individuals + n_unc_ind
                , n_con_ind
                , popsize - n_unc_ind
                );
            if (rank != UINT_MAX) {
                crowding (ctx, all_individuals + n_unc_ind, n_con_ind, rank);
            }
            qsort \
                ( all_individuals + n_unc_ind
                , n_con_ind
                , sizeof (all_individuals [0])
                , nondom_cmp
                );
        }
    }
    /* Now we can put all better individuals sorted first into the
     * places of the individuals sorted last, if the source is in
     * PGA_NEWPOP and the target is in PGA_OLDPOP.
     */

    constrained = all_individuals;
    unconstrained = all_individuals + popsize + numreplace;
    while (constrained - all_individuals < popsize) {
        while ((*constrained)->pop == ctx->ga.oldpop) {
            constrained++;
        }
        while ((*(unconstrained-1))->pop == ctx->ga.newpop) {
            unconstrained--;
        }
        if (constrained - all_individuals < popsize) {
            assert (constrained < unconstrained);
            INDCopyIndividual (*constrained, *(unconstrained-1));
            constrained++;
            unconstrained--;
        }
    }
    assert (constrained - all_individuals <= popsize + numreplace);
    assert (constrained - all_individuals >= popsize);

    /* Exchange old/newpop, will be done again by PGAUpdateGeneration,
     * we just make sure that the result of above replacements is now
     * newpop.
     */
    temp           = ctx->ga.oldpop;
    ctx->ga.oldpop = ctx->ga.newpop;
    ctx->ga.newpop = temp;
    PGADebugExited("PGA_NSGA_II_Replacement");
}

