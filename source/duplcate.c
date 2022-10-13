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
*     FILE: duplicate.c: This file contains the routines that have to do with
*                        testing for duplicate strings
*
*     Authors: David M. Levine, Philip L. Hallstrom, David M. Noelle,
*              Brian P. Walenz
*****************************************************************************/

#include "pgapack.h"

/*U****************************************************************************
  PGADuplicate - determines if a specified string is a duplicate of one
  already in an existing population

  Category: Generation

  Inputs:
    ctx  - context variable
    p    - string index
    pop1 - symbolic constant of the population containing string p
    pop2 - symbolic constant of the (possibly partial) population containing
           strings to compare string p against

  Outputs:
    Returns PGA_TRUE if PGAGetNoDuplicatesFlag() returns PGA_TRUE and
    string p in population pop1 is a duplicate of at least one strings
    already inserted into population pop2.  Otherwise returns PGA_FALSE

  Example:
    Check the current to-be-inserted string if it is a copy of any of
    the strings in PGA_NEWPOP. Note that the check relies on all
    individuals in PGA_NEWPOP to also be inserted into the duplicate
    hash, see PGAHashIndividual.

    PGAContext *ctx;
    int p;
    :

    while (PGADuplicate (ctx, p, PGA_NEWNEW, PGA_NEWPOP)) {
        PGAChange (ctx, p, PGA_NEWPOP);
    }

****************************************************************************U*/
int PGADuplicate (PGAContext *ctx, int p, int pop1, int pop2)
{
    int RetVal = PGA_FALSE;
    PGADebugEntered ("PGADuplicate");

    if (ctx->ga.NoDuplicates) {
        int p2, fp;
        size_t idx = PGAIndividualHashIndex (ctx, p, pop1);
        PGAIndividual *ind = NULL;
    
        assert (pop2 == PGA_NEWPOP);
        for (ind = ctx->scratch.hashed [idx]; ind; ind = ind->next_hash) {
            p2 = ind - ind->pop;

            if (ctx->fops.Duplicate) {
                int p2f = p2 + 1;
                fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
                if ((*ctx->fops.Duplicate)(&ctx, &fp, &pop1, &p2f, &pop2)) {
                    return PGA_TRUE;
                }
            } else {
                if ((*ctx->cops.Duplicate)(ctx, p, pop1, p2, pop2)) {
                    return PGA_TRUE;
                }
            }
        }
    }
    PGADebugExited ("PGADuplicate");
    return RetVal;
}


/*U****************************************************************************
  PGAChange - Repeatedly apply mutation to a string (with an increasing
  mutation rate) until one or more mutations have occurred.  This routine is
  usually used with PGADuplicate to modify a duplicate string.  It is not
  intended to replace PGAMutation

  Category: Generation

  Inputs:
    ctx  - context variable
    p    - string index
    pop  - symbolic constant of the population containing string p

  Outputs:
    Mutates string p in population pop via side effect.

  Example:
    Check the current to-be-inserted string if it is a copy of any of
    the strings in PGA_NEWPOP. Note that the check relies on all
    individuals in PGA_NEWPOP to also be inserted into the duplicate
    hash, see PGAHashIndividual.

    PGAContext *ctx;
    int p;
    :

    while (PGADuplicate (ctx, p, PGA_NEWNEW, PGA_NEWPOP)) {
        PGAChange (ctx, p, PGA_NEWPOP);
    }


****************************************************************************U*/
void PGAChange (PGAContext *ctx, int p, int pop)
{
    int    changed = PGA_FALSE;
    int    fp, nflips;
    double mr;

    PGADebugEntered ("PGAChange");

    mr = ctx->ga.MutationProb;
    if (mr == 0) {
        mr = 1.0 / ctx->ga.StringLen;
    }

    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGAChange", " mr = "
        , PGA_DOUBLE, (void *) &mr
        );

    while (!changed) {
	if (ctx->fops.Mutation) {
	    fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
            nflips = (*ctx->fops.Mutation)(&ctx, &fp, &pop, &mr);
	} else {
	    nflips = (*ctx->cops.Mutation)( ctx, p, pop, mr );
	}

        if (nflips > 0) {
            changed = PGA_TRUE;
        } else {
            if (mr >= 1) {
                break;
            }
            mr = 1.1 * mr;
            if (mr > 1) {
                mr = 1;
            }
        }
    }

    if (changed) {
	PGASetEvaluationUpToDateFlag (ctx, p, pop, PGA_FALSE);
    } else {
	PGAError (ctx, "Could not change string:", PGA_WARNING, PGA_VOID, NULL);
	PGAPrintString (ctx, stderr, p, pop);
    }

    PGADebugExited ("PGAChange");
}

/*U****************************************************************************
  PGAHashIndividual - Compute Hash for this individual and insert into
  hash-table

  Category: Generation

  Inputs:
     ctx  - context variable
     p    - string index
     pop  - symbolic constant of the population containing string p

  Outputs:
     Computes hash of given individual and inserts it into hash table.

****************************************************************************U*/
void PGAHashIndividual (PGAContext *ctx, int p, int pop)
{
    if (ctx->ga.NoDuplicates) {
        PGAIndividual *ind = PGAGetIndividual (ctx, p, pop);
        size_t idx = PGAIndividualHashIndex (ctx, p, pop);

        if (ctx->scratch.hashed [idx]) {
            assert (ind->pop == ctx->scratch.hashed [idx]->pop);
            ind->next_hash = ctx->scratch.hashed [idx];
        }
        ctx->scratch.hashed [idx] = ind;
    }
}

/*U****************************************************************************
  PGASetNoDuplicatesFlag - A boolean flag to indicate if duplicate strings are
  allowed in the population. Valid choices are PGA_TRUE and PGA_FALSE.  The
  default is PGA_FALSE -- allow duplicates.

  Category: Generation

  Inputs:
    ctx  - context variable
    flag - PGA_TRUE or PGA_FALSE

  Outputs:
    None

  Example:
    Set the NoDuplicates flag to require that all strings are unique.

    PGAContext *ctx;
    :
    PGASetNoDuplicatesFlag (ctx,PGA_TRUE);

****************************************************************************U*/
void PGASetNoDuplicatesFlag (PGAContext *ctx, int no_dup)
{
    PGADebugEntered ("PGASetNoDuplicatesFlag");

    switch (no_dup) {
        case PGA_TRUE:
        case PGA_FALSE:
            ctx->ga.NoDuplicates = no_dup;
            break;
        default:
            PGAError
                ( ctx, "PGASetNoDuplicatesFlag: Invalid value of no_dup:"
                , PGA_FATAL, PGA_INT, (void *) &no_dup
                );
            break;
    }

    PGADebugExited ("PGASetNoDuplicatesFlag");
}

/*U***************************************************************************
  PGAGetNoDuplicatesFlag - Returns PGA_TRUE if duplicates are not allowed,
  else returns PGA_FALSE.

  Category: Generation

  Inputs:
    ctx - context variable

  Outputs:
    The value of the NoDuplicates flag.

  Example:
    PGAContext *ctx;
    int nodups;
    :
    nodups = PGAGetNoDuplicatesFlag (ctx);
    switch (nodups) {
    case PGA_TRUE:
        printf ("Duplicate strings not allowed in population\n");
        break;
    case PGA_FALSE:
        printf ("Duplicate strings allowed in population\n");
        break;
    }

***************************************************************************U*/
int PGAGetNoDuplicatesFlag (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetNoDuplicatesFlag");

    PGAFailIfNotSetUp ("PGAGetNoDuplicatesFlag");

    PGADebugExited ("PGAGetNoDuplicatesFlag");

    return (ctx->ga.NoDuplicates);
}
