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
 *    \file
 * This file contains the routines specific to the character datatype.
 * \authors Authors:
 *          David M. Levine, Philip L. Hallstrom, David M. Noelle,
 *          Brian P. Walenz, Ralf Schlatterbeck
 *****************************************************************************/

#include <pgapack.h>

/*!****************************************************************************
    \brief Sets the value of an allele in a string of the character data
           type.
    \ingroup allele
    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \param   i    allele index
    \param   val  character value to set the allele to
    \return  The allele is changed by side-effect.

    \rst

    Example
    -------

    Copies the alleles from member ``p`` in :c:macro:`PGA_OLDPOP` to
    member ``q`` in :c:macro:`PGA_NEWPOP`.
    Assumes the strings are of the same length.

    .. code-block:: c

        PGAContext *ctx;
        int p, q, i;
        int l;

        ...
        l = PGAGetStringLength (ctx);
        for (i=0; i<l; i++) {
            char a = PGAGetCharacterAllele (ctx, p, PGA_OLDPOP, i);
            PGASetCharacterAllele (ctx, q, PGA_NEWPOP, i, a);
        }

    \endrst

******************************************************************************/
void PGASetCharacterAllele (PGAContext *ctx, int p, int pop, int i, char val)
{
    PGAIndividual *ind;

    PGADebugEntered  ("PGASetCharacterAllele");
    PGACheckDataType ("PGASetCharacterAllele", PGA_DATATYPE_CHARACTER);

    ind = PGAGetIndividual (ctx, p, pop);
    ((PGACharacter *)ind->chrom)[i] = val;

    PGADebugExited ("PGASetCharacterAllele");
}

/*!****************************************************************************
    \brief Return the value of character allele in a string of the
           character data type.
    \ingroup allele
    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \param   i    allele index
    \return  The value of allele i in string p

    \rst

    Example
    -------

    Copies the alleles from member ``p`` in :c:macro:`PGA_OLDPOP` to
    member ``q`` in :c:macro:`PGA_NEWPOP`.
    Assumes the strings are of the same length.

    .. code-block:: c

        PGAContext *ctx;
        int p, q, i;
        int l;

        ...
        l = PGAGetStringLength (ctx);
        for (i=0; i<l; i++) {
            char a = PGAGetCharacterAllele (ctx, p, PGA_OLDPOP, i);
            PGASetCharacterAllele (ctx, q, PGA_NEWPOP, i, a);
        }

    \endrst

******************************************************************************/
char PGAGetCharacterAllele (PGAContext *ctx, int p, int pop, int i)
{
     PGAIndividual *ind;

    PGADebugEntered("PGAGetCharacterAllele");
     PGACheckDataType("PGAGetCharacterAllele", PGA_DATATYPE_CHARACTER);

     ind = PGAGetIndividual ( ctx, p, pop );

    PGADebugExited("PGAGetCharacterAllele");

     return (((PGACharacter *)ind->chrom)[i]);
}


/*!****************************************************************************
    \brief Sets a flag to specify whether the character strings will be
           exclusively lowercase, exclusively uppercase, or a mixure of
           both cases.
    \ingroup init

    \param   ctx    context variable
    \param   value  symbolic constant specifying which case
    \return  None

    \rst

    Description
    -----------

    Legal flags are :c:macro:`PGA_CINIT_UPPER`,
    :c:macro:`PGA_CINIT_LOWER`, and
    :c:macro:`PGA_CINIT_MIXED`.  Default is :c:macro:`PGA_CINIT_LOWER`.
    See :ref:`group:const-randinit` for the constants and section
    :ref:`sec:initialization` of the user guide for details.

    Example
    -------

    Set program to generate exclusively uppercase letters.

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetCharacterInitType (ctx, PGA_CINIT_UPPER);

    \endrst

******************************************************************************/
void PGASetCharacterInitType (PGAContext *ctx, int value)
{
    PGADebugEntered  ("PGASetCharacterInitType");
    PGACheckDataType ("PGASetCharacterInitType", PGA_DATATYPE_CHARACTER);

    switch (value)
    {
    case PGA_CINIT_UPPER:
    case PGA_CINIT_LOWER:
    case PGA_CINIT_MIXED:
         ctx->init.CharacterType = value;
         break;
    default:
         PGAError
            ( ctx, "PGASetCharacterInitType: Invalid case type:"
            , PGA_FATAL, PGA_INT, (void *)&value
            );
         break;
    }

    PGADebugExited ("PGASetCharacterInitType");
}

/*!****************************************************************************
    \brief Allocate memory for a string of type character.
    \ingroup internal
    \param   ctx       context variable
    \param   p         string index
    \param   pop       symbolic constant of the population string p is in
    \param   initflag  A true/false flag used in conjunction with
                       ctx->ga.RandomInit to initialize the string
                       either randomly or set to zero
    \return  Member p in population pop is allocated and initialized.

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the create
    string user function for the char datatype by default.

    Example
    -------

    Allocates memory and assigns the address of the allocated memory to
    the string field ``ind->chrom`` of the individual.  Additionally, the
    string is initialized to zero.

    .. code-block:: c

        PGAContext *ctx;
        int p;

        ...
        PGACharacterCreateString (ctx, p, PGA_NEWPOP, PGA_FALSE);

    \endrst

******************************************************************************/
void PGACharacterCreateString (PGAContext *ctx, int p, int pop, int initflag)
{
    int i, fp;
    PGACharacter *c;
    PGAIndividual *new = PGAGetIndividual(ctx, p, pop);

    PGADebugEntered ("PGACharacterCreateString");

    new->chrom = (void *)malloc (ctx->ga.StringLen * sizeof(PGACharacter));
    if (new->chrom == NULL) {
        PGAError
            ( ctx, "PGACharacterCreateString: No room to allocate new->chrom"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    c = (PGACharacter *)new->chrom;
    if (initflag) {
        if (ctx->fops.InitString) {
            fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
            (*ctx->fops.InitString)(&ctx, &fp, &pop);
        } else {
            (*ctx->cops.InitString)(ctx, p, pop);
        }
    } else {
        for (i=0; i<ctx->ga.StringLen; i++) {
            c[i] = 0;
        }
    }

    PGADebugExited("PGACharacterCreateString");
}

/*!****************************************************************************
    \brief Randomly mutates a character-valued gene with a specified
           probability.
    \ingroup internal

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population string p is in
    \param   mr   probability of mutating an character-valued gene
    \return  Returns the number of mutations

    \rst

    Description
    -----------

    This routine is called from :c:func:`PGAMutate`.

    Note that this function is set in :c:func:`PGASetUp` as the mutation
    user function for the char datatype by default.


    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p;
       int NumMutations;

       ...
       NumMutations = PGACharacterMutation (ctx, p, PGA_NEWPOP, 0.01);

    \endrst
******************************************************************************/
int PGACharacterMutation (PGAContext *ctx, int p, int pop, double mr)
{
    PGACharacter *c;
    int i, j;
    int count = 0;

    PGADebugEntered("PGACharacterMutation");

    c = (PGACharacter *)PGAGetIndividual(ctx, p, pop)->chrom;
    for (i=0; i<ctx->ga.StringLen; i++) {
        /* randomly choose an allele */
        if (PGARandomFlip(ctx, mr)) {
             switch (ctx->init.CharacterType) {
             case PGA_CINIT_LOWER:
                  c [i] = PGARandomInterval (ctx, 'a', 'z');
                  break;
             case PGA_CINIT_UPPER:
                  c [i] = PGARandomInterval (ctx, 'A', 'Z');
                  break;
             case PGA_CINIT_MIXED:
                  j = PGARandomInterval (ctx, 0, 51);
                  if (j < 26) {
                       c [i] = 'A' + j;
                  } else {
                       c [i] = 'a' + j - 26;
                  }
                  break;
             }
             count++;
        }
    }

    PGADebugExited ("PGACharacterMutation");
    return (count);
}

/*!****************************************************************************
    \brief Perform one-point crossover on two parent strings producing
           two children via side-effect.
    \ingroup internal

    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing
                   string p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population to contain
                   string c1 and c2
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the char datatype when selecting
    one-point crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``,
    producing children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGACharacterOneptCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGACharacterOneptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGACharacter *parent1, *parent2, *child1, *child2;
    int i, xsite;

    PGADebugEntered ("PGACharacterOneptCrossover");

    parent1 = (PGACharacter *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent2 = (PGACharacter *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child1  = (PGACharacter *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child2  = (PGACharacter *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    xsite = PGARandomInterval (ctx, 1,ctx->ga.StringLen-1);

    for(i=0; i<xsite; i++) {
        child1 [i] = parent1 [i];
        child2 [i] = parent2 [i];
    }

    for(i=xsite; i<ctx->ga.StringLen; i++) {
        child1 [i] = parent2 [i];
        child2 [i] = parent1 [i];
    }

    PGADebugExited ("PGACharacterOneptCrossover");
}

/*!****************************************************************************
    \brief Perform two-point crossover on two parent strings producing
           two children via side-effect.
    \ingroup internal

    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing
                   string p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population to contain
                   string c1 and c2
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the char datatype when selecting
    two-point crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGACharacterTwoptCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGACharacterTwoptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGACharacter *parent1, *parent2, *child1, *child2;
    int i, temp, xsite1, xsite2;

    PGADebugEntered ("PGACharacterTwoptCrossover");

    parent1 = (PGACharacter *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent2 = (PGACharacter *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child1  = (PGACharacter *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child2  = (PGACharacter *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    /* pick two cross sites such that xsite2 > xsite1 */
    xsite1 = PGARandomInterval (ctx, 1,ctx->ga.StringLen-1);
    xsite2 = xsite1;
    while (xsite2 == xsite1) {
        xsite2 = PGARandomInterval (ctx, 1,ctx->ga.StringLen-1);
    }
    if (xsite1 > xsite2) {
        temp   = xsite1;
        xsite1 = xsite2;
        xsite2 = temp;
    }
    for (i=0; i<xsite1; i++) {
        child1 [i] = parent1 [i];
        child2 [i] = parent2 [i];
    }
    for (i=xsite1; i<xsite2; i++) {
        child1 [i] = parent2 [i];
        child2 [i] = parent1 [i];
    }
    for (i=xsite2; i<ctx->ga.StringLen; i++) {
        child1 [i] = parent1 [i];
        child2 [i] = parent2 [i];
    }

    PGADebugExited ("PGACharacterTwoptCrossover");
}


/*!****************************************************************************
    \brief Perform uniform crossover on two parent strings producing two
           children via side-effect
    \ingroup internal

    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing
                   string p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population to contain
                   string c1 and c2
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the char datatype when selecting
    uniform crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGACharacterUniformCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGACharacterUniformCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGACharacter *parent1, *parent2, *child1, *child2;
    int i;

    PGADebugEntered ("PGACharacterUniformCrossover");

    parent1 = (PGACharacter *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    parent2 = (PGACharacter *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    child1  = (PGACharacter *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    child2  = (PGACharacter *)PGAGetIndividual (ctx, c2, pop2)->chrom;

    for (i=0; i<ctx->ga.StringLen; i++) {
        if (parent1 [i] == parent2 [i]) {
             child1 [i] = parent1 [i];
             child2 [i] = parent2 [i];
        } else if (PGARandomFlip (ctx, ctx->ga.UniformCrossProb)) {
             child1 [i] = parent1 [i];
             child2 [i] = parent2 [i];
        } else {
             child1 [i] = parent2 [i];
             child2 [i] = parent1 [i];
        }
    }

    PGADebugExited ("PGACharacterUniformCrossover");
}

/*!****************************************************************************
    \brief Write a character-valued string to a file.
    \ingroup internal

    \param   ctx  context variable
    \param   fp   file pointer to file to write the string to
    \param   p    index of the string to write out
    \param   pop  symbolic constant of the population string p is in
    \return  None

    \rst

    Example
    -------

    Write string s to stdout.

    .. code-block:: c

       PGAContext *ctx;
       int p;

       ...
       PGACharacterPrintString (ctx, stdout, p, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGACharacterPrintString (PGAContext *ctx, FILE *fp, int p, int pop)
{
    PGACharacter *c;
    int           i, pos, len;

    PGADebugEntered ("PGACharacterPrintString");

    c = (PGACharacter *)PGAGetIndividual (ctx, p, pop)->chrom;
    len = PGAGetStringLength (ctx);

    pos = 0;
    while (len > 0) {
        fprintf (fp, "#%5d: [", pos);
        for (i=0; i<50 && len>0; i++,len--,c++) {
            fputc (*c, fp);
        }
        pos+=50;
        fprintf (fp, "]\n");
    }
    fprintf (fp, "\n");

    PGADebugExited ("PGACharacterPrintString");
}

/*!****************************************************************************
    \brief Copy one character-valued string to another, assumes the
           strings are of the same length.
    \ingroup internal

    \param   ctx   context variable
    \param   p1    string to copy
    \param   pop1  symbolic constant of population containing string p1
    \param   p2    string to copy p1 to
    \param   pop2  symbolic constant of population containing string p2
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the copy
    string user function for the char datatype by default.

    Example
    -------

    Copy character string ``x`` to ``y`` (both are implicitly assumed to
    be the same length)

    .. code-block:: c

       PGAContext *ctx;
       int x, y;

       ...
       PGACharacterCopyString (ctx, x, PGA_OLDPOP, y, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGACharacterCopyString
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    void *source, *dest;
    int len;

    PGADebugEntered ("PGACharacterCopyString");

    source = PGAGetIndividual (ctx, p1, pop1)->chrom;
    dest   = PGAGetIndividual (ctx, p2, pop2)->chrom;
    len    = PGAGetStringLength (ctx);
    memcpy (dest, source, len * sizeof (PGACharacter));

    PGADebugExited ("PGACharacterCopyString");
}

/*!****************************************************************************
    \brief Return true if string p1 in pop1 is a duplicate of string p2
           in pop2, else returns false, assumes the strings are the same
           length.
    \ingroup internal

    \param   ctx   context variable
    \param   p1    string index of the first string to compare
    \param   pop1  symbolic constant of the population string p1 is in
    \param   p2    string index of the second string to compare
    \param   pop2  symbolic constant of the population string p2 is in
    \return  Returns true if strings are duplicates

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    duplicate checking user function for the char datatype by default.

    Example
    -------

    Compare string ``x`` with ``y`` to see if they are duplicates

    .. code-block:: c

       PGAContext *ctx;
       int x, y;

       ...
       if (PGACharacterDuplicate (ctx, x, PGA_NEWPOP, y, PGA_NEWPOP)) {
           printf ("strings are duplicates\n");
       }

    \endrst

******************************************************************************/
int PGACharacterDuplicate (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    void *a, *b;
    int len;

    PGADebugEntered ("PGACharacterDuplicate");

    a = PGAGetIndividual (ctx, p1, pop1)->chrom;
    b = PGAGetIndividual (ctx, p2, pop2)->chrom;
    len = PGAGetStringLength (ctx);

    PGADebugExited ("PGACharacterDuplicate");

    return (!memcmp (a, b, len * sizeof (PGACharacter)));
}

/*!****************************************************************************
    \brief Return hash value of given gene.
    \ingroup internal

    \param   ctx   context variable
    \param   p     string index of the string to hash
    \param   pop   symbolic constant of the population string p is in
    \return  Hash value for string

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    hash user function for the char datatype by default.

    \endrst

******************************************************************************/
PGAHash PGACharacterHash (PGAContext *ctx, int p, int pop)
{
    void *a = PGAGetIndividual(ctx, p, pop)->chrom;
    PGAHash hash = PGAUtilHash
        (a, sizeof (PGACharacter) * ctx->ga.StringLen, PGA_INITIAL_HASH);
    return hash;
}

/*!****************************************************************************
    \brief Randomly initialize a string of type character.
    \ingroup internal
    \param   ctx  context variable
    \param   p    index of string to randomly initialize
    \param   pop  symbolic constant of the population string p is in
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    init string user function for the char datatype by default.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p;

       ...
       PGACharacterInitString (ctx, p, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGACharacterInitString (PGAContext *ctx, int p, int pop)
{
    int len, i, j;
    PGACharacter *c;

    PGADebugEntered ("PGACharacterInitString");

    len = ctx->ga.StringLen;
    c = (PGACharacter *)PGAGetIndividual (ctx, p, pop)->chrom;
    switch (ctx->init.CharacterType) {
    case PGA_CINIT_LOWER:
        for (i = 0; i < len; i++) {
            c [i] = PGARandomInterval (ctx, 'a', 'z');
        }
        break;
    case PGA_CINIT_UPPER:
        for (i = 0; i < len; i++) {
            c [i] = PGARandomInterval (ctx, 'A', 'Z');
        }
        break;
    case PGA_CINIT_MIXED:
        for (i = 0; i < len; i++) {
            j = PGARandomInterval (ctx, 0, 51);
            if (j < 26) {
                c [i] = 'A' + j;
            } else {
                c [i] = 'a' + j - 26;
            }
        }
        break;
    }
    PGADebugExited ("PGACharacterInitString");
}

/*!****************************************************************************
    \brief Build an MPI datatype for a character string.
    \ingroup internal
    \param    ctx   context variable
    \param    p     index of the string to build a datatype from
    \param    pop   symbolic constant of the population string p is in
    \return   MPI_Datatype

    \rst

    Description
    -----------

    Called only by MPI routines.  Not for user consumption.

    \endrst

******************************************************************************/
MPI_Datatype PGACharacterBuildDatatype (PGAContext *ctx, int p, int pop)
{
    int idx = 0;
    /* Number of elements in each block (array of integer) */
    int            counts [PGA_MPI_HEADER_ELEMENTS + 1];
    /* byte displacement of each block (array of integer) */
    MPI_Aint       displs [PGA_MPI_HEADER_ELEMENTS + 1];
    /* type of elements in each block (array of handles to datatype objects) */
    MPI_Datatype   types [PGA_MPI_HEADER_ELEMENTS + 1];
    MPI_Datatype   individualtype; /* new datatype (handle) */
    PGAIndividual *traveller;      /* address of individual in question */

    PGADebugEntered ("PGACharacterBuildDatatype");

    traveller = PGAGetIndividual (ctx, p, pop);

    idx = PGABuildDatatypeHeader (ctx, p, pop, counts, displs, types);

    MPI_Get_address (traveller->chrom, &displs [idx]);
    counts [idx] = ctx->ga.StringLen;
    types  [idx] = MPI_CHAR;
    idx++;

    MPI_Type_create_struct (idx, counts, displs, types, &individualtype);
    MPI_Type_commit (&individualtype);

    PGADebugExited ("PGACharacterBuildDatatype");

    return individualtype;
}

/*!****************************************************************************
    \brief Compute genetic difference of two strings.
    \ingroup internal

    \param   ctx    context variable
    \param   p1     first string index
    \param   pop1   symbolic constant of the population the first string is in
    \param   p2     second string index
    \param   pop2   symbolic constant of the population the second string is in
    \return  genetic distance of the two strings

    \rst

    Description
    -----------

    Return sum of the absolute values of the differences of each allele.
    Internal function.  Use :c:func:`PGAUserFunctionGeneDistance`.

    \endrst

******************************************************************************/
double PGACharacterGeneDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGACharacter *c1 = (PGACharacter *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGACharacter *c2 = (PGACharacter *)PGAGetIndividual (ctx, p2, pop2)->chrom;
    int ret = 0;
    int i;

    PGADebugEntered ("PGACharacterGeneDistance");
    for (i=0; i<ctx->ga.StringLen; i++) {
        if (c1 [i] != c2 [i]) {
            ret++;
        }
    }
    PGADebugExited ("PGACharacterGeneDistance");
    return ret;
}

