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
* \file
* This file contains routines specific to the binary datatype.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

/*!***************************************************************************
 *  \defgroup allele Allele Values
 *  \brief Functions for getting and setting allele values
 *****************************************************************************/
/*!***************************************************************************
 *  \defgroup init Initialization
 *  \brief Functions used to change initialization
 *****************************************************************************/
/*!***************************************************************************
 *  \defgroup query Query parameters
 *  \brief Functions to query PGA parameters
 *****************************************************************************/
/*!***************************************************************************
 *  \defgroup explicit Functions for explicit usage
 *  \brief These functions are needed when explicitly programming the GA.
 *  \rsts
 *  See chapter :ref:`chp:explicit` for more details.
 *  \endrst
 *****************************************************************************/
/*!***************************************************************************
 *  \defgroup internal Functions for internal usage
 *  \brief These functions are only used internally
 *****************************************************************************/
/*!***************************************************************************
 *  \defgroup reporting Reporting of algorithm progress
 *  \brief Reporting, output printing
 *****************************************************************************/
/*!***************************************************************************
 *  \defgroup standard-api Standard API
 *  \brief The standard API
 *****************************************************************************/


/*!****************************************************************************
    \brief Sets a binary allele to the specified value.
    \ingroup allele

    \param ctx context variable
    \param p   string index
    \param pop symbolic constant of the population the string is in
    \param i   allele index
    \param val binary value (either 1 or 0) to set the allele to
    \return The allele is changed by side-effect

    \rst

    Example
    -------

    Copies the alleles from member ``p`` in :c:macro:`PGA_OLDPOP` to member
    ``q`` in :c:macro:`PGA_NEWPOP`.

    .. code-block:: c

       PGAContext *ctx;
       int p, q, i;
       int l;

       ...
       l = PGAGetStringLength (ctx);
       for (i=0 i<l; i++) {
           int a = PGAGetBinaryAllele (ctx, p, PGA_OLDPOP, i);
           PGASetBinaryAllele (ctx, q, PGA_NEWPOP, i, a);
       }

    \endrst

******************************************************************************/
void PGASetBinaryAllele (PGAContext *ctx, int p, int pop, int i, int val)
{
    int windex;        /* index of the computer word allele i is in      */
    int bix;           /* bit position in word chrom[windex] of allele i */
    PGAIndividual *ind;
    PGABinary *chrom;

    PGADebugEntered  ("PGASetBinaryAllele");
    PGACheckDataType ("PGASetBinaryAllele", PGA_DATATYPE_BINARY);

    INDEX (windex, bix, i, WL);
    ind = PGAGetIndividual (ctx, p, pop);
    chrom = (PGABinary *)ind->chrom;
    if (val == 0) {
        UNSET (bix, chrom[windex]);
    } else {
        SET (bix, chrom[windex]);
    }

    PGADebugExited ("PGASetBinaryAllele");
}

/*!****************************************************************************
    \brief Return the value of a (binary) allele.
    \ingroup allele
    \param  ctx  context variable
    \param  p    string index
    \param  pop  symbolic constant of the population the string is in
    \param  i    allele index
    \return The value of the ith allele of string p in population pop

    \rst

    Description
    -----------
    
    Applies to the binary data type, set with parameter
    :c:macro:`PGA_DATATYPE_BINARY` of :c:func:`PGACreate`.

    Example
    -------

    Copies the alleles from member ``p`` in :c:macro:`PGA_OLDPOP` to
    member ``q`` :c:macro:`PGA_NEWPOP`.
    Assumes the strings are of the same length.

    .. code-block:: c

        PGAContext *ctx;
        int p, q, i;

        ...
        for (i=PGAGetStringLength (ctx)-1; i>=0; i--) {
            int a = PGAGetBinaryAllele (ctx, p, PGA_OLDPOP, i);
            PGASetBinaryAllele (ctx, q, PGA_NEWPOP, i, a);
        }

    \endrst

******************************************************************************/
int PGAGetBinaryAllele (PGAContext *ctx, int p, int pop, int i)
{

    int windex;        /* index of the computer word allele i is in      */
    int bix;           /* bit position in word chrom[windex] of allele i */
    PGAIndividual *ind;
    PGABinary *chrom;

    PGADebugEntered  ("PGAGetBinaryAllele");
    PGACheckDataType ("PGAGetBinaryAllele", PGA_DATATYPE_BINARY);

    INDEX (windex,bix,i,WL);
    ind = PGAGetIndividual (ctx, p, pop);
    chrom = (PGABinary *)ind->chrom;

    PGADebugExited ("PGAGetBinaryAllele");
    return (BIT (bix, chrom [windex]) != 0);
}

/*!****************************************************************************
    \brief Specify the probability of initializing an allele to "1"
           for the binary data type.
    \ingroup init
    \param   ctx  context variable
    \param   p    the binary initialization probability
    \return  None

    \rst

    Description
    -----------

    The default value is 0.5.
    This is used during string creation of a
    :c:macro:`PGA_DATATYPE_BINARY` string.

    Example
    -------

    Set approximately 1 percent of all binary alleles to "1" when randomly
    initializing the population.

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetBinaryInitProb (ctx, 0.01);

    \endrst

******************************************************************************/
void PGASetBinaryInitProb (PGAContext *ctx, double p)
{
    PGADebugEntered  ("PGASetBinaryInitProb");
    PGAFailIfSetUp   ("PGASetBinaryInitProb");
    PGACheckDataType ("PGASetBinaryInitProb", PGA_DATATYPE_BINARY);

    if ((p <= 1.0) && (p >= 0.0)) {
        ctx->init.BinaryProbability = p;
    } else {
        PGAError
            ( ctx, "PGASetBinaryInitProb: Invalid value of p:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &p
            );
    }
    PGADebugExited ("PGASetBinaryInitProb");
}

/*!***************************************************************************
    \brief Return the probability that an allele will be randomly
           initialized to "1" for data type binary.

    \ingroup query

    \param  ctx  context variable
    \return The probability that a bit will be randomly initialized to one

    \rst

    Description
    -----------

    This is used during string creation of a
    :c:macro:`PGA_DATATYPE_BINARY` string.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double prob;

       ...
       prob = PGAGetBinaryInitProb (ctx);

    \endrst

*****************************************************************************/
double PGAGetBinaryInitProb (PGAContext *ctx)
{
    PGADebugEntered("PGAGetBinaryInitProb");
    PGAFailIfNotSetUp("PGAGetBinaryInitProb");
    PGACheckDataType("PGAGetBinaryInitProb", PGA_DATATYPE_BINARY);

    PGADebugExited("PGAGetBinaryInitProb");
    return(ctx->init.BinaryProbability);
}


/*!****************************************************************************
    \brief Allocate a binary string for member p of population pop.
    \ingroup internal
    \param   ctx       context variable
    \param   p         string index
    \param   pop       symbolic constant of the population string p is in
    \param   initflag  a flag, if set, randomly initialize, else clear alleles
    \return  Member p in population pop is allocated and initialized.

    \rst

    Description
    -----------
    If initflag is :c:macro:`PGA_TRUE`, randomly initialize all alleles,
    otherwise clear all alleles. Applies to
    :c:macro:`PGA_DATATYPE_BINARY` strings.

    Note that this function is set in :c:func:`PGASetUp` as the create
    string user function for the binary datatype by default.

    Example
    -------

    Allocates and clears alleles for all strings in
    :c:macro:`PGA_NEWPOP`.

    .. code-block:: c

        PGAContext *ctx;
        int p;

        ...
        for (p=PGAGetPopSize(ctx)-1; p>=0; p--) {
            PGABinaryCreateString (ctx, p, PGA_NEWPOP, PGA_FALSE);
        }

    \endrst
******************************************************************************/
void PGABinaryCreateString (PGAContext *ctx, int p, int pop, int initflag)
{
    int i, fp;
    PGABinary *s;
    PGAIndividual *new = PGAGetIndividual(ctx, p, pop);

    PGADebugEntered ("PGABinaryCreateString");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGABinaryCreateString"
        , "initflag = ", PGA_INT, (void *) &initflag
        );

    new->chrom = (void *)malloc(ctx->ga.tw * sizeof(PGABinary));
    if (new->chrom == NULL)
        PGAError
            ( ctx, "PGABinaryCreateString: No room to allocate new->chrom"
            , PGA_FATAL, PGA_VOID, NULL
            );

    s = (PGABinary *)new->chrom;
    if (initflag) {
        if (ctx->fops.InitString) {
            fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
            (*ctx->fops.InitString)(&ctx, &fp, &pop);
        } else {
            (*ctx->cops.InitString)(ctx, p, pop);
        }
    } else {
        for (i=0; i<ctx->ga.tw; i++) {
            s[i] = 0;
        }
    }

    PGADebugExited("PGABinaryCreateString");
}

/*!****************************************************************************
    \brief Randomly mutates a bit with a specified probability.
    \ingroup internal

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant for the population string p is in
    \param   mr   probability of mutating (toggling) a bit
    \return  Return the number of mutations

    \rst

    Description
    -----------

    This routine is called from :c:func:`PGAMutate`.

    Note that this function is set in :c:func:`PGASetUp` as the mutation
    user function for the binary datatype by default.

    Example
    -------

    Mutates string ``p`` in population :c:macro:`PGA_NEWPOP` with a
    probability of 0.001 for each bit.

    .. code-block:: c

       PGAContext *ctx;
       int p;

       ...
       PGABinaryMutation (ctx, p, PGA_NEWPOP, .001);

    \endrst

******************************************************************************/
int PGABinaryMutation (PGAContext *ctx, int p, int pop, double mr)
{
    int i,wi;
    int count = 0;
    PGABinary *c;

    PGADebugEntered ("PGABinaryMutation");

    c = (PGABinary *)PGAGetIndividual (ctx, p, pop)->chrom;
    for (wi=0; wi<ctx->ga.fw; wi++) {
        for (i=0; i<(int)WL; ++i) {
            if (PGARandomFlip (ctx, mr)) {
                TOGGLE (i,c[wi]);
                count++;
            }
        }
    }

    /* clean up the partial word if eb > 0 */
    if (ctx->ga.eb > 0) {
        for (i=0; i<ctx->ga.eb; ++i) {
            if (PGARandomFlip (ctx, mr)) {
                TOGGLE (i,c[ctx->ga.fw]);
                count++;
            }
        }
    }

    PGADebugExited ("PGABinaryMutation");
    return (count);

}

/*!****************************************************************************
    \brief Performs one-point crossover on two parent strings to create
           two children via side-effect.
    \ingroup internal

    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population containing c1 and c2
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the binary datatype when selecting
    one-point crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGABinaryOneptCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGABinaryOneptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGABinary *parent1 = (PGABinary *)PGAGetIndividual(ctx, p1, pop1)->chrom;
    PGABinary *parent2 = (PGABinary *)PGAGetIndividual(ctx, p2, pop1)->chrom;
    PGABinary *child1  = (PGABinary *)PGAGetIndividual(ctx, c1, pop2)->chrom;
    PGABinary *child2  = (PGABinary *)PGAGetIndividual(ctx, c2, pop2)->chrom;

    /*
      If the bits are numbered from 0 as follows:

      b   b   b   b   b   b   b   b          b  b
      0   1   2   3   4   5   6   7         30 31

      Then if the cross site is bit 5 (which is the sixth bit by our
      numbering scheme) we would get

      o   o   o   o   o   n   n   n          n  n
      0   1   2   3   4   5   6   7         30 31

      where o indicates the original bit and n is a new bit from the crossover
      operator.
    */

    PGABinary mask;
    int windex;   /* index of the word the crossover bit position is in */
    int bix;      /* bit position to perform crossover (mod WL)         */
    int i;
    int xsite;

    PGADebugEntered("PGABinaryOneptCrossover");

    xsite = PGARandomInterval(ctx, 1,ctx->ga.StringLen-1);

    INDEX(windex,bix,xsite,WL);

    for(i=0;i<windex;i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }

    mask = ~0;
    mask = mask >> bix;

    child1[windex] = (~mask & parent1[windex])|(mask & parent2[windex]);
    child2[windex] = (~mask & parent2[windex])|(mask & parent1[windex]);

    for(i=windex+1;i<ctx->ga.tw;i++) {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }

    PGADebugExited("PGABinaryOneptCrossover");
}


/*!****************************************************************************
    \brief Perform two-point crossover on two parent strings producing
           two children via side-effect
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
    crossover user function for the binary datatype when selecting
    two-point crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGABinaryTwoptCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGABinaryTwoptCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGABinary *parent1 = (PGABinary *)PGAGetIndividual(ctx, p1, pop1)->chrom;
    PGABinary *parent2 = (PGABinary *)PGAGetIndividual(ctx, p2, pop1)->chrom;
    PGABinary *child1  = (PGABinary *)PGAGetIndividual(ctx, c1, pop2)->chrom;
    PGABinary *child2  = (PGABinary *)PGAGetIndividual(ctx, c2, pop2)->chrom;

    PGABinary mask, mask1, mask2;
    int windex1, windex2;
    int bix1, bix2;
    int i;
    int xsite1, xsite2;
    int temp;

    PGADebugEntered("PGABinaryTwoptCrossover");

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

    INDEX(windex1,bix1,xsite1,WL);
    INDEX(windex2,bix2,xsite2,WL);

    if ( windex1 == windex2 ) {     /* both cross sites in the same word */

        for(i=0;i<windex1;i++) {
            child1[i] = parent1[i];
            child2[i] = parent2[i];
        }

        mask1 = ~0;
        if (bix1 == 0)
             mask1 = 0;
        else
             mask1 = mask1 << (WL-bix1);
        mask2 = ~0;
        mask2 = mask2 >> bix2;
        mask  = mask1 | mask2;

        child1[windex1] = (mask & parent1[windex1])|(~mask & parent2[windex1]);
        child2[windex1] = (mask & parent2[windex1])|(~mask & parent1[windex1]);

        for(i=windex1+1;i<ctx->ga.tw;i++) {
            child1[i] = parent1[i];
            child2[i] = parent2[i];
        }
    }
    else {                          /* cross sites in different words */

        for(i=0;i<windex1;i++) {
            child1[i] = parent1[i];
            child2[i] = parent2[i];
        }

        mask = ~0;
        mask = mask >> bix1;

        child1[windex1] = (~mask & parent1[windex1])|(mask & parent2[windex1]);
        child2[windex1] = (~mask & parent2[windex1])|(mask & parent1[windex1]);

        for(i=windex1+1; i<windex2; i++) {
            child1[i] = parent2[i];
            child2[i] = parent1[i];
        }

        mask = ~0;
        mask = mask >> bix2;

        child1[windex2] = (mask & parent1[windex2])|(~mask & parent2[windex2]);
        child2[windex2] = (mask & parent2[windex2])|(~mask & parent1[windex2]);

        for(i=windex2+1; i<ctx->ga.tw; i++) {
            child1[i] = parent1[i];
            child2[i] = parent2[i];
        }
    }

    PGADebugExited("PGABinaryTwoptCrossover");
}


/*!****************************************************************************
    \brief Perform uniform crossover on two parent strings producing two
           children via side-effect.
    \ingroup internal

    \param   ctx   context variable
    \param   p1    the first parent string
    \param   p2    the second parent string
    \param   pop1  symbolic constant of the population containing string
                   p1 and p2
    \param   c1    the first child string
    \param   c2    the second child string
    \param   pop2  symbolic constant of the population to contain string
                   c1 and c2
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    crossover user function for the binary datatype when selecting
    uniform crossover.

    Example
    -------

    Performs crossover on the two parent strings ``m`` and ``d``, producing
    children ``s`` and ``b``.

    .. code-block:: c

       PGAContext *ctx;
       int m, d, s, b;

       ...
       PGABinaryUniformCrossover (ctx, m, d, PGA_OLDPOP, s, b, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGABinaryUniformCrossover
    (PGAContext *ctx, int p1, int p2, int pop1, int c1, int c2, int pop2)
{
    PGABinary *parent1 = (PGABinary *)PGAGetIndividual (ctx, p1, pop1)->chrom;
    PGABinary *parent2 = (PGABinary *)PGAGetIndividual (ctx, p2, pop1)->chrom;
    PGABinary *child1  = (PGABinary *)PGAGetIndividual (ctx, c1, pop2)->chrom;
    PGABinary *child2  = (PGABinary *)PGAGetIndividual (ctx, c2, pop2)->chrom;
    PGABinary mask;
    int j,wi;

    PGADebugEntered ("PGABinaryUniformCrossover");

    for (wi=0; wi<ctx->ga.tw; wi++) {
        if (parent1[wi] == parent2[wi]) {
            child1[wi] = parent1[wi];
            child2[wi] = parent2[wi];
        } else {
            mask = 0;
            for (j=0; j<(int)WL; j++) {
                if (PGARandomFlip (ctx, ctx->ga.UniformCrossProb)) {
                    SET(j,mask);
                }
            }
            child1[wi] = (mask & parent1 [wi]) | (~mask & parent2 [wi]);
            child2[wi] = (mask & parent2 [wi]) | (~mask & parent1 [wi]);
        }
    }

    PGADebugExited ("PGABinaryUniformCrossover");
}

/*!****************************************************************************
    \brief Write a bit string to a file.
    \ingroup internal

    \param   ctx    context variable
    \param   fp     file to write the bit string to
    \param   chrom  pointer to the bit string to write
    \param   nb     number of bits to write out
    \return  None

    \rst

    Description
    -----------

    Puts the binary
    representation of the bit string pointed to by ``chrom`` into a character
    string and writes that out. Assumes the maximum length of string to
    print is ``WL``, and that all bits are in the same word.

    Internal function.  Use :c:func:`PGABinaryPrintString` to print a
    binary string.
    \endrst

******************************************************************************/
static
void PGABinaryPrint (PGAContext *ctx, FILE *fp, PGABinary *chrom, int nb)
{
     char *s, string[WL+1];
     PGABinary mask;
     int i;

     mask = ((PGABinary)1)<<(WL-1);
     s = string;
     for(i=0; i<nb; mask>>=1,i++)              /* mask each bit and set the  */
          *s++ = (mask&(*chrom)?'1':'0');      /* appropriate character      */
     *s=0;                                     /* string terminator          */
     fprintf(fp, "%s", string);                /* print out character string */
}


/*!****************************************************************************
    \brief Write a bit string to a file.
    \ingroup internal
    \param   ctx  context variable
    \param   fp   file pointer to file to write bit string to
    \param   p    index of the string to write out
    \param   pop  symbolic constant of the population string p is in
    \return  None

    \rst

    Example
    -------

    Write string ``s`` to stdout.

    .. code-block:: c

       PGAContext *ctx;
       int s;

       ...
       PGABinaryPrintString (ctx, stdout, s, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGABinaryPrintString (PGAContext *ctx, FILE *fp, int p, int pop)
{
    PGABinary *c = (PGABinary *)PGAGetIndividual (ctx, p, pop)->chrom;
    int wsize = 64 / WL;
    int i, j;
    int done = 0;

    PGADebugEntered ("PGABinaryPrintString");

    for (i=0; i<ctx->ga.fw; i+=wsize) {
         fprintf (fp,"[ ");
         for (j=0; j<wsize && i+j<ctx->ga.fw; j++) {
             PGABinaryPrint (ctx, fp, (c+i+j), WL);
         }
         if (j < wsize && ctx->ga.eb > 0) {
             PGABinaryPrint (ctx, fp, (c+ctx->ga.fw), ctx->ga.eb);
             done = 1;
         }
         fprintf (fp," ]\n");
    }
    if (ctx->ga.eb > 0 && !done) {
         fprintf (fp,"[ ");
         PGABinaryPrint (ctx, fp, (c+ctx->ga.fw), ctx->ga.eb);
         fprintf (fp," ]\n");
    }

    PGADebugExited ("PGABinaryPrintString");
}

/*!****************************************************************************
    \brief Copy one bit string to another.
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
    string user function for the binary datatype by default.

    Example
    -------

    Copy bit string ``x`` to ``y`` (both are implicitly assumed to have
    the same length).

    .. code-block:: c

       PGAContext *ctx;
       int x, y

       ...
       PGABinaryCopyString (ctx, x, PGA_OLDPOP, y, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGABinaryCopyString (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
    PGABinary *source = (PGABinary *)PGAGetIndividual(ctx, p1, pop1)->chrom;
    PGABinary *dest   = (PGABinary *)PGAGetIndividual(ctx, p2, pop2)->chrom;
    int i;

    PGADebugEntered("PGABinaryCopyString");

    for (i = ctx->ga.tw-1; i>=0; i--)
        dest[i] = source[i];

    PGADebugExited("PGABinaryCopyString");
}

/*!****************************************************************************
    \brief Return true if bit string a is a duplicate of bit string b,
           else returns false.
    \ingroup internal

    \param   ctx   context variable
    \param   p1    string index of the first string to compare
    \param   pop1  symbolic constant of the population string p1 is in
    \param   p2    string index of the second string to compare
    \param   pop2  symbolic constant of the population string p2 is in
    \return  Return true/false if strings are duplicates

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    duplicate checking user function for the binary datatype by default.

    Example
    -------

    Compare bit string ``x`` with ``y`` and print a message if they are
    the same.

    .. code-block:: c

       PGAContext *ctx;
       int x, y;

       ...
       if (PGABinaryDuplicate (ctx, x, PGA_NEWPOP, y, PGA_NEWPOP)) {
           printf ("strings are duplicates\n");
       }

    \endrst

******************************************************************************/
int PGABinaryDuplicate (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
     PGABinary *a = (PGABinary *)PGAGetIndividual (ctx, p1, pop1)->chrom;
     PGABinary *b = (PGABinary *)PGAGetIndividual (ctx, p2, pop2)->chrom;
     int wi;

     PGADebugEntered ("PGABinaryDuplicate");

     for (wi=0; wi<ctx->ga.tw; wi++) {
        if (a [wi] != b [wi]) {
            break;
        }
     }

     PGADebugExited ("PGABinaryDuplicate");

     return wi == ctx->ga.tw ? PGA_TRUE : PGA_FALSE;
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
    hash user function for the binary datatype by default.

    \endrst

******************************************************************************/
PGAHash PGABinaryHash (PGAContext *ctx, int p, int pop)
{
     void *a = PGAGetIndividual (ctx, p, pop)->chrom;
     PGAHash hash = PGAUtilHash
        (a, sizeof (PGABinary) * ctx->ga.tw, PGA_INITIAL_HASH);
     return hash;
}

/*!****************************************************************************
    \brief Randomly initialize a string of the binary data type.
    \ingroup internal

    \param   ctx  context variable
    \param   p    index of string to randomly initialize
    \param   pop  symbolic constant of the population string p is in
    \return  None

    \rst

    Description
    -----------

    Note that this function is set in :c:func:`PGASetUp` as the
    init string user function for the binary datatype by default.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p;

       ...
       PGABinaryInitString (ctx, p, PGA_NEWPOP);

    \endrst

******************************************************************************/
void PGABinaryInitString (PGAContext *ctx, int p, int pop)
{
    PGABinary *c = (PGABinary *)PGAGetIndividual (ctx, p, pop)->chrom;
    int i;
    int windex;        /* index of the computer word allele i is in      */
    int bix;           /* binary position in word chrom[windex] of allele i */

    PGADebugEntered ("PGABinaryInitString");

    for (i = 0; i < ctx->ga.tw; i++) {
        c[i] = 0;
    }
    for (i = 0; i < ctx->ga.StringLen; i++) {
        INDEX(windex,bix,i,WL);
        if (PGARandomFlip(ctx, ctx->init.BinaryProbability)) {
            SET(bix, c[windex]);
        }
    }

    PGADebugExited ("PGABinaryInitString");
}

/*!****************************************************************************
    \brief Build an MPI datatype for a binary string datatype.
    \ingroup internal

    \param    ctx   context variable
    \param    p     index of the string to build a datatype from
    \param    pop   symbolic constant of the population string p is in
    \return   MPI_Datatype.

    \rst

    Description
    -----------

    Called only by MPI routines.  Not for user consumption.

    \endrst

******************************************************************************/
MPI_Datatype PGABinaryBuildDatatype (PGAContext *ctx, int p, int pop)
{
    int            idx = 0;
    /* Number of elements in each block (array of integer) */
    int            counts [PGA_MPI_HEADER_ELEMENTS + 1];
    /* byte displacement of each block (array of integer) */
    MPI_Aint       displs [PGA_MPI_HEADER_ELEMENTS + 1];
    /* type of elements in each block (array of handles to datatype objects) */
    MPI_Datatype   types  [PGA_MPI_HEADER_ELEMENTS + 1];
    MPI_Datatype   individualtype; /* new datatype (handle) */
    PGAIndividual *traveller;      /* address of individual in question */

    PGADebugEntered ("PGABinaryBuildDatatype");

    traveller = PGAGetIndividual (ctx, p, pop);
    idx = PGABuildDatatypeHeader (ctx, p, pop, counts, displs, types);

    MPI_Get_address (traveller->chrom, &displs [idx]);
    counts [idx] = ctx->ga.tw;
    types  [idx] = MPI_UNSIGNED_LONG;
    idx++;

    MPI_Type_create_struct (idx, counts, displs, types, &individualtype);
    MPI_Type_commit (&individualtype);

    PGADebugExited ("PGABinaryBuildDatatype");

    return individualtype;
}


/*!****************************************************************************
    \brief Return the Hamming distance between two strings.
    \ingroup internal

    \param   ctx  context variable
    \param   s1   the first string to compare
    \param   s2   the second string to compare
    \return  The Hamming distance between two strings

    \rst

    Description
    -----------

    Note that this function is used in :c:func:`PGABinaryGeneDistance`
    for implementing a generic genetic distance function. It is used in
    the default setting of the :c:macro:`PGA_USERFUNCTION_GEN_DISTANCE`
    user function for the binary data type.


    Example
    -------

    Return the Hamming distance between bit strings ``x`` and ``y``.

    .. code-block:: c

       PGAContext *ctx;
       PGABinary *x, *y;
       int d;

       ...
       d = PGABinaryHammingDistance (ctx, x, y);

    \endrst

******************************************************************************/
int PGABinaryHammingDistance (PGAContext *ctx, PGABinary *s1, PGABinary *s2)
{
    int        j, wi, distance;
    PGABinary  t1, t2, mask;

    PGADebugEntered ("PGABinaryHammingDistance");

    distance = 0;
    for (wi=0; wi<ctx->ga.tw; wi++) { /* step through each word in string   */
        if (s1[wi] != s2[wi]) {     /* if equal, no bits are different      */
            mask = 1;
            for (j=0; j<(int)WL; ++j) {  /* not equal, compare all bits     */
                /* Build bit mask in position j. Mask bit from each         */
                /* string into t1 and t2 and test if bits are the same      */
                t1 = s1 [wi] & mask;
                t2 = s2 [wi] & mask;
                if (t1 != t2) {
                    distance++;
                }
                mask <<= 1;          /* shift mask 1 position */
            }
        }
    }

    PGADebugExited ("PGABinaryHammingDistance");

    return (distance);
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

    For binary genes this is the Hamming distance.  It is used in
    the default setting of the :c:macro:`PGA_USERFUNCTION_GEN_DISTANCE`
    user function for the binary data type.
    Internal function.  Use the gene distance user function
    :c:func:`PGAUserFunctionGeneDistance`.

    \endrst

******************************************************************************/
double PGABinaryGeneDistance
    (PGAContext *ctx, int p1, int pop1, int p2, int pop2)
{
     PGABinary *c1 = (PGABinary *)PGAGetIndividual (ctx, p1, pop1)->chrom;
     PGABinary *c2 = (PGABinary *)PGAGetIndividual (ctx, p2, pop2)->chrom;

     PGADebugEntered("PGABinaryGeneDistance");
     PGADebugExited("PGABinaryGeneDistance");
     return (double)PGABinaryHammingDistance (ctx, c1, c2);
}
