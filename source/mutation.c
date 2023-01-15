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
* This file contains the data structure neutral mutation routines.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

/*!****************************************************************************
    \brief Perform mutation on a string.
    \ingroup explicit

    \param  ctx  context variable
    \param  p    index of string to mutate
    \param  pop  symbolic constant of the population containing p
    \return The number of mutations performed.  Member p in population pop is
            mutated by side-effect

    \rst

    Description
    -----------

    The type of mutation depends on the data type.  Refer to section
    :ref:`sec:mutation` in the user guide for data-specific examples.

    Example
    -------

    Mutate the best string in the population, until 10 or more mutations
    have occured.

    .. code-block:: c

        PGAContext *ctx;
        int p, count = 0;

        ...
        p = PGAGetBestIndex (ctx, PGA_NEWPOP);
        while (count < 10) {
            count += PGAMutate (ctx, p, PGA_NEWPOP);
        }

    \endrst

******************************************************************************/
int PGAMutate (PGAContext *ctx, int p, int pop)
{
    double mr;
    int count;
    int fp;
    PGADebugEntered ("PGAMutate");

    mr    = ctx->ga.MutationProb;
    if (ctx->fops.Mutation) {
        fp = ((p == PGA_TEMP1) || (p == PGA_TEMP2)) ? p : p+1;
        count = (*ctx->fops.Mutation)(&ctx, &fp, &pop, &mr);
    } else {
        count = (*ctx->cops.Mutation)( ctx, p, pop, mr );
    }

    if (count > 0) {
        PGASetEvaluationUpToDateFlag (ctx, p, pop, PGA_FALSE);
    }

    PGADebugExited ("PGAMutate");

    return count;
}

/*!****************************************************************************
    \brief Set type of mutation to use.
    \ingroup init
    \param  ctx            context variable
    \param  mutation_type  symbolic constant to specify the mutation type
    \return None

    \rst

    Description
    -----------

    Only effects integer- and real-valued strings.
    Binary-valued strings are always complemented.
    In character-valued strings, one alphabetic character is replaced with
    another chosen uniformly randomly.  The alphabetic characters will be lower,
    upper, or mixed case depending on how the strings were initialized.

    Valid choices are :c:macro:`PGA_MUTATION_CONSTANT` (Real/Integer),
    :c:macro:`PGA_MUTATION_RANGE` (Real/Integer),
    :c:macro:`PGA_MUTATION_UNIFORM` (Real),
    :c:macro:`PGA_MUTATION_GAUSSIAN` (Real),
    :c:macro:`PGA_MUTATION_PERMUTE` (Integer),
    :c:macro:`PGA_MUTATION_DE` (Real), and
    :c:macro:`PGA_MUTATION_POLY` (Real/Integer).
    The default for integer-valued strings conforms to how the strings
    were initialized.  The default for real-valued strings is
    :c:macro:`PGA_MUTATION_GAUSSIAN`.
    See :ref:`group:const-mutation` for the constants and section
    :ref:`sec:mutation` in the user guide for more details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMutationType (ctx, PGA_MUTATION_UNIFORM);

    \endrst

******************************************************************************/
void PGASetMutationType (PGAContext *ctx, int mutation_type)
{
    PGADebugEntered ("PGASetMutationType");

    switch (mutation_type) {
    case PGA_MUTATION_CONSTANT:
    case PGA_MUTATION_RANGE:
    case PGA_MUTATION_UNIFORM:
    case PGA_MUTATION_GAUSSIAN:
    case PGA_MUTATION_PERMUTE:
    case PGA_MUTATION_DE:
    case PGA_MUTATION_POLY:
         ctx->ga.MutationType = mutation_type;
         break;
    default:
         PGAError
            ( ctx, "PGASetMutationType: Invalid value of mutation_type:"
            , PGA_FATAL, PGA_INT, (void *) &mutation_type
            );
         break;
    }

    PGADebugExited ("PGASetMutationType");
}

/*!***************************************************************************
    \brief Return the type of mutation used.
    \ingroup query
    \param  ctx  context variable
    \return Returns the integer corresponding to the symbolic constant
            used to specify the type of mutation specified

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int mutatetype;

       ...
       mutatetype = PGAGetMutationType (ctx);
       switch (mutatetype) {
       case PGA_MUTATION_CONSTANT:
           printf ("Mutation Type = PGA_MUTATION_CONSTANT\n");
           break;
       case PGA_MUTATION_RANGE:
           printf ("Mutation Type = PGA_MUTATION_RANGE\n");
           break;
       case PGA_MUTATION_UNIFORM:
           printf ("Mutation Type = PGA_MUTATION_UNIFORM\n");
           break;
       case PGA_MUTATION_GAUSSIAN:
           printf ("Mutation Type = PGA_MUTATION_GAUSSIAN\n");
           break;
       case PGA_MUTATION_PERMUTE:
           printf ("Mutation Type = PGA_MUTATION_PERMUTE\n");
           break;
       case PGA_MUTATION_DE:
           printf ("Mutation Type = PGA_MUTATION_DE\n");
           break;
       case PGA_MUTATION_POLY:
           printf ("Mutation Type = PGA_MUTATION_POLY\n");
           break;
       }

    \endrst

*****************************************************************************/
int PGAGetMutationType (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetMutationType");
    PGAFailIfNotSetUp ("PGAGetMutationType");
    PGADebugExited ("PGAGetMutationType");
    return ctx->ga.MutationType;
}

/*!****************************************************************************
    \brief Set multiplier to mutate strings of data type real with.
    \ingroup init

    \param  ctx  context variable
    \param  val  the mutation value to use for Real mutation
    \return None

    \rst

    Description
    -----------

    The use of this value depends on the type of mutation being used.
    The default value is 0.1 unless the mutation type is
    :c:macro:`PGA_MUTATION_CONSTANT` in which case the default is 0.01.
    See section :ref:`sec:mutation` in the user guide for more details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMutationRealValue (ctx, 50.0);

    \endrst

******************************************************************************/
void PGASetMutationRealValue (PGAContext *ctx, double val)
{
    PGADebugEntered ("PGASetMutationRealValue");

    if (val < 0.0) {
        PGAError
            ( ctx, "PGASetMutationRealValue: Invalid value of val:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    } else {
        ctx->ga.MutateRealValue = val;
    }

    PGADebugExited ("PGASetMutationRealValue");
}

/*!***************************************************************************
    \brief Return the value of the multiplier used to mutate
           strings of data type real with.
    \ingroup query

    \param  ctx  context variable
    \return The value of the multiplier

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double val;

       ...
       val = PGAGetMutationRealValue (ctx);

    \endrst

*****************************************************************************/
double PGAGetMutationRealValue (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetMutationRealValue");
    PGAFailIfNotSetUp ("PGAGetMutationRealValue");

    PGADebugExited ("PGAGetMutationRealValue");

    return ctx->ga.MutateRealValue;
}

/*!****************************************************************************
    \brief Set multiplier to mutate data type integer strings with.
    \ingroup init

    \param  ctx  context variable
    \param  val  the mutation value to use for Integer mutation
    \return None

    \rst

    Description
    -----------

    The use of this value depends on the type of mutation being used.
    The default value is 1.  See section :ref:`sec:mutation` of the user
    guide for more details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMutationIntegerValue (ctx, 5);

    \endrst

******************************************************************************/
void PGASetMutationIntegerValue (PGAContext *ctx, int val)
{
    PGADebugEntered ("PGASetMutationIntegerValue");

    if (val < 0.0) {
        PGAError
            ( ctx, "PGASetMutationIntegerValue: Invalid value of val:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    } else {
        ctx->ga.MutateIntegerValue = val;
    }

    PGADebugExited ("PGASetMutationIntegerValue");
}


/*!***************************************************************************
    \brief Return the value of the multiplier used to mutate
           data type integer strings with.
    \ingroup query

    \param  ctx  context variable
    \return The value of the multiplier

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;
        int ival;

        ...
        ival = PGAGetMutationIntegerValue (ctx);

    \endrst

*****************************************************************************/
int PGAGetMutationIntegerValue (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetMutationIntegerValue");
    PGAFailIfNotSetUp ("PGAGetMutationIntegerValue");

    PGADebugExited ("PGAGetMutationIntegerValue");

    return ctx->ga.MutateIntegerValue;
}

/*!****************************************************************************
    \brief If this flag is set to true, then for Integer and Real
           strings whenever a gene is mutated, if it underflows
           (overflows) the lower (upper) bound it is reset to the lower
           (upper) bound.
    \ingroup init

    \param  ctx   context variable
    \param  val   flag to indicate if mutation is bounded
    \return None

    \rst

    Description
    -----------

    When this is enabled, all allele values remain within the range the
    integer strings were initialized on.  If this flag is
    :c:macro:`PGA_FALSE` (the default), the alleles may take any values.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMutationBoundedFlag (ctx, PGA_TRUE);

    \endrst

******************************************************************************/
void PGASetMutationBoundedFlag (PGAContext *ctx, int val)
{
    PGADebugEntered ("PGASetMutationBoundedFlag");

    switch (val) {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.MutateBoundedFlag = val;
         break;
    default:
         PGAError
            ( ctx, "PGASetMutationBoundedFlag: Invalid value:"
            , PGA_FATAL, PGA_INT, (void *) &val
            );
         break;
    }

    PGADebugExited ("PGASetMutationBoundedFlag");
}


/*!****************************************************************************
    \brief Return flag to indicate whether mutated integer strings
           remain in the initialization range.

    \ingroup query

    \param  ctx  context variable
    \return true if restricted to the given range

    \rst

    Description
    -----------

    The initialization range is specified when initialized with
    :c:func:`PGASetIntegerInitRange`.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int val;

       ...
       val = PGAGetMutationBoundedFlag (ctx);

    \endrst

******************************************************************************/
int PGAGetMutationBoundedFlag (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetMutationBoundedFlag");
    PGAFailIfNotSetUp ("PGAGetMutationBoundedFlag");
    PGADebugExited    ("PGAGetMutationBoundedFlag");
    return ctx->ga.MutateBoundedFlag;
}

/*!****************************************************************************
    \brief If this flag is set to true, then for Integer and Real
           strings whenever a gene is mutated, if it underflows
           (overflows) the lower (upper) bound it is reset to a random
           value between the old value and the violated bound.
    \ingroup init

    \param  ctx   context variable
    \param  val   either PGA_TRUE or PGA_FALSE
    \return None

    \rst

    Description
    -----------

    When this flag is set, all allele values remain within the range the
    strings were initialized on by bouncing values that violate the
    bound back from the boundary.  If this flag is :c:macro:`PGA_FALSE`
    (the default), the alleles may take any values. See also
    :c:func:`PGASetMutationBoundedFlag`.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMutationBounceBackFlag (ctx, PGA_TRUE);

    \endrst

******************************************************************************/
void PGASetMutationBounceBackFlag (PGAContext *ctx, int val)
{
    PGADebugEntered ("PGASetMutationBounceBackFlag");

    switch (val) {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.MutateBounceFlag = val;
         break;
    default:
         PGAError
            ( ctx, "PGASetMutationBounceBackFlag: Invalid value:"
            , PGA_FATAL, PGA_INT, (void *) &val
            );
         break;
    }

    PGADebugExited ("PGASetMutationBounceBackFlag");
}


/*!****************************************************************************
    \brief Return flag to indicate whether mutated strings remain within
           the initialization range by bouncing them from the boundary.
    \ingroup query

    \param  ctx  context variable
    \return true if restricted to the given range

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int val;

       ...
       val = PGAGetMutationBounceBackFlag (ctx);

    \endrst

******************************************************************************/
int PGAGetMutationBounceBackFlag (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetMutationBounceBackFlag");
    PGAFailIfNotSetUp ("PGAGetMutationBounceBackFlag");
    PGADebugExited    ("PGAGetMutationBounceBackFlag");
    return ctx->ga.MutateBounceFlag;
}



/*!****************************************************************************
    \brief Specify the probability that a given allele will be mutated.
    \ingroup init
    \param  ctx           context variable
    \param  mutation_prob the mutation probability
    \return None

    \rst

    Description
    -----------

    If this is called without calling :c:func:`PGASetMutationType`, the
    default mutation type is :c:macro:`PGA_MUTATION_CONSTANT`. The default
    probability is the reciprocal of the string length.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMutationProb (ctx, 0.001);

    \endrst

******************************************************************************/
void PGASetMutationProb (PGAContext *ctx, double mutation_prob)
{
    PGADebugEntered ("PGASetMutationProb");

    if ((mutation_prob < 0.0) || (mutation_prob > 1.0)) {
        PGAError
            ( ctx, "PGASetMutationProb: Invalid value of mutation_prob:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &mutation_prob
            );
    } else {
        ctx->ga.MutationProb = mutation_prob;
    }

    PGADebugExited ("PGASetMutationProb");
}

/*!***************************************************************************
    \brief Return the probability of mutation.
    \ingroup query
    \param  ctx  context variable
    \return The mutation probability

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double pm;

       ...
       pm = PGAGetMutationProb (ctx);

    \endrst

*****************************************************************************/
double PGAGetMutationProb (PGAContext *ctx)
{
    PGADebugEntered   ("PGAGetMutationProb");
    PGAFailIfNotSetUp ("PGAGetMutationProb");
    PGADebugExited    ("PGAGetMutationProb");
    return ctx->ga.MutationProb;
}

/*!****************************************************************************
    \brief Set Eta for polynomial mutation.
    \ingroup init
    \param  ctx  context variable
    \param  eta  the polynomial mutation eta
    \return None

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMutationPolyEta (ctx, 200);

    \endrst

******************************************************************************/
void PGASetMutationPolyEta (PGAContext *ctx, double eta)
{
    if (eta < 0.0) {
        PGAError
            ( ctx, "PGASetMutationPolyEta: Invalid value of eta:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &eta
            );
    }
    ctx->ga.MutatePolyEta = eta;
}

/*!***************************************************************************
    \brief Return the Eta for polynomial mutation.
    \ingroup query
    \param  ctx  context variable
    \return The eta for polynomial mutation

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double eta;

       ...
       eta = PGAGetMutationPolyEta (ctx);

    \endrst

*****************************************************************************/
double PGAGetMutationPolyEta (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetMutationPolyEta");
    return ctx->ga.MutatePolyEta;
}

/*!****************************************************************************
    \brief Specify the constant for polynomial mutation.
    \ingroup init
    \param  ctx  context variable
    \param  v    the polynomial mutation constant
    \return None

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetMutationPolyValue (ctx, 2.5);

    \endrst

******************************************************************************/
void PGASetMutationPolyValue (PGAContext *ctx, double v)
{
    if (v < 0.0) {
        PGAError
            ( ctx, "PGASetMutationPolyConstant: Invalid value:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &v
            );
    }
    ctx->ga.MutatePolyValue = v;
}

/*!***************************************************************************
    \brief Return the value for polynomial mutation.
    \ingroup query
    \param  ctx  context variable
    \return The value for polynomial mutation

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double v;

       ...
       v = PGAGetMutationPolyValue (ctx);

    \endrst

*****************************************************************************/
double PGAGetMutationPolyValue (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetMutationPolyValue");
    return ctx->ga.MutatePolyValue;
}


/*!****************************************************************************
    \brief Set the variant used for Differential Evolution.
    \ingroup init
    \param  ctx     context variable
    \param  variant symbolic constant for variant
    \return None

    \rst

    Description
    -----------

    Only used if the mutation type is Differential Evolution
    :c:macro:`PGA_MUTATION_DE`. The possible variants are
    :c:macro:`PGA_DE_VARIANT_RAND`, :c:macro:`PGA_DE_VARIANT_BEST`, and
    :c:macro:`PGA_DE_VARIANT_EITHER_OR`. See
    :ref:`group:const-de-variant` for the constants and section
    :ref:`sec:mutation` in the user guide for details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDEVariant (ctx, PGA_DE_VARIANT_BEST);

    \endrst

******************************************************************************/
void PGASetDEVariant (PGAContext *ctx, int variant)
{
    PGAFailIfSetUp ("PGASetDEVariant");
    switch (variant) {
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

/*!***************************************************************************
    \brief Return the variant of Differential Evolution
    \ingroup query
    \param  ctx  context variable
    \return Returns the integer corresponding to the symbolic constant
            used to specify the variant of differential evolution

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int variant;

       ...
       variant = PGAGetDEVariant (ctx);
       switch (variant) {
       case PGA_DE_VARIANT_RAND:
           printf ("DE Variant = PGA_DE_VARIANT_RAND\n");
           break;
       case PGA_DE_VARIANT_BEST:
           printf ("DE Variant = PGA_DE_VARIANT_BEST\n");
           break;
       case PGA_DE_VARIANT_EITHER_OR:
           printf ("DE Variant = PGA_DE_VARIANT_EITHER_OR\n");
           break;
       }

    \endrst

*****************************************************************************/
int PGAGetDEVariant (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDEVariant");
    return ctx->ga.DEVariant;
}

/*!****************************************************************************
    \brief Set the scale factor F for Differential Evolution
    \ingroup init
    \param  ctx  context variable
    \param  val  the scale factor
    \return None

    \rst

    Description
    -----------

    The default for the scale factor :math:`F` is 0.9. For details see section
    :ref:`sec:mutation` in the user guide.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDEScaleFactor (ctx, 0.75);

    \endrst

******************************************************************************/
void PGASetDEScaleFactor (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 2.0) {
        PGAError
            ( ctx, "PGASetDEScaleFactor: Invalid value of F:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    } else {
        ctx->ga.DEScaleFactor = val;
    }
}

/*!***************************************************************************
    \brief Return the value of the scale factor F for Differential Evolution.
    \ingroup query

    \param  ctx  context variable
    \return The value of the scale factor

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double val;

       ...
       val = PGAGetDEScaleFactor (ctx);

    \endrst

*****************************************************************************/
double PGAGetDEScaleFactor (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDEScaleFactor");
    return ctx->ga.DEScaleFactor;
}

/*!****************************************************************************
    \brief Set the auxiliary factor K for Differential Evolution.
    \ingroup init

    \param  ctx  context variable
    \param  val  the auxiliary factor
    \return None

    \rst

    Description
    -----------

    The default for the aux factor :math:`K` of Differential Evolution is
    :math:`0.5 * (F + 1)` where :math:`F` is the Differential Evolution
    scale factor, see :c:func:`PGASetDEScaleFactor`.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDEAuxFactor (ctx, 0.75);

    \endrst

******************************************************************************/
void PGASetDEAuxFactor (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 2.0) {
        PGAError
            ( ctx, "PGASetDEAuxFactor: Invalid value of K:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    } else {
        ctx->ga.DEAuxFactor = val;
    }
}

/*!***************************************************************************
    \brief Return the value of the auxiliary factor K for Differential
           Evolution.
    \ingroup query

    \param  ctx  context variable
    \return The value of the auxiliary factor

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double val;

       ...
       val = PGAGetDEAuxFactor (ctx);

    \endrst

*****************************************************************************/
double PGAGetDEAuxFactor (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDEAuxFactor");
    return ctx->ga.DEAuxFactor;
}

/*!****************************************************************************
    \brief Set the crossover probability for Differential Evolution.
    \ingroup init

    \param  ctx  context variable
    \param  val  the crossover probability
    \return None

    \rst

    Description
    -----------

    The default for the Differential Evolution crossover probability is
    0.9.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDECrossoverProb (ctx, 0.75);

    \endrst

******************************************************************************/
void PGASetDECrossoverProb (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 1.0) {
        PGAError
            ( ctx, "PGASetDECrossoverProb: Invalid value of crossover:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    } else {
        ctx->ga.DECrossoverProb = val;
    }
}

/*!***************************************************************************
    \brief Return the crossover probability for Differential Evolution.
    \ingroup query

    \param  ctx  context variable
    \return The value of the Differential Evolution crossover probability

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double val;

       ...
       val = PGAGetDECrossoverProb (ctx);

    \endrst

*****************************************************************************/
double PGAGetDECrossoverProb (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDECrossoverProb");
    return ctx->ga.DECrossoverProb;
}

/*!****************************************************************************
    \brief Set the jitter for Differential Evolution.
    \ingroup init

    \param  ctx  context variable
    \param  val  the jitter for Differential Evolution
    \return None

    \rst

    Description
    -----------

    By default jitter is turned off (the value is 0 by default).
    Very small amounts (on the order of 0.001) have been
    recommended for some problems like digital filter design.
    See section :ref:`sec:mutation` in the user guide for details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDEJitter (ctx, 0.001);

    \endrst

******************************************************************************/
void PGASetDEJitter (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 2.0) {
        PGAError
            ( ctx, "PGASetDEJitter: Invalid value of jitter:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    } else {
        ctx->ga.DEJitter = val;
    }
}

/*!***************************************************************************
    \brief Return the jitter for Differential Evolution.
    \ingroup query

    \param  ctx  context variable
    \return The value of the Differential Evolution jitter

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double val;

       ...
       val = PGAGetDEJitter (ctx);

    \endrst

*****************************************************************************/
double PGAGetDEJitter (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDEJitter");
    return ctx->ga.DEJitter;
}

/*!****************************************************************************
    \brief Set the either-or probability of Differential Evolution.
    \ingroup init

    \param  ctx  context variable
    \param  val  the either-or probability
    \return None

    \rst

    Description
    -----------

    The default for this probability is 0.5, it is only used
    when the either-or variant of Differential Evolution has
    been selected by calling :c:func:`PGASetDEVariant` with
    parameter ``PGA_DE_VARIANT_EITHER_OR``.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDEProbabilityEO (ctx, 0.75);

    \endrst

******************************************************************************/
void PGASetDEProbabilityEO (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 1.0) {
        PGAError
            ( ctx, "PGASetDEProbabilityEO: Invalid value of EO probabilty:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    } else {
        ctx->ga.DEProbabilityEO = val;
    }
}

/*!***************************************************************************
    \brief Return the probability of the either-or variant of
           Differential Evolution.
    \ingroup query

    \param  ctx  context variable
    \return The value of the Differential Evolution EO-Probability

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double val;

       ...
       val = PGAGetDEProbabilityEO (ctx);

    \endrst

*****************************************************************************/
double PGAGetDEProbabilityEO (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDEProbabilityEO");
    return ctx->ga.DEProbabilityEO;
}

/*!****************************************************************************
    \brief Set the number of differences for Differential Evolution.
    \ingroup init

    \param  ctx  context variable
    \param  val  the number of differences
    \return None

    \rst

    Description
    -----------

    Some variants of Differential Evolution can specify the number
    of differences that go into the new value of an allele. By
    default this number is 1.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDENumDiffs (ctx, 2);

    \endrst

******************************************************************************/
void PGASetDENumDiffs (PGAContext *ctx, int val)
{
    if (val < 1 || val > 2) {
        PGAError
            ( ctx, "PGASetDENumDiffs: Invalid value of num diffs:"
            , PGA_FATAL, PGA_INT, (void *) &val
            );
    } else {
        ctx->ga.DENumDiffs = val;
    }
}

/*!***************************************************************************
    \brief Return the number of differences for Differential Evolution.
    \ingroup query

    \param  ctx  context variable
    \return The value of the Differential Evolution number of differences

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int val;

       ...
       val = PGAGetDENumDiffs (ctx);

    \endrst

*****************************************************************************/
int PGAGetDENumDiffs (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDENumDiffs");
    return ctx->ga.DENumDiffs;
}

/*!****************************************************************************
    \brief Set the crossover type for Differential Evolution.
    \ingroup init

    \param  ctx  context variable
    \param  val  the crossover type
    \return None

    \rst

    Description
    -----------

    Differential Evolution supports the crossover types
    :c:macro:`PGA_DE_CROSSOVER_BIN` and :c:macro:`PGA_DE_CROSSOVER_EXP`,
    the default is :c:macro:`PGA_DE_CROSSOVER_BIN`. See
    :ref:`group:const-de-cross` for the constants and section
    :ref:`sec:mutation` of the user guide for details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDECrossoverType (ctx, PGA_DE_CROSSOVER_EXP);

    \endrst

******************************************************************************/
void PGASetDECrossoverType (PGAContext *ctx, int val)
{
    switch (val) {
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

/*!***************************************************************************
    \brief Return the Differential Evolution crossover type.
    \ingroup query
    \param  ctx  context variable
    \return The value of the Differential Evolution crossover type

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int val;

       ...
       val = PGAGetDECrossoverType (ctx);

    \endrst

*****************************************************************************/
int PGAGetDECrossoverType (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDECrossoverType");
    return ctx->ga.DECrossoverType;
}

/*!****************************************************************************
    \brief Set the Differential Evolution dither range (+/-).
    \ingroup init
    \param  ctx  context variable
    \param  val  the dither range
    \return None

    \rst

    Description
    -----------

    By default dither is turned off (the value is 0 by default).
    Quite large amounts (on the order of 0.5) are recommended
    for some problems like digital filter design.
    See section :ref:`sec:mutation` in the user guide for details.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDEDither (ctx, 0.5);

    \endrst

******************************************************************************/
void PGASetDEDither (PGAContext *ctx, double val)
{
    if (val < 0.0 || val > 1.0) {
        PGAError
            ( ctx, "PGASetDEDither: Invalid value of Dither:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    } else {
        ctx->ga.DEDither = val;
    }
}

/*!***************************************************************************
    \brief Return the Differential Evolution dither value.
    \ingroup query
    \param  ctx  context variable
    \return The value of the Differential Evolution dither

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double val;

       ...
       val = PGAGetDEDither (ctx);

    \endrst

*****************************************************************************/
double PGAGetDEDither (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetDEDither");
    return ctx->ga.DEDither;
}

/*!****************************************************************************
    \brief Set if Differential Evolution dither is per individual.
    \ingroup init
    \param  ctx   context variable
    \param  val   boolean flag
    \return None

    \rst

    Description
    -----------

    If this is set to :c:macro:`PGA_TRUE`, then for Differential
    Evolution if the dither value is non-zero we produce a new random
    value to add to the scale factor :math:`F` *for each individual*.
    Otherwise if the flag is not set (:c:macro:`PGA_FALSE`), then we
    produce a new value in each generation, the same value for *all*
    individuals.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetDEDitherPerIndividual (ctx, PGA_TRUE);

    \endrst

******************************************************************************/
void PGASetDEDitherPerIndividual (PGAContext *ctx, int val)
{
    switch (val) {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.DEDitherPerIndividual = val;
         break;
    default:
         PGAError
            ( ctx, "PGASetDEDitherPerIndividual: Invalid value:"
            , PGA_FATAL, PGA_INT, (void *) &val
            );
         break;
    }
}


/*!****************************************************************************
    \brief Return boolean flag to indicate whether the dither is applied
           anew for each individual or if the value is re-used for all
           individuals in one generation.
    \ingroup query
    \param  ctx  context variable
    \return flag to indicate if dither is applied for each individual

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int val;

       ...
       val = PGAGetDEDitherPerIndividual (ctx);

    \endrst

******************************************************************************/
int PGAGetDEDitherPerIndividual (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetMutationBounceBackFlag");
    return ctx->ga.DEDitherPerIndividual;
}
