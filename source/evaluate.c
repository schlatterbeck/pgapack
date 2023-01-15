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
* This file contains routines specific to the evaluation of the strings.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/
/*!***************************************************************************
 *  \defgroup evaluation Evaluation
 *  \brief Functions used during evaluation and fitness computation
 *****************************************************************************/

#include "pgapack.h"

/*!****************************************************************************
    \brief Maps the value v defined on [a,b] to r defined on [l,u].
    \ingroup internal
    \param  ctx   context variable
    \param  v     value from original interval (usually the decoded bit string)
    \param  a     lower bound of integer interval
    \param  b     upper bound of integer interval
    \param  l     lower bound of real interval
    \param  u     upper bound of real interval
    \return Scaled value of v defined on [l,u]

    \rst

    Description
    -----------

    In the context of PGAPack :math:`[a,b]` is the discrete interval
    :math:`[0,2^{nbits}-1]` (i.e., the number of bits in a binary
    string) and :math:`[l,u]` represent the range of possible values of
    the real number :math:`r`.

    Example
    -------

    Map a five bit (that is, an integer with a range of [0, 31]) integer
    5 to a real in the range [0, 3.14].

    .. code-block:: c

        PGAContext *ctx;
        double x;

        ...
        x = PGAMapIntegerToReal (ctx, 5, 0, 31, 0.0, 3.14);

    \endrst

******************************************************************************/
static double PGAMapIntegerToReal
    (PGAContext *ctx, int v, int a, int b, double l, double u)
{
    double retval = ((v-a) * (u-l) / (b-a) + l);
    /* This may exceed the upper bound due to imprecision */
    if (retval > u) {
        return u;
    }
    return retval;
}

/*!****************************************************************************
    \brief Map the value r defined on [l,u] to v defined on [a,b].
    \ingroup internal

    \param   ctx       context variable
    \param   r         real value defined on [l,u]
    \param   l         lower bound of real interval
    \param   u         upper bound of real interval
    \param   a         lower bound of integer interval
    \param   b         upper bound of integer interval
    \return Scaled value of r defined on [a,b]

    \rst

    Description
    -----------

    In the context of PGAPack :math:`[a,b]` is the discrete interval
    :math:`[0,2^{nbits}-1]`
    (i.e., the number of bits in a binary string) and :math:`[l,u]`
    represent the range of possible values of the real number :math:`r`.

    Example
    -------

    Map the value r on the interval [0, 3.14] to a five bit integer v.

    .. code-block:: c

      PGAContext *ctx;
      double r;
      int v;

      ...
      v = PGAMapRealToInteger (ctx, r, 0.0, 3.14, 0, 31);

    \endrst

******************************************************************************/
static int PGAMapRealToInteger
    (PGAContext *ctx, double r, double l, double u, int a, int b)
{
    return PGARound (ctx, (b - a) * (r - l) / (u - l) + a);
}

/*!****************************************************************************
    \brief Set the evaluation function value for a string to a specified
           value.
    \ingroup evaluation

    \param   ctx   context variable
    \param   p     string index
    \param   pop   symbolic constant of the population string p is in
    \param   val   the (user) evaluation value to assign to string p
    \param   aux   Auxiliary evaluations
    \return  Sets the evaluation function value of string p and the
             EvalUpToDate flag via side effect

    \rst

    Description
    -----------

    This uses a trick for implementing optional arguments in C. The real
    function to use is without leading underscore. There is a
    macro that makes the last argument optional.
    Also sets the evaluation up to date flag to :c:macro:`PGA_TRUE`.

    Example
    -------

    Set the evaluation function value of string ``p`` in population
    :c:macro:`PGA_NEWPOP` to 123.456.

    .. code-block:: c

       PGAContext *ctx;
       int p;
       double aux [...];

       ...
       PGASetEvaluation (ctx, p, PGA_NEWPOP, 123.456);

    or

    .. code-block:: c

       PGASetEvaluation (ctx, p, PGA_NEWPOP, 123.456, aux);

    \endrst

******************************************************************************/
void _PGASetEvaluation
    (PGAContext *ctx, int p, int pop, double val, const double *aux)
{
    PGAIndividual *ind;

    PGADebugEntered ("PGASetEvaluation");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGASetEvaluation"
        , "p = ", PGA_INT, (void *) &p
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGASetEvaluation"
        , "pop = ", PGA_INT, (void *) &pop
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGASetEvaluation"
        , "val = ", PGA_DOUBLE, (void *) &val
        );

    ind               = PGAGetIndividual (ctx, p, pop);
    ind->evalue       = val;
    ind->evaluptodate = PGA_TRUE;
    /* In serial evaluation it will occur that the aux array is modified
     * in place, so we do not need to copy if it's the same pointer
     */
    if (ctx->ga.NumAuxEval) {
        if (aux != NULL && aux != ind->auxeval) {
            memcpy (ind->auxeval, aux, sizeof (double) * ctx->ga.NumAuxEval);
        }
        ind->auxtotalok = PGA_FALSE;
    }

    PGADebugExited ("PGASetEvaluation");
}

/*!***************************************************************************
    \brief Return the evaluation function value for string p in
           population pop and optionally a pointer to the auxiliary
           evaluations.
    \ingroup evaluation

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \param   aux  Pointer to Auxiliary evaluations
    \return  The evaluation function value for string p in population pop

    \rst

    Description
    -----------

    This uses a trick for implementing optional arguments in C. The real
    function to use is without leading underscore. There is a
    macro that makes the last argument optional.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p;
       double eval;

       ...
       eval = PGAGetEvaluation (ctx, p, PGA_NEWPOP);

    or

    .. code-block:: c

       const double *p;
       ...
       eval = PGAGetEvaluation (ctx, p, PGA_NEWPOP, &p);

    \endrst

*****************************************************************************/
double _PGAGetEvaluation (PGAContext *ctx, int p, int pop, const double **aux)
{
    PGAIndividual *ind;

    PGADebugEntered ("PGAGetEvaluation");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGAGetEvaluation"
        , "p = ", PGA_INT, (void *) &p
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGAGetEvaluation"
        , "pop = ", PGA_INT, (void *) &pop
        );

    ind = PGAGetIndividual (ctx, p, pop);

#ifndef OPTIMIZE
    if (ind->evaluptodate != PGA_TRUE) {
        PGAError
            ( ctx, "Evaluation not up to date.  Returning old evaluation."
            , PGA_WARNING, PGA_VOID, NULL
            );
    }
#endif

    if (aux) {
        *aux = ind->auxeval;
    }
    PGADebugExited ("PGAGetEvaluation");
    return ind->evalue;
}

/*!***************************************************************************
    \brief Return the auxiliary evaluation for string p in population
           pop.
    \ingroup evaluation

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \return  The evaluation function value for string p in population pop

    \rst

    Description
    -----------

    This is mostly used internally: the evaluation function will get a
    pointer to this anyway, this is the only point where the aux
    evaluations should be modified.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int p;
       double *aux;

       ...
       aux = PGAGetAuxEvaluation (ctx, p, PGA_NEWPOP);

    \endrst

*****************************************************************************/
double *PGAGetAuxEvaluation (PGAContext *ctx, int p, int pop)
{
    PGAIndividual *ind;

    PGADebugEntered ("PGAGetAuxEvaluation");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGAGetAuxEvaluation"
        , "p = ", PGA_INT, (void *) &p
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGAGetAuxEvaluation"
        , "pop = ", PGA_INT, (void *) &pop
        );

    ind = PGAGetIndividual (ctx, p, pop);

#ifndef OPTIMIZE
    if (ind->evaluptodate != PGA_TRUE) {
        PGAError
            ( ctx, "Evaluation not up to date.  Returning old evaluation."
            , PGA_WARNING, PGA_VOID, NULL
            );
    }
#endif

    PGADebugExited ("PGAGetAuxEvaluation");
    return ind->auxeval;
}

/*!****************************************************************************
    \brief Set the flag to indicate whether the evaluate function value is
           up-to-date or not.
    \ingroup evaluation
    \param  ctx    context variable
    \param  p      string index
    \param  pop    symbolic constant of the population string p is in
    \param  status boolean for whether up-to-date
    \return Sets the flag associated with the evaluation
            function value of string p via side effect

    \rst

    Description
    -----------

    Note that this flag is always set to :c:macro:`PGA_TRUE` when
    :c:func:`_PGASetEvaluation` is called.

    Example
    -------

    Set the evaluation function flag for string ``p`` in population
    :c:macro:`PGA_NEWPOP` to :c:macro:`PGA_FALSE` (as might happen
    after, for example, calling a hill-climbing routine that modified
    this string).

    .. code-block:: c

        PGAContext *ctx;
        int p;

        ...
        PGASetEvaluationUpToDateFlag (ctx, p, PGA_NEWPOP, PGA_FALSE);

    \endrst

******************************************************************************/
void PGASetEvaluationUpToDateFlag (PGAContext *ctx, int p, int pop, int status)
{
    PGAIndividual *ind;

    PGADebugEntered ("PGASetEvaluationUpToDateFlag");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGASetEvaluationUpToDateFlag"
        , "p = ", PGA_INT, (void *) &p
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGASetEvaluationUpToDateFlag"
        , "pop = ", PGA_INT, (void *) &pop
        );

    ind = PGAGetIndividual (ctx, p, pop);

    switch (status) {
    case PGA_TRUE:
    case PGA_FALSE:
        ind->evaluptodate = status;
        /* Invalidate cached auxtotal */
        if (!status) {
            ind->auxtotalok = PGA_FALSE;
        }
        break;
    default:
        PGAError
            ( ctx, "PGASetEvaluationUpToDateFlag: Invalid value of status:"
            , PGA_FATAL, PGA_INT, (void *) &status
            );
        break;
    }

    PGADebugExited ("PGASetEvaluationUpToDateFlag");
}

/*!***************************************************************************
    \brief Return boolean to indicate whether the evaluate function
           value is up to date
    \ingroup evaluation

    \param   ctx  context variable
    \param   p    string index
    \param   pop  symbolic constant of the population the string is in
    \return  Return true if the evaluate function value is up to date

    \rst

    Example
    -------

    .. code-block:: c

        PGAContext *ctx;

        ...
        if (PGAGetEvaluationUpToDateFlag (ctx)) {
            printf ("Evaluation function value current\n");
        } else {
            printf ("Evaluation function value out-of-date\n");
        }

    \endrst

*****************************************************************************/
int PGAGetEvaluationUpToDateFlag (PGAContext *ctx, int p, int pop)
{
    PGAIndividual *ind;

    PGADebugEntered ("PGAGetEvaluationUpToDateFlag");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGAGetEvaluationUpToDateFlag"
        , "p = ", PGA_INT, (void *) &p
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGAGetEvaluationUpToDateFlag"
        , "p = ", PGA_INT, (void *) &pop
        );

    ind = PGAGetIndividual (ctx, p, pop);

    PGADebugExited ("PGAGetEvaluationUpToDateFlag");
    return ind->evaluptodate;
}

/*!****************************************************************************
    \brief Interpret a binary string as encoding a real value and
           return the real value it represents.
    \ingroup allele

    \param    ctx    context variable
    \param    p      string index
    \param    pop    symbolic constant of the population the string is in
    \param    start  starting bit position in the binary representation
    \param    end    ending bit position in the binary representation
    \param    lower  lower bound of the interval the real number is defined on
    \param    upper  upper bound of the interval the real number is defined on
    \return  The real value encoded by the binary string

    \rst

    Example
    -------

    Decode a real value from the string ``p`` in population
    :c:macro:`PGA_NEWPOP`.  The value to decode lies on the interval
    :math:`[-10,20]` and is represented using the 20 bits in bit positions
    10--29.

    .. code-block:: c

        double x;

        ...
        x = PGAGetRealFromBinary (ctx, p, PGA_NEWPOP, 10, 29, -10.0, 20.0);

    \endrst

******************************************************************************/
double PGAGetRealFromBinary
    ( PGAContext *ctx
    , int p, int pop, int start, int end
    , double lower, double upper
    )
{
    int length, sum;
    double value;

    PGADebugEntered  ("PGAGetRealFromBinary");
    PGACheckDataType ("PGAGetRealFromBinary", PGA_DATATYPE_BINARY);

    length = end - start + 1;

    if (start < 0) {
        PGAError
            ( ctx, "PGAGetRealFromBinary: start less than 0:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (end >= PGAGetStringLength (ctx)) {
        PGAError
            ( ctx, "PGAGetRealFromBinary: end greater than string length:"
            , PGA_FATAL, PGA_INT, (void *) &end
            );
    }
    if (start >= end) {
        PGAError
            ( ctx, "PGAGetRealFromBinary: start exceeds end:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (lower >= upper) {
        PGAError
            ( ctx, "PGAGetRealFromBinary: lower exceeds upper:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &lower
            );
    }

    sum = PGAGetIntegerFromBinary (ctx, p, pop, start, end);
    value = PGAMapIntegerToReal
        ( ctx, sum, 0
        , (length == sizeof (unsigned) * 8 - 1) ? INT_MAX : (1u << length) - 1
        , lower, upper
        );

    PGADebugExited ("PGAGetRealFromBinary");

     return (value);
}

/*!****************************************************************************
    \brief Interpret a binary reflected Gray code sequence in a binary
           string as encoding a real value and return the real value it
           represents.
    \ingroup allele

    \param  ctx    context variable
    \param  p      string index
    \param  pop    symbolic constant of the population the string is in
    \param  start  starting bit position in the binary representation
    \param  end    ending bit position in the binary representation
    \param  lower  lower bound of the interval the real number is defined on
    \param  upper  upper bound of the interval the real number is defined on
    \return The real value encoded by the binary reflected Gray code sequence

    \rst

    Example
    -------

    Decode a real value from the string ``p`` in population
    :c:macro:`PGA_NEWPOP`.  The value to decode lies on the interval
    :math:`[-10,20]` and is represented using the 20 bits in bit
    positions 10--29.

    .. code-block:: c

        double x;

        ...
        x = PGAGetRealFromGrayCode (ctx, p, PGA_NEWPOP, 10, 29, -10.0, 20.0);

    \endrst

******************************************************************************/
double PGAGetRealFromGrayCode
    ( PGAContext *ctx
    , int p, int pop, int start, int end
    , double lower, double upper
    )
{
    int length, sum;
    double value;

    PGADebugEntered  ("PGAGetRealFromGrayCode");
    PGACheckDataType ("PGAGetRealFromGrayCode", PGA_DATATYPE_BINARY);

    length = end - start + 1;

    if (start < 0) {
        PGAError
            ( ctx, "PGAGetRealFromGrayCode: start less than 0:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (end >= PGAGetStringLength (ctx)) {
        PGAError
            ( ctx, "PGAGetRealFromGrayCode: end greater than string length:"
            , PGA_FATAL, PGA_INT, (void *) &end
            );
    }
    if (start >= end) {
        PGAError
            ( ctx, "PGAGetRealFromGrayCode: start exceeds end:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (lower >= upper) {
        PGAError
            ( ctx, "PGAGetRealFromGrayCode: lower exceeds upper:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &lower
            );
    }

    sum = PGAGetIntegerFromGrayCode (ctx, p, pop, start, end);
    value = PGAMapIntegerToReal
        ( ctx, sum, 0
        , (length == sizeof (unsigned) * 8 - 1) ? INT_MAX : (1u << length) - 1
        , lower, upper
        );

    PGADebugExited ("PGAGetRealFromGrayCode");

    return value;
}

/*!****************************************************************************
    \brief Encode a real value as a binary string.
    \ingroup allele

    \param ctx    context variable
    \param p      string index
    \param pop    symbolic constant of the population the string is in
    \param start  starting bit position in p to encode val in
    \param end    ending bit position in p to encode val in
    \param low    lower bound of the interval the val is defined on
    \param high   upper bound of the interval the val is defined on
    \param val    the real number to be represented as a binary string
    \return The string is modified by side-effect

    \rst

    Example
    -------

    Encode 3.14 from the interval :math:`[0,10]` in 30 bits in bit
    positions 0--29 in string ``p`` in population :c:macro:`PGA_NEWPOP`.

    .. code-block:: c

        PGAContext *ctx;
        int p;

        ...
        PGAEncodeRealAsBinary (ctx, p, PGA_NEWPOP, 0, 29, 0.0, 10.0, 3.14);

    \endrst

******************************************************************************/
void PGAEncodeRealAsBinary
    ( PGAContext *ctx
    , int p, int pop, int start, int end
    , double low, double high, double val
    )
{
    int length, d;

    PGADebugEntered  ("PGAEncodeRealAsBinary");
    PGACheckDataType ("PGAEncodeRealAsBinary", PGA_DATATYPE_BINARY);

    length = end - start + 1;
    if (start < 0) {
        PGAError
            ( ctx, "PGAEncodeRealAsBinary: start less than 0:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (end >= PGAGetStringLength (ctx)) {
        PGAError
            ( ctx, "PGAEncodeRealAsBinary: end greater than string length:"
            , PGA_FATAL, PGA_INT, (void *) &end
            );
    }
    if (start >= end) {
        PGAError
            ( ctx, "PGAEncodeRealAsBinary: start exceeds end:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (low >= high) {
        PGAError
            ( ctx, "PGAEncodeRealAsBinary: low exceeds high:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &low
            );
    }
    if (val < low || val > high) {
        PGAError
            ( ctx, "PGAEncodeRealAsBinary: val outside of bounds:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    }

    d = PGAMapRealToInteger
        ( ctx, val, low, high, 0
        , (length == sizeof (unsigned) * 8 - 1) ? INT_MAX : (1u << length) - 1
        );
    PGAEncodeIntegerAsBinary (ctx, p, pop, start, end, d);

    PGADebugExited ("PGAEncodeRealAsBinary");
}

/*!****************************************************************************
    \brief Encode a real value as a binary reflected Gray code sequence.
    \ingroup allele
    \param  ctx    context variable
    \param  p      string index
    \param  pop    symbolic constant of the population the string is in
    \param  start  starting bit position in p to encode val in
    \param  end    ending bit position in p to encode val in
    \param  low    lower bound of the interval the val is defined on
    \param  high   upper bound of the interval the val is defined on
    \param  val    the real number to be represented as a binary string
    \return The string is modified by side-effect

    \rst

    Example
    -------

    Encode 3.14 from the interval :math:`[0,10]` in 30 bits in bit
    positions 0--29 in string ``p`` in population :c:macro:`PGA_NEWPOP`
    as a binary reflected Gray code sequence.

    .. code-block:: c

        PGAContext *ctx;
        int p;

        ...
        PGAEncodeRealAsGrayCode (ctx, p, PGA_NEWPOP, 0, 29, 0.0, 10.0, 3.14);

    \endrst

******************************************************************************/
void PGAEncodeRealAsGrayCode
    ( PGAContext *ctx
    , int p, int pop, int start, int end
    , double low, double high, double val
    )
{
    int length, d;

    PGADebugEntered  ("PGAEncodeRealAsGrayCode");
    PGACheckDataType ("PGAEncodeRealAsGrayCode", PGA_DATATYPE_BINARY);

    length = end - start + 1;
    if (start < 0) {
        PGAError
            ( ctx, "PGAEncodeRealAsGrayCode: start less than 0:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (end >= PGAGetStringLength (ctx)) {
        PGAError
            (ctx, "PGAEncodeRealAsGrayCode: end greater than string length:"
            , PGA_FATAL, PGA_INT, (void *) &end
            );
    }
    if (start >= end) {
        PGAError
            ( ctx, "PGAEncodeRealAsGrayCode: start exceeds end:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (low >= high) {
        PGAError
            ( ctx, "PGAEncodeRealAsGrayCode: low exceeds high:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &low
            );
    }
    if (val < low || val > high) {
        PGAError
            ( ctx, "PGAEncodeRealAsGrayCode: val outside of bounds:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &val
            );
    }

    d = PGAMapRealToInteger
        ( ctx, val, low, high, 0
        , (length == sizeof (unsigned) * 8 - 1) ? INT_MAX : (1u << length) - 1
        );
    PGAEncodeIntegerAsGrayCode (ctx, p, pop, start, end, d);

    PGADebugExited ("PGAEncodeRealAsGrayCode");
}


/*!****************************************************************************
    \brief Interpret a binary string as encoding an integer value and
           return the integer value it represents.
    \ingroup allele

    \param  ctx    context variable
    \param  p      string index
    \param  pop    symbolic constant of the population the string is in
    \param  start  starting bit position in the binary representation
    \param  end    ending bit position in the binary representation
    \return The integer value encoded by the binary string

    \rst

    Example
    -------

    Get an integer ``j`` from bits 10--29 of string ``p`` in population
    :c:macro:`PGA_NEWPOP`.

    .. code-block:: c

        PGAContext *ctx;
        int j, p;

        ...
        j = PGAGetIntegerFromBinary (ctx, p, PGA_NEWPOP, 10, 29);

    \endrst

******************************************************************************/
unsigned int PGAGetIntegerFromBinary
    (PGAContext *ctx, int p, int pop, int start, int end)
{
    size_t length;
    int i;
    unsigned int val, power2;

    PGADebugEntered  ("PGAGetIntegerFromBinary");
    PGACheckDataType ("PGAGetIntegerFromBinary", PGA_DATATYPE_BINARY);

    length = end - start + 1;
    if (length > sizeof (int) * 8 - 1) {
        PGAError
            ( ctx
            , "PGAGetIntegerFromBinary: "
              "length of bit string exceeds sizeof type int:"
            , PGA_FATAL, PGA_INT, (void *) &length
            );
    }
    if (start < 0) {
        PGAError
            ( ctx, "PGAGetIntegerFromBinary: start less than 0:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (end >= PGAGetStringLength (ctx)) {
        PGAError
            ( ctx, "PGAGetIntegerFromBinary: end greater than string length:"
            , PGA_FATAL, PGA_INT, (void *) &end
            );
    }
    if (start >= end) {
        PGAError
            ( ctx, "PGAGetIntegerFromBinary: start exceeds end:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    val = 0;
    power2 = 1u << (length - 1);
    for (i=start; i<=end; i++) {
        if (PGAGetBinaryAllele (ctx, p, pop, i)) {
            val += power2;
        }
        power2 >>= 1;
    }

    PGADebugExited ("PGAGetIntegerFromBinary");

    return (val);
}

/*!****************************************************************************
    \brief Interpret a binary reflected Gray code sequence as encoding
           an integer value and return the integer value it represents.
    \ingroup allele

    \param  ctx    context variable
    \param  p      string index
    \param  pop    symbolic constant of the population the string is in
    \param  start  starting bit position in the binary representation
    \param  end    ending bit position in the binary representation
    \return The integer value encoded by the binary reflected Gray code sequence

    \rst

    Example
    -------

    Get an integer ``j`` from bits 10--29 of string ``p`` in population
    :c:macro:`PGA_NEWPOP`.  The string is encoded in Gray code.

    .. code-block:: c

        PGAContext *ctx;
        int j, p;

        ...
        j = PGAGetIntegerFromGrayCode (ctx, p, PGA_NEWPOP, 10, 29);

    \endrst

******************************************************************************/
unsigned int PGAGetIntegerFromGrayCode
    (PGAContext *ctx, int p, int pop, int start, int end)
{
    size_t length, i;
    unsigned int val;
    int *BitString;
    unsigned power2;

    PGADebugEntered  ("PGAGetIntegerFromGrayCode");
    PGACheckDataType ("PGAGetIntegerFromGrayCode", PGA_DATATYPE_BINARY);

    length = end - start + 1;
    if (length > sizeof (int) * 8 - 1) {
        PGAError
            ( ctx
            , "PGAGetIntegerFromGrayCode: "
              "length of binary string exceeds size of type int:"
            , PGA_FATAL, PGA_INT, (void *) &length
            );
    }
    if (start < 0) {
        PGAError
            ( ctx, "PGAGetIntegerFromGrayCode: start less than 0:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (end >= PGAGetStringLength (ctx)) {
        PGAError
            ( ctx, "PGAGetIntegerFromGrayCode: end greater than string length:"
            , PGA_FATAL, PGA_INT, (void *) &end
            );
    }
    if (start >= end) {
        PGAError
            ( ctx, "PGAGetIntegerFromGrayCode: start exceeds end:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }

    BitString = malloc (length * sizeof (int));
    if (!BitString) {
        PGAError
            ( ctx, "PGAGetIntegerFromGrayCode: No room for BitString"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    BitString[0] = PGAGetBinaryAllele (ctx, p, pop, start);

    for (i=1; i<length; i++) {
        BitString [i] =
            BitString [i-1] ^ PGAGetBinaryAllele (ctx, p, pop, start + i);
    }
    val = 0;
    power2 = 1u << (length - 1);
    for (i=0; i<length; i++) {
        if (BitString[i]) {
            val += power2;
        }
        power2 >>= 1;
    }
    free (BitString);

    PGADebugExited ("PGAGetIntegerFromGrayCode");
    return val;
}

/*!****************************************************************************
    \brief Encode an integer value as a binary string.
    \ingroup allele

    \param  ctx    context variable
    \param  p      string index
    \param  pop    symbolic constant of the population the string is in
    \param  start  starting bit position in p to encode val in
    \param  end    ending bit position in p to encode val in
    \param  val    the integer value to be represented as a binary string
    \return The string is modified by side-effect

    \rst

    Example
    -------

    Encode an integer 7 in 20 bits in bit positions 0--19 in string
    ``p`` in population :c:macro:`PGA_NEWPOP`.

    .. code-block:: c

        PGAContext *ctx;
        int p;

        ...
        PGAEncodeIntegerAsBinary (ctx, p, PGA_NEWPOP, 0, 19, 7);

    \endrst

******************************************************************************/
void PGAEncodeIntegerAsBinary
    (PGAContext *ctx, int p, int pop, int start, int end, unsigned int val)
{
    size_t length, i;
    unsigned power2;

    PGADebugEntered  ("PGAEncodeIntegerAsBinary");
    PGACheckDataType ("PGAEncodeIntegerAsBinary", PGA_DATATYPE_BINARY);

    length = end - start + 1;

    if (length > sizeof (int) * 8 - 1) {
        PGAError
            ( ctx
            , "PGAEncodeIntegerAsBinary: "
              "length of bit string exceeds size of type int:"
            , PGA_FATAL, PGA_INT, (void *) &length
            );
    }
    if (start < 0) {
        PGAError
            ( ctx, "PGAEncodeIntegerAsBinary: start less than 0:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (end >= PGAGetStringLength (ctx)) {
        PGAError
            ( ctx, "PGAEncodeIntegerAsBinary: end greater than string length:"
            , PGA_FATAL, PGA_INT, (void *) &end
            );
    }
    if (start >= end) {
        PGAError
            ( ctx, "PGAEncodeIntegerAsBinary: start exceeds end:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if ((val > (1u << length) - 1) && (length != sizeof (int) * 8) - 1) {
        PGAError
            ( ctx
            , "PGAEncodeIntegerAsBinary: Integer too big for string length:"
            , PGA_FATAL, PGA_INT, (void *) &val
            );
    }

    power2 = 1u << (length - 1);
    for (i=0; i<length; i++) {
        if (val >= power2) {
            PGASetBinaryAllele (ctx, p, pop, start + i, 1);
            val -= power2;
        } else {
            PGASetBinaryAllele (ctx, p, pop, start + i, 0);
        }
        power2 >>= 1;
    }

    PGADebugExited ("PGAEncodeIntegerAsBinary");
}

/*!****************************************************************************
    \brief Encode a real value as a binary reflected Gray code sequence.
    \ingroup allele
    \param  ctx    context variable
    \param  p      string index
    \param  pop    symbolic constant of the population the string is in
    \param  start  starting bit position in p to encode val in
    \param  end    ending bit position in p to encode val in
    \param  val    the integer value to be represented as a binary reflected
                   Gray code sequence
    \return The string is modified by side-effect

    \rst

    Example
    -------

    Encode an integer 7 in 20 bits in bit positions 0--19 in string
    ``p`` in population :c:macro:`PGA_NEWPOP` using Gray code.

    .. code-block:: c

        PGAContext *ctx;
        int p;

        ...
        PGAEncodeIntegerAsGrayCode (ctx, p, PGA_NEWPOP, 0, 19, 7);

    \endrst

******************************************************************************/
void PGAEncodeIntegerAsGrayCode
    (PGAContext *ctx, int p, int pop, int start, int end, unsigned int val)
{
    size_t length, i;
    int *bit;
    unsigned power2;

    PGADebugEntered  ("PGAEncodeIntegerAsGrayCode");
    PGACheckDataType ("PGAEncodeIntegerAsGrayCode", PGA_DATATYPE_BINARY);

    length = end - start + 1;

    if (length > sizeof (int) * 8 - 1) {
        PGAError
            ( ctx
            , "PGAEncodeIntegerAsGrayCode: "
              "length of bit string exceeds size of type int:"
            , PGA_FATAL, PGA_INT, (void *) &length
            );
    }
    if (start < 0) {
        PGAError
            ( ctx, "PGAEncodeIntegerAsGrayCode: start less than 0:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if (end >= PGAGetStringLength (ctx)) {
        PGAError
            ( ctx
            , "PGAEncodeIntegerAsGrayCode: end greater than string length:"
            , PGA_FATAL, PGA_INT, (void *) &end
            );
    }
    if (start >= end) {
        PGAError
            ( ctx, "PGAEncodeIntegerAsGrayCode: start exceeds end:"
            , PGA_FATAL, PGA_INT, (void *) &start
            );
    }
    if ((val > (1u << length) - 1) && (length != sizeof (int) * 8 - 1)) {
        PGAErrorPrintf
            ( ctx, PGA_FATAL
            , "PGAEncodeIntegerAsGrayCode: Integer too big "
              "for string length: %u"
            , val
            );
    }
    bit = malloc (length * sizeof (int));
    if (bit == NULL) {
        PGAError
            ( ctx, "PGAEncodeIntegerAsGrayCode: No room to allocate bit"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }
    power2 = 1u << (length - 1);
    for (i=0; i<length; i++) {
        if (val >= power2) {
            bit [i] = 1;
            val -= power2;
        } else {
            bit [i] = 0;
        }
        power2 >>= 1;
    }
    PGASetBinaryAllele (ctx, p, pop, start, bit[0]);
    for (i=1; i<length; i++) {
        PGASetBinaryAllele (ctx, p, pop, start + i, bit [i-1] ^ bit [i]);
    }
    free (bit);

    PGADebugExited ("PGAEncodeIntegerAsGrayCode");
}
