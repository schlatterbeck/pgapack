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
 *  \file
 *  This file contains the data structure neutral crossover routines.
 *  \authors Authors:
 *           David M. Levine, Philip L. Hallstrom, David M. Noelle,
 *           Brian P. Walenz, Ralf Schlatterbeck
 *****************************************************************************/

#include "pgapack.h"

/*!***************************************************************************
 *  \defgroup operators Operators
 *  \brief Genetic Operators used during genetic search
 *****************************************************************************/

/*!***************************************************************************
    \brief Perform crossover on two parent strings to create two child strings
           (via side-effect).
    \ingroup explicit

    \param ctx  the context variable
    \param p1   the first parent string
    \param p2   the second parent string
    \param pop1 symbolic constant of the population containing string
                p1 and p2
    \param c1   the first child string
    \param c2   the second child string
    \param pop2 symbolic constant of the population to contain string
                c1 and c2
    \return c1 and c2 in pop2 are children of p1 and p2 in pop1.  p1 and
            p2 are not modified.

    \rst

    Description
    -----------

    The type of crossover performed is either the default or that
    specified by :c:func:`PGASetCrossoverType`.


    Example
    -------

    Perform crossover on the two parent strings ``mom`` and ``dad`` in
    population :c:macro:`PGA_OLDPOP`, and insert the child strings,
    ``child1`` and ``child1``, in population :c:macro:`PGA_NEWPOP`.

    .. code-block:: c

      PGAContext *ctx;
      int mom, dad, child1, child2;

      ...
      PGACrossover (ctx, mom, dad, PGA_OLDPOP, child1, child2, PGA_NEWPOP);

    \endrst

*****************************************************************************/

void PGACrossover ( PGAContext *ctx, int p1, int p2, int pop1,
                    int c1, int c2, int pop2 )
{
    int fp1, fp2, fc1, fc2;

    PGADebugEntered ("PGACrossover");
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGACrossover", " p1 = "
        , PGA_INT, (void *) &p1
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGACrossover", " p2 = "
        , PGA_INT, (void *) &p2
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGACrossover", " c1 = "
        , PGA_INT, (void *) &c1
        );
    PGADebugPrint
        ( ctx, PGA_DEBUG_PRINTVAR, "PGACrossover", " c2 = "
        , PGA_INT, (void *) &c2
        );

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

/*!***************************************************************************
    \brief Return the type of crossover selected.
    \ingroup query
    \param  ctx context variable
    \return Return the integer corresponding to the symbolic constant
            used to specify the crossover type.

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int crosstype;

       ...
       crosstype = PGAGetCrossoverType (ctx);
       switch (crosstype) {
       case PGA_CROSSOVER_ONEPT:
           printf ("Crossover Type = PGA_CROSSOVER_ONEPT\n");
           break;
       case PGA_CROSSOVER_TWOPT:
           printf ("Crossover Type = PGA_CROSSOVER_TWOPT\n");
           break;
       case PGA_CROSSOVER_UNIFORM:
           printf ("Crossover Type = PGA_CROSSOVER_UNIFORM\n");
           break;
       case PGA_CROSSOVER_SBX:
           printf ("Crossover Type = PGA_CROSSOVER_SBX\n");
           break;
       case PGA_CROSSOVER_EDGE:
           printf ("Crossover Type = PGA_CROSSOVER_EDGE\n");
           break;
       }

    \endrst

*****************************************************************************/

int PGAGetCrossoverType (PGAContext *ctx)
{
    PGADebugEntered("PGAGetCrossoverType");

    PGAFailIfNotSetUp("PGAGetCrossoverType");

    PGADebugExited("PGAGetCrossoverType");

    return(ctx->ga.CrossoverType);
}

/*!***************************************************************************
    \brief Return the crossover probability.
    \ingroup query
    \param  ctx context variable
    \return The crossover probability

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double pc;

       ...
       pc = PGAGetCrossoverProb (ctx);

    \endrst

*****************************************************************************/
double PGAGetCrossoverProb (PGAContext *ctx)
{
    PGADebugEntered("PGAGetCrossoverProb");

    PGAFailIfNotSetUp("PGAGetCrossoverProb");

    PGADebugExited("PGAGetCrossoverProb");

    return(ctx->ga.CrossoverProb);
}

/*!***************************************************************************
    \brief Return the probability of an allele being selected from the
           first child string in uniform crossover.
    \ingroup query
    \param  ctx context variable
    \return The uniform crossover probability

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;
      double pu;

      ...
      pu = PGAGetUniformCrossoverProb (ctx);

    \endrst

*****************************************************************************/
double PGAGetUniformCrossoverProb (PGAContext *ctx)
{
    PGADebugEntered("PGAGetUniformCrossoverProb");

    PGAFailIfNotSetUp("PGAGetUniformCrossoverProb");

    PGADebugExited("PGAGetUniformCrossoverProb");

    return(ctx->ga.UniformCrossProb);
}

/*!****************************************************************************
    \brief Specify the type of crossover to use.
    \ingroup init
    \param ctx            context variable
    \param crossover_type symbolic constant to specify crossover type
    \return None

    \rst

    Description
    -----------

    Valid choices are :c:macro:`PGA_CROSSOVER_ONEPT`,
    :c:macro:`PGA_CROSSOVER_TWOPT`, or :c:macro:`PGA_CROSSOVER_UNIFORM`,
    :c:macro:`PGA_CROSSOVER_SBX`, and :c:macro:`PGA_CROSSOVER_EDGE`
    for one-point, two-point, uniform, simulated binary (SBX), and edge
    crossover, respectively.  The default is :c:macro:`PGA_CROSSOVER_TWOPT`.
    Edge crossover is only defined for integer genes and SBX is only
    defined for integer and real genes. See :ref:`group:const-crossover`
    for the constants and section :ref:`sec:crossover` in the user guide
    for details.

    Example
    -------

    Use uniform crossover when crossingover strings.

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGASetCrossoverType (ctx, PGA_CROSSOVER_UNIFORM);

    \endrst

******************************************************************************/
void PGASetCrossoverType (PGAContext *ctx, int crossover_type)
{

    PGADebugEntered("PGASetCrossoverType");

    switch (crossover_type) {
        case PGA_CROSSOVER_ONEPT:
        case PGA_CROSSOVER_TWOPT:
        case PGA_CROSSOVER_UNIFORM:
        case PGA_CROSSOVER_SBX:
        case PGA_CROSSOVER_EDGE:
            ctx->ga.CrossoverType = crossover_type;
            break;
        default:
            PGAError( ctx,
                     "PGASetCrossoverType: Invalid value of crossover_type:",
                      PGA_FATAL, PGA_INT, (void *) &crossover_type );
    };

    PGADebugExited("PGASetCrossoverType");
}


/*!****************************************************************************
    \brief Set Probability that a selected string will undergo crossover.
    \ingroup init
    \param  ctx  context variable
    \param  p    the crossover probability
    \return  None

    \rst

    Description
    -----------

    The default is 0.85.

    Example
    -------

    Make crossover happen infrequently.

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetCrossoverProb (ctx, 0.001);

    \endrst

******************************************************************************/
void PGASetCrossoverProb (PGAContext *ctx, double p)
{
    PGADebugEntered ("PGASetCrossoverProb");

    if ((p < 0.0) || (p > 1.0)) {
        PGAError
            ( ctx, "PGASetCrossoverProb: Invalid value of p:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &p
            );
    } else {
        ctx->ga.CrossoverProb = p;
    }

    PGADebugExited ("PGASetCrossoverProb");
}

/*!****************************************************************************
    \brief Set probability used in uniform crossover to specify that an
           allele value value be selected from a particular parent.
    \ingroup init
    \param  ctx  context variable
    \param  p    the crossover probability
    \return  None

    \rst

    Description
    -----------

    The default is 0.6. The crossover type must have been set to
    :c:macro:`PGA_CROSSOVER_UNIFORM` with :c:func:`PGASetCrossoverType`
    for this function call to have any effect.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetUniformCrossoverProb (ctx, 0.9);

    \endrst

******************************************************************************/
void PGASetUniformCrossoverProb (PGAContext *ctx, double p)
{
    PGADebugEntered ("PGASetUniformCrossoverProb");

    if ((p < 0.0) || (p > 1.0)) {
        PGAError
            ( ctx, "PGASetUniformCrossoverProb: Invalid value of p:"
            , PGA_FATAL, PGA_DOUBLE, (void *) &p
            );
    } else {
        ctx->ga.UniformCrossProb = p;
    }

    PGADebugExited ("PGASetUniformCrossoverProb");
}

/*!****************************************************************************
    \brief If this flag is set to true, then for Integer and Real
           strings with simulated binary crossover (SBX) crossed over
           values that exceed the bounds are confined to the bounds by
           setting them to the boundary.
    \ingroup init

    \param  ctx   context variable
    \param  flag  to indicate if strings should be constrained to boundary
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetCrossoverBoundedFlag (ctx, PGA_TRUE);

    \endrst

******************************************************************************/
void PGASetCrossoverBoundedFlag (PGAContext *ctx, int flag)
{
    switch (flag)
    {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.CrossBoundedFlag = flag;
         break;
    default:
         PGAError
            ( ctx, "PGASetCrossoverBoundedFlag: Invalid value:"
            , PGA_FATAL, PGA_INT, (void *) &flag
            );
         break;
    }
}

/*!****************************************************************************
    \brief Return boolean value to indicate whether crossed
           over strings remain in the range specified.
    \ingroup query

    \param  ctx  context variable
    \returns boolean value indicating if strings remain in range

    \rst

    Description
    -----------

    Return :c:macro:`PGA_TRUE` if restricted to the given range,
    otherwise :c:macro:`PGA_FALSE`.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int val;

       ...
       val = PGAGetCrossoverBoundedFlag (ctx);

    \endrst

******************************************************************************/
int PGAGetCrossoverBoundedFlag (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetCrossoverBoundedFlag");
    return (ctx->ga.CrossBoundedFlag);
}

/*!****************************************************************************
    \brief If this flag is set to true, then for Integer and Real
           strings with simulated binary crossover (SBX) crossed over
           values that exceed the bounds are confined to the bounds by
           bouncing them back to a random value between the boundary and
           the neares parent.
    \ingroup init
    \param  ctx   context variable
    \param  flag  to indicate whether out-of-range values should be bounced
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetCrossoverBounceBackFlag (ctx, PGA_TRUE);

    \endrst

******************************************************************************/
void PGASetCrossoverBounceBackFlag (PGAContext *ctx, int flag)
{
    switch (flag)
    {
    case PGA_TRUE:
    case PGA_FALSE:
         ctx->ga.CrossBounceFlag = flag;
         break;
    default:
         PGAError(ctx, "PGASetCrossoverBounceBackFlag: Invalid value:",
                  PGA_FATAL, PGA_INT, (void *) &flag);
         break;
    }
}

/*!****************************************************************************
    \brief Return boolean to indicate whether crossed over strings are
           bounced back when exceeding the bounds.
    \ingroup query
    \param  ctx  context variable
    \return flag indicating whether out-of-range values are bounced back


    \rst

    Description
    -----------

    Return :c:macro:`PGA_TRUE` if restricted to the given range by
    bouncing out-of-range values back from the boundary, otherwise
    :c:macro:`PGA_FALSE`.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int val;

       ...
       val = PGAGetCrossoverBounceBackFlag (ctx);

    \endrst

******************************************************************************/
int PGAGetCrossoverBounceBackFlag (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetCrossoverBounceBackFlag");
    return (ctx->ga.CrossBounceFlag);
}

/*!****************************************************************************
    \brief Set the eta parameter for simulated binary crossover (SBX).
    \ingroup init
    \param  ctx   context variable
    \param  eta   eta >= 0
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGASetCrossoverSBXEta (ctx, 10);

    \endrst

******************************************************************************/
void PGASetCrossoverSBXEta (PGAContext *ctx, double eta)
{
    if (eta < 0) {
        PGAError(ctx, "PGASetCrossoverSBXEta: Invalid value:",
                 PGA_FATAL, PGA_DOUBLE, (void *) &eta);
    }
    ctx->ga.CrossSBXEta = eta;
}

/*!****************************************************************************
    \brief Return simulated binary crossover (SBX) eta value.
    \ingroup query
    \param  ctx  context variable
    \return  The SBX eta value

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double eta;

       ...
       eta = PGAGetCrossoverSBXEta (ctx);

    \endrst

******************************************************************************/
double PGAGetCrossoverSBXEta (PGAContext *ctx)
{
    PGAFailIfNotSetUp ("PGAGetCrossoverSBXEta");
    return (ctx->ga.CrossSBXEta);
}

/*!****************************************************************************
    \brief Compute random number for simulated binary crossover (SBX)
           polynomial distribution only once per string/individual.
    \ingroup init
    \param  ctx   context variable
    \param  val   flag indicating if random number is computed once per string
    \return None

    \rst

    Description
    -----------

    If set to :c:macro:`PGA_TRUE` all alleles will use the same value
    which means that the resulting string will point into the same
    direction as the vector between both parents. The default is
    :c:macro:`PGA_FALSE` indicating that a new random number is used for
    each string.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetCrossoverSBXOncePerString (ctx, PGA_TRUE);
    \endrst

******************************************************************************/
void PGASetCrossoverSBXOncePerString (PGAContext *ctx, int val)
{
    if (val != PGA_TRUE && val != PGA_FALSE) {
        PGAError(ctx, "PGASetCrossoverSBXOncePerString: Invalid value:",
                 PGA_FATAL, PGA_INT, (void *) &val);
    }
    ctx->ga.CrossSBXOnce = val;
}

/*!****************************************************************************
    \brief Return simulated binary crossover (SBX) setting if random number
           for SBX polynomial distribution is computed once per string.
    \ingroup query
    \param  ctx  context variable
    \return The SBX once-per-string value

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int r;

       ...
       r = PGAGetCrossoverSBXOncePerString (ctx);

    \endrst

******************************************************************************/
int PGAGetCrossoverSBXOncePerString (PGAContext *ctx)
{
    PGAFailIfNotSetUp("PGAGetCrossoverSBXOncePerString");
    return (ctx->ga.CrossSBXOnce);
}


/*!****************************************************************************
    \brief Cross over two parent alleles with simulated binary crossover (SBX).
    \ingroup internal
    \param   ctx context variable
    \param   p1  (double) Allele of first string
    \param   p2  (double) Allele of second string
    \param   u   Random value between 0 and 1
    \param   c1  pointer to new first child allele
    \param   c2  pointer to new second child allele
    \return  None

    \rst

    Description
    -----------

    This uses double for both parent alleles but is re-used in
    both, integer and real SBX crossover in
    :c:func:`PGAIntegerSBXCrossover` and :c:func:`PGARealSBXCrossover`,
    respectively. The probability is used to
    compute the new alleles from the polynomial distribution.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double p1, p2, u;
       double c1, c2;
       double result;

       ...
       u = PGARandom01 (ctx, 0);
       result = PGACrossoverSBX (ctx, p1, p2, u, &c1, &c2);

    \endrst

******************************************************************************/
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
