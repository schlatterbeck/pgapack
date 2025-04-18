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
* This file contains routines to generate randomness.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/
/*!***************************************************************************
 *  \defgroup random Functions for randomness
 *****************************************************************************/

#include "pgapack.h"

/*!****************************************************************************
    \brief Flip a biased coin and return true if the coin is a "winner".
    \ingroup random

    \param   ctx  context variable
    \param   p    biased probability (.5 is a fair coin)
    \return  true if coin flip is a winner

    \rst

    Example
    -------

    To return :c:macro:`PGA_TRUE` approximately seventy percent of the
    time, use

    .. code-block:: c

       PGAContext *ctx;
       int bit;

       ...
       bit = PGARandomFlip (ctx, 0.7)

    \endrst

******************************************************************************/
int PGARandomFlip (PGAContext *ctx, double p)
{
    PGADebugEntered ("PGARandomFlip");
    PGADebugExited  ("PGARandomFlip");

    return (PGARandom01 (ctx, 0) < p) ? PGA_TRUE : PGA_FALSE;
}


/*!****************************************************************************
    \brief Return a uniform random number on the specified interval.
    \ingroup random

    \param   ctx    context variable
    \param   start  starting (integer) value of the interval
    \param   end    ending   (integer) value of the interval
    \return  A uniformly distributed random number in the interval [start, end]

    \rst

    Example
    -------

    Generate a value uniformly random from the interval :math:`[0,99]`

    .. code-block:: c

       PGAContext *ctx;
       int r;

       ...
       r = PGARandomInterval (ctx, 0, 99);

    \endrst

******************************************************************************/
int PGARandomInterval (PGAContext *ctx, int start, int end)
{
    PGADebugEntered ("PGARandomInterval");
    PGADebugExited  ("PGARandomInterval");

    return (int)floor (PGARandom01 (ctx, 0) * (double)(end-start+1)) + start;
}

/*!****************************************************************************
    \brief Generate a uniform random number on the interval [0,1).
    \ingroup random

    \param   ctx      context variable
    \param   newseed  either 0 to get the next random number, or nonzero
                      to reseed
    \return  A random number on the interval [0,1)

    \rst

    Description
    -----------

    If the second argument is 0 it returns the next random number in the
    sequence.  Otherwise, the second argument is used as a new seed for the
    population.

    This is a C language implementation of the universal random number
    generator proposed by George Marsaglia, Arif Zaman, and Wai Wan Tsang
    [MZT90]_ and translated from F. James' version [Jam90]_.

    This algorithm is a combination of a lagged Fibonacci and arithmetic
    sequence (F. James) generator with period of :math:`2^{144}`.  It
    provides 32-bit floating point numbers in the range from zero to
    one.  It is claimed to be portable and provides bit-identical
    results on all machines with at least 24-bit mantissas.

    :c:func:`PGARandom01` should be initialized with a 32-bit integer
    seed such that :math:`0 \le seed \le 900,000,000`.
    Each of these 900,000,000 values gives rise to an independent
    sequence of :math:`\approx 10^{30}`.

    Example
    -------

    To get the next random number use

    .. code-block:: c

       PGAContext *ctx;
       double r;

       ...
       r = PGARandom01 (ctx, 0);

    \endrst

******************************************************************************/
double PGARandom01 (PGAContext *ctx, int newseed)
{

    /* initialization variables */
    int ij, kl, i, j, k, l, m, ii, jj;
    float s, t;
    PGARandomState *st = ctx->randstate;

    PGADebugEntered ("PGARandom01");

    /* initialization */

    if (newseed != 0) {

        st->seed = newseed % 900000000;
        ij   = st->seed / 30082;
        kl   = st->seed - 30082 * ij;
        i    = ((ij/177) % 177) + 2;
        j    = ( ij      % 177) + 2;
        k    = ((kl/169) % 178) + 1;
        l    = ( kl      % 169);

        for (ii=0; ii<97; ii++) {

            s = 0.0;
            t = 0.5;

            for (jj=0; jj<24; jj++) {

                m = (((i*j) % 179) * k) % 179;
                i = j;
                j = k;
                k = m;
                l = ((53*l) + 1) % 169;
                if (((l*m) % 64) >= 32) {
                    s += t;
                }
                t *= .5;
            }

            st->u [ii] = s;
        }

        st->c   = 362436.  /16777216.;
        st->cd  = 7654321. /16777216.;
        st->cm  = 16777213./16777216.;
        st->i96 = 96;
        st->j96 = 32;
    }

    /* random number generation */
    st->uni = st->u [st->i96] - st->u [st->j96];
    if (st->uni < 0.) {
        st->uni += 1.0;
    }
    st->u [st->i96] = st->uni;
    st->i96--;
    if (st->i96 < 0) {
        st->i96 = 96;
    }
    st->j96--;
    if (st->j96 < 0) {
        st->j96 = 96;
    }
    st->c   -= st->cd;
    if (st->c < 0.) {
        st->c += st->cm;
    }
    st->uni -= st->c;
    if (st->uni < 0.) {
        st->uni += 1.0;
    }

    PGADebugExited ("PGARandom01");
    return (double) st->uni;
}

/*!****************************************************************************
    \brief Return a uniform random number on the interval [start,end]
    \ingroup random
    \param   ctx    context variable
    \param   start  starting (double) value of the interval
    \param   end    ending   (double) value of the interval
    \return  A random number on the interval [start,end]

    \rst

    Example
    -------

    Generate a uniform random number on the interval :math:`[-0.5, 1.5]`

    .. code-block:: c

       PGAContext *ctx;
       double r;

       ...
       r = PGARandomUniform (ctx, -0.5, 1.5);

    \endrst

******************************************************************************/
double PGARandomUniform (PGAContext *ctx, double start, double end)
{
    double val, r;

    PGADebugEntered ("PGARandomUniform");

    r = PGARandom01 (ctx, 0);
    val = (end-start) * r + start;

    PGADebugExited ("PGARandomUniform");

    return val;
}


/*!****************************************************************************
    \brief Return an approximation to a Gaussian random number
    \ingroup random
    \param   ctx    context variable
    \param   mean   the mean of the Gaussian distribution
    \param   sigma  the standard deviation of the Gaussian distribution
    \return  A random number selected from a Gaussian distribution with
             given mean and standard deviation

    \rst

    Example
    -------

    To generate a Gaussian random number with mean 0.0 and standard
    deviation 1.0 use

    .. code-block:: c

       PGAContext *ctx;
       double r;

       ...
       r = PGARandomGaussian (ctx, 0.0, 1.0);

    \endrst

******************************************************************************/
double PGARandomGaussian (PGAContext *ctx, double mean, double sigma)
{
    int i;
    double sum = 0.;

    PGADebugEntered ("PGARandomGaussian");

    for (i=11; i>=0; i--) {
        sum += PGARandom01 (ctx, 0);
    }

    PGADebugExited ("PGARandomGaussian");

    return (sum - 6.0) * sigma + mean;
}

/*!***************************************************************************
    \brief Return the integer the random number generator was seeded with.
    \ingroup query
    \param   ctx  context variable
    \return  The seed for the random number generator

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int seed;

       ...
       seed = PGAGetRandomSeed (ctx);

    \endrst

*****************************************************************************/
int PGAGetRandomSeed (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetRandomSeed");

    PGADebugExited  ("PGAGetRandomSeed");

    /* When initializing random number generator, seed is modulo 900000000 */
    /* Larger numbers today happen when using time(0) to init seed */
    return ctx->init.RandomSeed % 900000000;
}

/*!****************************************************************************
    \brief Set a seed for the random number generator.
    \ingroup init
    \param   ctx   context variable
    \param   seed  seed  for the random number generator
    \return  None

    \rst

    Description
    -----------

    The default is to initialize the seed randomly (from the time).
    Specifying a seed explicitly allows for reproducibility of runs.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetRandomSeed (ctx, 1);

    \endrst

******************************************************************************/
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
#define MAX_PROCESSORS 2048
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

void PGASetRandomSeed (PGAContext *ctx, int seed)
{

    PGADebugEntered ("PGASetRandomSeed");
    PGAFailIfSetUp  ("PGASetRandomSeed");

    if ((seed < 1) || (seed + MAX_PROCESSORS > 900000000)) {
        PGAError
            ( ctx, "PGASetRandomSeed: Invalid value of seed:"
            , PGA_FATAL, PGA_INT, (void *) &seed
            );
    } else {
        ctx->init.RandomSeed = seed;
    }

    PGADebugExited ("PGASetRandomSeed");
}

/*!****************************************************************************
    \brief Set flag that a deterministic random seed should be used in
    parallel processes
    \ingroup init
    \param   ctx   context variable
    \return  None

    \rst

    Description
    -----------

    By default each parallel MPI process uses its own random number
    generator. So when using random numbers in evaluation or
    hillclimbing the random numbers we get are not reproduceable because
    it is a stochastic process which parallel process evaluates which
    individual. When this function is called we set a flag that uses a
    second random number generator which is deterministically re-seeded
    from the rank-0 process for each evaluation. The net effect is that
    hillclimbing runs (or evaluations that depend on random numbers) can
    be made deterministic (when using the PGAPack random number
    generator) when running again with the same random seed.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;

       ...
       PGASetRandomDeterministic (ctx, PGA_TRUE);

    \endrst

******************************************************************************/
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
#define MAX_PROCESSORS 2048
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

void PGASetRandomDeterministic (PGAContext *ctx, int flag)
{

    PGAFailIfSetUp  ("PGASetRandomDeterministic");

    if (flag) {
        ctx->init.RandomDeterministic = PGA_TRUE;
    } else {
        ctx->init.RandomDeterministic = PGA_FALSE;
    }
}


#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
#define CUTOFF 13
static int sample_a1 (PGASampleState *state);
static int sample_a2 (PGASampleState *state);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Init random sampling of k out of n without replacement.
    \ingroup random
    \param   ctx   context variable
    \param   state pointer to PGASampleState, needs to be allocated by caller
    \param   k     k of the k out of n
    \param   n     n of the k out of n
    \return  None

    \rst

    Description
    -----------

    Algorithm from [Vit87]_.  This is implemented as an iterator, i.e.,
    this init method that initializes ``n`` and ``k`` and a function
    :c:func:`PGARandomNextSample` that returns the next sample index.
    The algorithm guarantees that sample indexes are returned sorted in
    index order.
    Note that the internal implementation variable ``CUTOFF`` is
    arbitrary and no measurements were performed -- modern CPUs can
    probably iterate a lot when not accessing memory, so this can
    probably be set a lot higher.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       PGASampleState state;
       int s;

       ...
       PGARandomSampleInit (ctx, &state, 3, 6);
       s = PGARandomNextSample (&state);

    \endrst

******************************************************************************/

void PGARandomSampleInit (PGAContext *ctx, PGASampleState *state, int k, int n)
{
    PGADebugEntered ("PGARandomSampleInit");
    if (k <= 0) {
        PGAError
            ( ctx, "PGARandomSampleInit: Invalid value of k:"
            , PGA_FATAL, PGA_INT, (void *) &k
            );
    }
    if (n <= 0) {
        PGAError
            ( ctx, "PGARandomSampleInit: Invalid value of n:"
            , PGA_FATAL, PGA_INT, (void *) &n
            );
    }
    if (k > n) {
        PGAError
            ( ctx, "PGARandomSampleInit: Invalid value of k:"
            , PGA_FATAL, PGA_INT, (void *) &k
            );
    }
    memset (state, 0, sizeof (*state));
    state->n   = n;
    state->k   = k;
    state->idx = 0;
    state->ctx = ctx;
    PGADebugExited ("PGARandomSampleInit");
}

/*!****************************************************************************
    \brief Get next sample index for k out of n.
    \ingroup random

    \param   state pointer to PGASampleState, needs to be allocated by caller
                   and initialized by \ref PGARandomSampleInit
    \return  next sample

    \rst

    Description
    -----------

    See :c:func:`PGARandomSampleInit` for documentation and example.

    \endrst
******************************************************************************/

int PGARandomNextSample (PGASampleState *state)
{
    int ret = 0;
    PGAContext *ctx = state->ctx; /* Needed for debug below */
    PGADebugEntered ("PGARandomNextSample");
    if (state->k <= 0) {
        PGAError
            ( state->ctx, "PGARandomNextSample: Invalid value of k:"
            , PGA_FATAL, PGA_INT, (void *) &(state->k)
            );
    }
    if (state->k > 1 && state->n - state->k > CUTOFF) {
        ret = sample_a2 (state);
    } else {
        ret = sample_a1 (state);
    }
    PGADebugExited ("PGARandomNextSample");
    return ret;
}

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

/* This is algorithm A1 from paper above */
static int sample_a1 (PGASampleState *state)
{
    int n = state->n;
    int k = state->k;
    double v = PGARandom01 (state->ctx, 0);
    double top;
    double quot;
    int s = 0;

    if (k == 1) {
        s = (int)(n * v);
    } else {
        top  = n - k;
        quot = top / n;
        while (quot > v) {
            s++;
            top--;
            n--;
            quot = quot * top / n;
        }
    }
    state->idx += s + 1;
    n--;
    k--;
    state->k = k;
    state->n = n;
    return state->idx - 1;
}

/* This is algorithm A2 from paper above but the case check is done in
 * function PGARandomNextSample above
 */
static int sample_a2 (PGASampleState *state)
{
    int n = state->n;
    int k = state->k;
    int qu1 = n - k + 1;
    double vprime;
    double x = 0.0;
    double u;
    double y1, y2;
    int t, top, bottom, limit;
    int s = qu1;

    while (1) {
        while (s >= qu1) {
            vprime = exp (log (PGARandom01 (state->ctx, 0)) / k);
            x = n * (1.0 - vprime);
            s = (int)x;
        }
        u = PGARandom01 (state->ctx, 0);
        y1 = exp (log ((u * n) / qu1) / (k - 1.0));
        vprime = y1 * (-x / (n + 1.0)) * ((double)qu1 / (qu1 - s));
        if (vprime <= 1.0) {
            break;
        }
        y2 = 1.0;
        top = n - 1;
        if (k - 1 > s) {
            bottom = n - k;
            limit  = n - s;
        } else {
            bottom = n - s - 1;
            limit  = qu1;
        }
        for (t = n - 1; t >= limit; t--) {
            y2 = y2 * top / bottom;
            top--;
            bottom--;
        }
        if (n / (n - x) >= y1 * exp (log (y2) / (k - 1))) {
            break;
        }
    }
    state->idx += s + 1;
    n -= s + 1;
    k--;
    state->k = k;
    state->n = n;
    return state->idx - 1;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
