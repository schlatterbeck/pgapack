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
* This file contains system routines such as errors and exits.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
char PGAProgram [100];    /* Holds argv[0] for PGAUsage() call */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Report error messages.
    \ingroup reporting
    \param   ctx       context variable
    \param   msg       the error message to print
    \param   level     indicate the error's severity
    \param   datatype  the data type of the following argument
    \param   data      the address of the data to be written out, cast as a void
                       pointer
    \return  None

    \rst

    Description
    -----------

    Prints out the message supplied, and the value of a piece of data.
    Terminates if :c:macro:`PGA_FATAL`. See :ref:`group:const-printflags`
    for the ``level`` constants and :ref:`subsec:other` in the user
    guide for more information.


    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int         val;

       ...
       PGAError
        ( ctx, "Some Non Fatal Error: val = "
        , PGA_WARNING, PGA_INT, (void *) &val
        );

       PGAError (ctx, "A Fatal Error!", PGA_FATAL, PGA_VOID, NULL);

    \endrst

******************************************************************************/
void PGAError (PGAContext *ctx, char *msg, int level, int datatype, void *data)
{

    PGADebugEntered ("PGAError");

    switch (datatype) {
      case PGA_INT:
        fprintf (stderr, "%s %d\n", msg, *(int *)    data);
        break;
      case PGA_DOUBLE:
        fprintf (stderr, "%s %f\n", msg, *(double *) data);
        break;
      case PGA_CHAR:
        fprintf (stderr, "%s %s\n", msg,  (char *)   data);
        break;
      case PGA_VOID:
        fprintf (stderr, "%s\n", msg);
        break;
    }
    if (level == PGA_FATAL) {
        fprintf (stderr, "PGAError: Fatal\n");
        PGADestroy (ctx);
        exit (-1);
    }
    PGADebugExited ("PGAError");
}

/*!****************************************************************************
    \brief Report error messages using printf-like interface.
    \ingroup reporting
    \param   ctx       context variable
    \param   level     indicate the error's severity
    \param   fmt       the printf-style format to print
    \param   ...       additional parameters depending on format string in msg
    \return  None

    \rst

    Description
    -----------

    Prints out the message supplied, and optional parameters.
    Uses a printf-style interface. Terminates if :c:macro:`PGA_FATAL`.
    Note that compared to :c:func:`PGAError` ``msg`` and ``level`` are
    exchanged, seems more natural.  See :ref:`group:const-printflags`
    for the ``level`` constants and :ref:`subsec:other` in the user
    guide for more information.

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int         val;

       ...
       PGAErrorPrintf (ctx, PGA_WARNING, "Some Non Fatal Error: val = %d", val);
       PGAErrorPrintf (ctx, PGA_FATAL, "A Fatal Error!");

    \endrst

******************************************************************************/
void PGAErrorPrintf (PGAContext *ctx, int level, char *fmt, ...)
{
    va_list args;
    va_start (args, fmt);
    vfprintf (stderr, fmt, args);
    va_end (args);
    if (level == PGA_FATAL) {
        fprintf (stderr, "\nPGAError: Fatal\n");
        PGADestroy (ctx);
        exit (-1);
    }
}

/*!****************************************************************************
    \brief Deallocate memory for this instance of PGAPack, if this
           context initialized MPI, finalize MPI as well.
    \ingroup standard-api
    \param   ctx    context variable
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

      PGAContext *ctx;

      ...
      PGADestroy (ctx);

    \endrst

******************************************************************************/
void PGADestroy (PGAContext *ctx)
{
    int i;

    PGADebugEntered ("PGADestroy");

    /*  These are allocated by PGASetUp.  Free then only if PGASetUp
     *  was called.
     */
    if (ctx->sys.SetUpCalled) {
        /*  Free the population...fly little birdies!  You're FREE!!!  */
        for (i=0; i<ctx->ga.PopSize + 2; i++) {
            ctx->cops.ChromFree (ctx->ga.oldpop + i);
            ctx->cops.ChromFree (ctx->ga.newpop + i);
        }
        free (ctx->ga.oldpop);
        free (ctx->ga.newpop);

        /*  Free the scratch space.  */
        free (ctx->scratch.intscratch);
        if (ctx->scratch.permute != NULL) {
            free (ctx->scratch.permute);
        }
        free (ctx->scratch.dblscratch);
        free (ctx->ga.selected);
        free (ctx->ga.sorted);
        /* Need to close output file if we opened it */
        if (  ctx->ga.OutFileName != NULL
           && PGAGetRank (ctx, MPI_COMM_WORLD) == 0
           )
        {
            int err = fclose (ctx->ga.OutputFile);
            if (err != 0) {
                PGAErrorPrintf
                    ( ctx, PGA_WARNING
                    , "Warning: Closing output file returned: %s"
                    , strerror (errno)
                    );
            }
        }
    }

    /*  These are allocated by PGACreate  */
    if (ctx->ga.datatype == PGA_DATATYPE_REAL) {
        free (ctx->init.RealMax);
        free (ctx->init.RealMin);
    } else if (ctx->ga.datatype == PGA_DATATYPE_INTEGER) {
        free (ctx->init.IntegerMax);
        free (ctx->init.IntegerMin);
    }

    /*  We want to finalize MPI only if it was not started for us (as
     *  fortran would do) AND it is actually running.  It would not be
     *  running if, for example, -pgahelp is specified on the command
     *  line.
     */
    MPI_Initialized (&i);
    if ((ctx->par.MPIAlreadyInit == PGA_FALSE) && i) {
        MPI_Finalize ();
    }

    /*  We really should perform a PGADebugPrint here, but we can't;
     *  we've already deallocated most of the stuff we need!!
     */
    free (ctx);
}

/*!***************************************************************************
    \brief Return the largest integer of the current machine
    \ingroup query
    \param   ctx  context variable
    \return  The largest integer of the given machine

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int intmax;

       ...
       intmax = PGAGetMaxMachineIntValue (ctx);

    \endrst

*****************************************************************************/
int PGAGetMaxMachineIntValue (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetMaxMachineIntValue");

    PGADebugExited  ("PGAGetMaxMachineIntValue");

    return ctx->sys.PGAMaxInt;
}

/*!***************************************************************************
    \brief Return the smallest integer of the current machine
    \ingroup query
    \param   ctx  context variable
    \return  The smallest integer of the given machine

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       int intmin;

       ...
       intmin = PGAGetMinMachineIntValue (ctx);

    \endrst

*****************************************************************************/
int PGAGetMinMachineIntValue (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetMinMachineIntValue");

    PGADebugExited  ("PGAGetMinMachineIntValue");

    return ctx->sys.PGAMinInt;
}

/*!***************************************************************************
    \brief Return the largest double of the current machine
    \ingroup query
    \param   ctx  context variable
    \return  The largest double of the given machine

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double big;

       ...
       big = PGAGetMaxMachineDoubleValue (ctx);

    \endrst

*****************************************************************************/
double PGAGetMaxMachineDoubleValue (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetMaxMachineDoubleValue");

    PGADebugExited  ("PGAGetMaxMachineDoubleValue");

    return ctx->sys.PGAMaxDouble;
}

/*!***************************************************************************
    \brief Return the smallest double of the current machine.
    \ingroup query
    \param   ctx  context variable
    \return  The smallest double of the given machine

    \rst

    Example
    -------

    .. code-block:: c

       PGAContext *ctx;
       double small;

       ...
       small = PGAGetMinMachineDoubleValue (ctx);

    \endrst

*****************************************************************************/
double PGAGetMinMachineDoubleValue (PGAContext *ctx)
{
    PGADebugEntered ("PGAGetMinMachineDoubleValue");

    PGADebugExited  ("PGAGetMinMachineDoubleValue");

    return ctx->sys.PGAMinDouble;
}


/*!****************************************************************************
    \brief Print list of available parameters and quit
    \ingroup reporting
    \param   ctx  context variable
    \return  print list of available parameters

    \rst

    Description
    -----------

    Print the usage info out if MPI isn't running (thus, only one process
    is probably running), or if we actually are the rank-0 process.

    Example
    -------

    .. code-block:: c

       PGAContext ctx;

       ...
       PGAUsage (ctx);

    \endrst

******************************************************************************/
void PGAUsage (PGAContext *ctx)
{
    if (!ctx->par.MPIAlreadyInit || (PGAGetRank(ctx, MPI_COMM_WORLD) == 0)) {
        PGAPrintVersionNumber (ctx);
        printf ("PGAPack usage: %s [pga options]\n", PGAProgram);
        printf ("Valid PGAPack options:\n");
        printf ("\t-pgahelp          \tget this message\n");
        printf ("\t-pgahelp debug    \tlist of debug options\n");
        printf ("\t-pgadbg <option>  \tset debug option\n");
        printf ("\t-pgadebug <option>\tset debug option\n");
        printf ("\t-pgaversion       \tprint current PGAPack version number\n");
        printf ("\n");
    }
    PGADestroy (ctx);
    exit (-1);
}

/*!****************************************************************************
    \brief Print PGAPack version number
    \ingroup notimplemented
    \param   ctx  context variable
    \return  print PGAPack version number

    \rst

    Description
    -----------

    This currently prints a meaningless hard-coded string.

    Example
    -------

    .. code-block:: c

       PGAContext ctx;

       ...
       PGAPrintVersionNumber (ctx);

    \endrst

******************************************************************************/
void PGAPrintVersionNumber (PGAContext *ctx)
{
    if (!ctx->par.MPIAlreadyInit || (PGAGetRank(ctx, MPI_COMM_WORLD) == 0)) {
#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
#ifdef FAKE_MPI
#define PRINT1  "Sequential"
#else
#define PRINT1  "Parallel"
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

        printf
            ( "\nPGAPack version 1.0: (%s, %s)\n\n"
            , (OPTIMIZE) ? "Optimized"  : "Debug"
            , PRINT1
            );
    }
}
