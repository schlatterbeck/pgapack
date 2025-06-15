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
* This file contains routines needed to parse the command line.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/
/*!***************************************************************************
 *  \defgroup cmdline Functions for command line parsing
 *****************************************************************************/

#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
extern char PGAProgram[100];
#define bad_arg(a)    ( ((a)==NULL) || ((*(a)) == '-') )
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#if OPTIMIZE==0
/*!****************************************************************************
    \brief Routine to parse debug command line options, and set the
           appropriate debug level.
    \ingroup internal

    \param   ctx  context variable
    \param   st   debug command line options
    \return  None

    \rst

    Description
    -----------

    Internal function.  Called only by :c:func:`PGAReadCmdLine`.
    Uses :c:func:`PGASetDebugLevel` for setting the debug level.

    \endrst

******************************************************************************/
static void PGAParseDebugArg (PGAContext *ctx, char *st)
{
    size_t        index, length = strlen (st);
    int           num2index = 0, num1index = 0, num1 = 0, num2 = 0, x;
    char          range = 0, num1ch[4], num2ch[4];

    for (index=0; index<length; index++) {
        if (!isdigit (st [index]) && st [index] != ',' && st [index] != '-') {
            PGAFatalPrintf
                (ctx, "PGASetDebugLevel: Invalid Debug Value: %s", st);
        }
        if (st [index] == '-') {
            range = 1;
            num1ch [num1index] = '\0';
            num1 = atoi (num1ch);
            if (num1 < 0 || num1 > PGA_DEBUG_MAXFLAGS) {
                PGAFatalPrintf
                    ( ctx
                    , "PGASetDebugLevel: Lower Limit Out of Range: %d"
                    , num1
                    );
            }
            num1index = 0;
        } else {
            if (isdigit(st[index])) {
                if (range) {
                    num2ch [num2index++] = st [index];
                } else {
                    num1ch [num1index++] = st [index];
                }
            }
            if (st[index] == ',' || index == length) {
                if (range) {
                    num2ch [num2index] = '\0';
                    num2 = atoi (num2ch);
                    if (num2 < 0 || num2 > PGA_DEBUG_MAXFLAGS) {
                        PGAFatalPrintf
                            ( ctx
                            , "PGASetDebugLevel: Upper Limit Out of Range: %d"
                            , num2
                            );
                    }
                    if (num1 <= num2) {
                        for (x = num1; x <= num2; x++) {
                            if (x == 212) {
                                fprintf (stdout, "%s %s\n", num1ch, num2ch);
                            }
                            PGASetDebugLevel (ctx, x);
                        }
                    } else {
                        PGAFatalPrintf
                            ( ctx
                            , "PGASetDebugLevel: Lower Limit Exceeds Upper: %d"
                            , num1
                            );
                    }
                    num2index = 0;
                    range = 0;
                } else {
                    num1ch [num1index] = '\0';
                    num1 = atoi (num1ch);
                    if (num1 < 0 || num1 > PGA_DEBUG_MAXFLAGS) {
                        PGAFatalPrintf
                            ( ctx
                            , "PGASetDebugLevel: Debug Number Out of Range: %d"
                            , num1
                            );
                    }
                    if (num1 == 212) {
                        fprintf (stdout, "%s\n", num1ch);
                    }

                    PGASetDebugLevel (ctx, num1);
                    num1index = 0;
                }
            }
        }
    }
}
#endif

/*!****************************************************************************
    \brief Code to strip arguments out of command list.
    \ingroup internal

    \param   argc  address of the count of the number of command line arguments
    \param   argv  array of command line arguments
    \param   c     current command-line index
    \param   num   number of args to strip
    \return  None

    \rst

    Description
    -----------

    Internal function.  Called only by :c:func:`PGAReadCmdLine`.

    \endrst

******************************************************************************/
static void PGAStripArgs (char **argv, int *argc, int *c, int num)
{
    char **a;
    int i;

    /* Strip out the argument. */
    for (a = argv, i = (*c); i <= *argc; i++, a++)
        *a = (*(a + num));
    (*argc) -= num;
}

/*!****************************************************************************
    \brief Code that looks at the arguments, recognizes any that are for
           PGAPack, uses the arguments, and removes them from the
           command line args.
    \ingroup cmdline

    \param   ctx   context variable
    \param   argc  address of the count of the number of command line arguments
    \param   argv  array of command line arguments
    \return  None

    \rst

    Example
    -------

    .. code-block:: c

       int main (int argc, char **argv)
       {
           PGAContext *ctx;

           PGAReadCmdLine (ctx, &argc, argv);
       }

    \endrst

******************************************************************************/
void PGAReadCmdLine (PGAContext *ctx, int *argc, char **argv)
{
    int c;
    char *s, **a;

    /* Put name of called program (according to the args) into PGAProgram */
    s = strrchr (*argv, '/');
    if (s) {
        strncpy (PGAProgram, s + 1, sizeof (PGAProgram) - 1);
    } else {
        strncpy (PGAProgram, *argv, sizeof (PGAProgram) - 1);
    }

    /* Set all command line flags (except procgroup) to their defaults */

    /* Move to last argument, so that we can go backwards. */
    a = &argv [*argc - 1];

    /*
     * Loop backwards through arguments, catching the ones that start with
     * '-'.  Backwards is more efficient when you are stripping things out.
     */
    for (c = (*argc); c > 1; c--, a--) {
        if (**a != '-') {
            continue;
        }

        if (!strcmp (*a, "-pgadbg") || !strcmp (*a, "-pgadebug")) {
            if bad_arg (a[1]) {
                PGAUsage (ctx);
            }
#if OPTIMIZE==0
            PGAParseDebugArg (ctx, a[1]);
#endif
            PGAStripArgs (a, argc, &c, 2);
            continue;
        }

        if (!strcmp (*a, "-pgaversion")) {
            PGAStripArgs (a, argc, &c, 1);
            PGAPrintVersionNumber (ctx);
            PGADestroy (ctx);
            exit (-1);
        }

        if (!strcmp (*a, "-pgahelp")) {
            if (a[1] == NULL) {
                 PGAUsage(ctx);
            } else {
                if (!strcmp (a[1], "debug")) {
                     PGAPrintDebugOptions (ctx);
                } else {
                    fprintf (stderr, "Invalid option following -pgahelp.\n");
                }
            }
        }
    }
}
