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
* This file implements setting of user functions.
* \authors Authors:
*          David M. Levine, Philip L. Hallstrom, David M. Noelle,
*          Brian P. Walenz, Ralf Schlatterbeck
*****************************************************************************/

#include "pgapack.h"

/*!****************************************************************************
    \brief Specify the name of a user-written function to provide a
           specific GA capability (e.g., crossover, mutation, etc.).
    \ingroup init
    \param   ctx       context variable
    \param   constant  symbolic constant of the user function to set
    \param   f         name of the function to use
    \return  None

    \rst

    Description
    -----------

    This function *must* be used when using a non-native
    datatype and must be called once for each of:

    - :c:macro:`PGA_USERFUNCTION_CREATESTRING`     -- String creation
    - :c:macro:`PGA_USERFUNCTION_MUTATION`         -- Mutation
    - :c:macro:`PGA_USERFUNCTION_CROSSOVER`        -- Crossover
    - :c:macro:`PGA_USERFUNCTION_PRINTSTRING`      -- String Output
    - :c:macro:`PGA_USERFUNCTION_COPYSTRING`       -- Duplication
    - :c:macro:`PGA_USERFUNCTION_DUPLICATE`        -- Duplicate Checking
    - :c:macro:`PGA_USERFUNCTION_GEN_DISTANCE`     -- Genetic Distance
    - :c:macro:`PGA_USERFUNCTION_INITSTRING`       -- Initialization
    - :c:macro:`PGA_USERFUNCTION_BUILDDATATYPE`    -- MPI Datatype creation
    - :c:macro:`PGA_USERFUNCTION_STOPCOND`         -- Stopping conditions
    - :c:macro:`PGA_USERFUNCTION_ENDOFGEN`         --
      Auxiliary functions at the end of each generation
    - :c:macro:`PGA_USERFUNCTION_PRE_EVAL`         --
      Auxiliary functions before evaluation but after crossover and
      mutation
    - :c:macro:`PGA_USERFUNCTION_HASH`             -- Hashing of genes
    - :c:macro:`PGA_USERFUNCTION_SERIALIZE`        -- Serialize userdefined gene
    - :c:macro:`PGA_USERFUNCTION_DESERIALIZE`      --
      Deserialize userdefined gene
    - :c:macro:`PGA_USERFUNCTION_SERIALIZE_FREE`   -- Free serialized version
    - :c:macro:`PGA_USERFUNCTION_CHROM_FREE`       -- Free chromosome

    It *may* be called when using a native datatype to replace the built-in
    functions PGAPack has for that datatype (For example, if the Integer data
    type is used for a traveling salesperson problem, the user may want to
    provide their own custom crossover operator).  See
    :ref:`group:const-ufun` for the constants and chapters
    :ref:`chp:custom1` and :ref:`chp:new-data` in the user guide and the
    examples in the examples directory for more details.

    Example
    -------

    .. code-block:: c

       void MyStringInit (PGAContext *, void *);
       PGAContext *ctx;

       ...
       PGASetUserFunction (ctx, PGA_USERFUNCTION_INITSTRING, MyStringInit);

    \endrst

******************************************************************************/
void PGASetUserFunction (PGAContext *ctx, int constant, void *f)
{
    PGADebugEntered ("PGASetUserFunction");

    if (f == NULL) {
        PGAError
            ( ctx, "PGASetUserFunction: Invalid function"
            , PGA_FATAL, PGA_VOID, NULL
            );
    }

    switch (constant) {
      case PGA_USERFUNCTION_CREATESTRING:
        if (ctx->sys.UserFortran) {
            PGAError
                ( ctx
                , "PGASetUserFunction: Cannot call "
                  "PGA_USERFUNCTION_CREATESTRING from Fortran."
                , PGA_FATAL, PGA_VOID, NULL
                );
        } else {
            ctx->cops.CreateString = (void(*)(PGAContext *, int, int, int))f;
        }
        break;
      case PGA_USERFUNCTION_SERIALIZE:
        if (ctx->sys.UserFortran) {
            PGAError
                ( ctx
                , "PGASetUserFunction: Cannot call "
                  "PGA_USERFUNCTION_SERIALIZE from Fortran."
                , PGA_FATAL, PGA_VOID, NULL
                );
        } else {
            ctx->cops.Serialize =
                (size_t(*)(PGAContext *, int, int, const void **))f;
        }
        break;
      case PGA_USERFUNCTION_DESERIALIZE:
        if (ctx->sys.UserFortran) {
            PGAError
                ( ctx
                , "PGASetUserFunction: Cannot call "
                  "PGA_USERFUNCTION_DESERIALIZE from Fortran."
                , PGA_FATAL, PGA_VOID, NULL
                );
        } else {
            ctx->cops.Deserialize =
                (void(*)(PGAContext *, int, int, const void *, size_t))f;
        }
        break;
      case PGA_USERFUNCTION_SERIALIZE_FREE:
        if (ctx->sys.UserFortran) {
            PGAError
                ( ctx
                , "PGASetUserFunction: Cannot call "
                  "PGA_USERFUNCTION_SERIALIZE_FREE from Fortran."
                , PGA_FATAL, PGA_VOID, NULL
                );
        } else {
            ctx->cops.SerializeFree = (void(*)(void *))f;
        }
        break;
      case PGA_USERFUNCTION_CHROM_FREE:
        if (ctx->sys.UserFortran) {
            PGAError
                ( ctx
                , "PGASetUserFunction: Cannot call "
                  "PGA_USERFUNCTION_CHROM_FREE from Fortran."
                , PGA_FATAL, PGA_VOID, NULL
                );
        } else {
            ctx->cops.ChromFree = (void(*)(PGAIndividual *))f;
        }
        break;
      case PGA_USERFUNCTION_MUTATION:
        if (ctx->sys.UserFortran) {
            ctx->fops.Mutation = (int(*)(void *, void *, void *, void *))f;
        } else {
            ctx->cops.Mutation = (int(*)(PGAContext *, int, int, double))f;
        }
        break;
      case PGA_USERFUNCTION_CROSSOVER:
        if (ctx->sys.UserFortran) {
            ctx->fops.Crossover =
                (void(*)( void *, void *, void *, void *, void *, void *
                        , void *))f;
        } else {
            ctx->cops.Crossover =
                (void(*)(PGAContext *, int, int, int, int, int, int))f;
        }
        break;
      case PGA_USERFUNCTION_PRINTSTRING:
        if (ctx->sys.UserFortran) {
            ctx->fops.PrintString = (void(*)(void *, void *, void *, void *))f;
        } else {
            ctx->cops.PrintString =  (void(*)(PGAContext *, FILE *, int, int))f;
        }
        break;
      case PGA_USERFUNCTION_COPYSTRING:
        if (ctx->sys.UserFortran) {
            PGAError
                ( ctx
                , "PGASetUserFunction: Cannot call "
                  "PGA_USERFUNCTION_COPYSTRING from Fortran."
                , PGA_FATAL, PGA_VOID, NULL
                );
        } else {
            ctx->cops.CopyString = (void(*)(PGAContext *, int, int, int, int))f;
        }
        break;
      case PGA_USERFUNCTION_DUPLICATE:
        if (ctx->sys.UserFortran) {
            ctx->fops.Duplicate =
                (int(*)(void *, void *, void *, void *, void *))f;
        } else {
            ctx->cops.Duplicate = (int(*)(PGAContext *, int, int, int, int))f;
        }
        break;
      case PGA_USERFUNCTION_INITSTRING:
        if (ctx->sys.UserFortran) {
            ctx->fops.InitString = (void(*)(void *, void *, void *))f;
        } else {
            ctx->cops.InitString = (void(*)(PGAContext *, int, int))f;
        }
        break;
      case PGA_USERFUNCTION_BUILDDATATYPE:
        if (ctx->sys.UserFortran) {
            PGAError
                ( ctx
                , "PGASetUserFunction: Cannot call "
                  "PGA_USERFUNCTION_BUILDDATATYPE from Fortran."
                , PGA_FATAL, PGA_VOID, NULL
                );
        } else {
            ctx->cops.BuildDatatype =
                (MPI_Datatype(*)(PGAContext *, int, int))f;
        }
        break;
      case PGA_USERFUNCTION_STOPCOND:
        if (ctx->sys.UserFortran) {
            ctx->fops.StopCond = (int(*)(void *))f;
        } else {
            ctx->cops.StopCond = (int(*)(PGAContext *))f;
        }
        break;
      case PGA_USERFUNCTION_ENDOFGEN:
        if (ctx->sys.UserFortran) {
            ctx->fops.EndOfGen = (void(*)(void *))f;
        } else {
            ctx->cops.EndOfGen = (void(*)(PGAContext *))f;
        }
        break;
      case PGA_USERFUNCTION_PRE_EVAL:
        if (ctx->sys.UserFortran) {
            ctx->fops.PreEval = (void(*)(void *, void *))f;
        } else {
            ctx->cops.PreEval = (void(*)(PGAContext *, int))f;
        }
        break;
      case PGA_USERFUNCTION_GEN_DISTANCE:
        if (ctx->sys.UserFortran) {
            ctx->fops.GeneDistance =
                (double(*)(void *, void *, void *, void *, void *))f;
        } else {
            ctx->cops.GeneDistance =
                (double(*)(PGAContext *, int, int, int, int))f;
        }
        break;
      case PGA_USERFUNCTION_HASH:
        if (ctx->sys.UserFortran) {
            ctx->fops.Hash = (PGAHash(*)(void *, void *, void *))f;
        } else {
            ctx->cops.Hash = (PGAHash(*)(PGAContext *, int, int))f;
        }
        break;
      default:
        PGAError
            ( ctx, "PGASetUserFunction: Invalid constant:"
            , PGA_FATAL, PGA_INT, (void *) &constant
            );
        break;
    }

    PGADebugExited ("PGASetUserFunction");
}
