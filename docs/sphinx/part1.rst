.. _part:started:

Getting Started
---------------

.. _chp:intro:

Introduction
++++++++++++

PGAPack is a parallel genetic algorithm library that is intended to
provide most capabilities desired in a genetic algorithm package, in an
integrated, seamless, and portable manner. Key features of PGAPack are
as follows:

-  Ability to be called from Fortran or C.

-  Executable on uniprocessors, multiprocessors, multicomputers, and
   workstation networks.

-  Binary-, integer-, real-, and character-valued native data types.

-  Object-oriented data structure neutral design.

-  Parameterized population replacement.

-  Multiple choices for selection, crossover, and mutation operators.

-  Easy integration of hill-climbing heuristics.

-  Easy-to-use interface for novice and application users.

-  Multiple levels of access for expert users.

-  Full extensibility to support custom operators and new data types.

-  Extensive debugging facilities.

-  Large set of example problems.

.. _chp:examples:

Examples
++++++++

This chapter presents some simple PGAPack programs. The problem chosen
is the Maxbit problem. The objective is to maximize the number of 1-bits
in a string.

Section :ref:`sec:simple-example` presents a simple PGAPack program
in C whose structure is sufficient to solve many problems.
Section :ref:`sec:simple-example-fortran` presents this same program
in Fortran. Section :ref:`sec:default-values` shows how to change
default values in PGAPack.
Section :ref:`sec:parallel-simple-example` contains an example that
shows how keyboard input may be read in an MPI environment. Finally,
Section :ref:`sec:execute` shows how to compile, link, and execute a
PGAPack program. These and other examples may be found in the
``./examples/c`` and ``./examples/fortran`` directories.

.. _sec:simple-example:

Maxbit Problem in C
~~~~~~~~~~~~~~~~~~~

.. _example:simple-main:

.. code-block:: c
   :caption: PGAPack C Program for the Maxbit Example

   #include "pgapack.h"
   double evaluate (PGAContext *ctx, int p, int pop, double *aux);

   int main(int argc, char **argv)
   {
       PGAContext *ctx; 
       ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
       PGASetUp        (ctx);
       PGARun          (ctx, evaluate);
       PGADestroy      (ctx);
       return;
   }

   double evaluate (PGAContext *ctx, int p, int pop, double *aux)
   {
       int i, nbits, stringlen;

       stringlen = PGAGetStringLength (ctx);
       nbits     = 0;
       for (i=0; i<stringlen; i++)
           if (PGAGetBinaryAllele (ctx, p, pop, i))
               nbits++;
       return (double) nbits;
   }

Figure :ref:`example:simple-main` shows a minimal
program and evaluation function in C for the Maxbit problem. All
PGAPack C programs *must* include the header file ``pgapack.h``. The
:c:func:`PGACreate` call is always the first function called in a
PGAPack program. It initializes the context variable, ``ctx``. The
parameters to :c:func:`PGACreate` are the arguments to the program (given by
``argc`` and ``argv``), the data type selected
(:c:macro:`PGA_DATATYPE_BINARY`), the string length (``100``), and the
direction of optimization (:c:macro:`PGA_MAXIMIZE`). The
:c:func:`PGASetUp` call initializes all parameters and function pointers
not explicitly set by the user to default values.

:c:func:`PGARun` executes the genetic algorithm. Its second argument is the
name of a user-defined function (``evaluate``) that will be called to
evaluate the strings. :c:func:`PGADestroy` releases all memory allocated by
PGAPack. Note that all PGAPack functions take the context variable as an
argument (except :c:func:`PGACreate`, which creates the context variable).

The ``evaluate`` function must be written by the user, must return a
``double``, and must follow the exact calling sequence shown. An
evaluation function may return more values than just the return value in
the array pointed to by ``aux``. This can be used for evaluating
constraints for constrained problems or for multi-objective
optimization. Often the number of auxiliary return values is zero and
the ``aux`` argument is ignored. For details, see
section :ref:`sec:evaluation`. :c:func:`PGAGetStringLength` returns the
string length. :c:func:`PGAGetBinaryAllele` returns the value of the ``i``\ th
bit of string ``p`` in population ``pop``.

.. _sec:simple-example-fortran:

Maxbit Problem in Fortran
~~~~~~~~~~~~~~~~~~~~~~~~~

.. _example:maxbit-fortran:

.. code-block:: fortran
   :caption: PGAPack Fortran Program for the Maxbit Example

         include "pgapackf.h"
         external evaluate
         integer ctx
         ctx = PGACreate  (PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE)
         call  PGASetUp   (ctx)
         call  PGARun     (ctx, evaluate)
         call  PGADestroy (ctx)
         stop
         end

         double precision function evaluate (ctx, p, pop)
         include "pgapackf.h"
         integer ctx, p, pop, i, bit, nbits, stringlen
         stringlen = PGAGetStringLength(ctx)
         nbits     = 0
         do i=1, stringlen
            bit = PGAGetBinaryAllele(ctx, p, pop, i)
            if (bit .eq. 1) then 
               nbits = nbits + 1
            endif
         enddo
         evaluate = dble(nbits)
         return
         end

The Fortran Maxbit problem in Figure :ref:`example:maxbit-fortran`
is similar to the C version above. The Fortran
include file is ``pgapackf.h`` and should be included in every Fortran
function or subroutine that makes PGAPack calls [#]_. Since Fortran
provides no standard mechanism for specifying command line arguments,
these are omitted from the :c:func:`PGACreate` function call. The context
variable, ``ctx``, is declared ``integer`` in Fortran.

The evaluation function ``evaluate`` must contain exactly the calling
sequence shown and must return a ``double precision`` value. The ``aux``
value is optional in fortran because fortran does less strict
type-checking than C and with standard calling conventions the argument
can be omitted if not needed. Note that ``evaluate`` is declared in an
``external`` statement in the program unit in which it is used as an
actual argument. This is a requirement of the Fortran language. In
Fortran, the range of allele values is ``1:stringlen``, rather than
``0:stringlen-1`` as in C.

.. _sec:default-values:

Specifying Nondefault Values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _example:soph-main:

.. code-block:: c
   :caption: Specifying Nondefault values

   #include "pgapack.h"
   double evaluate (PGAContext *ctx, int p, int pop, double *aux);

   int main(int argc, char **argv)
   {
       PGAContext *ctx; 
       ctx = PGACreate     (&argc, argv, PGA_DATATYPE_BINARY, 100, PGA_MAXIMIZE);
       PGASetPopSize       (ctx, 500);
       PGASetFitnessType   (ctx, PGA_FITNESS_RANKING);
       PGASetCrossoverType (ctx, PGA_CROSSOVER_UNIFORM);
       PGASetUp            (ctx);
       PGARun              (ctx, evaluate);
       PGADestroy          (ctx);
       return;
   }

PGAPack offers a wide range of choices for parameter values, operators,
and algorithmic choices. These will be set to default values in
:c:func:`PGASetUp` if the user does not explicitly set a value for them. A
nondefault value may be set by using the :ref:`PGASet <group:init>`
family of calls *after* :c:func:`PGACreate` has been called, but
*before* :c:func:`PGASetUp` has been called.

In Figure :ref:`example:soph-main` the :ref:`PGASet <group:init>`
calls change the default values for population size, fitness
calculation, and crossover type. :c:func:`PGASetPopSize` changes the
population size to 500. :c:func:`PGASetFitnessType` specifies that the fitness
values be determined by using a ranking procedure rather than by direct
use of the evaluation function values. :c:func:`PGASetCrossoverType` specifies
that uniform crossover, rather than the default of two-point crossover
is to be used. Most :ref:`PGASet <group:init>` calls are discussed in
Chapter :ref:`chp:functionality`.

.. _sec:differential-evolution:

Differential Evolution
~~~~~~~~~~~~~~~~~~~~~~

Differential Evolution (DE) is an evolutionary algorithm (EA) invented
by Price and Storn in the 1990’s [SP95]_, [SP97]_, [PSL05]_. It
is used with floating-point genes and uses differences of individuals
(floating-point vectors) which are added to another vector to form a
*donor* vector which is then crossed-over with an existing individual.
The algorithm is described in more detail in
section :ref:`sec:mutation`. Since in PGAPack the DE
algorithm is implemented in a mutation strategy, typically for DE a
strategy with only mutation is selected, see :c:func:`PGASetMixingType` with
option :c:macro:`PGA_MIX_MUTATE_ONLY` in
section :ref:`sec:population-replacement`.

DE applies selection pressure during population replacement: A
newly-mutated string replaces its parent if it has the same or a better
fitness. There is no selection mechanism during the selection phase like
in other EAs. To emulate this (non-) selection, PGAPack introduces a new
selection type, linear selection, enabled with the parameter
:c:macro:`PGA_SELECT_LINEAR` of :c:func:`PGASetPopReplaceType`
which just returns all individuals in
sequence and is no selection operator in the genetic-algorithm sense
because no selection pressure is applied. More details of the selection
operator for DE are given in section :ref:`sec:selection`.

.. _example:de-settings:

.. code-block:: c
   :caption: Specifying Nondefault Values for Differential Evolution

       PGASetPopSize (ctx, 30);
       PGASetNumReplaceValue (ctx, 30);

       PGASetSelectType (ctx, PGA_SELECT_LINEAR);
       PGASetPopReplaceType (ctx, PGA_POPREPL_PAIRWISE_BEST);
       PGASetMixingType (ctx, PGA_MIX_MUTATE_ONLY);
       PGASetMutationType (ctx, PGA_MUTATION_DE);
       PGASetMutationBounceBackFlag (ctx, PGA_TRUE);


For the population replacement strategy the pairwise-best replacement
type is introduced for DE, enabled with the parameter
:c:macro:`PGA_POPREPL_PAIRWISE_BEST` of :c:func:`PGASetPopReplaceType`,
which can also be used in other EA variants
due to the modular nature of PGAPack, more details are given in
section :ref:`sec:population-replacement`. Typical settings for
Differential Evolution are given in
figure :ref:`example:de-settings`. In that example
the population size and the number of individuals replaced in each
generation are set to the same value (with :c:func:`PGASetPopSize` and
:c:func:`PGASetNumReplaceValue`, respectively):
Since DE’s replacement strategy
replaces individuals only if they are better than an existing individual
this strategy is elitist and it makes sense to apply DE to all
individuals in each generation.

.. _sec:parallel-simple-example:

Parallel I/O
~~~~~~~~~~~~

.. _example:parallel-simple-main:

.. code-block:: c
   :caption: PGAPack Maxbit Example in C with I/O

   #include "pgapack.h"
   double evaluate (PGAContext *ctx, int p, int pop, double *aux);

   int main (int argc, char **argv)
   {
        PGAContext *ctx;
        int myid, len, maxiter;

        MPI_Init (&argc, &argv);
        MPI_Comm_rank (MPI_COMM_WORLD, &myid);
        if (myid == 0) {                  /* Process 0 has a dialog */
            printf ("String length? ");   /* with the user and      */
            scanf ("%d", &len);           /* broadcasts the user's  */
            printf ("Max iterations? ");  /* parameters to all      */
            scanf ("%d", &maxiter);       /* other processes        */
        }
        MPI_Bcast (&len,     1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast (&maxiter, 1, MPI_INT, 0, MPI_COMM_WORLD);

        ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, len, PGA_MAXIMIZE);
        PGASetMaxGAIterValue (ctx, maxiter);
        PGASetUp (ctx);
        PGARun (ctx, evaluate);
        PGADestroy (ctx);

        MPI_Finalize ();
        return 0;
   }

The examples in
Figures :ref:`example:parallel-simple-main` (C) and
:ref:`example:parallel-simple-main-f77`
(Fortran) read values for the two parameters ``len`` (string length) and
``maxiter`` (maximum number of GA iterations) from standard input. These
examples will work correctly with either a sequential or parallel
version of PGAPack. However, the explicit use of MPI calls for I/O is
necessary *only* if a parallel version of PGAPack is used, and parameter
values are read from standard input. The purpose is to be sure that each
process receives a copy of the input values. See
Appendix :ref:`app:par-background` for further details.

``MPI_Init(&argc, &argv)`` is always the first function called in any
MPI program. Each process executes
``MPI_Comm_rank(MPI_COMM_WORLD, &myid)`` to determine its unique rank in
the communicator [#]_ ``MPI_COMM_WORLD``. The logic used in this program
is to have process ``0`` read and write from/to standard input/output
and broadcast (using ``MPI_Bcast``) the parameters to the other
processes. The PGAPack function calls are similar to those in the
previous examples. If the user called ``MPI_Init``, the user must also
call ``MPI_Finalize`` before exiting.

We elaborate here on the ``MPI_Bcast`` function because of its practical
value in the model of parallel I/O shown. For more detailed discussion
of MPI concepts and functions, the user should consult
[MPI94]_, [GLS94]_ or a newer MPI reference [MPI21]_.

The C binding for ``MPI_Bcast`` is

.. code-block:: c

   int MPI_Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm)

and the Fortran binding

.. code-block:: fortran

   MPI_BCAST(buffer, count, datatype, root, comm, ierror)
   <type> buffer(*)
   integer count, datatype, root, comm, ierror

``MPI_Bcast`` will result in every process in communicator ``comm``
receiving a copy of the contents of ``buf``/``buffer``. The other
parameters are the number of items (``count``), the datatype
(``datatype``), which may be one of ``MPI_DOUBLE``, ``MPI_INT``,
``MPI_CHAR``, ``MPI_UNSIGNED``, or ``MPI_LONG``; the rank of the process
with the original copy (``root``); the MPI communicator (``comm``); and,
for Fortran, a variable to handle an error return code (``ierror``).

.. _example:parallel-simple-main-f77:

.. code-block:: fortran
   :caption: PGAPack Maxbit Example in Fortran with I/O

         include 'pgapackf.h'
         include 'mpif.h'

         double precision evaluate
         external         evaluate

         integer(8) ctx
         integer myid, len, maxiter, ierror

         call MPI_Init(ierror)
         call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierror)

         if (myid .eq. 0) then
            print *, 'String length?'
            read  *, len
            print *, 'Max iterations?'
            read  *, maxiter
         endif
         call MPI_Bcast(len,     1, MPI_INT, 0, MPI_COMM_WORLD, ierror)
         call MPI_Bcast(maxiter, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)

         ctx = PGACreate(PGA_DATATYPE_BINARY, len, PGA_MAXIMIZE)
         call PGASetMaxGAIterValue(ctx, maxiter)
         call PGASetUp(ctx)
         call PGARun(ctx, evaluate)
         call PGADestroy(ctx)

         call MPI_Finalize(ierror)

         stop
         end

.. _sec:execute:

Compiling, Linking, and Execution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When PGAPack was installed, the makefiles in the ``./examples/c`` and
``./examples/fortran`` directories were correctly configured for the
machine PGAPack was installed on using the version of MPI specified (if
any). To run your own programs, it is best to copy the appropriate
makefile (C or Fortran) to your directory and modify it to use your
source code files. The makefile will compile your source code files,
link in the PGAPack library (and MPI library if a parallel version of
PGAPack was built), and build your executable.

How you execute your program will depend on whether a sequential or
parallel version of PGAPack was built, the MPI implementation used and
the machine you are running on. If a sequential version of PGAPack was
built (i.e., one where the user did not supply a version of MPI), the
executable ``maxbit`` can be executed on a uniprocessor Unix system by
typing ``maxbit``. If the ``MPICH`` implementation of MPI was used, it
may be executed (using four processes) by ``mpirun maxbit -np 4``.
Appendix :ref:`chp:start-up` contains some examples.

.. [#]
   Since not all Fortran compilers support the ``-I`` mechanism for
   specifying the include file search path, you may need to copy or set
   up a symbolic link to ``pgapackf.h`` from the directory you are
   compiling a Fortran program in.
.. [#] See Appendix :ref:`app:par-background`
