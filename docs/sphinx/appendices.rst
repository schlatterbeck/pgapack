.. _part:appendices:

Appendixes
==========

.. _app:default-values:

Default Values
--------------

.. _tab:default-values:

.. table:: PGAPack Default Values

 +----------------------+--------------------+-----------------------------------------+
 | Population size      | 100                | :c:func:`PGASetPopSize`                 |
 +----------------------+--------------------+-----------------------------------------+
 | Maximum iterations   | 1000               | :c:func:`PGASetMaxGAIterValue`          |
 +----------------------+--------------------+-----------------------------------------+
 | Maximum no change    | 100                | :c:func:`PGASetMaxNoChangeValue`        |
 | iters                |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Max. population      | 95                 | :c:func:`PGASetMaxSimilarityValue`      |
 | homogeneity before   |                    |                                         |
 | stopping             |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Number of new        | 10% of population  | :c:func:`PGASetNumReplaceValue`         |
 | strings to generate  | size (10)          |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Probability of       | 0.85               | :c:func:`PGASetCrossoverProb`           |
 | crossover            |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Uniform crossover    | 0.6                | :c:func:`PGASetUniformCrossoverProb`    |
 | bias                 |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Mutation probability | 1/L                | :c:func:`PGASetMutationProb`            |
 +----------------------+--------------------+-----------------------------------------+
 | Real mutation        | 0.1 or 0.01        | :c:func:`PGASetMutationRealValue`       |
 | constant             |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Integer mutation     | 1                  | :c:func:`PGASetMutationIntegerValue`    |
 | constant             |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Mutation range       | ``PGA_FALSE``      | :c:func:`PGASetMutationBoundedFlag`     |
 | bounded              |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Probabilistic binary | 0.6                | :c:func:`PGASetPTournamentProb`         |
 | tournament parameter |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Use restart operator | ``PGA_FALSE``      | :c:func:`PGASetRestartFlag`             |
 +----------------------+--------------------+-----------------------------------------+
 | Restart frequency    | 50                 | :c:func:`PGASetRestartFrequencyValue`   |
 +----------------------+--------------------+-----------------------------------------+
 | Restart allele       | 0.5                | :c:func:`PGASetRestartAlleleChangeProb` |
 | mutation rate        |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Allow no duplicate   | ``PGA_FALSE``      | :c:func:`PGASetNoDuplicatesFlag`        |
 | strings              |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Multiplier for       | 1.01               | :c:func:`PGASetFitnessCmaxValue`        |
 | minimization         |                    |                                         |
 | problems             |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Parameter MAX in     | 1.2                | :c:func:`PGASetMaxFitnessRank`          |
 | fitness by ranking   |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Frequency of         | 10                 | :c:func:`PGASetPrintFrequencyValue`     |
 | statistics printing  |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Print strings        | ``PGA_FALSE``      | :c:func:`PGASetPrintOptions`            |
 +----------------------+--------------------+-----------------------------------------+
 | Print offline        | ``PGA_FALSE``      | :c:func:`PGASetPrintOptions`            |
 | statistics           |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Print online         | ``PGA_FALSE``      | :c:func:`PGASetPrintOptions`            |
 | statistics           |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Print best string    | ``PGA_FALSE``      | :c:func:`PGASetPrintOptions`            |
 +----------------------+--------------------+-----------------------------------------+
 | Print worst string   | ``PGA_FALSE``      | :c:func:`PGASetPrintOptions`            |
 +----------------------+--------------------+-----------------------------------------+
 | Print genetic        | ``PGA_FALSE``      | :c:func:`PGASetPrintOptions`            |
 | distance             |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Randomly initialize  | ``PGA_TRUE``       | :c:func:`PGASetRandomInitFlag`          |
 | population           |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Probability of       | 0.5                | :c:func:`PGASetBinaryInitProb`          |
 | initializing a bit   |                    |                                         |
 | to one               |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | How to initialize    | Range              | :c:func:`PGASetRealInitRange`           |
 | real strings         |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Real initialization  | :math:`[0,1]`      | :c:func:`PGASetRealInitRange`           |
 | range                |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | How to initialize    | Permutation        | :c:func:`PGASetIntegerInitPermute`      |
 | integer strings      |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Integer              | :math:`[0,L-1]`    | :c:func:`PGASetIntegerInitPermute`      |
 | initialization range |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Seed random number   | ``PGA_TRUE``       | :c:func:`PGASetRandomSeed`              |
 | with clock           |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | Default MPI          | ``MPI_COMM_WORLD`` | :c:func:`PGASetCommunicator`            |
 | communicator         |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | DE scale factor      | 0.9                | :c:func:`PGASetDEScaleFactor`           |
 | :math:`F`            |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | DE auxiliary factor  | :math:`0.5 \cdot   | :c:func:`PGASetDEAuxFactor`             |
 | :math:`K`            | (F + 1)`           | `                                       |
 +----------------------+--------------------+-----------------------------------------+
 | DE Crossover prob    | 0.9                | :c:func:`PGASetDECrossoverProb`         |
 | :math:`Cr`           |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | DE Dither            | 0.0                | :c:func:`PGASetDEDither`                |
 +----------------------+--------------------+-----------------------------------------+
 | DE Jitter            | 0.0                | :c:func:`PGASetDEJitter`                |
 +----------------------+--------------------+-----------------------------------------+
 | DE Either/Or         | 0.5                | :c:func:`PGASetDEProbabilityEO`         |
 | Probability          |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+
 | DE Number of         | 1                  | :c:func:`PGASetDENumDiffs`              |
 | differences          |                    |                                         |
 +----------------------+--------------------+-----------------------------------------+

.. _chp:function-bindings:

Function Bindings
-----------------

.. _app:symbolic-constants:

Symbolic Constants
~~~~~~~~~~~~~~~~~~

PGAPack defines many symbolic constants that are used as arguments to
PGAPack functions. These constants are the same for both Fortran and C.
The constants are listed in section :ref:`sec:constant-definitions`
These constants are the same for both Fortran and C. The constant groups
also define the default values.

.. _app:bindings-c:

C Bindings
~~~~~~~~~~

See :ref:`sec:function-group-standard` for functions of the top-level
API using :c:func:`PGARun`, for :ref:`chp:explicit` you want to consult
the :ref:`sec:function-group-explicit` and if looking at the internal
implementation the function calls are documented in
:ref:`sec:function-group-internal`.

.. _app:bindings-fortran:

Fortran 77 Bindings
~~~~~~~~~~~~~~~~~~~

Use the rules defined in Chapter :ref:`chp:fortran` (and the
machine-specific idiosyncrasies noted in Appendix :ref:`chp:start-up`)
to determine the Fortran bindings.

.. _app:par-background:

Parallelism Background
----------------------

Parallel Computer Taxonomy
~~~~~~~~~~~~~~~~~~~~~~~~~~

Traditionally, parallel computers are classified according to Flynn’s
taxonomy [Fly72]_. Flynn’s classification distinguishes
parallel computers according to the number of instruction streams and
data operands being computed on simultaneously.

Flynn’s single-instruction single-data (SISD) model is the traditional
sequential computer. A single program counter fetches instructions from
memory. The instructions are executed on *scalar* operands. There is no
parallelism in this model.

In the single-instruction multiple-data (SIMD) model there is again a
single program counter fetching instructions from memory. However, now
the operands of the instructions can be one of two types: either scalar
or array. If the instruction calls for execution involving only scalar
operands, it is executed by the control processor (i.e., the central
processing unit fetching instructions from memory). If, on the other
hand, the instruction calls for execution using array operands, it is
broadcast to the *array* of processing elements. The processing elements
are separate computing devices that rely upon the control processor to
determine the instructions they will execute.

In a multiple-instruction multiple-data (MIMD) computer there exist
multiple processors each of which has its own program counter.
Processors execute independently of each other according to whatever
instruction the program counter points to next. MIMD computers are
usually further subdivided according to whether the processors share
memory or each has its own memory.

In a shared-memory MIMD computer both the program’s instructions and the
part of the program’s data to be shared exist within a single shared
memory. Additionally, some data may be private to a processor and not be
globally accessible by other processors. The processors execute
asynchronously of each other. Communication and synchronization between
the processors are handled by having them each read or write a
shared-memory location.

A distributed-memory MIMD computer consists of multiple “nodes.” A node
consists of a processor, its own memory, a network interface, and
sometimes a local disk. The program instructions and data reside in the
node’s memory. The nodes are connected via some type of network that
allows them to communicate with each other. Parallelism is achieved by
having each processor compute simultaneously on the data in its own
memory. Communication and synchronization are handled by passing of
messages (a destination node address and the local data to be sent) over
the interconnection network.

Processes vs. Processors
~~~~~~~~~~~~~~~~~~~~~~~~

We distinguish the two terms process and processor. A *process* is a
software abstraction with a unique address space that can be scheduled
by the operating system. A *processor* is the physical computer hardware
on which computations take place.

On MIMD parallel computers, usually one process executes on each
processor (although this is not required). On a uniprocessor, multiple
processes timeshare the single processor.

Message-Passing Programming Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the message-passing programming model multiple processes communicate
by passing messages—transferring data from the address space of one
process into the address space of another process. When a process needs
data stored in the memory of another process, the data must be sent from
the process that “owns” it, to the memory of the process that needs it.

The message-passing programming model is currently one of the most
popular. One reason for the popularity is portability. Message passing
is the natural programming model on distributed-memory MIMD computers.
Each process is naturally mapped to one of the machine’s nodes. A
similar implementation is common on a workstation network where one
process runs on each workstation. On a shared-memory MIMD computer
multiple processes can emulate message passing by communicating only via
logical message queues—areas of shared memory partitioned by process. On
a uniprocessor the multiple processes that timeshare the physical
processor can also emulate the idea of logical message queues for their
communication.

One example of the message-passing programming model is the
controller/responder model. In this model a *controller* process
distributed work (computation to be performed) to the responder
processes.  The responders perform the work and return the result to the
controller. In many implementations the controller plays a bookkeeping
role only and does not perform any computation.

Parallel Genetic Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using the term “parallel genetic algorithm” it is important to
distinguish between parallel models, their (parallel or sequential)
implementation, and the computer hardware.

Models
^^^^^^

A sequential GA model (more accurately called a *global* model) has a
single population and no restrictions (partitioning) upon which strings
recombine with which. The sequential GA is the traditional GA model
given in the literature. In a parallel GA model there are either
multiple populations (an island model) or a partitioning of a single
population (often called a fine-grained model).

Implementations
^^^^^^^^^^^^^^^

Both parallel and sequential GA models can have parallel or sequential
implementations. A sequential implementation of the global model is the
traditional approach discussed in the GA literature. One process,
running on a uniprocessor (PCs and workstations), performs all the
calculations. In a parallel implementation of the global model the steps
of the GA (some or all of selection, crossover, mutation, and fitness
calculation) are executed simultaneously by multiple processes running
on a parallel computer or workstation network.

In a sequential implementation of a parallel GA model, multiple
processes, each responsible for a subpopulation or partition of the full
population, time share the processor of the uniprocessor they are
running on. In a parallel implementation of a parallel GA model, the
multiple processes each run on a unique processor of a parallel computer
or workstation network.

MPI
~~~

MPI (Message Passing Interface) is a *specification* of a
message-passing library for parallel computers and workstation
networks—it defines a set of functions and their behavior. The actual
*implementation* of this interface is left up to vendors and researchers
to develop. At the time of this writing several implementations of MPI,
both proprietary and freely available, exist. MPI was designed by a
large group of parallel computer vendors, computer researchers, and
application developers as a standard for message passing.

Communicators
^^^^^^^^^^^^^

Almost all MPI functions require a *communicator*. If MPI routines are
called directly, the user must supply a communicator. Also, if any of
PGAPack’s parallel routines, other than :c:func:`PGARun`, are used, the user
must supply a communicator as well.

A communicator combines the notions of context and group. A *context* is
an extension of the notion of a “tag” used in many other message-passing
systems to identify a message. Contexts differ from tags in that they
are allocated by the system, not the user, and that no wild-card
matching among contexts is allowed. A *group* contains :math:`n`
processes whose *rank* is an integer between :math:`0,\ldots,n-1`.
Processes may belong to more than one group and have a unique rank
within each.

Any MPI implementation will always supply the default communicator
``MPI_COMM_WORLD``. This communicator contains all processes that were
created when MPI was initialized. For most users this is all they have
to know about communicators. Simply supply ``MPI_COMM_WORLD`` whenever a
communicator is required as an argument. For more sophisticated use,
users are referred to [MPI94]_, [GLS94]_, [MPI21]_.

Parallel I/O
^^^^^^^^^^^^

The issue of parallel I/O is independent of PGAPack. However, since we
assume many PGAPack users will wish to do I/O we address this topic. A
primary consideration has to do with whether one or all processors do
I/O. Consider the following two code fragments, keeping in mind that
they are being executed simultaneously by *multiple* processes:

.. code-block:: c

   ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 30, PGA_MINIMIZE) 

and

.. code-block:: c

   int len;
   scanf ("%d",&len);
   ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, len, PGA_MINIMIZE);

In the first case, all processes will receive the value of 30 for the
string length since it is a constant. In the second case, however, the
value of the string length is determined at run time. Whether one or all
processes execute the ``scanf`` function is undefined in MPI and depends
on the particular parallel computing environment. In PGAPack we require
that all processes have the same values for all fields in the context
variable. Since some of these fields may be set by using values
specified at run time, we suggest that your I/O that reads in
PGAPack parameters be done as follows:

.. code-block:: c

   #include "pgapack.h"
   double evaluate (PGAContext *ctx, int p, int pop, double *aux);

   int main( int argc, char **argv )
   {
        PGAContext *ctx;
        int myid, len;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        if (myid == 0) {                        /* Process 0 has a dialog */
            printf("String length? ");          /* with the user and      */
            scanf("%d", &len);                  /* broadcasts the user's  */
        }
        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

        ctx = PGACreate(&argc, argv, PGA_DATATYPE_BINARY, len, PGA_MAXIMIZE);
        PGASetUp(ctx);
        PGARun(ctx, evaluate);
        PGADestroy(ctx);

        MPI_Finalize();
        return(0);
   }

The key point is that *only* process 0 (as determined by
``MPI_Comm_rank``) performs I/O and the value of ``len`` is then
broadcast (using ``MPI_Bcast``) to the other processes.

.. _chp:start-up:

Machine Idiosyncrasies
----------------------

Data Type Sizes
~~~~~~~~~~~~~~~

PGAPack is written entirely in ANSI C. However, because it is callable
from Fortran, and no standards exist for interlanguage communication,
problems may arise. These have to do with a lack of consistency in the
size of data types between the two languages.

On all machines we have tested, an ``integer`` declaration in Fortran is
the same size as an ``int`` declaration in C and everything works
properly. For floating-point numbers, however, we have found at least
one inconsistency. The requirement is for the Fortran floating-point
number to be the same size as a C ``double``. On most machines a Fortran
``double precision`` declaration is the equivalent size.

Since Fortran does not support pointers, an ``integer`` variable is used
to hold the address of the context variable (and possibly MPI
communicator addresses as well). Therefore, a Fortran ``integer`` must
be “large enough” to hold an address on the machine. For all 32-bit
address space machines we have tested this is the case. On machines with
a 64-bit address space, however, this may not be true. Therefore we use
constructs in Fortran to select an integer data type that is large
enough, see chapter :ref:`chp:fortran` for details.

Startup
~~~~~~~

The MPI standard provides for *source code* portability. However, the
MPI standard does *not* specify how an MPI program shall be started or
how the number of processes in the computation is specified. These will
vary according to the computer being used and the choice of MPI
implementation. This section used to have documentation about a lot of
machines that no longer exist today. We refer you to the documentation
of OpenMPI [OMPI23]_ or MPICH [MPIC23]_ or the documentation of whatever
MPI implementation you are using.

.. _chp:problems:

Common Problems
---------------

This collects some problems seen over the years, some may be specific to
MPI versions or variants that are no longer in use, since it is hard to
know what is still relevant all information has been left in.

-  When reading input value to be used as parameters in
   :ref:`PGASet <group:init>` calls, the :ref:`PGAset <group:init>`
   calls themselves may not be executed until *after*
   :c:func:`PGACreate` has been called.

-  In C, when reading input parameters which are of type ``double``, the
   ``scanf`` conversion specification should be of the form ``%lf``,
   *not* ``%f`` which is appropriate for ``float``\ s.

-  An infinite loop can occur if the number of permutations of the bit
   string is less than the population size. For example, for a
   binary-valued string of length four, there are :math:`2^4 = 16`
   possibilities. If the population size is greater than 16, and
   duplicate strings are not allowed in the population, an infinite loop
   will occur.

-  Erroneous results can occur if the name of a user’s function
   conflicts with a library function used by PGAPack. For example, if a
   program defined its own ``ceil`` function, this would conflict with
   the C math library function of the same name.

-  All floating point constants and variables used in PGAPack are of
   type ``double``. Particularly from Fortran, the user should be
   careful to make sure that they pass a ``double precision`` constant
   or variable.

-  :c:func:`PGACreate` removes command line arguments. One consequence is that
   if :c:func:`PGACreate` is called twice in the same program (unusual, but
   legal), the second :c:func:`PGACreate` call will *not* receive the
   command-line arguments.

-  If one includes ``mpi.h`` (or ``mpif.h``) when it should not be,
   errors will result, as well as warnings about redefining macros and
   typedefs. This usually happens when a sequential version of
   PGAPack is used (with “fake” MPI stub routines and definitions) and
   the user’s program explicitly includes “real” ``mpi.h`` or ``mpif.h``
   header files.

-  If one fails to include ``mpi.h`` (or ``mpif.h``) when it should be
   (such as calling MPI functions directly) errors may result. Since
   ``pgapack.h`` includes ``mpi.h`` this should not happen in C. The
   Fortran include file, ``pgapackf.h``, however, does *not* include
   ``mpif.h``. The user must explicitly include it in every subroutine
   and function that makes MPI calls. Not including ``mpif.h`` could
   result in any of several different errors, including

   -  syntax errors when compiling (for example, ``MPI_COMM_WORLD``
      being undefined)

   -  general errors in the computed results

   -  the program crashing when it calls the undefined subroutine
      ``MPI_Init``

   -  general MPI errors such as:

      ::

             0 - Error in MPI_COMM_RANK : Invalid communicator
             [0] Aborting program!

   We have also seen the following error from not including ``bmpif.h``
   in the main program:

   .. code-block:: none

      PGACreate: Invalid value of datatype: 0
      PGAError: Fatal

-  If the ``ch_p4`` device in ``MPICH`` is used to run on workstations
   one must have a correct processor group file (``procgroup``). The
   error message

   .. code-block:: none

      (ptera-36%)a.out
      p0_18429:  p4_error: open error on procgroup file (procgroup): 0
      (ptera-37%)

   may occur if the processor group file is not specified correctly. See
   the ``MPICH`` users guide for more details.

-  A common error with the ``procgroup`` file when using the ``ch_p4``
   device in ``MPICH`` is to have an incorrect path to the executable.

-  When compiling the ``examples`` directory we have seen “multiply
   defined” error messages. For example:

   .. code-block:: none

      Making C examples
        Compiling classic
      ld: /usr/local/mpi/lib/sun4/ch_p4/libmpi.a(initialize.o): _MPI_Initialized: multiply defined
      collect2: ld returned 2 exit status

   We have seen this error occur when a sequential version of
   PGAPack was built and the library (``./lib/arch/libpgag.a`` or
   ``./lib/arch/libpgaO.a``) was not deleted before attempting to build
   a new, parallel version of PGAPack. The “fake” MPI stub routines are
   in the sequential library and have name conflicts when a “real” MPI
   library is referenced. The solution is to delete the old ``.a`` file
   and rerun ``make install``. The ``Makefile`` target ``clobber`` takes
   care of deleting all exiting libraries::

     make clobber
