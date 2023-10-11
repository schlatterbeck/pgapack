.. |--| unicode:: U+2013   .. en dash

.. |examples/c/namefull.c| replace:: ``examples/c/namefull.c``
.. |examples/c/udtstr.c| replace:: ``examples/c/udtstr.c``
.. |examples/fortran/namefull.f| replace:: ``examples/fortran/namefull.f``
.. |examples/fortran/constr.f| replace:: ``examples/fortran/constr.f``
.. |examples/nsgaii/crowdingplot| replace:: ``examples/nsgaii/crowdingplot``

PGAPack
+++++++

PGAPack is a general-purpose, data-structure-neutral, parallel genetic
algorithm library originally developed by David Levine at Argonne
National Laboratory. It has libraries for C and Fortran. There are
companion projects:

- PGAPy_, a Python wrapper for PGApack https://github.com/schlatterbeck/pgapy
- `debian-pgapack`_, a project for building debian packages (for all MPI
  libraries supported by Debian):
  https://github.com/schlatterbeck/debian-pgapack

Documentation is on `Read the Docs`_.

Updates
=======

2nd update Oct 2023:

- Implement Negative Assortative Mating, a mating restriction that tries
  to mate individuals with large genetic distance.

Oct 2023:

- Add Differential Evolution with integer genes
- Fix feature interaction bug with population replacement RTR or
  population replacement pairwise best and the no duplicates flag.

Apr 2023:

- Add MPI_Abort to the fake mpi wrapper
- Add missing MPI_Finalize() prototype to fakempi include

3rd update Jan 2023:

- Now the user-guide is converted to sphinx. The old LaTeX version can
  still be built but will vanish eventually
- You can access the documentation on `Read the Docs`_ now.

2nd update Jan 2023:

- Generalize mean hamming distance reporting to mean genetic distance
  reporting. This now works with *all* data types not just binary. This
  uses the already-existing genetic distance user function which is
  implemented for all builtin data types and can be overridden to use
  euclidian distance instead of the default manhattan distance for
  integer and real data types.
- Deprecated PGAHammingDistance in favor of PGAGeneDistance. There is a
  backward-compatible define for C but not for Fortran. It's doubful
  anybody has ever used this in custom code (it was mainly used in
  built-in reporting of the hamming distance of binary strings).

Update Jan 2023:

- Add Sphinx documentation
- Generate manual pages from Sphinx docs, this has a lot of bug-fixes in
  the manpages (e.g. wrong documentation of return value) and documents
  some functions that did previously not have a manual page
- Default for PGASetSumConstraintsFlag is now PGA_TRUE, works more
  reliable, this is only relevant when using one of the NSGA multi
  objective algorithms

Update Dec 2022:

- Bug fixes discovered during implementation of a regression test for
  the python wrapper
- Output of string printing and result printing can now be redirected to
  a file

Update Oct 2022:

- Use hashing for comparing individuals when the NoDuplicates flag is
  turned on. Previously the comparison operations were O(n²) in the
  population size. Now the effort is linear. The downside is that when
  you have a custom comparison function with PGA_USERFUNCTION_DUPLICATE
  you also need a hash function with PGA_USERFUNCTION_HASH. This always
  applies when you have user-defined datatypes (and want to use the
  NoDuplicates flag). Examples for C are in |examples/c/namefull.c|_ and
  in |examples/c/udtstr.c|_ (for a user defined datatype) and for Fortran
  in |examples/fortran/namefull.f|_ (there are no user-defined datatypes
  in Fortran). The good news is that there is a utility function
  ``PGAUtilHash`` that you can use when implementing a custom hashing
  function.
- Factor MPI serialization: When serializing an Individual, pgapack needs
  certain fields to be sent together with the gene, this is now in its
  own function. This function can also be used for MPI serialization
  when using user-defined datatypes with MPI: In that case a user
  function PGA_USERFUNCTION_BUILDDATATYPE has to be written. The new
  code substantially reduces the boilerplate code for writing a function
  for building an MPI datatype and will probably need no updates when
  the pgapack internal information changes. An example is in
  |examples/c/udtstr.c|_ and in the user guide.
- Add a new serialization API for MPI serialization that can be used
  instead of PGA_USERFUNCTION_BUILDDATATYPE. This is especially useful
  when the user-defined datatype is variable length. We send the length
  of the serialization in a first MPI message before sending the
  (variable length) individual. Since we're not using multicast, this
  works fine for transferring variable-size information with MPI.
  This new API will be used in the companion-project PGAPy_ for user
  defined datatypes in python.
- Bug-fix in multi-objective optimization: When evaluations are exactly
  equal the ranking would not correctly compute the dominance relation
- Bug-fix in multi-objective optimization: The crowding metric was not
  properly initialized resulting sometimes in different optimization
  paths when compiled with/without optimization (-O2 and -O3 in gcc)
- Fix feature interaction between multi-objective optimization and the
  NoDuplicates flag: When combining two populations in the multi
  objective optimization algorithms (NSGA-II and NSGA-III) where both
  populations contain instances of the same indidivual, duplicates would
  result.
- The script for plotting the pareto front for 3-dimensional problems
  used to be in ``examples/nsgaiii/crowdingplot3`` (this was already a
  symbolic link in the latest releases) and is now gone, use the ``-3``
  option for |examples/nsgaii/crowdingplot|_.

The bug-fixes in multi-objective optization will result in different
optimization paths being taken compared to previous versions (because of
different sorting).

Update Aug 2022:

- Add a crossover method for permutations (e.g. traveling salesman)
- Add Epsilon-Constrained optimization, see `blogpost on epsilon
  constrained optimization`_

Update Mar 2022:

- Attempt to get everything compiled with visual studio compiler. This
  compiler is stuck in the 1990s of the last millenium because it does
  not support dynamically sized arrays on the stack. This is part of the
  C99 standard. The workarounds involve some ugly macros.
- Bug-Fix in the genetic distance function PGARealGeneDistance which is
  used for RTR population replacement. This converted the distances to
  int which is wrong.

Update Jan 2022:

- Now the tournament size can be a floating-point value implementing
  fine-grained tournaments (the fractional part is used to add an
  additional tournament participant probabilistically). See the
  Selection chapter in the user guide and the citations on the topic.
  For Fortran this could mean changes to the constant passed to
  PGASetTournamentSize.
- Implement simulated binary crossover (SBX) and polynomial mutation,
  see user guide.
- NSGA-III for many-objective optimization is now implemented
- There is a small plotting-utility ``examples/nsgaiii/crowdingplot3``
  written in python that can plot three function values in a 3D-plot.
  It can directly use the output of an optimization, e.g.::

    examples/nsgaiii/crowdingplot3 test/nsgaiii_optimize_13_1.data

Second Update December 2021:

- Now the multiobjective optimization algorithm NSGA-II (Nondominated
  Sorting Genetic Algorithm) by Deb et. al. is implemented. Like for
  constrained optimization this uses multiple objective functions.
- There are examples from the original paper (see README.rst) in the
  directory ``examples/nsgaii``, both with and without constraints.
- Note that multiobjective optimization is considered experimental:
  There are interaction with other parts of the API of the library,
  e.g., functions dealing with the *best* evaluation like
  ``PGAGetBestIndex`` currently no longer have a valid semantic
  interpretation with multiobjective optimization, they sort by
  nondominance-rank now. And reporting has been rewritten to provide a
  meaningful output, in particular the optimization result prints all
  non-dominated solutions.
- A Fortran example with constraints *and* multi-objective optimization
  can be found in |examples/fortran/constr.f|_
- There is a small plotting-utility |examples/nsgaii/crowdingplot|_
  written in python that can plot one function value (in the objective
  space) against a second function value, similar to the graphics in the
  NSGA-II paper.
- You also want to check the next section for news.

First Update December 2021:

- If you're upgrading: The signature of your evaluation function has
  changed, it has grown a new parameter at the end. If you're not using
  constrained optimization you will only have to change your objective
  function to add this parameter, it is unused in that configuration.
  In Fortran you can get away without any changes.
- This release probably changes the path an optimization takes because we
  use a new (stable) sort for sorting populations during copying of
  individuals for elitist algorithms. This can result in different
  individuals being copied (which have the same evaluation but might have
  different genetic material).
- Add auxiliary evaluations, currently only used for constrained
  optimization from a paper by Deb, 2000 (see user guide for citation).
  To find out about the new feature see the user guide, section 4.9
  "String Evaluation and Fitness". You may also want to look at the
  examples in examples/deb.
- Fixes for Fortran on 64-bit machines: The context variable is a
  pointer that didn't fit into a 4-byte integer on these machines
  resulting in a core-dump.
- Regression tests that use the alreay-coded examples as tests, this
  includes the Fortran examples.
  You can run them with "make test". Or, e.g., "make MPI=openmpi test"
  The default for MPI is to run with 4 processors and use the machine
  file .mpi-${MPI}-machinefile in your home directory (${MPI} is replaced
  by the mpi implementation given to the make command, openmpi in this
  example).
- New examples for constrained optimization using all the examples from
  Deb 2000.
- Tested MPI on a multiprocessor machine (a bunch of Orange-Pi computers
  acting as a (slow :-) multiprocessor). Works fine with Debian's
  OpenMPI and MPICH MPI implementations. Does not work for me with LAM,
  there is a debian bug-report #1000446.

Updated September 2020:

- Add Differential Evolution (DE) as a new Mutation Strategy
- Add more options to fully emulate Differential Evolution
- Update Docs for DE

Updated May 2020:

- Add Tournament Selection *without* replacement as an option
- Add Truncation Selection
- Update Documentation and manual pages

Updated March 2020:

- Add restricted tournament replacement, see updated user guide for
  details and references
- Fix some compiler warnings
- Implement Tournament Selection with more than 2 individuals, new
  parameter settable with ``PGASetTournamentSize``, the default is the old
  default of 2.

Updated Sept 2017: new installation instructions, availability:

- Bug fixes in MPI code: Now compiles against all MPI implementations
  shipped with Debian Linux (openmpi, mpich, lam).
- Bug fix in ``PGAChange`` that did not call ``PGASetEvaluationUpToDateFlag``:
  This would result in occasional wrong evaluation of individuals,
  noteably the evaluation went *down* even with an elitist strategy.
- Bug fix for restart with an integer gene: According to the user guide
  this should use ``PGA_MUTATION_CONSTANT`` but tried to use
  ``PGA_MUTATION_UNIFORM`` which is undefined for integer genes
- Fixes to the user guide with new documentation, the old original
  postscript is still available. Notably documentation bugs reported via
  the debian project were fixed. The user guide can be built from source
  again (after probably a *very* long time).
- Make Fortran compile again

Updated March 2008:

- PGAPack has also been built successfully against LAM/MPI and Open MPI.

Copyright
=========

See the file COPYING for Copyright and disclaimer information.

Introduction
============

PGAPack is a general-purpose, data-structure-neutral, parallel genetic
algorithm library developed at Argonne National Laboratory.
Key features are:

- Callable from Fortran or C.
- Runs on uniprocessors, parallel computers, and workstation networks.
- Binary-, integer-, real-, and character-valued native data types.
- Object-oriented data structure neutral design.
- Parameterized population replacement.
- Multiple choices for selection, crossover, and mutation operators.
- An implementation of Differential Evolution
- Optimization with constraints
- Epsilon-constrained optimization
- Multi-objective optimization with NSGA-II
- Many-objective optimization with NSGA-III
- Easy integration of hill-climbing heuristics.
- Easy-to-use interface for novice and application users.
- Fully extensible to support custom operators and new data types.
- Extensive debugging facilities.
- A large set of example problems.
- It is released under the MPICH2 license (also used by the MPICH2 MPI
  implementation from Argonne National Laboratory).
- A separate package with Python bindings PGAPy_


Availability
============

PGAPack is freely available.

The latest version can be obtained from github at
https://github.com/schlatterbeck/pgapack

The distribution contains all source code, installation instructions,
users guide, and a collection of examples in C and Fortran.

Older versions of the distribution are still available by anonymous ftp
from ftp://ftp.mcs.anl.gov/pub/pgapack

Note that the github project contains all older releases in the git
repo.


Computational Environment
=========================

PGAPack is written in ANSI C and uses the MPI message passing interface
and should run on most uniprocessors, parallel computers, and workstation
networks.  PGAPack has been tested on the workstations and parallel computers
specified by the ARCH_TYPE variable below.

Documentation
=============

* Documentation is now on `Read the Docs`_.
* The PGAPack users guide which used to be in LaTeX is now converted to
  Sphinx with cross-links to a reference documentation.
* The old LaTeX version is still available in the directory ``docs`` but
  no longer built by default. The ancient original documentation is
  still preserved as ``docs/user_guide-orig.ps`` for historical reasons.
  It is not recommended for a reference.
* Man pages for PGAPack functions are in the ``./man`` directory. They
  are created automatically from the Sphinx documentation in
  ``docs/sphinx`` using some postprocessing from the manual page export
  of Sphinx. But the man-pages are still checked into git and only
  rebuilt when something changes. The reason is that the manpages should
  be easily installable.
* For building the man page sources a Sphinx setup is needed, see below in
  `Building the documentation`_.
* Installation instructions are in this ``README.rst`` file.
* Example problems are in the ``./examples`` directory.

Building the documentation
--------------------------

To build the Sphinx documentation you should install into a `Sphinx
virtual environment`_: This uses a Python virtual environment and
installs Sphinx and all the necessary addons into this environment.
In addition to Sphinx proper you also need the additional packages in
``docs/sphinx/requirements.txt``. You can install with::

 pip install -r docs/sphinx/requirements.txt

But be sure that you have activated the virtual environment before
issuing this command, otherwise you install into the global python
interpreter or your user configuration.

You also need install ``doxygen``, ``latexmk``, ``texlive-latex-extra``,
``inkscape`` for pdf file generation, on a Debian-based system (applies
also to Ubuntu) you can achieve this with::
  
  sudo apt install doxygen latexmk texlive-latex-extra inkscape

After this you can change to ``docs/sphinx`` directory and build the
html documentation with::

 make html

Alternatively you can build manual pages with the target ``fixedman``
and a pdf file with the target ``latexpdf``. The default if no target is
given is to build all three. The latter can also be achieved by::

 make documentation

from the top-level. Note that you need to have the sphinx virtual
environment activated for this to work. This is also the reason why the
documentation is no longer built by default with the default make target
from the top-level Makefile.

Currently the Sphinx documentation uses some hacks by modifying
subprograms in memory while building the documentation. The Python
community calls this `monkey patching`_. This is because exhale
hard-codes some of the section headings in the documentation and I did
not want to have 'Classes' when the code is in C which doesn't have
classes. And I like the functions in the function groups sorted by name
which originally was not supported by breathe but a patch from me has
been accepted and I expect this to be available in a future version.
In short this means that you may be unable to build the documentation
when a new version comes along. Please open a bug report on github if
this occurs to you.


Installation Requirements
=========================

To compile you must have an ANSI C compiler that includes a full
implementation of the Standard C library and related header files.  To build a
*parallel* version of PGAPack you must provide an implementation of MPI
(Message Passing Interface) for the parallel computer or workstation network
you are running on.

Most of our testing and development was done using MPICH, a freely available
implementation of MPI.  MPICH runs on many parallel computers and
workstation networks and is publicly available and free.  The complete
distribution is available by anonymous ftp from ftp://ftp.mcs.anl.gov.
Take the file ``mpich.tar.gz`` from the directory ``pub/mpi``.  Additional
information about MPICH is avaliable on the World Wide Web at
http://www.mcs.anl.gov/mpi. Note that MPI today is shipped with some
Linux distributions, noteably Debian Linux.

In addition to MPICH, the current installation was compiled successfully
with openmpi and lam.

Installation Instructions
=========================

When installing PGAPack you make two choices: whether to build a sequential
(the default) or parallel version, and whether to build a debug or optimized
(the default) version.  In broad outline, the
installation steps are as follows.

1.  Check out from github
2.  Run ::

      make MPI=$MPIVERSION

    replacing ``$MPIVERSION`` with either ``serial``, ``openmpi``,
    ``mpich``, or ``lam``.  If this doesn't work, you can specify
    ``MPI_LIB`` and/or ``MPI_INCLUDE`` in addition.
    The original targets of the old configure were preserved for
    historical reasons, so you may want to build with::

      make ARCH_TYPE=$ARCHITECTURE

    replacing ``$ARCHITECTURE`` with one of the following:

    ============== ================================================
    Architecture   Description
    ============== ================================================
    sun4           for Sun SparcStations workstations,
    next           for NeXT workstations,
    rs600          for IBM RS6000 workstations,
    irix           for Silicon Graphics workstations,
    hpux           for Hewlett Packard workstations,
    alpha          for DEC Alpha workstations,
    linux          for machines running Linux,
    freebsd        for machines running FreeBSD,
    generic        for generic 32-bit machines,
    powerchallenge for the Silicon Graphics Power Challenge Array,
    challenge      for the Silicon Graphics Challenge,
    t3d            for the Cray T3D,
    sp2            for the IBM SP2,
    paragon        for the Intel Paragon, or
    exemplar       for the Convex  Exemplar.
    ============== ================================================

    The full make options are ``ARCH_TYPE``, ``CC``,
    ``CFLAGS``, ``FC``, ``FFLAGS``, ``DEBUG``, ``MPI_INC``, ``MPI_LIB``

    In addition it is now possible to *add* C-compiler options with
    ADD_CFLAGS and Fortran compiler options with ADD_FFLAGS. The latter
    may be needed with Gnu Fortran compilers prior to major version 10
    because of a `bug in constant declarations`_. Use::

        make MPI=$MPIVERSION ADD_FFLAGS=-fno-range-check

    All parameters are optional and do the following:

    =========== =============================================================
    Parameter   Description
    =========== =============================================================
    CC          The name of the ANSI C compiler, cc by default.
    CPPFLAGS    C Preprocessor flags (later appended to ``CFLAGS``)
    CFLAGS      Options passed to the C compiler including necessary
                options for include file location.
    ADD_CFLAGS  Additional options passed to C compiler.
                This is easier to use than FFLAGS because no knowledge
                of include directives is necessary.
    DEBUG       If specified, enables the debugging features
                and compiles the source code with the ``-g`` flag.
    FC          The name of the Fortran 77 compiler, f77 by default.
                (The Fortran compiler is used only to compile the Fortran
                examples in the ``./examples/`` directory.)
    FFLAGS      Options passed to the Fortran compiler including
                necessary options for include file location.
    ADD_FFLAGS  Additional options passed to the Fortran compiler.
                This is easier to use than FFLAGS because no knowledge
                of include directives is necessary.
    INCLUDES    Include options (usually ``-I directory``) but see the
                ``MPI_INC`` below
    LDFLAGS     Linker options
    ADD_LDFLAGS Additional linker options (in addition to to the
                defaults computed for the current architecture)
    LIBS        Additional libraries, note that you probably have to
                include the math library with ``-lm``
    MPI         Specify one of the known MPI types, one of ``openmpi``,
                ``mpich``, ``lam``, or ``serial``
                (for a non-MPI implementation)
    MPI_INC     The Include-Option where MPI include files are located.
    MPI_LIB     The Linker options for the MPI library, can also be the
                library file to link.
    OPT         The optimization option your compiler understands
    SHAREDLIBS  If set to something different from ``yes`` will not build
                shared libraries
    =========== =============================================================

    If the ``MPI`` or ``MPI_INC``, ``MPI_LIB`` options are specified, a
    parallel version of PGAPack will be built, unless you explicitly
    specify ``MPI=serial``.
    If these flags are not specified, a rudimentary check for a default
    MPI installation is done. If no MPI installation is found, a sequential
    version of PGAPack will be built.

    Note that older versions required to set the ``WL`` (word length)
    preprocessor define. This is no longer required, unless you have a
    very unusual machine where the C-expression::

      sizeof(unsigned long) * 8

    is not the number of bits in an unsigned long (e.g. if you have a
    different size of character).

3.  Add PGAPack's man pages to your man page path::

      setenv MANPATH "$MANPATH"":/home/pgapack/man"

4.  Execute a simple test problem

    Sequential version::

        C:        ``/usr/local/pga/examples/c/maxbit``
        Fortran:  ``/usr/local/pga/examples/fortran/maxbit``

    Parallel version::

        C:        ``mpirun -np 4 /usr/local/pga/examples/c/maxbit``
        Fortran:  ``mpirun -np 4 /usr/local/pga/examples/fortran/maxbit``

    If a parallel version of PGAPack was used, the actual commands to execute
    a parallel program depend on the particular MPI implementation and
    parallel computer.  If the MPICH implementation was used the ``mpirun``
    command can be used to execute a parallel program on most systems.

Compiling without Fortran
-------------------------

Note that Fortran is used only for the Fortran examples in
``examples/fortran`` and ``examples/mgh``. But these are also used in
the tests. If you can live without all test tests passing you can simply
override the ``FC`` (Fortran Compiler) Makefile variable like so::

    make MPI=serial FC=

This will set the Fortran compiler to an empty string and no attempt to
compile fortran code is made. Of course you may chose a different
setting for the MPI variable (e.g. ``MPI=openmpi``).
If you add the ``test`` target::

    make MPI=serial FC= test

Only the tests that do not need a Fortran compiler are run.


Using OpenMPI (Debian, Ubuntu Linux)
====================================

1. Install openmpi::

    sudo apt install libopenmpi-dev

2. Run::

    make MPI=openmpi

3. Execute a simple test problem in examples/c folder:

   - Sequential version::

        ./maxbit

   - Parallel version::

        mpirun -np 4 ./maxbit

   If you want Open MPI to default to the number of hardware threads
   instead of the number of processor cores, use the ``--use-hwthread-cpus``
   option::

        mpirun --use-hwthread-cpus ./maxbit

   Don't be surprised when the parallel version actually runs *slower*
   than the sequential version *on this problem*: The parallel version
   needs additional communication overhead which results in faster
   execution only when the execution time of the evaluation is large
   compared to the communication overhead.

Structure of the Distribution Directory
=======================================

============= ============================================================
File/Dir      Description
============= ============================================================
CHANGES       Changes new to this release of PGAPack.
COPYING       Copyright and disclaimer information.
README.rst    This file.
Makefile      Makefile to build everything
docs          Directory containing documentation. This builds the manual
              from LaTeX sources
examples      A directory containing C and Fortran examples.
include       The PGAPack include directory.
lib           The directory the library will be installed in.
man           The directory containing the PGAPack man pages.
source        The source code for the PGAPack system.
test          A directory containing programs to verify the installation.
              This now runs all the examples including the Fortran
              examples. With no Fortran compiler only the C-Tests are run.
============= ============================================================


Contributions
=============

PGAPack was written to be extensible in two ways: adding new operators that
work with existing data types, and defining new data types.  Enhancements of
either type that you wish to share are welcome for possible inclusion in
future versions of PGAPack.


Acknowledgment
==============

Users of PGAPack are asked to acknowledge its use in any document referencing
work based on the program, such as published research.  Also, please supply
to us a copy of any published research referencing work based on the software.

History
=======

David Levine is the principal author of pgagpack and wrote most of the code
during the mid-1990s. Dirk Eddelbuettel became its Debian maintainer in 2008,
organised a relicensing by Argonne National Laboratories under the MPICH2
license and was the effective upstream maintainer until 2017.

In 2017 maintenance (and some development) was taken over be Ralf
Schlatterbeck, who maintains the github project at
https://github.com/schlatterbeck/pgapack

This repository contains the original 1996, 2008, and 2009 releases as
distributed by Argonne National Laboratories as the first commits. It
then has changes from the google code project (now archived by google at
https://code.google.com/archive/p/pgapack/source) which later became the
git repo of Dirk Eddelbuettel at https://github.com/eddelbuettel/pgapack
Note that the changes by Allan Clark in that repository that introduced
a new automake/autoconf configuration is currently on the autoconf
branch |--| it did not work to build against different variants of MPI
implementations (or against the serial version without MPI). There are
currently no plans to incorporate automake again |--| computer
architectures have become more similar in recent years so that the effort
of maintaining a working automake environment seems not justified.

.. _PGAPy: https://github.com/schlatterbeck/pgapy
.. _`blogpost on epsilon constrained optimization`:
    https://blog.runtux.com/posts/2022/08/29/
.. _`debian-pgapack`: https://github.com/schlatterbeck/debian-pgapack
.. _`examples/c/namefull.c`:
    https://github.com/schlatterbeck/pgapack/blob/master/examples/c/namefull.c
.. _`examples/fortran/namefull.f`:
    https://github.com/schlatterbeck/pgapack/blob/master/examples/fortran/namefull.f
.. _`examples/fortran/constr.f`:
    https://github.com/schlatterbeck/pgapack/blob/master/examples/fortran/constr.f
.. _`examples/c/udtstr.c`:
    https://github.com/schlatterbeck/pgapack/blob/master/examples/c/udtstr.c
.. _`examples/nsgaii/crowdingplot`:
    https://github.com/schlatterbeck/pgapack/blob/master/examples/nsgaii/crowdingplot
.. _`bug in constant declarations`: https://godbolt.org/z/ahMrv4r1E
.. _`Read the Docs`: https://pgapack.readthedocs.io/en/latest/
.. _`Sphinx virtual environment`:
    https://www.sphinx-doc.org/en/master/usage/installation.html#using-virtual-environments
.. _`monkey patching`: https://en.wikipedia.org/wiki/Monkey_patch
