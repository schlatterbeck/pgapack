.. |--| unicode:: U+2013   .. en dash

PGAPack
+++++++

PGAPack is a general-purpose, data-structure-neutral, parallel genetic
algorithm library originally developed by David Levine at Argonne
National Laboratory. It has libraries for C and Fortran. There are
companion projects:

- PGAPy, a Python wrapper for PGApack https://github.com/schlatterbeck/pgapy
- debian-pgapack, a project for building debian packages (for all MPI
  libraries supported by Debian):
  https://github.com/schlatterbeck/debian-pgapack

Updates 
=======

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
- Make fortran compile again

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
- Easy integration of hill-climbing heuristics.
- Easy-to-use interface for novice and application users.
- Fully extensible to support custom operators and new data types.
- Extensive debugging facilities.
- A large set of example problems.
- It is released under the MPICH2 license (also used by the MPICH2 MPI
  implementation from Argonne National Laboratory).


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

* The PGAPack users guide is in ``./docs/user_guide.pdf`` after building the
  guide from sources (see ``Makefile``). The old original version was
  preserved as ``docs/user_guide-orig.ps`` |--| it is recommended to use the
  latest version that had some fixes and documentation updates for newer
  features.
* Man pages for PGAPack functions are in the ``./man`` directory.
* Installation instructions are in this ``README.rst`` file and the
  users guide.
* Example problems are in the ``./examples`` directory.


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

    where all parameters are optional and do the following:

    =========== =============================================================
    Parameter   Description
    =========== =============================================================
    CC          The name of the ANSI C compiler, cc by default.
    CPPFLAGS    C Preprocessor flags (later appended to ``CFLAGS``)
    CFLAGS      Options passed to the C compiler.
    DEBUG       If specified, enables the debugging features
                and compiles the source code with the ``-g`` flag.
    FC          The name of the Fortran 77 compiler, f77 by default.
                (The Fortran compiler is used only to compile the Fortran
                examples in the ``./examples/`` directory.)
    FFLAGS      Options passed to the Fortran compiler.
    INCLUDES    Include options (usually ``-I directory``) but see the
                ``MPI_INC`` below
    LDFLAGS     Linker options
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



PGAPack on PCs 
===================

PGAPack has not been ported to MS-DOS, Windows 3.1, Windows 95, or Apple OS.
As mentioned earlier, however, PGAPack is written in ANSI standard C and
should compile in these environments.  Be aware, however, that PGAPack's
random number generator, PGARandom01, assumes certain machine characteristics
for ints and floats that may not correspond with what your PC and/or compiler
support, resulting in erroneous values.


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
              Note that the verification is known to fail.
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
