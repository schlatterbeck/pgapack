## About

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
- It is released under the MPICH2 license (also used by the MPICH2 MPI implementation from Argonne National Laboratory).

## History

David Levine is the principal author of pgagpack and wrote most of the code
during the mid-1990s. Dirk Eddelbuettel became its Debian maintainer in 2008,
organised a relicensing by Argonne National Laboratories under the MPICH2
license and is currently also the effective upstream maintainer.

This repository contains the original 1996, 2008, and 2009 releases as
distributed by Argonne National Laboratories as the first commits. It
then has changes from the google code project (now archived by google at
https://code.google.com/archive/p/pgapack/source) which later became the
git repo of Dirk Eddelbuettel at https://github.com/eddelbuettel/pgapack
Note that the changes by Allan Clark in that repository that introduced
a new automake/autoconf configuration is currently on the autoconf
branch -- it did not work to build against different variants of MPI
implementations (or against the serial version without MPI). It is
intended to eventually merge this and implement configuration options
that can easily configure for different MPI implementations.
