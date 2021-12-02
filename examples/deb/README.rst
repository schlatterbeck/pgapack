Constraint Handling Example
===========================

This directory contains C source code of constrained optimization
problems from Kalyanmoy Deb's constraint handling method [1]_.

The second problem with 38 constraints was first published by Hock and
Schittkowski [2]_ and is number 85 in their publication.
The best value documented by Schittkowski is -1.9051338 for which the
input values are given (see comments at the end of deb2.c). But he claims
that his optimizer found a solution with value -1.9051553 (which is
better and also found by the differential evolution implemenation in this
directory). Deb gives the evaluation value -1.914595 which is even lower,
but it produces constraint violations, probably because of a typo in one
of the tables, see also comments in deb2.c.

_[1]: Kalyanmoy Deb. An efficient constraint handling method
      for genetic algorithms. Computer Methods in Applied Mechanics and
      Engineering, 186(2–4):311–338, June 2000.
_[2]: W. Hock and K. Schittkowski. Test examples for nonlinear
      programming codes. Lecture Notes in Economics and Mathematical
      Systems, 187, 1981.
