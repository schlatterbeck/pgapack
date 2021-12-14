Constraint Handling Example
===========================

This directory contains C source code of constrained optimization
problems from Kalyanmoy Deb's constraint handling method [1]_.

The second problem with 38 constraints was first published by Hock and
Schittkowski [2]_ and is number 85 in their publication.
The best value documented by Schittkowski is -1.9051338 for which the
input values are given (see comments at the end of deb2.c). He claims
that his optimizer found a solution with value -1.9051553 (which is
better and also found by the differential evolution implemenation in this
directory). Deb gives the evaluation value -1.914595 which is even lower,
but it produces constraint violations, probably because of a typo in one
of the tables, see also comments in deb2.c.

The last problem (number 9) is the welded beam design from Deb [1]_ which
is introduced in the beginning of the paper and revisited at the end.

Note that it was tried to solve the given problems using Differential
Evolution [3]_, [4]_, [5]_ with a lower population size as in Deb [1]_
and often less generations -- although some tough problems require a
larger number of generations, but due to the smaller population size with
still equal or smaller number of function evaluations. Consider problem 1
as an example: It has a population size of only 4 (compared to 20 in [1]_)
but 150 instead of 50 generations. The number of function evaluations is
still only a little more than half compared to the paper [1]_.  The
solutions found (with only a hard-coded random number initialisation of
1, i.e.  only a single experiment) are generally better than the ones
achieved in the original paper [1]_.

An exception is problem 7 which has 3 equality constraints. The problem
is solved with one order of magnitude less precision (10e-2 for epsilon
instead of 10e-3 in [1]_) with only a population size of 4 (!) and 2000
generations. With a smaller epsilon it converges to a local minimum that
is one order of magnitude worse than the best solution. I've not tried if
with different random number initialization the problem could be solved
with a smaller epsilon.

.. [1] Kalyanmoy Deb. An efficient constraint handling method
       for genetic algorithms. Computer Methods in Applied Mechanics and
       Engineering, 186(2–4):311–338, June 2000.
.. [2] W. Hock and K. Schittkowski. Test examples for nonlinear
       programming codes. Lecture Notes in Economics and Mathematical
       Systems, 187, 1981.
.. [3] Rainer Storn and Kenneth Price. Differential evolution – a simple
       and efficient adaptive scheme for global optimization over
       continuous spaces. Technical Report TR-95-012, International
       Computer Science Institute (ICSI), March 1995.
.. [4] Rainer Storn and Kenneth Price. Differential evolution – a simple
       and efficient heuristic for global optimization over continuous spaces.
       Journal of Global Optimization, 11(4):341–359, December 1997.
.. [5] Kenneth V. Price, Rainer M. Storn, and Jouni A. Lampinen.
       Differential Evolution: A Practical Approach to Global
       Optimization.  Springer, Berlin, Heidelberg, 2005.
