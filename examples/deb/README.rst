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
is solved with one order of magnitude less precision (1e-2 for epsilon
instead of 1e-3 in [1]_) with only a population size of 20 (!) and 2000
generations. In the paper [1]_ Deb uses a population size of 50 and 7000
generations. With a smaller epsilon it converges to a local minimum that
is one order of magnitude worse than the best solution. I've not tried if
with different random number initialization the problem could be solved
with a smaller epsilon.

Here the Epsilon Constraint method [6]_ comes into play: This allows
violation of constraints up to an epsilon bound. The epsilon bound is
decreased during the optimization run and at some generation G the
bound is set to 0, in our case G=1500. This has the effect that while
searching for solutions without constraint violations, the objective
function is optimized, too. With this method, again with a population
size of 20 and 2000 generations, the search finds a solution with a
bound of 1e-6 (instead of 1e-3 in the original paper [1]) and a solution
that is only a little better than the optimal solution cited (the
optimal solution doesn't allow violations of bounds, we violate the
bounds by 1e-6, so we find a slightly better solution).

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
.. [6] Tetsuyuki Takahama and Setsuko Sakai. Constrained optimization by
       the ε constrained differential evolution with an archive and
       gradient-based mutation. In IEEE Congress on Evolutionary
       Computation (CEC), Barcelona, Spain, July 2010.
