Sequencing Examples
===================

This directory contains examples of using sequence-preserving crossover
and mutation operators. It is a simple traveling salesperson problem
with 5 cities. We test the Edge Recombination crossover, first
introduced in 1989 [1]_ and later improved [3]_ as the crossover
operator.  In addition a permutation mutation operator is used.  These
operators both preserve the semantics that the used integer genes form a
permutation. The Edge Recombination crossover allows to fix certain
edges (which may not be modified by the crossover operation). This
feature is also implemented and available as a command-line switch.

Note that there are *a lot* of sequence-preserving crossover variants in
the literature. This particular operand works well for traveling
salesperson but not so well for other problems. Different
sequence-preserving crossover operators perform differently on different
problems. Also note that a pure genetic algorithm (GA) implementation
(without coupling the GA with a good heuristic for the problem) will
probably perform less satisfactory than a good heuristic. One of the
earliest heuristics for finding good (but typically not optimal)
solutions to traveling salesperson problems is the 1973 paper by Lin and
Kernighan [5]_ which performs better than the GA given in this example
(for larger problems). One of the reasons for this is that the GA sees
only the result from the whole tour in the evaluation result. Goldberg
has called this the "blind traveling salesman" [6]_, p. 170.

.. [1] Darrell Whitley, Tim Starkweather, and D’Ann Fuquay. Scheduling
       problems and traveling salesman: The genetic edge recombination
       operator. In Schaffer [2]_, pages 133–140.
.. [2] J. David Schaffer, editor. *Proceedings of the Third International
       Conference on Genetic Algorithms (ICGA).* Morgan Kaufmann, June 1989.
.. [3] Darrel Whitley, Timothy Starkweather, and Daniel Shaner.
       The traveling salesman and sequence scheduling: Quality solutions
       using genetic edge recombination. In Davis [4]_, chapter 22,
       pages 350–372.
.. [4] Lawrence Davis, editor. *Handbook of Genetic Algorithms.* Van
       Nostrand Reinhold, 1991.
.. [5] S. Lin and B. W. Kernighan. An effective heuristic algorithm for
       the traveling-salesman problem. Operations Research,
       21(2):498–516, March 1973.
.. [6] David E. Goldberg. Genetic Algorithms in Search, Optimization &
       Machine Learning. Addison Wesley, October 1989.
