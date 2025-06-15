Multiobjective Example
======================

This directory contains C source code of multiobjective problems
from the NSGA-II paper by Kalyanmoy Deb et. al. [1]_.

Note that it was tried to solve the given problems using Differential
Evolution [2]_, [3]_, [4]_ -- except for the population replacement
scheme which uses NSGA-II.

Some of the problems were given more generations than in the paper. Note
that the ZDT4 example has a typo in the NSGA-II paper [1]_ and should
remove the factor g(x) from f2(x), see the original paper by Zitzler,
Deb, and Thiele [5]_ p.178 formula (10). Even with more generations, the
current implementation fails to find a solution similar to the NSGA-II
paper [1]_: It is better in the lower values of f1 but worse in the
higher values. It's currently unsure if the test problem is correctly
implemented or if there are other errors.

There is a small section where a rotated problem is tried. This is
implemented in ``rotated.c`` and can be called with index 13 of the
optimize program. Since the paper does not mention the rotation, we
rotate around the plane produced by the angle bisector of ``(x [0], x [1])``
and ``(x [2], x [3])``. The rotation angle is 45°.

We give it a *lot* more generations than in the paper
(2500 instead of 500) and achieve quite good results. With 5000
generations the solution is almost as good as with 500 generations of
the unrotated problem. Note that we force ``yi`` to be between -0.3 and
0.3 as in the paper. This avoids the multiple pareto fronts that would
occur with a larger range as in the ZDT4 example.

The multi-objective problems *with* constraints are also implemented.
Note that in the "water" problem there is also a small typo on the
formula for g6 which should have a ``/`` (division instead of
multiplication like for all the other constraint functions) after the
first constant [6]_. This typo has not much effect on the results.

There is a small python script ``crowdingplot`` which can be used to
do a scatter-plot of one objective against another (in objective space).
By default the first variable is plotted against the second but the
variables to plot can be given with the -x and -y options. This
reproduces some of the results in the paper [1]_. Just pipe the output
of the optimization run through the script or use the ``nsgaii*.data``
files in the test directory as input.

The water problem (f-index 12 when calling optimize) has strong correlations
between objectives 0 and 2, we can see this when we plot the two
objectives using::

    crowdingplot -x=0 -y=2 ../../test/nsgaii_optimize_12.data

This means that getting a better evaluation for objective 0 also gets a
better evaluation of objective 2 and vice-versa, i.e. the objectives are
not contradictory. Furthermore objective 3 is redundant [7]_ with the
other objectives and can be left out when plotting objectives. So the
5-dimensional (the problem has 5 objectives and 7 constraints) pareto
front can be reduced to a 3-dimensional front and can be plotted with
crowdingplot and yield quite similar results (when scaling to the same
visible range)::

 crowdingplot -3 -x=0 -y=1 -z=4 ../../test/nsgaii_optimize_12.data
 crowdingplot -3 -x=2 -y=1 -z=4 ../../test/nsgaii_optimize_12.data

This means we can reduce the number of objectives not just for plotting,
but during optimization, too [7]_.

Details on the "water" benchmark including graphics can be found on my
blog [8]_.

.. [1] Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, and T. Meyarivan.
       A fast and elitist multiobjective genetic algorithm: NSGA-II.
       IEEE Transactions on Evolutionary Computation, 6(2):182–197,
       April 2002.
.. [2] Rainer Storn and Kenneth Price. Differential evolution – a simple
       and efficient adaptive scheme for global optimization over
       continuous spaces. Technical Report TR-95-012, International
       Computer Science Institute (ICSI), March 1995.
.. [3] Rainer Storn and Kenneth Price. Differential evolution – a simple
       and efficient heuristic for global optimization over continuous spaces.
       Journal of Global Optimization, 11(4):341–359, December 1997.
.. [4] Kenneth V. Price, Rainer M. Storn, and Jouni A. Lampinen.
       Differential Evolution: A Practical Approach to Global
       Optimization.  Springer, Berlin, Heidelberg, 2005.
.. [5] Eckart Zitzler, Kalyanmoy Deb, and Lothar Thiele. Comparison of
       multiobjective evolutionary algorithms: Empirical results.
       Evolutionary Computation, 8(2):173–195, 2000.
.. [6] Tapabrata Ray, Kang Tai, and Kin Chye Seow. Multiobjective design
       optimization by an evolutionary algorithm. Engineering Optimization,
       33(4):399–424, 2001.
.. [7] Hemant Kumar Singh, Amitay Isaacs, and Tapabrata Ray.  A pareto
       corner search evolutionary algorithm and dimensionality reduction
       in many-objective optimization problems. IEEE Transactions on
       Evolutionary Computation, 15(4):539–556, August 2011.
.. [8] Ralf Schlatterbeck. `Water: A multi-objective benchmark problem`_.
       Blog post, Open Source Consulting, June 2025.

.. _`Water: A multi-objective benchmark problem`:
    https://blog.runtux.com/posts/2025/06/15/
