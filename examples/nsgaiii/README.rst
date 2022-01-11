Multiobjective Example
======================

This directory contains C source code of multiobjective problems
from the NSGA-III papers by Kalyanmoy Deb and Himanshu Jain [1]_ [2]_.
The problems are described in another paper [4]_.

I've also used implementation notes from [3]_ and implemented the idea
of reference directions [5]_. Note that we're using only 3 objectives to
be able to plot the results (and because we're using it for regression
testing which shouldn't use too much time).

The multi-objective problems *with* constraints [2]_ are also implemented.

There is a small python script ``crowdingplot3`` which can be used to
do a 3D-scatter-plot of three objective (in objective space).

.. [1] Kalyanmoy Deb and Himanshu Jain. An evolutionary many-objective
       optimization algorithm using reference-point-based nondominated
       sorting approach, part I: Solving problems with box constraints.
       IEEE Transactions on Evolutionary Computation, 18(4):577–601,
       August 2014.
.. [2] Himanshu Jain and Kalyanmoy Deb. An evolutionary many-objective
       optimization algorithm using reference-point-based nondominated
       sorting approach, part II: Handling constraints and extending to
       an adaptive approach. IEEE Transactions on Evolutionary Computation,
       18(4):602–622, August 2014.
.. [3] Julian Blank, Kalyanmoy Deb, and Proteek Chandan Roy. Investigating
       the normalization procedure of NSGA-III. In Kalyanmoy Deb, Erik
       Goodman, Carlos A. Coello Coello, Kathrin Klamroth, Kaisa Miettinen,
       Sanaz Mostaghim, and Patrick Reed, editors, Evolutionary
       Multi-Criterion Optimization, 10th International Conference (EMO),
       volume 11411 of Lecture Notes in Computer Science, pages 229–240.
       Springer, East Lansing, MI, USA, March 2019.
.. [4] Kalyanmoy Deb, Lothar Thiele, Marco Laumanns, and Eckart Zitzler.
       Scalable test problems for evolutionary multiobjective optimization.
       In Ajith Abraham, Lakhmi Jain, and Robert Goldberg, editors,
       Evolutionary Multiobjective Optimization, Theoretical Advances
       and Applications, Advanced Information and Knowledge Processing,
       pages 105–145. Springer, 2005.
.. [5] Yash Vesikar, Kalyanmoy Deb, and Julian Blank. Reference point
       based NSGA-III for preferred solutions. In Suresh Sundaram, editor,
       IEEE Symposium Series on Computational Intelligence (SSCI),
       pages 1587–1594. Bengaluru, India, November 2018.
