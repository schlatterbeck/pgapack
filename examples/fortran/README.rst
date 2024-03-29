C Examples
==========

Seven example FORTRAN programs come with PGAPack.  What follows is a listing 
of the names of the programs and a brief description of their function.


+-----------+------------------------------------------------------------------+
| classic.f | Optimizes Griewank's, Rastrigin's or Schwefel's test function.   |
|           | All use a real valued datatype.  The problems are:               |
|           |                                                                  |
|           |   1.  Griewank                                                   |
|           |   2.  Rastrigin                                                  |
|           |   3.  Schwefel                                                   |
+-----------+------------------------------------------------------------------+
| dejong.f  | This program optimizes one of the five functions from            |
|           | the DeJong test suite with a binary string as the chromosome.    |
|           | Graycoding is an option.                                         |
+-----------+------------------------------------------------------------------+
| example.f | An example from the User's Guide.                                |
+-----------+------------------------------------------------------------------+
| maxbit.f  | A very simple example that maximizes the number of 1's in a      |
|           | binary string.  Does no I/O.                                     |
+-----------+------------------------------------------------------------------+
| maxchar.f | Maximizes the number of z's in a character string.  Does I/O.    |
+-----------+------------------------------------------------------------------+
| maxint.f  | Maximizes the sum of all alleles in the gene.                    |
+-----------+------------------------------------------------------------------+
| name.f    | Evolves a character string to match a hard-coded string.  Does   |
|           | no I/O, but shows how to use userfunctions InitString, Mutation, |
|           | and StopCond.  Also uses an elitist model of evolution with no   |
|           | crossover.                                                       |
+-----------+------------------------------------------------------------------+
| namefull.c| Same idea as name.c, but uses all native user defined functions. |
+-----------+------------------------------------------------------------------+
