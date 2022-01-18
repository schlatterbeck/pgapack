#include <stdio.h>
#include "pgapack.h"

void p_gcd (int a, int b)
{
    printf ("gcd (%d, %d): %d\n", a, b, LIN_gcd (a, b));
}

void p_binom (int a, int b)
{
    printf ("binom (%d, %d): %d\n", a, b, LIN_binom (a, b));
}

void p_vec (int dim, int n, void *v)
{
    double (*vec)[dim] = v;
    int i, d;
    for (i=0; i<n; i++) {
        for (d=0; d<dim; d++) {
            printf ("%e ", vec [i][d]);
        }
        printf ("\n");
    }
    printf ("\n");
}

void p_dasdennis (int dim, int npart)
{
    void *p = NULL;
    int n = LIN_dasdennis (dim, npart, &p, 0, 1, NULL);
    printf ("dasdennis (%d, %d):\n", dim, npart);
    p_vec (dim, n, p);
}

int main ()
{
    int row, col, i;
    double m1 [][3] = {{0.8, 0.5, 0.5}, {0.1, 0.3, 0.9}, {0.1, 0.3, 0.9}};
    double m2 [][3] = {{1.0, 0.2, 0.0}, {0.4, 0.1, 0.4}, {0.1, 0.0, 1.0}};
    double m3 [3][3];
    double m4 [][3] = {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}};
    double v [3], v1 [3];
    double v2 [] = {1, 2, 3};
    void *vec = NULL;
    int vecl = 0;
    double dir [3] = {1, 1, 1};
    double dir1 [3] = {3, 4, 5};
    double dir2 [3] = {7, 6, 2};
    //LIN_print_matrix (3, m1);
    //LIN_print_matrix (3, m2);
    for (i=0; i<3; i++) {
        for (col=0; col<3; col++) {
            for (row=0; row<2; row++) {
                m3 [col][row] = m2 [row+1][col] - m2 [0][col];
            }
        }
        for (col=0; col<3; col++) {
            m3 [col][2] = col == i ? -1 : 0;
            v [col] = -m2 [0][col];
        }
        //LIN_print_matrix (3, m3);
        //LIN_print_vector (3, v);
        LIN_solve (3, m3, v);
        //LIN_print_matrix (3, m3);
        LIN_print_vector (3, v);
    }
    LIN_solve (3, m4, v2);
    LIN_print_vector (3, v2);

    p_gcd (66, 7);
    p_gcd (77, 7);
    p_gcd (0x1C8CFC00, 2*3*4*5*7);
    p_gcd (35, 42);
    p_binom (3+7-1, 7);
    p_binom (3+6-1, 6);
    p_binom (3+5-1, 5);
    p_binom (3+4-1, 4);
    p_binom (4+4-1, 4);
    p_binom (4+3-1, 3);
    p_dasdennis (2, 1);
    p_dasdennis (2, 2);
    p_dasdennis (3, 1);
    p_dasdennis (3, 2);
    p_dasdennis (3, 4);
    p_dasdennis (3, 5);
    p_dasdennis (4, 1);
    p_dasdennis (4, 2);
    p_dasdennis (4, 3);
    p_dasdennis (4, 4);
    vecl = LIN_dasdennis (3, 2, &vec, vecl, 1, NULL);
    vecl = LIN_dasdennis (3, 1, &vec, vecl, 0.5, dir);
    p_vec (3, vecl, vec);
    free (vec);
    vec  = NULL;
    vecl = 0;
    vecl = LIN_dasdennis (3, 2, &vec, vecl, 1, NULL);
    vecl = LIN_dasdennis (3, 7, &vec, vecl, 0.1, dir1);
    vecl = LIN_dasdennis (3, 5, &vec, vecl, 0.1, dir2);
    p_vec (3, vecl, vec);
    v  [0] = v  [1] = v  [2] = 1;
    v1 [0] = v1 [1] = v1 [2] = 4;
    printf ( "Distance (1,1,1), (4,4,4) = %e\n"
           , LIN_euclidian_distance (3, v, v1)
           );
    v [0] = v [1] = v [2] = 1.0 / sqrt (3);
    printf ("norm2 (1/sqrt(3),1/sqrt(3),1/sqrt(3)) = %e\n", LIN_2norm (3, v));

    return 0;
}
