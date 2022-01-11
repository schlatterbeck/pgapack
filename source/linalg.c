/*
COPYRIGHT

The following is a notice of limited availability of the code, and disclaimer
which must be included in the prologue of the code and in all source listings
of the code.

(C) COPYRIGHT 2022 Dr. Ralf Schlatterbeck Open Source Consulting

Permission is hereby granted to use, reproduce, prepare derivative works, and
to redistribute to others. This software was authored by:

Ralf Schlatterbeck
Open Source Consulting
*/

/*****************************************************************************
*     File: linalg.c: Linear Algebra Algorithms
*
*     Authors: Ralf Schlatterbeck
*****************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include "pgapack.h"
#define LIN_ERROR_SINGULAR 1

/*
 * LIN_solve:
 * Solve a linear matrix equation, or system of linear scalar equations.
 * Linear matrix equation ax = b
 * Where a is an n * n matrix and b is a vector of length n.
 * Returns the solution in b, a and b are modified in-place.
 * Returns 0 if no error, a positive error-code otherwise.
 */
int LIN_solve (int n, void *a, double *b)
{
    int col, row, row2;
    int rowidx [n];
    double (*m) [n] = a;
    double r [n];

    for (row=0; row<n; row++) {
        rowidx [row] = row;
    }

    for (row=0; row<n-1; row++) {
        for (row2=row; row2<n; row2++) {
            if (m [rowidx [row2]][rowidx [row]] != 0) {
                int tmp = rowidx [rowidx [row]];
                rowidx [rowidx [row]]  = rowidx [rowidx [row2]];
                rowidx [rowidx [row2]] = tmp;
                break;
            }
            if (row2 == n) {
                return LIN_ERROR_SINGULAR;
            }
        }
        for (row2=row+1; row2<n; row2++) {
            double c = m [rowidx [row2]][row] / m [rowidx [row]][row];
            for(col=row+1; col<n; col++) {
                m [rowidx [row2]][col] -= c * m [rowidx [row]][col];
            }
            m [rowidx [row2]][row] = 0;
            b [rowidx [row2]] -= c * b [rowidx [row]];
        }
    }
    for (row=n-1; row>=0; row--) {
        double result = 0;
        for (col=row+1; col<n; col++) {
            result += m [rowidx [row]][col] * b [rowidx [col]];
        }
        b [rowidx [row]] = (b [rowidx [row]] - result) / m [rowidx [row]][row];
        r [row] = b [rowidx [row]];
        if (isnan (r [row])) {
            return LIN_ERROR_SINGULAR;
        }
    }
    memcpy (b, r, sizeof (double) * n);
    return 0;
}

void LIN_print_matrix (int n, void *a)
{
    double (*m) [n] = a;
    int row, col;
    for (row=0; row<n; row++) {
        for (col=0; col<n; col++) {
            printf ("%e ", m [row][col]);
        }
        printf ("\n");
    }
    printf ("\n");
}

void LIN_print_vector (int n, double *v)
{
    int col;
    for (col=0; col<n; col++) {
        printf ("%e ", v [col]);
    }
    printf ("\n\n");
}

double LIN_dot (int dim, double *v1, double *v2)
{
    int i;
    double ret = 0;
    for (i=0; i<dim; i++) {
        ret += v1 [i] * v2 [i];
    }
    return ret;
}

double LIN_euclidian_distance (int dim, double *v1, double *v2)
{
    int i;
    double ret = 0;
    for (i=0; i<dim; i++) {
        ret += pow (v1 [i] - v2 [i], 2);
    }
    return sqrt (ret);
}

double LIN_2norm (int dim, double *v)
{
    return sqrt (LIN_dot (dim, v, v));
}

/*
 * LIN_gcd: Greates common divisor of two integers
 */
int LIN_gcd (int a, int b)
{
    int m;
    if (b > a) {
        return LIN_gcd (b, a);
    }
    m = a % b;
    if (!m) {
        return b;
    }
    return LIN_gcd (b, m);
}

/*
 *                                           ( a )
 * LIN_binom: Compute binom of two integers (     )
 *                                           ( b )
 * Will return -1 on overflow
 */
int LIN_binom (int a, int b)
{
    int i, j;
    int numer [b];
    int denom [b];
    int idxn = 0, idxd = 0;
    unsigned int r;
    assert (a > b);
    assert (b >= 1);
    if (b > a / 2) {
        b = a - b;
    }
    for (i=1; i<=b; i++) {
        int n = a - b + i;
        int d = i;
        for (j=0; j<idxn; j++) {
            int g = LIN_gcd (numer [j], d);
            if (g > 1) {
                numer [j] /= g;
                d /= g;
            }
            if (d == 1) {
                break;
            }
        }
        if (d > 1) {
            denom [idxd++] = d;
        }
        for (j=0; j<idxd; j++) {
            int g = LIN_gcd (n, denom [j]);
            if (g > 1) {
                denom [j] /= g;
                n /= g;
            }
            if (n == 1) {
                break;
            }
        }
        if (n > 1) {
            numer [idxn++] = n;
        }
    }
    for (i=0; i<idxd; i++) {
        assert (denom [i] == 1);
    }
    r = 1;
    for (i=0; i<idxn; i++) {
        int m = r * numer [i];
        if (m < r || m > INT_MAX) {
            return -1;
        }
        r = m;
    }
    return r;
}

/*
 * LIN_normalize_to_refplane
 * Normalize a vector with dimension dim to the n-dimensional reference
 * hyperplane. And, yes, hyperplane is probably a misnomer, it's a
 * 3-dimensional tetraeder for dimension 4.
 */
void LIN_normalize_to_refplane (int dim, double *v)
{
    int j;
    double sq = sqrt (dim);
    double norm = 0;

    /* Scalar product of v with the middle of the reference plane */
    for (j=0; j<dim; j++) {
        norm += v [j] * (1.0 / sq);
    }
    for (j=0; j<dim; j++) {
        v [j] = v [j] / norm / sq;
    }
}

/* Static recursive function for LIN_dasdennis, see below */
static int dasdennis (int dim, int npart, int depth, int sum, void *p)
{
    double (*vec)[dim] = p;
    int n = npart - sum + 1;
    int i, offset = 0;

    if (depth == dim - 1) {
        vec [0][depth] = 1.0 - (double)sum / npart;
        return 1;
    }
    for (i=0; i<n; i++) {
        double v = (double)i / npart;
        if (i && depth) {
            memcpy (vec [offset], vec [0], sizeof (double) * dim);
        }
        vec [offset][depth] = v;
        offset += dasdennis (dim, npart, depth + 1, sum + i, vec [offset]);
    }
    return offset;
}

/* Static function for LIN_dasdennis scaling of points
 * Note that when scaling we translate the scaled-down points from the
 * centroid of the shifted points to the reference direction on the
 * reference plane.
 */
void dasdennisscale (int dim, int npoints, double scale, double *dir, void *v)
{
    int i, j;
    double dir_normed [dim];
    double centroid [dim];
    double (*vec)[dim] = v;
    assert (scale > 0);
    assert (scale < 1);
    for (i=0; i<npoints; i++) {
        for (j=0; j<dim; j++) {
            vec [i][j] *= scale;
        }
    }
    for (j=0; j<dim; j++) {
        dir_normed [j] = dir [j];
    }
    LIN_normalize_to_refplane (dim, dir_normed);
    for (j=0; j<dim; j++) {
        centroid [j] = dir_normed [j] * scale;
    }

    /* Now shift points back to reference plane */
    for (i=0; i<npoints; i++) {
        for (j=0; j<dim; j++) {
            vec [i][j] += dir_normed [j] - centroid [j];
        }
    }
}

/*
 * Compute Das & Dennis points
 * This is the case where the memory is already allocated
 * For details see below.
 */
void LIN_dasdennis_allocated
    (int dim, int npart, double scale, double *dir, int npoints, void *mem)
{
    dasdennis (dim, npart, 0, 0, mem);
    if (scale != 1 && dir != NULL) {
        dasdennisscale (dim, npoints, scale, dir, mem);
    }
}

/*
 * Compute Das & Dennis points
 * This gets the dimension of the problem, the number of partitions,
 * a list of existing points to extend, and the number of existing
 * points. It will re-alloc the exiting array pointer pointed to by
 * result (this must be a NULL pointer if no pre-existing points are
 * given) and return the new number of points. Note that if there are no
 * pre-existing points, the pointer pointed to by result must be NULL
 * and exiting must be 0.
 * Optionally the points can be scaled (with 0 < scale <= 1) and shifted
 * in the direction of a given point back onto the reference hyperplane.
 * This is not done if dir == NULL or scale == 1.
 * Will return -1 on error. A previously allocated result will be
 * de-allocated in case of error.
 */
int LIN_dasdennis
    (int dim, int npart, void *result, int nexist, double scale, double *dir)
{
    void **r = result;
    double (*vec)[dim];
    int npoints = LIN_binom (dim + npart - 1, npart);
    void *new = NULL;
    assert (  (*r == NULL && nexist == 0)
           || (*r != NULL && nexist != 0)
           );
    if (npoints < 0) {
        goto err;
    }
    if (nexist) {
        new = realloc (*r, (nexist + npoints) * sizeof (double) * dim);
    } else {
        new = malloc (sizeof (double) * dim * npoints);
    }
    if (new == NULL) {
        goto err;
    }
    vec = *r = new;
    LIN_dasdennis_allocated (dim, npart, scale, dir, npoints, vec + nexist);
    return nexist + npoints;
err:
    if (nexist) {
        free (*r);
        *r = NULL;
    }
    return -1;
}
