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

/*!***************************************************************************
* \file
* This file contains linear algebra algorithms.
* \authors Author: Ralf Schlatterbeck
*****************************************************************************/
/*!***************************************************************************
 *  \defgroup linalg Linear Algebra Algorithms
 *****************************************************************************/

#include <stdlib.h>
#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
#define LIN_ERROR_SINGULAR 1
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!***************************************************************
    \brief Solve a linear matrix equation, or system of linear scalar equations.
    \ingroup linalg
    \param n  size of the matrix
    \param a  n * n matrix
    \param b  vector of lenth n
    \return 0 if no error, a positive error-code otherwise, returns the
            solution in b, a and b are modified in-place

    \rst
    Description
    -----------

    Linear matrix equation :math:`ax = b`
    \endrst
******************************************************************/
int LIN_solve (int n, void *a, double *b)
{
    int col, row, row2;
    DECLARE_DYNARRAY (int, rowidx, n);
    DECLARE_DYNPTR (double, m, n) = a;
    DECLARE_DYNARRAY (double, r, n);

    for (row=0; row<n; row++) {
        rowidx [row] = row;
    }

    for (row=0; row<n-1; row++) {
        for (row2=row; row2<n; row2++) {
            if (DEREF2_DYNPTR (m, n, rowidx [row2], rowidx [row]) != 0) {
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
            /* m [rowidx [row2]][row] / m [rowidx [row]][row]; */
            double c = DEREF2_DYNPTR (m, n, rowidx [row2], row)
                     / DEREF2_DYNPTR (m, n, rowidx [row], row);
            for(col=row+1; col<n; col++) {
                /* m [rowidx [row2]][col] -= c * m [rowidx [row]][col]; */
                DEREF2_DYNPTR (m, n, rowidx [row2], col)
                    -= c * DEREF2_DYNPTR (m, n, rowidx [row], col);
            }
            /* m [rowidx [row2]][row] = 0; */
            DEREF2_DYNPTR (m, n, rowidx [row2], row) = 0;
            b [rowidx [row2]] -= c * b [rowidx [row]];
        }
    }
    for (row=n-1; row>=0; row--) {
        double result = 0;
        for (col=row+1; col<n; col++) {
            /* result += m [rowidx [row]][col] * b [rowidx [col]]; */
            result += DEREF2_DYNPTR (m, n, rowidx [row], col)
                    * b [rowidx [col]];
        }
        /* ... / m [rowidx [row]][row] */
        b [rowidx [row]] = (b [rowidx [row]] - result)
                         / DEREF2_DYNPTR (m, n, rowidx [row], row);
        r [row] = b [rowidx [row]];
        if (isnan (r [row])) {
            return LIN_ERROR_SINGULAR;
        }
    }
    memcpy (b, r, sizeof (double) * n);
    return 0;
}

/*!***************************************************************
    \brief Print matrix.
    \ingroup linalg
    \param n  size of the matrix
    \param a  n * n matrix
    \return None

    \rst
    Description
    -----------

    Print matrix :math:`a`. Mainly used for testing.
    \endrst
******************************************************************/
void LIN_print_matrix (int n, void *a)
{
    DECLARE_DYNPTR (double, m, n) = a;
    int row, col;
    for (row=0; row<n; row++) {
        for (col=0; col<n; col++) {
            printf ("%e ", DEREF2_DYNPTR (m, n, row, col));
        }
        printf ("\n");
    }
    printf ("\n");
}

/*!***************************************************************
    \brief Print vector.
    \ingroup linalg
    \param n  size of the vector
    \param v  vector of length n
    \return None

    \rst
    Description
    -----------

    Print vector :math:`v`. Mainly used for testing.
    \endrst
******************************************************************/
void LIN_print_vector (int n, double *v)
{
    int col;
    for (col=0; col<n; col++) {
        printf ("%e ", v [col]);
    }
    printf ("\n\n");
}

/*!***************************************************************
    \brief Dot product of two vectors.
    \ingroup linalg
    \param dim  size of the vectors
    \param v1   first vector of length n
    \param v2   second vector of length n
    \return The dot product of v1 and v2

******************************************************************/
double LIN_dot (int dim, double *v1, double *v2)
{
    int i;
    double ret = 0;
    for (i=0; i<dim; i++) {
        ret += v1 [i] * v2 [i];
    }
    return ret;
}

/*!***************************************************************
    \brief Euclidian distance of two vectors.
    \ingroup linalg
    \param dim  size of the vectors
    \param v1   first vector of length n
    \param v2   second vector of length n
    \return The euclidian distance of v1 and v2

******************************************************************/
double LIN_euclidian_distance (int dim, double *v1, double *v2)
{
    int i;
    double ret = 0;
    for (i=0; i<dim; i++) {
        ret += pow (v1 [i] - v2 [i], 2);
    }
    return sqrt (ret);
}

/*!***************************************************************
    \brief Euclidian norm (2-norm) of a vector.
    \ingroup linalg
    \param dim  size of the vectors
    \param v    vector of length n
    \return The euclidian norm of v

******************************************************************/
double LIN_2norm (int dim, double *v)
{
    return sqrt (LIN_dot (dim, v, v));
}

/*!***************************************************************
    \brief Greatest common divisor of two integers
    \ingroup linalg
    \param a  first integer
    \param b  second integer
    \return The greatest common divisor (gcd) of a and b

******************************************************************/
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
 * Will return 0 on overflow
 */
/*!***************************************************************
    \brief Compute binomial coefficient of two integers
    \ingroup linalg
    \param a  first integer
    \param b  second integer
    \return Binomial coefficient

    \rst
    Description
    -----------

    Computes :math:`\binom{a}{b}`
    \endrst
******************************************************************/
size_t LIN_binom (int a, int b)
{
    int i, j;
    DECLARE_DYNARRAY (size_t, numer, b);
    DECLARE_DYNARRAY (size_t, denom, b);
    int idxn = 0, idxd = 0;
    size_t r;
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
        unsigned int m = r * numer [i];
        if (m < r || m > INT_MAX) {
            return 0;
        }
        r = m;
    }
    return r;
}

/*!***************************************************************
    \brief Normalize a vector with dimension dim to the n-dimensional
           reference hyperplane.
    \ingroup linalg
    \param dim  dimension
    \param v    vector of dimension dim
    \return The vector v is modified in-place

    \rst
    Description
    -----------

    Hyperplane is probably a misnomer, it's a 3-dimensional tetraeder
    for dimension 4. It *is* a plane for the 3-dimensional case going
    through the point 1 on each axis.
    \endrst
******************************************************************/
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

/*!***************************************************************
    \brief Static recursive function for \ref LIN_dasdennis.
    \ingroup internal
    \param dim   dimension
    \param npart Number of partitions
    \param depth Recursion depth
    \param sum   Number of points so far
    \param p     Pointer to created points, memory must be already allocated
    \return The number of points allocated

******************************************************************/
static int dasdennis (int dim, int npart, int depth, int sum, void *p)
{
    DECLARE_DYNPTR(double, vec, dim) = p;
    int n = npart - sum + 1;
    int i, offset = 0;

    if (depth == dim - 1) {
        DEREF2_DYNPTR (vec, dim, 0, depth) = 1.0 - (double)sum / npart;
        return 1;
    }
    for (i=0; i<n; i++) {
        double v = (double)i / npart;
        if (i && depth) {
            memcpy
                ( DEREF1_DYNPTR (vec, dim, offset)
                , DEREF1_DYNPTR (vec, dim, 0)
                , sizeof (double) * dim
                );
        }
        DEREF2_DYNPTR (vec, dim, offset, depth) = v;
        offset += dasdennis
            (dim, npart, depth + 1, sum + i, DEREF1_DYNPTR (vec, dim, offset));
    }
    return offset;
}

/*!***************************************************************
    \brief Static function for \ref LIN_dasdennis scaling of points
    \ingroup internal
    \param dim     dimension
    \param npoints Number of points
    \param scale   Scaling factor
    \param dir     Direction vector
    \param v       Pointer to points, points must have been created
    \return None

    \rst
    Description
    -----------

    When scaling we translate the scaled-down points from the
    centroid of the shifted points to the reference direction on the
    reference plane.
    \endrst
******************************************************************/
STATIC
void dasdennisscale (int dim, int npoints, double scale, double *dir, void *v)
{
    int i, j;
    DECLARE_DYNARRAY (double, dir_normed, dim);
    DECLARE_DYNARRAY (double, centroid, dim);
    DECLARE_DYNPTR (double, vec, dim) = v;
    assert (scale > 0);
    assert (scale < 1);
    for (i=0; i<npoints; i++) {
        for (j=0; j<dim; j++) {
            DEREF2_DYNPTR (vec, dim, i, j) *= scale;
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
            DEREF2_DYNPTR (vec, dim, i, j) += dir_normed [j] - centroid [j];
        }
    }
}

/*!***************************************************************
    \brief Compute Das & Dennis points to allocated memory.
    \ingroup linalg
    \param dim     dimension
    \param npart   Number of partitions
    \param scale   Scaling factor
    \param dir     Direction vector
    \param npoints Number of points
    \param mem     Memory for points, memory must have been allocated
    \return None

    \rst

    Description
    -----------

    This is the case where the memory is already allocated.
    For more details see :c:func:`LIN_dasdennis`.

    \endrst
******************************************************************/
void LIN_dasdennis_allocated
    (int dim, int npart, double scale, double *dir, int npoints, void *mem)
{
    dasdennis (dim, npart, 0, 0, mem);
    if (scale != 1 && dir != NULL) {
        dasdennisscale (dim, npoints, scale, dir, mem);
    }
}

/*!***************************************************************
    \brief Allocate memory and compute Das & Dennis points.
    \ingroup standard-api
    \param dim     dimension
    \param npart   Number of partitions
    \param result  List of existing points to extend
    \param nexist  Number of existing points
    \param scale   Scaling factor
    \param dir     Direction vector
    \return Number of points allocated, -1 on error

    \rst
    Description
    -----------

    It will re-alloc the exiting array pointer pointed to by
    result (this must be a ``NULL`` pointer if no pre-existing points are
    given) and return the new number of points. Note that if there are no
    pre-existing points, the pointer pointed to by ``result`` must be
    ``NULL`` and ``nexist`` must be 0.
    Optionally the points can be scaled (with 0 < scale <= 1) and shifted
    in the direction of a given point back onto the reference hyperplane.
    This is not done if ``dir == NULL`` or ``scale == 1``.
    A previously allocated result will be de-allocated in case of error.
    \endrst
******************************************************************/
int LIN_dasdennis
    (int dim, int npart, void *result, int nexist, double scale, double *dir)
{
    void **r = result;
    DECLARE_DYNPTR(double, vec, dim);
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
    LIN_dasdennis_allocated
        (dim, npart, scale, dir, npoints, DEREF1_DYNPTR (vec, dim, nexist));
    return nexist + npoints;
err:
    if (nexist) {
        free (*r);
        *r = NULL;
    }
    return -1;
}
