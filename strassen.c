/* Implementation of fast matrix multiplication algorithm using Strassen method
 * 2016 by Anton Leshchenko
 * */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "strassen.h"

int strassen_threshold = 64;

static void def()
{
    fprintf(stderr, "Using default threshold: %d\n", strassen_threshold);
}

int main(int argc, char *argv[])
{
    int ret = EXIT_SUCCESS;
    if (argc < 2) {
        fprintf(stderr, "Usage: %s N <filename>, where N is size of input matrices\n", argv[0]);
        return EXIT_FAILURE;
    }
    size_t n = atoi(argv[1]);
    if (argc > 2) {
        int thr = atoi(argv[2]);
        if (thr == next_pow2(thr))
            strassen_threshold = thr;
        else
            def();

    }
    else
        def();
    DYN_ARRAY(a, n);
    DYN_ARRAY(b, n);
    DYN_ARRAY(mul, n);
    DYN_ARRAY(res, n);
    if (input_read(n, a, b)) {
        fprintf(stderr, "Error reading input data!\n");
        ret = EXIT_FAILURE;
        goto end;
    }

    strassen_mul(n, a, b, mul);
    sum(n, a, mul, res);
    print_matrix(n, res);
end:
    free(a);
    free(b);
    free(mul);
    free(res);
    return ret;
}

static void print_matrix(size_t n, data *matrix)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j)
                printf("\t%lf", matrix[i*n+j]);
            else
                printf("%lf", matrix[i*n+j]);
        }
        printf("\n");
    }
}


/* Matrix operations: multiply, sum, subtract;
 * builtins are used to tell compiler that
 * input data is aligned, so that SIMD operations
 * can be safely applied
 * Restrict assumes memory areas do not overlap
 *	A*B */
static void
mul_ikj(const size_t n, const data *restrict a,
        const data *restrict b, data *restrict res)
{
    int i, j, k;

    data *x = __builtin_assume_aligned(a, 16);
    data *y = __builtin_assume_aligned(b, 16);
    data *z = __builtin_assume_aligned(res, 16);

    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) { /* Optimize memory access for better cache
				     hits */
            const data aink = x[i*n + k];
            for (j = 0; j < n; j++)
                z[i*n + j] += aink * y[k*n + j];
        }
    }
}
/*	A + B */
static void
sum(const size_t n, const data *restrict a,
    const data *restrict b, data *restrict res)
{
    int i;

    data *x = __builtin_assume_aligned(a, 16);
    data *y = __builtin_assume_aligned(b, 16);
    data *z = __builtin_assume_aligned(res ,16);

    for (i = 0; i < n*n; i++) {
        z[i] = x[i] + y[i];
    }
}
/*	A - B */
static void
subtr(const size_t n, const data *restrict a,
      const data *restrict b, data *restrict res)
{
    int i;

    data *x = __builtin_assume_aligned(a, 16);
    data *y = __builtin_assume_aligned(b, 16);
    data *z = __builtin_assume_aligned(res ,16);

    for (i = 0; i < n*n; i++) {
        z[i] = x[i] - y[i];
    }
}
/* Actual strassen multiply function */
static void
strassen_r(const size_t n, const data *restrict a, const data *restrict b,
           data *restrict res)
{
    if (n <= strassen_threshold) {
        mul_ikj(n, a, b, res);
        return;
    }
    size_t half = n / 2;
    DYN_ARRAY(a_res, half);
    DYN_ARRAY(b_res, half);

    DYN_ARRAY(a11, half);
    DYN_ARRAY(a12, half);
    DYN_ARRAY(a21, half);
    DYN_ARRAY(a22, half);

    DYN_ARRAY(b11, half);
    DYN_ARRAY(b12, half);
    DYN_ARRAY(b21, half);
    DYN_ARRAY(b22, half);

    DYN_ARRAY(c11, half);
    DYN_ARRAY(c12, half);
    DYN_ARRAY(c21, half);
    DYN_ARRAY(c22, half);

    DYN_ARRAY(p1, half);
    DYN_ARRAY(p2, half);
    DYN_ARRAY(p3, half);
    DYN_ARRAY(p4, half);
    DYN_ARRAY(p5, half);
    DYN_ARRAY(p6, half);
    DYN_ARRAY(p7, half);

    /* Fill four sub-matrices */
    for (int i = 0; i < half; i++) {
        for (int j = 0; j < half; j++) {
            a11[i*half + j] = a[i*n + j];
            a12[i*half + j] = a[i*n + j + half];
            a21[i*half + j] = a[(i + half)*n + j];
            a22[i*half + j] = a[(i + half)*n + j + half];

            b11[i*half + j] = b[i*n + j];
            b12[i*half + j] = b[i*n + j + half];
            b21[i*half + j] = b[(i + half)*n + j];
            b22[i*half + j] = b[(i + half)*n + j + half];
        }
    }
    /* p1 = (a11+a22)*(b11+b22) */
    sum(half, a11, a22, a_res);
    sum(half, b11, b22, b_res);
    strassen_r(half, a_res, b_res, p1);

    /* p2 = (a21+a22)*b11 */
    sum(half, a21, a22, a_res);
    strassen_r(half, a_res, b11, p2);

    /* p3 = a11*(b12 - b22) */
    subtr(half, b12, b22, b_res);
    strassen_r(half, a11, b_res, p3);

    /* p4 = a22*(b21 - b11) */
    subtr(half, b21, b11, b_res);
    strassen_r(half, a22, b_res, p4);

    /* p5 = (a11 + a12)*b22 */
    sum(half, a11, a12, a_res);
    strassen_r(half, a_res, b22, p5);

    /* p6 = (a21 - a11)*(b11 + b12) */
    subtr(half, a21, a11, a_res);
    sum(half, b11, b12, b_res);
    strassen_r(half, a_res, b_res, p6);

    /* p7 = (a12 - a22)*(b21 + b22) */
    subtr(half, a12, a22, a_res);
    sum(half, b21, b22, b_res);
    strassen_r(half, a_res, b_res, p7);



    /* c11 = p1 + p4 - p5 + p7 */
    sum(half, p1, p4, a_res);
    sum(half, a_res, p7, b_res);
    subtr(half, b_res, p5, c11);
    /* c12 = p3 + p5 */
    sum(half, p3, p5, c12);
    /* c21 = p2 + p4 */
    sum(half, p2, p4, c21);
    /* c22 = p1 + p3 - p2 + p6 */
    sum(half, p1, p3, a_res);
    sum(half, a_res, p6, b_res);
    subtr(half, b_res, p2, c22);

    data *ac11 = __builtin_assume_aligned(c11, 16);
    data *ac12 = __builtin_assume_aligned(c12, 16);
    data *ac21 = __builtin_assume_aligned(c21, 16);
    data *ac22 = __builtin_assume_aligned(c22, 16);
    data *ares = __builtin_assume_aligned(res, 16);

    for (int i = 0; i < half; i++) {
        for (int j = 0; j < half; j++) {
            ares[i*n + j] = ac11[i*half + j];
            ares[i*n + j + half] = ac12[i*half + j];
            ares[(i + half)*n + j] = ac21[i*half + j];
            ares[(i + half)*n + j + half] = ac22[i*half + j];
        }
    }
    free(a_res);
    free(b_res);
    free(a11);
    free(a12);
    free(a21);
    free(a22);

    free(b11);
    free(b12);
    free(b21);
    free(b22);

    free(c11);
    free(c12);
    free(c21);
    free(c22);

    free(p1);
    free(p2);
    free(p3);
    free(p4);
    free(p5);
    free(p6);
    free(p7);
}

static unsigned int next_pow2(int n)
{
    return pow(2, ceil(log(n)/log(2)));

}
/* Helper function for cases when n is not power of 2 */
static void
strassen_mul(const size_t n, const data *restrict a, const data *restrict b,
             data *restrict res)
{
    size_t m = next_pow2(n);
    DYN_ARRAY(a_tmp, m);
    DYN_ARRAY(b_tmp, m);
    DYN_ARRAY(res_tmp, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a_tmp[i*m + j] = a[i*n + j];
            b_tmp[i*m + j] = b[i*n + j];
        }
    }

    strassen_r(m, a_tmp, b_tmp, res_tmp);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            res[i*n+j] = res_tmp[i*m+j];
    free(a_tmp);
    free(b_tmp);
    free(res_tmp);
}
/* Returns number of successfully read items */
static size_t matrix_read(const size_t n, data *m)
{
    int res = 0;
    for (int i = 0; i < n*n; i++) {
        res = scanf("%lf", &m[i]);
        if (!res)
            return 0;
    }
    return n*n;
}
/* Reads two matrices from stdin */
static int input_read(const int n, data *a, data *b)
{
    size_t res = matrix_read(n, a);
    if (res != n*n) {
        fprintf(stderr, "Error reading matrix a: read %zu items, expected %d\n", res, n*n);
        return EXIT_FAILURE;
    }
    res = matrix_read(n, b);
    if (res != n*n) {
        fprintf(stderr, "Error reading matrix b: read %zu items\n", res);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
