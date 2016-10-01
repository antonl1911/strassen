#ifndef STRASSEN_H
#define STRASSEN_H

#define align(x) __attribute__((align(x)))
typedef double data;
/* Allocate flat memory chunk, to help compiler optimize loops */
#define DYN_ARRAY(name, size) data *name = calloc(size*size, sizeof *name); \
                                    if (!name) \
                                        abort();

static void
strassen_mul(const size_t n, const data *restrict a, const data *restrict b,
		  data *restrict res);
static unsigned int next_pow2(int n);
static void print_matrix(size_t n, data *matrix);
static void sum(const size_t n, const data *restrict a,
    const data *restrict b, data *restrict res);
static int input_read(const int n, data *a, data *b);


#endif
