#ifndef MPG_H
#define MPG_H

#ifndef __GMP_H__
#include <gmp.h>
#endif /* __GMP_H__ */

typedef struct {
        mpz_t real;
        mpz_t imaginary;
} _mpg_struct;

typedef _mpg_struct mpg_t[1];
typedef _mpg_struct *mpg_ptr;

/* Primary decomposition of a Gaussian integers: any Gaussian integer
 * can be written as the form i^(unitexp) * (1 + i)^(evenexp) *
 * primary, where primary is a primary Gaussian integer,
 * i.e. congruent to 1 modulo 2 + 2i. */
typedef struct {
        unsigned int unitexp;
        unsigned long int evenexp;
        mpg_t primary;
} _mpg_pdec_struct;

typedef _mpg_pdec_struct mpg_pdec_t[1];
typedef _mpg_pdec_struct *mpg_pdec_ptr;



#endif
