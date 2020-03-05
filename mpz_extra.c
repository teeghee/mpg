#include <gmp.h>

void
mpz_nearest_z (mpz_t rop, const mpz_t n, const mpz_t d)
{
        /* Find the integer rop closest to the quotient n/d. */
        mpz_t frem, bound, tmp;
        mpz_inits(frem, bound, tmp, (mpz_ptr) NULL);

        mpz_fdiv_qr(rop, frem, n, d);
        mpz_abs(tmp, d);
        mpz_fdiv_q_ui(bound, tmp, 2);

        if (mpz_cmp(frem, bound) >= 0) {
                mpz_add_ui(rop, rop, 1);
        }

        mpz_clears(frem, bound, tmp, (mpz_ptr) NULL);
}

/* int main(void) */
/* { */
/*         mpz_t nearest, n, d; */
/*         mpz_inits(nearest, n, d, (mpz_ptr) NULL); */

/*         mpz_set_si(n, -771327); */
/*         mpz_set_ui(d, 6); */

/*         mpz_nearest_z(nearest, n, d); */
/*         gmp_printf("%Zd\n", nearest); */

/*         mpz_clears(nearest, n, d, (mpz_ptr) NULL); */
/* } */
