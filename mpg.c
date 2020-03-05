#include <stdio.h>
#include <stdarg.h>
#include <gmp.h>

#include "mpg.h"
#include "mpz_extra.h"


void
mpg_init (mpg_t x)
{
        /* Initialize x, and set its value to 0. */
        mpz_init(x->real);
        mpz_init(x->imaginary);
}

void
mpg_inits (mpg_ptr x, ...)
{
        va_list ap;
        va_start(ap, x);

        while (x != NULL) {
                mpg_init(x);
                x = va_arg(ap, mpg_ptr);
        }

        va_end(ap);
}

void
mpg_pdec_init (mpg_pdec_t dec)
{
        mpg_init(dec->primary);
}


void
mpg_pdec_inits (mpg_pdec_ptr x, ...)
{
        va_list ap;
        va_start(ap, x);

        while (x != NULL) {
                mpg_pdec_init(x);
                x = va_arg(ap, mpg_pdec_ptr);
        }

        va_end(ap);
}

void
mpg_clear (mpg_t x)
{
        /* Free the space occupied by x.  Call this function for all
         * 'mpg_t' variables when you are done with them. */
        mpz_clear(x->real);
        mpz_clear(x->imaginary);
}

void
mpg_clears (mpg_ptr x, ...)
{
        va_list ap;
        va_start(ap, x);

        while (x != NULL) {
                mpg_clear(x);
                x = va_arg(ap, mpg_ptr);
        }

        va_end(ap);
}

void
mpg_pdec_clear (mpg_pdec_t dec)
{
        mpg_clear(dec->primary);
}

void
mpg_pdec_clears (mpg_pdec_ptr x, ...)
{
        va_list ap;
        va_start(ap, x);

        while (x != NULL) {
                mpg_pdec_clear(x);
                x = va_arg(ap, mpg_pdec_ptr);
        }

        va_end(ap);
}

void
mpg_set (mpg_t rop, const mpg_t op)
{
        /* Set the value of rop to op. */
        mpz_set(rop->real, op->real);
        mpz_set(rop->imaginary, op->imaginary);
}

void
mpg_set_si (mpg_t rop, const long int a, const long int b)
{
        /* Set the value of rop to a + bi */
        mpz_set_si(rop->real, a);
        mpz_set_si(rop->imaginary, b);
}

void
mpg_init_set (mpg_t rop, const mpg_t op)
{
        mpg_init(rop);
        mpg_set(rop, op);
}

void
mpg_swap (mpg_t x, mpg_t y)
{
        mpg_t tmp;
        mpg_init_set(tmp, x);
        mpg_set(x, y);
        mpg_set(y, tmp);
        mpg_clear(tmp);
}

void
mpg_add (mpg_t rop, const mpg_t op1, const mpg_t op2)
{
        mpz_add(rop->real, op1->real, op2->real);
        mpz_add(rop->imaginary, op1->imaginary, op2->imaginary);
}

void
mpg_sub (mpg_t rop, const mpg_t op1, const mpg_t op2)
{
        mpz_sub(rop->real, op1->real, op2->real);
        mpz_sub(rop->imaginary, op1->imaginary, op2->imaginary);
}

void
mpg_mul (mpg_t rop, const mpg_t op1, const mpg_t op2)
{
        mpg_t tmp;
        mpg_init(tmp);

        mpz_mul(tmp->real, op1->real, op2->real);
        mpz_submul(tmp->real, op1->imaginary, op2->imaginary);
        mpz_mul(tmp->imaginary, op1->real, op2->imaginary);
        mpz_addmul(tmp->imaginary, op1->imaginary, op2->real);

        mpg_set(rop, tmp);

        mpg_clear(tmp);
}

void
mpg_mul_si (mpg_t rop, const mpg_t op1,
            const long int op2r, const long int op2i)
{
        mpg_t tmp;
        mpg_init(tmp);
        mpg_set_si(tmp, op2r, op2i);

        mpg_mul(rop, op1, tmp);

        mpg_clear(tmp);
}

void
mpg_submul (mpg_t rop, const mpg_t op1, const mpg_t op2)
{
        mpg_t tmp;
        mpg_init(tmp);
        mpg_mul(tmp, op1, op2);
        mpg_sub(rop, rop, tmp);
        mpg_clear(tmp);
}

#define mpg_zero_p(x) ((mpz_sgn((x)->real) == 0 && mpz_sgn((x)->imaginary) == 0) ? 1 : 0)

void
mpg_print (const mpg_t x)
{
        int sgnreal = mpz_sgn(x->real);
        int sgnimag = mpz_sgn(x->imaginary);

        if (sgnreal == 0 && sgnimag == 0) {
                gmp_printf("0\n");
        } else if (sgnreal != 0) {
                gmp_printf("%Zd", x->real);
                if (sgnimag == 1) {
                        gmp_printf("+");
                }
        }

        if (sgnimag != 0) {
                gmp_printf("%ZdI\n", x->imaginary);
        }
}

void
mpg_pdec_print (const mpg_pdec_t dec)
{
        gmp_printf("I^%u * (1 + I)^%lu * ", dec->unitexp, dec->evenexp);
        mpg_print(dec->primary);
}

int
mpg_equal_p (const mpg_t x, const mpg_t y)
{
        return (mpz_cmp(x->real, y->real) == 0 &&
                mpz_cmp(x->imaginary, y->imaginary) == 0);
}

void
mpg_norm (mpz_t norm, const mpg_t op)
{
        /* Set x = Norm(y) */
        mpz_mul(norm, op->real, op->real);
        mpz_addmul(norm, op->imaginary, op->imaginary);
}

int
mpg_norm_cmp (const mpg_t x, const mpg_t y)
{
        int result;
        mpz_t xnorm, ynorm;
        mpz_inits(xnorm, ynorm, (mpz_ptr) NULL);
        mpg_norm(xnorm, x);
        mpg_norm(ynorm, y);
        result = mpz_cmp(xnorm, ynorm);
        mpz_clears(xnorm, ynorm, (mpz_ptr) NULL);
        return result;
}


void
mpg_conj (mpg_t rop, const mpg_t op)
{
        /* Sets the value of rop to be the conjugate of op. */
        mpg_set(rop, op);
        mpz_neg(rop->imaginary, rop->imaginary);
}

void
mpg_div_qr (mpg_t q, mpg_t r, const mpg_t n, const mpg_t d)
{
        /* - q and r are the quotient and remainder of n divided by
             d. */
        /* - r satisfies 0 <= 2 * norm(r) <= norm(d). */
        mpg_t tmpq, tmpr;
        mpg_t tmpg;
        mpz_t tmpz, norm;
        mpg_inits(tmpq, tmpr, tmpg, (mpg_ptr) NULL);
        mpz_inits(tmpz, norm, (mpz_ptr) NULL);

        mpg_conj(tmpg, d);
        mpg_norm(norm, d);
        mpg_mul(tmpg, n, tmpg);

        mpz_nearest_z(tmpz, tmpg->real, norm);
        mpz_set(tmpq->real, tmpz);
        mpz_nearest_z(tmpz, tmpg->imaginary, norm);
        mpz_set(tmpq->imaginary, tmpz);

        mpg_set(tmpr, n);
        mpg_submul(tmpr, tmpq, d);

        mpg_set(q, tmpq);
        mpg_set(r, tmpr);

        mpg_clears(tmpg, tmpq, tmpr, (mpg_ptr) NULL);
        mpz_clears(tmpz, norm, (mpz_ptr) NULL);
}

void
mpg_div_q (mpg_t q, const mpg_t n, const mpg_t d)
{
        mpg_t tmpr;
        mpg_init(tmpr);
        mpg_div_qr(q, tmpr, n, d);
        mpg_clear(tmpr);
}

void
mpg_div_r (mpg_t r, const mpg_t n, const mpg_t d)
{
        mpg_t tmpq;
        mpg_init(tmpq);
        mpg_div_qr(tmpq, r, n, d);
        mpg_clear(tmpq);
}

#define mpg_even_p(x) ((mpz_even_p((x)->real) && mpz_even_p((x)->imaginary)) || (mpz_odd_p((x)->real) && mpz_odd_p((x)->imaginary)))

int
mpg_divisible_p (const mpg_t n, const mpg_t d)
{
        mpg_t rem;
        mpg_init(rem);
        mpg_div_r(rem, n, d);
        int result = mpg_zero_p(rem);
        mpg_clear(rem);
        return result;
}

int
mpg_congruent (const mpg_t x, const mpg_t y, const mpg_t modulus)
{
        mpg_t diff;
        mpg_init(diff);
        mpg_sub(diff, x, y);
        int result = mpg_divisible_p(diff, modulus);
        mpg_clear(diff);
        return result;
}

int
mpg_congruent_si (const mpg_t x, const long int y1, const long int y2,
                  const mpg_t modulus)
{
        int result;
        mpg_t y;
        mpg_init(y);

        mpg_set_si(y, y1, y2);
        result = mpg_congruent(x, y, modulus);

        mpg_clear(y);
        return result;
}

void
mpg_pow_ui (mpg_t rop, const mpg_t base, unsigned long int exp)
{
        if (exp == 0) {
                mpg_set_si(rop, 1, 0);
                return;
        }

        mpg_t tmp;
        mpg_init(tmp);
        mpg_set(rop, base);
        mpg_set_si(tmp, 1, 0);

        while (exp > 1) {
                if (exp & 1) {
                        // if exp is odd.
                        mpg_mul(tmp, rop, tmp);
                        mpg_mul(rop, rop, rop);
                        exp = (exp - 1) / 2;
                } else {
                        // if exp is even.
                        mpg_mul(rop, rop, rop);
                        exp = exp / 2;
                }
        }

        mpg_mul(rop, rop, tmp);
        mpg_clear(tmp);
}

void
mpg_primary_decompose (mpg_pdec_t rop, const mpg_t op)
{
        /* Decompose op into the product of a unit, a power of (1+i),
         * and a primary element. */
        mpg_t tmp;
        mpg_init(tmp);
        mpg_set_si(tmp, 1, 1);

        rop->evenexp = 0;
        mpg_set(rop->primary, op);
        while (mpg_even_p(rop->primary)) {
                (rop->evenexp)++;
                mpg_div_q(rop->primary, rop->primary, tmp);
        }

        mpg_set_si(tmp, -2, 2);
        for (unsigned int i = 0; i < 4; i++) {
                rop->unitexp = i;

                if (mpg_congruent_si(rop->primary, 1, 0, tmp)) {
                        break;
                }

                mpg_mul_si(rop->primary, rop->primary, 0, -1);
        }

        mpg_clear(tmp);
}

void
mpg_gcd(mpg_t g, const mpg_t x, const mpg_t y)
{
        if (mpg_zero_p(x) || mpg_zero_p(y)) {
                mpg_set_si(g, 0, 0);
                return;
        }

        mpg_pdec_t xdec, ydec;
        mpg_pdec_inits(xdec, ydec, (mpg_pdec_ptr) NULL);
        mpg_primary_decompose(xdec, x);
        mpg_primary_decompose(ydec, y);

        unsigned int expdiff;
        if (xdec->evenexp < ydec->evenexp) {
                expdiff = xdec->evenexp;
        } else {
                expdiff = ydec->evenexp;
        }
        mpg_t tmp;
        mpg_init(tmp);
        mpg_set_si(tmp, 1, 1);
        mpg_pow_ui(g, tmp, expdiff);

        mpg_t xnew, ynew;
        mpg_init(xnew);
        mpg_init(ynew);
        mpg_set(xnew, xdec->primary);
        mpg_set(ynew, ydec->primary);

        mpg_t z;
        mpg_init(z);
        mpg_pdec_t zdec;
        mpg_pdec_init(zdec);

        while (!mpg_equal_p(xnew, ynew)) {
                mpg_sub(z, xnew, ynew);
                mpg_primary_decompose(zdec, z);
                if (mpg_norm_cmp(xnew, ynew) > 0) {
                        mpg_set(xnew, zdec->primary);
                } else {
                        mpg_set(ynew, zdec->primary);
                }
        }

        mpg_mul(g, g, xnew);

        mpg_clears(tmp, xnew, ynew, z, (mpg_ptr) NULL);
        mpg_pdec_clears(xdec, ydec, zdec, (mpg_pdec_ptr) NULL);
}

void
mpg_quartic(mpg_t result, const mpg_t a, const mpg_t b)
{
        /* Computes the quartic residue symbol (a/b)_4.  Here d is
         * assumed to be prime to (1+i). */
        mpg_pdec_t adec, bdec;
        mpg_pdec_inits(adec, bdec, (mpg_pdec_ptr) NULL);
        mpg_primary_decompose(adec, a);
        mpg_primary_decompose(bdec, b);

        mpg_t tmp, part;
        mpg_inits(tmp, part, (mpg_ptr) NULL);
        mpg_set_si(tmp, 1, 0);
        mpg_sub(part, bdec->primary, tmp);
        mpg_set_si(tmp, 2, 2);
        mpg_div_q(part, part, tmp);

        unsigned long int m = mpz_fdiv_ui(part->real, 4);
        unsigned long int n = mpz_fdiv_ui(part->imaginary, 4);
        unsigned long int exp =
                ((n - m) * (adec->unitexp % 4) -
                 (n + (n + m) * (n + m)) * (adec->evenexp % 4)) % 4;

        mpg_t x, y;
        mpz_t xnorm, ynorm;
        mpg_inits(x, y, (mpg_ptr) NULL);
        mpz_inits(xnorm, ynorm, (mpz_ptr) NULL);
        mpg_set(x, adec->primary);
        mpg_set(y, bdec->primary);
        mpg_norm(xnorm, x);
        mpg_norm(ynorm, y);
        if (mpz_cmp(xnorm, ynorm) < 0) {
                mpg_swap(x, y);
                mpz_swap(xnorm, ynorm);
        }

        mpg_t z;
        mpg_pdec_t zdec;
        mpg_init(z);
        mpg_pdec_init(zdec);
        while (!mpg_equal_p(x, y)) {
                mpg_sub(z, x, y);
                mpg_primary_decompose(zdec, z);

                mpg_set_si(tmp, 1, 0);
                mpg_sub(part, y, tmp);
                mpg_set_si(tmp, 2, 2);
                mpg_div_q(part, part, tmp);
                m = mpz_fdiv_ui(part->real, 4);
                n = mpz_fdiv_ui(part->imaginary, 4);
                exp = (exp +
                       (n - m) * (zdec->unitexp % 4) -
                       (n + (n + m) * (n + m)) * (zdec->evenexp % 4)) % 4;

                mpg_set(x, zdec->primary);
                mpg_norm(xnorm, x);
                if (mpz_cmp(xnorm, ynorm) < 0) {
                        mpg_swap(x, y);
                        mpz_swap(xnorm, ynorm);
                }
        }

        mpg_set_si(tmp, 1, 0);
        if (!mpg_equal_p(x, tmp)) {
                mpg_set_si(result, 0, 0);
        } else {
                mpg_set_si(tmp, 0, 1);
                mpg_pow_ui(result, tmp, exp);
        }

        mpg_clears(tmp, part, x, y, z, (mpg_ptr) NULL);
        mpg_pdec_clears(adec, bdec, zdec, (mpg_pdec_ptr) NULL);
        mpz_clears(xnorm, ynorm, (mpz_ptr) NULL);
}

void
mpg_times_i (mpg_t rop, const mpg_t op)
{
        mpg_t tmp;
        mpg_init(tmp);
        mpg_set_si(tmp, 0, 1);
        mpg_mul(rop, op, tmp);
        mpg_clear(tmp);
}

int
mpg_associate_p (const mpg_t x, const mpg_t y)
{
        int result = 0;
        mpg_t tmp;
        mpg_init(tmp);
        mpg_set(tmp, y);
        for (int i = 0; i < 4; i++) {
                if (mpg_equal_p(x, tmp)) {
                        result = 1;
                        break;
                }
                mpg_times_i(tmp, tmp);
        }

        mpg_clear(tmp);
        return result;
}

int
mpg_probab_prime_p (const mpg_t n, int reps)
{
        int result;

        if (mpg_is_even(n)) {
                mpg_t tmp;
                mpg_init(tmp);
                mpg_set_si(tmp, 1, 1);
                if (mpg_associate_p(n, tmp)) {
                        result = 2;
                } else {
                        result = 0;
                }
                mpg_clear(tmp);
        } else {
                mpz_t norm;
                mpz_init(norm);
                mpg_norm(norm, n);
                result = mpz_probab_prime_p(norm, reps);
                if (result == 0 && mpz_perfect_square_p(norm)) {
                        mpz_sqrt(norm, norm);
                        if (mpz_congruent_ui_p(norm, 3, 4)) {
                                result = mpz_probab_prime_p(norm, reps);
                        }
                }
                mpz_clear(norm);
        }

        return result;
}

int main(void)
{

}
