#include "fp.h"

void add_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a + b;
	*c %= p;
}

void sub_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a - b;
	*c %= p;
}


void mul_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a * b;
	*c %= p;
}


void pow_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	mpz_powm(c->get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), p.get_mpz_t());
}


void inv_fp(const mpz_class& a, const mpz_class& p, mpz_class** c)
{
	if (a == 0) {
		cout << "zero inv" << endl;

	}
	else {
		pow_fp(a, p - 2, p, *c);
	}

}


void div_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	inv_fp(b, p, &c);
	mul_fp(a, *c, p, c);

}
