#include "fp.h"

void add_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p)
{
	*c = a + b;
	*c %= p;
}

void sub_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p)
{
	*c = a - b;
	*c %= p;
}


void mul_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p)
{
	*c = a * b;
	*c %= p;
}


void pow_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p)
{
	mpz_powm(c->get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), p.get_mpz_t());
}


void inv_fp(const mpz_class a, mpz_class** c, const mpz_class p)
{
	if (a == 0) {
		cout << "zero inv" << endl;

	}
	else {
		pow_fp(a, p - 2, *c, p);
	}

}


void div_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p)
{
	inv_fp(b, &c, p);
	mul_fp(a, *c, c, p);

}
