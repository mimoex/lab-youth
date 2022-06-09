#include "fp.h"

void add_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a + b;
	if (*c >= p) *c -= p;
	*c = *c % p;
}

void sub_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = (a - b);
	*c = (*c % p + p) % p;
}


void mul_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a * b;
	*c = *c%p;
}


void pow_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a;
	mpz_powm(c->get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), p.get_mpz_t());
}

//逆元
void inv_fp(const mpz_class& a, const mpz_class& p, mpz_class** c)
{
	if (a == 0) {
		cout << "zero inv" << endl;
		**c = 0;

	}
	else {
		pow_fp(a, p - 2, p, *c);
	}

}

//input  A, p; A∈[1,p-1]
//output A^-1 mod p
void algo5_inv(const mpz_class& A, const mpz_class& p, mpz_class** c)
{
	mpz_class u, v, r, s;


	u = p; v = A; r = 0; s = 1;

	if (A == 0) {
		cout << "zero inv" << endl;
	}
	else {

	
	while (u != 1 && v != 1) {
		while ((u & 1) == 0) {
			u = u / 2;
			if ((r & 1) == 0) {
				r = r / 2;
			}
			else {
				r = (r + p) / 2;
			}
		}
		while ((v & 1) == 0) {
			v = v / 2;
			if ((s & 1) == 0) {
				s = s / 2;
			}
			else {
				s = (s + p) / 2;
			}
		}
		if (u >= v) {
			u = u - v;
			r = r - s;
		}
		else {
			v = v - u;
			s = s - r;
		}
	}
	if (u == 1) {
		**c= r % p;
	}
	else {
		**c=s % p;
	}
	}
}


void div_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	algo5_inv(b, p, &c);
	mul_fp(a, *c, p, c);

}
