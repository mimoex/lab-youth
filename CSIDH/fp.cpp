#include "fp.h"

mpz_class mod(mpz_class a, mpz_class b) {
	return (a %= b) < 0 ? a + b : a;
}

void add_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a + b;
	if (*c >= p) *c -= p;
}

void sub_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	if (a >= b) {
		*c = a - b;
	}
	else {
		*c = a + p - b;
	}
}


void mul_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	*c = a * b;
	*c = *c % p;
}


void pow_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	//*c = a;
	mpz_powm(c->get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), p.get_mpz_t());
}

//逆元
void inv_fp(const mpz_class& a, const mpz_class& p, mpz_class* c)
{
	if (a == 0) {
		cout << "zero inv" << endl;
		*c = 0;

	}
	else {
		pow_fp(a, p - 2, p, c);
	}

}

//input  A, p; A∈[1,p-1]
//output A^-1 mod p
int algo5_inv(const mpz_class& A, const mpz_class& p, mpz_class* c)
{
	mpz_class u, v, r, s;


	u = p; v = A; r = 0; s = 1;

	if (A == 0) {
		cout << "zero inv" << endl;
		return -1;
	}
	else {


		while (u != 1 && v != 1) {
			while ((u & 1) == 0) {
				u >>= 1;
				if ((r & 1) == 0) {
					r >>=1;
				}
				else {
					r += p;
					r >>= 1;
				}
			}
			while ((v & 1) == 0) {
				v >>= 1;
				if ((s & 1) == 0) {
					s >>= 1;
				}
				else {
					s += p;
					s>>= 1;
				}
			}
			if (u >= v) {
				u -= v;
				r -= s;
			}
			else {
				v -= u;
				s -= r;
			}
		}
		if (u == 1) {
			*c =mod(r , p);
		}
		else {
			*c =mod(s , p);
		}
	}
	return 0;
}


void div_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c)
{
	algo5_inv(b, p, c);
	//inv_fp(b, p, c);
	mul_fp(a, *c, p, c);

}

//有限体pから乱数生成
mpz_class fp_random(const mpz_class& p)
{
	mpz_class x, cnt;
	size_t n = 512;

	random_device rnd;
	gmp_randclass r(gmp_randinit_default);
	r.seed(rnd());

	while (1) {
		x = r.get_z_bits(n);
		cnt++;

		if (x < p) {
			return x;
		}
	}
}