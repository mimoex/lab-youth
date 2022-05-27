#pragma once
#include <iostream>
#include <gmpxx.h>
#include <random>

using namespace std;

//有限体の加算	c=a+b%p
void add_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//有限体の減算	c=a-b%p
void sub_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//有限体の乗算	c=a*b%p
void mul_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//有限体のべき乗	c=a^b%p
void pow_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//有限体の逆元	c=a^-1
void inv_fp(const mpz_class a, mpz_class** c, const mpz_class p);
//有限体の除算	c=a/b%p
void div_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);

