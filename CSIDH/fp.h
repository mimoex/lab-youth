#pragma once
#include <iostream>
#include <gmpxx.h>
#include <random>

using namespace std;

//—LŒÀ‘Ì‚Ì‰ÁZ	c=a+b%p
void add_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//—LŒÀ‘Ì‚ÌŒ¸Z	c=a-b%p
void sub_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//—LŒÀ‘Ì‚ÌæZ	c=a*b%p
void mul_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//—LŒÀ‘Ì‚Ì‚×‚«æ	c=a^b%p
void pow_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//—LŒÀ‘Ì‚Ì‹tŒ³	c=a^-1
void inv_fp(const mpz_class a, mpz_class** c, const mpz_class p);
//—LŒÀ‘Ì‚ÌœZ	c=a/b%p
void div_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);

