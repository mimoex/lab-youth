#pragma once
#include <iostream>
#include <gmpxx.h>
#include <random>

using namespace std;

//�L���̂̉��Z	c=a+b%p
void add_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//�L���̂̌��Z	c=a-b%p
void sub_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//�L���̂̏�Z	c=a*b%p
void mul_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//�L���ׂ̂̂���	c=a^b%p
void pow_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);
//�L���̂̋t��	c=a^-1
void inv_fp(const mpz_class a, mpz_class** c, const mpz_class p);
//�L���̂̏��Z	c=a/b%p
void div_fp(const mpz_class a, const mpz_class b, mpz_class* c, const mpz_class p);

