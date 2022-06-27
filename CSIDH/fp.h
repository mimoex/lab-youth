#pragma once
#include <iostream>
#include <gmpxx.h>
#include <random>
#include <chrono>

using namespace std;

const size_t N = 74;	//log(11^74)/log(2) ≒ 255.9979

const int primes[N] = {	3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
						67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
						137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
						199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
						277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
						359, 367, 373, 587 };

struct Point {
	mpz_class x, y;
	bool inf = false;
};

//有限体の加算	c=a+b%p
void add_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体の減算	c=a-b%p
void sub_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体の乗算	c=a*b%p
void mul_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体のべき乗	c=a^b%p
void pow_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//有限体の逆元	c=a^-1
void inv_fp(const mpz_class& a, const mpz_class& p, mpz_class* c);
//有限体の除算	c=a/b%p
void div_fp(const mpz_class& a, const mpz_class& b, const mpz_class& p, mpz_class* c);
//p以下の乱数を生成
mpz_class random_fp(const mpz_class& p);


Point copy_point(Point& p);

//Weierstrass曲線の加算公式
Point ec_add(Point& p, Point& q, const mpz_class& a, const mpz_class& mod);
//Montgomery曲線の加算公式
Point montgomery_ec_add(Point& p, Point& q, const mpz_class& a, const mpz_class& b, const mpz_class& mod);

//Montgomery曲線のスカラー倍算(バイナリ法)
Point l_to_r_binary_mont(const Point& p, const mpz_class& a, const mpz_class& b, const mpz_class& mod, const mpz_class& n);

//楕円曲線上の点を生成
Point gen_point(const mpz_class& a, const mpz_class& b, const mpz_class& mod);
size_t check_point(const Point& p, const mpz_class& mod);