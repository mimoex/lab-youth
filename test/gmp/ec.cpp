#include <iostream>
#include <gmp.h>

typedef struct {
	mpz_t x;
	mpz_t y;
	bool isInf;
} Point;

/*---点の初期化---*/
Point init_point(Point p) {
	mpz_init(p.x);
	mpz_init(p.y);
	p.isInf = false;
	return p;
}


Point clean_point(Point p) {
	mpz_clear(p.x);
	mpz_clear(p.y);
	return p;
}


int hikaku_point(Point p, Point q) {
	if ((mpz_cmp(p.x, q.x) == 0) && (mpz_cmp(p.y, q.y) == 0)) return 0;
	else return 1;
}

int hikaku_point_neg(Point p, Point q) {
	mpz_neg(q.y, q.y);
	return hikaku_point(p, q);
}



/*モジュロ演算*/
void positive_modulo(mpz_t results, mpz_t a, mpz_t modulo) {
	mpz_t zero_value;
	mpz_init(zero_value);

	mpz_tdiv_r(results, a, modulo);
	if (mpz_cmp(results, zero_value) < 0) mpz_add(results, results, modulo);
}

/*---点の2倍算---*/
Point affine_doubling(Point p, mpz_t a, mpz_t modulo) {
	Point results{};
	results = init_point(results);

	if (p.isInf) return p;

	// y=0のとき
	mpz_t zero_value;
	mpz_init(zero_value);
	if (mpz_cmp(p.y, zero_value) == 0) {
		results.isInf = true;
		return results;
	}

	mpz_t lambda, mul, lh, rh;
	mpz_init(lambda); mpz_init(mul);
	mpz_init(lh); mpz_init(rh);

	/** Create (2y)^-1 % m */
	mpz_mul_2exp(lh, p.y, 1);
	// Assuming inverse exists (retval = 1)
	mpz_invert(lh, lh, modulo);

	// (3*(x1)^2)+a / 2*y1
	mpz_set_ui(mul, 3);
	mpz_mul(rh, p.x, p.x);
	positive_modulo(rh, rh, modulo); // x^2 % m
	mpz_mul(rh, rh, mul); // (3x^2)
	mpz_add(rh, rh, a); // (3x^2 + a)
	positive_modulo(rh, rh, modulo); // (3x^2 + a) % m

	mpz_mul(lambda, lh, rh);
	positive_modulo(lambda, lambda, modulo); // (3x^2 + a) . (2y)^-1 % m

	// x3 = lambda ^ 2 - 2 * x1
	mpz_mul(results.x, lambda, lambda);
	positive_modulo(results.x, results.x, modulo);
	mpz_sub(results.x, results.x, p.x);
	mpz_sub(results.x, results.x, p.x);
	positive_modulo(results.x, results.x, modulo);

	// y3=lambda(x1-x3) -y1
	mpz_sub(results.y, p.x, results.x);
	mpz_mul(results.y, results.y, lambda);
	mpz_sub(results.y, results.y, p.y);
	positive_modulo(results.y, results.y, modulo);

	return results;
}

Point affine_curve_addition(Point p, Point q, mpz_t a, mpz_t modulo) {
	//無限遠点
	if (q.isInf) return p;
	if (p.isInf) return q;

	Point results{};
	results = init_point(results);

	// P = -Qのとき
	if (hikaku_point_neg(p, q) == 0) {
		results.isInf = true;
		return results;
	}
	// P = Qのとき
	if (hikaku_point(p, q) == 0) return affine_doubling(p, a, modulo);

	mpz_t lambda, lh, rh;
	mpz_init(lambda);
	mpz_init(lh); mpz_init(rh);

	/* (q.x - p.x)の逆元mod */
	mpz_sub(lh, q.x, p.x);
	mpz_invert(lh, lh, modulo);

	mpz_sub(rh, q.y, p.y);
	/* λ= (q.y - p.y)*(q.x - p.x)逆元mod*/
	mpz_mul(lambda, lh, rh);
	positive_modulo(lambda, lambda, modulo);

	/* r.x = lambda^2 - p.x - q.x */
	mpz_mul(results.x, lambda, lambda);
	mpz_sub(results.x, results.x, p.x);
	mpz_sub(results.x, results.x, q.x);
	positive_modulo(results.x, results.x, modulo);

	/* r.y = lambda * (p.x - new_x) - p.y */
	mpz_sub(results.y, p.x, results.x);
	mpz_mul(results.y, results.y, lambda);
	mpz_sub(results.y, results.y, p.y);
	positive_modulo(results.y, results.y, modulo);

	return results;
}



// secp256k1(https://neuromancer.sk/std/secg/secp256k1)
const char* a_hex = "0000000000000000000000000000000000000000000000000000000000000000";
const char* b_hex = "0000000000000000000000000000000000000000000000000000000000000007";
const char* p_hex = "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f";
const char* gx_hex = "79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798";
const char* gy_hex = "483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8";
const char* order_hex = "ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551";

int main(void)
{
	mpz_t a, b, k, r, modulo;
	Point p{}, q{}, result_add[21], result_double[21];

	p = init_point(p);
	mpz_init(a);
	mpz_init(b);
	mpz_init(k);
	mpz_init(r); // order
	mpz_init(modulo);

	mpz_set_str(a, a_hex, 16);
	mpz_set_str(b, b_hex, 16);
	mpz_set_str(modulo, p_hex, 16);
	mpz_set_str(p.x, gx_hex, 16);
	mpz_set_str(p.y, gy_hex, 16);
	mpz_set_str(r, order_hex, 16);

	


	result_double[0] = p;
	result_add[0] = p;
	for (int i = 0; i < 20; i++) {
		result_double[i+1] = affine_doubling(result_double[i], a, modulo);
		gmp_printf("%dP   = (%Zd, %Zd)\n", (i+1)*2, result_double[i+1].x, result_double[i + 1].y);
		std::cout << result_double[i + 1].x<<" " << result_double[i + 1].y << std::endl;

		result_add[i + 1] = affine_curve_addition(result_add[i], p, a, modulo);
		gmp_printf("%dP+P   = (%Zd, %Zd)\n", i+1, result_add[i + 1].x, result_add[i + 1].y);

	}

	return 0;

}