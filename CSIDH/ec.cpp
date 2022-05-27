#include "fp.h"


Point ec_add(Point& p, Point& q, mpz_class& a, mpz_class& mod)
{
	if (p.inf) return q;
	if (q.inf) return p;

	mpz_class lambda, y_1;
	mpz_class lh, rh;
	Point result;

	if (p.x == q.x) {
		//P + (-P) = 0
		if (p.y == -q.y) {
			p.inf = true;
			return p;
		}
		//lambda=(3*x1*x1+a)/(2*y1)
		mul_fp(p.x, p.x, mod, &lambda);
		mul_fp(3, lambda, mod, &lambda);
		add_fp(lambda, a, mod, &lambda);

		mul_fp(2, p.y, mod, &y_1);

		div_fp(lambda, y_1, mod, &lambda);
	}
	else {
		//x1!=x2�̂Ƃ��C
		//lambda=(y2-y1)/(x2-x1)
		sub_fp(p.y, q.y, mod, &lh);
		sub_fp(p.x, q.x, mod, &rh);
		div_fp(lh, rh, mod, &lambda);
	}
	
	//P+Q��x���W
	mpz_class lambda2, x_temp;
	mul_fp(lambda, lambda, mod, &lambda2);
	add_fp(p.x, q.x, mod, &x_temp);
	sub_fp(lambda2, x_temp, mod, &result.x);

	//P+Q��y���W
	sub_fp(p.x, result.x, mod, &x_temp);
	mul_fp(lambda, x_temp, mod, &lambda2);
	sub_fp(lambda2, p.y, mod, &result.y);

	return result;
}
