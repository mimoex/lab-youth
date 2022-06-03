#include "fp.h"


Point ec_add(Point& p, Point& q, const mpz_class& a, const mpz_class& mod)
{
	if (p.inf) return q;
	if (q.inf) return p;

	mpz_class lambda_temp, y_1, lambda;
	mpz_class lh, rh;
	Point result;

	if (p.x == q.x) {
		//P + (-P) = 0
		if (p.y == -q.y) {
			p.inf = true;
			return p;
		}
		//lambda=(3*x1*x1+a)/(2*y1)
		mul_fp(p.x, p.x, mod, &lambda_temp);
		mul_fp(3, lambda_temp, mod, &lambda_temp);
		add_fp(lambda_temp, a, mod, &lambda_temp);

		add_fp(p.y, p.y, mod, &y_1);

		div_fp(lambda_temp, y_1, mod, &lambda);
	}
	else {
		//x1!=x2のとき，
		// cout << "Add" << endl;
		//lambda=(y2-y1)/(x2-x1)
		sub_fp(p.y, q.y, mod, &lh);
		sub_fp(p.x, q.x, mod, &rh);
		div_fp(lh, rh, mod, &lambda);
	}
	
	//P+Qのx座標
	mpz_class lambda2, x_temp;
	mul_fp(lambda, lambda, mod, &lambda2);
	add_fp(p.x, q.x, mod, &x_temp);
	sub_fp(lambda2, x_temp, mod, &result.x);

	//P+Qのy座標
	sub_fp(p.x, result.x, mod, &x_temp);
	
	mul_fp(lambda, x_temp, mod, &lambda2);
	sub_fp(lambda2, p.y, mod, &result.y);

	return result;
}



//Montgomery Curve
Point montgomery_ec_add(Point& p, Point& q, const mpz_class& a, const mpz_class& b, const mpz_class& mod)
{
	if (p.inf) { return q; cout << "p is inf" << endl; }
	if (q.inf) { return p; cout << "q is inf" << endl; }

	mpz_class x_temp, y_1, lambda;
	mpz_class lh, rh;
	Point result;


	//cout << "Px,Qx:" << p.x << "," << q.x << endl;
	if (p.x == q.x) {
		//P + (-P) = 0
		//cout << "Py,-Py:" << p.y << "," << -q.y << endl;
		if (p.y == -q.y) {
			p.inf = true;
			return p;
		}
		//x1!=x2のとき，
		//
		//x3 = b*(3*x1^2+2*a*x1+1)^2/(2*b*y1)^2-a-x1-x1

		//y3 = (2*x1+x1+a)*(3*x1^2+2*a*x1+1)/(2*b*y1)-b*(3*x1^2+2*a*x1+1)^3/(2*b*y1)^3-y1

	}
	else {
		//x1!=x2のとき，
		// 
		// x3 = b*(y2-y1)^2/(x2-x1)^2-a-x1-x2

		mpz_class llh, lrh, y_ltemp, subx2x1, suby2y1;

		sub_fp(q.x, p.x, mod, &subx2x1);	//x2-x1
		sub_fp(q.y, p.y, mod, &suby2y1);	//y2-y1

		mul_fp(suby2y1, suby2y1, mod, &lh);
		mul_fp(b, lh, mod, &lh);

		mul_fp(subx2x1, subx2x1, mod, &rh);

		div_fp(lh, rh, mod, &x_temp);

		sub_fp(x_temp, a, mod, &x_temp);
		sub_fp(x_temp, p.x, mod, &x_temp);
		sub_fp(x_temp, q.x, mod, &result.x);


		// y3 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-b*(y2-y1)^3/(x2-x1)^3-y1
		//
		//左辺
		mul_fp(2, p.x, mod, &llh);
		add_fp(llh, q.x, mod, &llh);
		add_fp(llh, a, mod, &llh);
		mul_fp(llh, suby2y1, mod, &y_ltemp);

		div_fp(y_ltemp, subx2x1, mod, &y_ltemp);

		//右辺
		pow_fp(suby2y1, 3, mod, &lh);
		mul_fp(b, lh, mod, &lh);

		pow_fp(subx2x1, 3, mod, &rh);

		div_fp(lh, rh, mod, &lrh);

		sub_fp(y_ltemp, lrh, mod, &y_ltemp);
		sub_fp(y_ltemp, q.y, mod, &result.y);
	}
	return result;
}




Point copy_point(Point &p) {
	Point result;
	result.x= p.x;
	result.y= p.y;
	result.inf = p.inf;
	return p;
}


mpz_class gen_prime()
{
	int i;
	mpz_class p, qq;
	qq = 4;
	for (i = 0; i < N; i++) {
		qq *= prime[i];
	}
	p = qq - 1;
	return p;
}

