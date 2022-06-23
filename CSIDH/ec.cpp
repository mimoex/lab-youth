#include "fp.h"

Point ec_add(Point& p, Point& q, const mpz_class& a, const mpz_class& mod)
{
	if (p.inf) return q;
	if (q.inf) return p;

	mpz_class lambda_temp, y_1, lambda;
	mpz_class lh, rh;
	Point result;


	//cout << "Px,Qx:" << p.x << "," << q.x << endl;
	if (p.x == q.x) {
		//P + (-P) = 0
		//cout << "Py,-Py:" << p.y << "," << -q.y << endl;
		if ((p.y + q.y == mod) || (p.y == q.y && p.y == 0)) {
			p.inf = true;
			return p;
		}
		//lambda=(3*x1*x1+a)/(2*y1)
		mul_fp(p.x, p.x, mod, &lambda_temp);
		mul_fp(3, lambda_temp, mod, &lambda_temp);
		add_fp(lambda_temp, a, mod, &lambda_temp);

		add_fp(p.y, p.y, mod, &y_1);

		div_fp(lambda_temp, y_1, mod, &lambda);
	} else {
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
	if (p.inf) return q;
	if (q.inf) return p;

	mpz_class x_temp, y_1, lambda;
	mpz_class lh, rh;
	Point result, temp;


	//cout << "Px,Qx:" << p.x << "," << q.x << endl;
	if (p.x == q.x) {
		//P + (-P) = 0
		//cout << "Py,-Py:" << p.y << "," << -q.y << endl;
		if ((p.y + q.y == mod) || (p.y == q.y && p.y == 0)) {
			p.inf = true;
			return p;
		}
		//x1=x2のとき，
		//
		//x3 = b*(3*x1^2+2*a*x1+1)^2/(2*b*y1)^2-a-x1-x1
		//A=3*x1, B=(A^2+2*a*x1+1), lh=b*B^2, C=2*b*y1, D=C^2

		mpz_class A, A2, B, B1, C, C1, D, D1;
		pow_fp(p.x, 2, mod, &A);
		mul_fp(3, A, mod, &A2);
		mul_fp(2, a, mod, &B);
		mul_fp(B, p.x, mod, &B1);
		add_fp(A2, B1, mod, &C);
		add_fp(C, 1, mod, &C1);

		mul_fp(2, b, mod, &D);
		mul_fp(D, p.y, mod, &D1);

		div_fp(C1, D1, mod, &lambda);

		//x
		mpz_class l2, a1,a2,a3,a4;
		pow_fp(lambda, 2, mod, &l2);
		mul_fp(b, l2, mod, &a1);
		mul_fp(2, p.x, mod, &a2);
		sub_fp(a1, a, mod, &a3);
		sub_fp(a3, a2, mod, &result.x);

		//y
		mpz_class temp10, temp11;
		sub_fp(p.x, result.x, mod, &temp10);
		mul_fp(lambda, temp10, mod, &temp11);
		sub_fp(temp11, p.y, mod, &result.y);

	}
	else {
		//x1!=x2のとき，
		// x3 = b*(y2-y1)^2/(x2-x1)^2-a-x1-x2

		sub_fp(q.y, p.y, mod, &lh);
		sub_fp(q.x, p.x, mod, &rh);
		div_fp(lh, rh, mod, &lambda);

		mpz_class lambda2, temp1;

		pow_fp(lambda, 2, mod, &lambda2);
		mul_fp(b, lambda2, mod, &temp1);
		sub_fp(temp1, a, mod, &result.x);
		sub_fp(result.x, p.x, mod, &result.x);
		sub_fp(result.x, q.x, mod, &result.x);

		// y3 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-b*(y2-y1)^3/(x2-x1)^3-y1
		//左辺
		mpz_class temp10, temp11;
		sub_fp(p.x, result.x, mod, &temp10);
		mul_fp(lambda, temp10, mod, &temp11);
		sub_fp(temp11, p.y, mod, &result.y);
	}
	result.inf = false;
	return result;
}


Point copy_point(Point &p) {
	Point result;
	result.x= p.x;
	result.y= p.y;
	result.inf = p.inf;
	return p;
}



/*** xを選んでからyを求める(未使用) ***/
mpz_class gen_y(const Point& p, const mpz_class& mod) {
	mpz_class y, x_cube, y_sqr;
	int test;
	pow_fp(p.x, 3, mod, &x_cube);
	add_fp(x_cube, p.x, mod, &y_sqr);
	mpz_sqrt(y.get_mpz_t(), y_sqr.get_mpz_t());
	return y;
}

/*** y^2=x^3+xのポイントを生成
	ベースポイントから乱数をスカラー倍としてポイントを生成***/
Point gen_point(const mpz_class& a, const mpz_class& b, const mpz_class& mod)
{
	/*** 乱数生成 ***/
	random_device rnd;
	gmp_randclass r(gmp_randinit_default);
	r.seed(rnd());

	size_t n = 512;
	mpz_class rand_num;

	Point p, gen_point;
	p.x = "1486724546773817160185706994400590876325221178963600688050988942095980519360894027838516528686056471052894965680483000185749997971778637728426028358414319";
	p.y = "3607689279296022073046436556648074498611954478353744067528678254744301698398578854112928964964392932579386282175432815139942539613739773902651999377835184";
	p.inf = false;
	rand_num = r.get_z_bits(n);

	gen_point = l_to_r_binary_mont(p, a, b, mod, rand_num);
	return gen_point;
}

/*** y^2=x^3+xのの点であれば0を返す ***/
size_t check_point(const Point& p, const mpz_class& mod)
{
	mpz_class test, x_cube, y_sqr, rh;
	pow_fp(p.x, 3, mod, &x_cube);
	add_fp(x_cube, p.x, mod, &rh);

	mul_fp(p.y, p.y, mod, &y_sqr);
	sub_fp(rh, y_sqr, mod, &test);

	if (test == 0)	return 0;
	else			return 1;
}
