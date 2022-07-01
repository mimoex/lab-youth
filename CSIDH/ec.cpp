#include "fp.h"

//Weierstrass Curve
Point ec_add(const Point& p, const Point& q, const mpz_class& a, const mpz_class& mod)
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
			result.x = 0; result.y = 0;
			result.inf = true;
			return result;
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
Point montgomery_ec_add(const Point& p, const Point& q, const mpz_class& a, const mpz_class& b, const mpz_class& mod)
{
	if (p.inf == 1) return q;
	if (q.inf == 1) return p;

	mpz_class x_temp, y_1, lambda;
	mpz_class lh, rh;
	Point result, temp;


	//cout << "Px,Qx:" << p.x << "," << q.x << endl;
	if (p.x == q.x) {
		//P + (-P) = 0
		//cout << "Py,-Py:" << p.y << "," << -q.y << endl;
		if ((p.y + q.y == mod) || (p.y == q.y && p.y == 0)) {
			result.x = 0; result.y = 0;
			result.inf = true;
			return result;
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
		mpz_class temp10, temp11,temp12, temp13,temp01, temp02, temp03, l3;
		add_fp(a2, p.x, mod, &temp12);
		add_fp(temp12, a, mod, &temp13);
		mul_fp(temp13, lambda, mod, &temp01);

		pow_fp(lambda, 3, mod, &l3);
		mul_fp(b, l3, mod, &temp02);
		
		sub_fp(temp01, temp02, mod, &temp10);
		sub_fp(temp10, p.y, mod, &result.y);

	}
	else {
		//x1!=x2のとき，
		// x3 = b*(y2-y1)^2/(x2-x1)^2-a-x1-x2

		sub_fp(q.y, p.y, mod, &lh);
		sub_fp(q.x, p.x, mod, &rh);
		div_fp(lh, rh, mod, &lambda);

		mpz_class lambda2, temp1,temp2, temp3;

		pow_fp(lambda, 2, mod, &lambda2);
		mul_fp(b, lambda2, mod, &temp1);
		add_fp(p.x, q.x, mod, &temp2);
		sub_fp(temp1, temp2, mod, &temp3);
		sub_fp(temp3, a, mod, &result.x);

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


//Tonelli-Shanks algorithm (参考：https://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm)
mpz_class TonelliShanks(const mpz_class & n, const mpz_class & mod)
{
	mpz_class s,r, q, result;
	s = 0;
	q = mod - 1;
	while ((q & 1) == 0) { q /= 2; ++s; }
	if (s == 1) {
		pow_fp(n, (mod + 1) / 4, mod, &r);
		if ((r * r) % mod == n) return r;
		return 0;
	}
	// Find the first quadratic non-residue z by brute-force search
	mpz_class z, c, t,tt, m, check2,b,b2, check3;
	z = 1;
	while ( check2 != mod - 1) pow_fp(++z, (mod - 1) / 2, mod, &check2);
	pow_fp(z, q, mod,&c);
	pow_fp(n, (q + 1) / 2, mod,&r);
	pow_fp(n, q, mod,&t);
	m = s;
	while (t != 1) {
		tt = t;
		int i = 0;
		while (tt != 1) {
			tt = (tt * tt) % mod;
			++i;
			if (i == m) return 0;
		}
		pow_fp(2, m - i - 1, mod - 1, &check3);
		pow_fp(c, check3, mod,&b);
		b2 = (b * b) % mod;
		r = (r * b) % mod;
		t = (t * b2) % mod;
		c = b2;
		m = i;
	}
	if ((r * r) % mod == n) return r;
	return 0;
}


mpz_class sqrt_mod(const mpz_class& n, const mpz_class& p)
{
	mpz_class result, check;
	pow_fp(n, (p - 1) / 2, p, &check);
	if (check == 1) {
		pow_fp(n, (p + 1) / 4, p, &result);
		return result;
	}
	else return 0;
	
}

/*** xを選んでからyを求める
	b*y^2=x^3+a*x^2+x ***/
Point gen_point_sqrt(const mpz_class& a, const mpz_class& b, const mpz_class& mod) {
	Point result;
	result.inf = false;
	size_t mask_bit = 32;
	mpz_class x, y, x_cube, x_sqr, y_temp, ax2, rh,div_rh, sqrt_result, div_x, div_a, testpow;
	result.x = random_fp(mod);
	int sqrt_check = 0;

	while (sqrt_check == 0) {
		pow_fp(result.x, 3, mod, &x_cube);

		pow_fp(result.x, 2, mod, &x_sqr);
		mul_fp(a, x_sqr, mod, &ax2);

		add_fp(x_cube, ax2, mod, &y_temp);
		add_fp(y_temp, result.x, mod, &rh);

		div_fp(rh, b, mod, &div_rh);
		
		//sqrt_result = TonelliShanks(div_rh, mod);
		sqrt_result = sqrt_mod(div_rh, mod);
		if (sqrt_result == 0) {	//not整数
			result.x++;
			//cout << "rand:" << result.x << endl;
		}
		else {
			sqrt_check = 1;
			//mpz_sqrt(result.y.get_mpz_t(), rh.get_mpz_t());
			result.y = sqrt_result;
			//cout <<"x,y :"<< result.x<<", " << result.y << endl;
		}
	}
	return result;
}

/*** y^2=x^3+xのポイントを生成
	ベースポイントから乱数をスカラー倍としてポイントを生成
	(未使用) 危険なので使わない
	2つのスカラーa,bからaPとbPを生成したときに，離散対数問題が解けてしまうおそれあり
	***/
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

/*** b*y^2=x^3+a*x^2+xの点であれば0を返す ***/
size_t check_point(const Point& p, const mpz_class a, const mpz_class b, const mpz_class& mod)
{
	mpz_class test, x_cube, x_sqr, y_sqr, rh, lh, ax2;
	pow_fp(p.x, 3, mod, &x_cube);

	pow_fp(p.x, 2, mod, &x_sqr);
	mul_fp(a, x_sqr, mod, &ax2);

	add_fp(x_cube, ax2, mod, &rh);
	add_fp(rh, p.x, mod, &rh);

	pow_fp(p.y, 2, mod, &y_sqr);
	mul_fp(b, y_sqr, mod, &lh);

	if (lh == rh)	return 0;
	else			return 1;
}