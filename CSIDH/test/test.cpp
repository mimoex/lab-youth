#include "fp.h"

//乗算チェック
void mul_check(void)
{
	mpz_class p = 7, a = 5, c = 8, check_mul;

	for (int i = 0; i < p; i++) {
		mul_fp(i, c, p, &check_mul);
		cout << i << "*" << c << "=" << check_mul << endl;
	}
}

//加算チェック
void ec_add_check()
{
	mpz_class a = 1, mod = 97;
	Point p, q, result;

	p.x = 31; p.y = 25;
	q.x = 96; q.y = 17;
	p.inf = false; q.inf = false; result.inf = false;
	result = ec_add(p, q, a, mod);

	cout << "[" << p.x << "," << p.y << "]+[" << q.x << "," << q.y << "]=[" << result.x << "," << result.y << "]" << endl;
}

//2倍算チェック
void ec_double_check()
{
	mpz_class a = 1, mod = 97;
	Point p, q, result;

	p.x = 66; p.y = 32;
	q.x = 66; q.y = 32;
	p.inf = false; q.inf = false; result.inf = false;
	result = ec_add(p, q, a, mod);

	cout << "[" << p.x << "," << p.y << "]+[" << q.x << "," << q.y << "]=[" << result.x << "," << result.y << "]" << endl;
}

//除算チェック
void div_check(void)
{
	mpz_class p = 7, a = 5, c, check_mul;

	for (int i = 0; i < p; i++) {
		div_fp(a, i, p, &c);
		mul_fp(i, c, p, &check_mul);
		cout << i << "*" << c << "=" << check_mul << endl;
	}
}



void CSIDH1()
{
	mpz_class a = 1, mod, n;
	Point p, result;

	//素数
	mod = gen_prime();
	cout << "p=" << mod << endl;

	//鍵
	int i;

	/*** 乱数生成 ***/
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr(0, 10);  //0から10までの乱数を生成 //論文では-5から5?

	size_t a_key[N], b_key[N]; //Aさんの鍵,Bさんの鍵

	for (i = 0; i < N; i++) {
		a_key[i] = distr(eng);
		b_key[i] = distr(eng);


	}
	cout << "Aさんの鍵:" << ends;
	for (i = 0; i < N; i++) {
		cout << a_key[i] << " " << ends;
	}
	cout << endl;

	cout << "Bさんの鍵:" << ends;
	for (i = 0; i < N; i++) {
		cout << b_key[i] << " " << ends;
	}
	cout << endl;


	//曲線定義
	p.x = 1; p.y = 0; p.inf = false;

	Point ea;
	ea = copy_point(p);


}


//テスト main関数
/*
int main(void)
{
	cout << "DH鍵共有" << endl;
	dh_key_exchange();

	cout << "\n\n除算チェック" << endl;
	div_check();

	cout << "\n\nEC add" << endl;
	CSIDH_1();

	gen_csidh_key();

	//mul_check();


	return 0;
}
*/