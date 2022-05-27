#include "fp.h"

//DH鍵共有
int dh_key_exchange(void)
{
	mpz_class p = 65537, g = 3;
	mpz_class a_sec, b_sec, a_pub, b_pub, s1,s2;

	gmp_randclass r(gmp_randinit_default);

	size_t n = 64;

	a_sec = r.get_z_bits(128);	//Aさんの秘密鍵生成(nビットの乱数)
	b_sec = r.get_z_bits(128);	//Bさんの秘密鍵生成(nビットの乱数)

	cout << "Aさんの秘密値:" << a_sec << endl;
	cout << "Bさんの秘密値:" << b_sec << endl;

	pow_fp(g, a_sec, &a_pub, p);	//Aさんの共有値
	pow_fp(g, b_sec, &b_pub, p);	//Bさんの共有値

	cout << "Aさんの共有値:" << a_pub << endl;
	cout << "Bさんの共有値:" << b_pub << endl;

	pow_fp(b_pub, a_sec, &s1, p);	//Aさんの処理
	pow_fp(a_pub, b_sec, &s2, p);	//Bさんの処理

	cout << "Aが求めた鍵:" << s1 << endl;
	cout << "Bが求めた鍵:" << s2 << endl;

	if (s1 - s2 == 0) {
		cout << "OK!" << endl;
		return 0;
	}
	else {
		cout << "NG" << endl;
		return -1;
	}
}

//除算チェック
void div_check(void)
{
	mpz_class p = 7, a = 5, c, check_mul;

	for (int i = 0; i < p; i++) {
		div_fp(a, i, &c, p);
		mul_fp(i, c, &check_mul, p);
		cout << i << "*" << c << "=" << check_mul << endl;
	}
}


int main(void)
{
	cout << "DH鍵共有" << endl;
	dh_key_exchange();

	cout << "\n\n除算チェック" << endl;
	div_check();


	return 0;
}




