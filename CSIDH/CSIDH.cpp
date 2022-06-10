#include "fp.h"

//鍵の生成
void gen_csidh_key()
{
	int i;

	/*** 乱数生成 ***/
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr(0, 10);  //0から10までの乱数を生成 //論文では-5から5?

	size_t a[N], b[N]; //Aさんの鍵,Bさんの鍵

	for (i = 0; i < N; i++) {
		a[i] = distr(eng);
		b[i] = distr(eng);


	}
	cout << "Aさんの鍵:" << ends;
	for (i = 0; i < N; i++) {
		cout << a[i] << " " << ends;
	}
	cout << endl;

	cout << "Bさんの鍵:" << ends;
	for (i = 0; i < N; i++) {
		cout << b[i] << " " << ends;
	}
	cout << endl;
}

//CSIDHに使う素数の生成
mpz_class gen_prime()
{
	int i;
	mpz_class p, qq;
	qq = 4;
	for (i = 0; i < N; i++) {
		qq *= prime[i];
	}
	p = qq - 1;
	cout << p << endl;
	string bit;
	bit = p.get_str(2);
	cout << bit.size() << endl;
	return p;
}