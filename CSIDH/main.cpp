#include "fp.h"

//DH�����L
int dh_key_exchange(void)
{
	mpz_class p = 65537, g = 3;
	mpz_class a_sec, b_sec, a_pub, b_pub, s1,s2;

	gmp_randclass r(gmp_randinit_default);

	size_t n = 64;

	a_sec = r.get_z_bits(128);	//A����̔閧������(n�r�b�g�̗���)
	b_sec = r.get_z_bits(128);	//B����̔閧������(n�r�b�g�̗���)

	cout << "A����̔閧�l:" << a_sec << endl;
	cout << "B����̔閧�l:" << b_sec << endl;

	pow_fp(g, a_sec, &a_pub, p);	//A����̋��L�l
	pow_fp(g, b_sec, &b_pub, p);	//B����̋��L�l

	cout << "A����̋��L�l:" << a_pub << endl;
	cout << "B����̋��L�l:" << b_pub << endl;

	pow_fp(b_pub, a_sec, &s1, p);	//A����̏���
	pow_fp(a_pub, b_sec, &s2, p);	//B����̏���

	cout << "A�����߂���:" << s1 << endl;
	cout << "B�����߂���:" << s2 << endl;

	if (s1 - s2 == 0) {
		cout << "OK!" << endl;
		return 0;
	}
	else {
		cout << "NG" << endl;
		return -1;
	}
}

//���Z�`�F�b�N
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
	cout << "DH�����L" << endl;
	dh_key_exchange();

	cout << "\n\n���Z�`�F�b�N" << endl;
	div_check();


	return 0;
}




