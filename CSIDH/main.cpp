#include "fp.h"

void mon_add_check()
{
	mpz_class a = 84, b= 1, mod = 251;
	Point p, q, result;

	//2倍
	//p.x = 173; p.y = 28;
	//q.x = 173; q.y = 28;

	//加算
	p.x = 173; p.y = 28;
	q.x = 22; q.y = 154;
	result=montgomery_ec_add(p,q,a,b,mod);
	cout << "[" << p.x << "," << p.y << "]+[" << q.x << "," << q.y << "]=[" << result.x << "," << result.y << "]" << endl;
}

void curve25519_add_check()
{
	mpz_class a = 486662, b = 1, mod;
	mod.set_str("7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed", 16);
	Point p, q, result;
	p.x.set_str("9", 16); p.y.set_str("20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9",16);
	result = p;
	for (int i = 1; i < 101; i++) {
		cout << i << "P=[" << result.x.get_str(16) << "," << result.y.get_str(16) << "]" << endl;
		result = montgomery_ec_add(p, result, a, b, mod);

	}
}


Point l_to_r_binary_mont(const Point& p, const mpz_class& a, const mpz_class& b, const mpz_class& mod, const mpz_class& n)
{
	string bit;
	size_t bit_size;
	Point result, temp_point;
	result.inf = true;
	temp_point = p;

	bit = n.get_str(2);
	bit_size = bit.size();

	//cout << "n=" << n << " = 0b" << bit << endl;

	result = p;

	for (int i = bit_size - 2; i >= 0; i--) {
		result = montgomery_ec_add(result, result, a, b, mod);	//double

		if (bit[bit_size - i - 1] == '1') {
			result = montgomery_ec_add(result, temp_point, a, b, mod);	//add
		}
	}
	return result;
}

void curve25519_binary_check()
{
	mpz_class a = 486662, b = 1, mod;
	mod.set_str("7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed", 16);
	Point p, q, result;
	p.x.set_str("9", 16); p.y.set_str("20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9", 16);
	result = p;
	for (int n = 1; n < 101; n++) {
		
		result = l_to_r_binary_mont(p, a, b, mod, n);
		if (result.inf) {
			cout << "Inf\n" << endl;
		}
		else {
			cout << n << "P=[" << result.x.get_str(16) << "," << result.y.get_str(16) << "]" << endl;
		}

	}
}

void x25519()
{
	size_t n = 256;
	random_device rnd;
	gmp_randclass r(gmp_randinit_default);
	r.seed(rnd());

	mpz_class mod, a = 486662, b = 1;

	mod.set_str("7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed", 16);
	cout << "p=" << mod << endl;

	Point p;
	p.x.set_str("9", 16);
	p.y.set_str("20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9", 16);

	//step1
	cout << "step1" << endl;
	mpz_class a_key;
	a_key=r.get_z_bits(n);
	cout << "a_key: " << a_key.get_str(16) << endl;

	//step2
	cout << "step2" << endl;
	Point GA;
	GA = l_to_r_binary_mont(p, a, b, mod, a_key);
	cout << "GA:" << GA.x.get_str(16) << ", " << GA.y << endl;

	//step3
	cout << "step3" << endl;
	Point GB;
	mpz_class b_key;
	b_key = r.get_z_bits(n);
	cout << "b_key: " << b_key.get_str(16) << endl;

	//step4
	cout << "step4" << endl;
	GB = l_to_r_binary_mont(p, a, b, mod, b_key);
	cout << GB.x.get_str(16) << ", " << GB.y << endl;

	//step5
	cout << "step5" << endl;
	Point KA;
	KA = l_to_r_binary_mont(GB, a, b, mod, a_key);
	cout << KA.x.get_str(16) << ", " << KA.y << endl;

	//step6
	cout << "step6" << endl;
	Point KB;
	KB = l_to_r_binary_mont(GA, a, b, mod, b_key);
	cout << KB.x.get_str(16) << ", " << KB.y << endl;

	//check
	Point check_share;
	check_share.x = KA.x - KB.x;
	check_share.y = KA.y - KB.y;

	if (check_share.x == 0 && check_share.y == 0) {
		cout << "ECDH OK!!!" << endl;
	}
	else {
		cout << "error" << endl;
	}

}


int main()
{
	x25519();

}
