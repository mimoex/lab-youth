#include "fp.h"

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


void ec_add_check()
{
	mpz_class a = 7, mod= 97;
	Point p, q, result;
	
	p.x = 24; p.y = 11;
	q.x = 96; q.y = 34;
	p.inf = false; q.inf = false; result.inf = false;
	result=ec_add(p, q, a, mod);

	cout << "[" << p.x <<","<< p.y<<"]+[" << q.x <<","<<q.y<<"]=["<< result.x<<","<<result.y<<"]" << endl;
}


Point binary_m(Point* p, const mpz_class* a, const mpz_class* mod, const mpz_class* n)
{
	string bit;
	size_t bit_size;
	Point result, temp_point;
	result.inf = true;
	temp_point = copy_point(*p);

	bit = n->get_str(2);
	bit_size = bit.size();

	cout << "n=" << *n <<" = 0b"<< bit << endl;

	result = *p;

	for (int i = bit_size-2; i >=0; i--) {
		result = ec_add(result, result, *a, *mod);	//doubling

		cout << "i: " << bit[i] << endl;
		if (bit[bit_size - i - 1] == '1') {
			result = ec_add(result, temp_point, *a, *mod);	//addition
		}

				
	}
	return result;
	
}

void binary_check()
{
	mpz_class a = 1, mod;
	mpz_class n = 3;
	Point p, result;

	//mod=gen_prime();
	mod = 97;
	cout << "p="<<mod << endl;

	p.x = 31; p.y = 72;
	p.inf = false; result.inf = false;
	for (n = 0; n < 10; n++) {
		result = binary_m(&p, &a, &mod, &n);
		cout << "nP=" << n << "[" << p.x << ", " << p.y << "] = [" << result.x << ", " << result.y << "]" << endl;

	}

}





int main()
{

	binary_check();

	//gen_csidh_key();

}
