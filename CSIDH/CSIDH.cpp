#include "fp.h"

//Veluの公式
size_t Velu(Point& P, mpz_class& a, mpz_class& b, const mpz_class& order, const mpz_class& mod)
{
	//if (check_point(P, mod) == 1) { cout << "違う" << endl; return 1; }

	Point Q, result;
	mpz_class tQ, tQtemp, uQ, uQtemp, wQ, wQtemp;
	mpz_class t = 0, w = 0;
	int i;

	Q.x = 0; Q.y = 0; Q.inf = true;

	for (i = 1; i < order; i++) {
		Q = montgomery_ec_add(P, Q, a, b, mod);

		//tQ = 3 * Q.x * Q.x + a;
		pow_fp(Q.x, 2, mod, &tQ);
		mul_fp(3, tQ, mod, &tQtemp);
		add_fp(tQtemp, b, mod, &tQ);
		
		//uQ = 2 * Q.y * Q.y;
		pow_fp(Q.y, 2, mod, &tQtemp);
		mul_fp(2, tQtemp, mod, &uQ);

		//wQ = uQ + tQ * Q.x;
		mul_fp(tQ, Q.x, mod, &wQtemp);
		add_fp(uQ, wQtemp, mod, &wQ);

		//t += tQ;
		add_fp(t, tQ, mod, &t);
		//w += wQ;
		add_fp(w, wQ, mod, &w);
	}
	//t *= 5;
	mul_fp(t, 5, mod, &t);
	//w *= 7;
	mul_fp(w, 7, mod, &w);
	//a -= t;
	sub_fp(b, t, mod, &b);
	//b -= w;
	sub_fp(a, w, mod, &a);
	return 0;
}

//鍵の生成
void gen_csidh_key()
{
	int i, j;

	/*** 乱数生成 ***/
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr(0, 2);  //0から10までの乱数を生成 //論文では-5から5?

	int a[N], b[N]; //Aさんの鍵,Bさんの鍵

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


	mpz_class mod;
	mod = "5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659";

	const mpz_class curve_a = 0, curve_b = 1;
	mpz_class bai, add_temp;;

	size_t check;

	/*** Aさんのstep1 ***/
	mpz_class A_curve_a = curve_a, A_curve_b = curve_b;

	Point PA;

	for (i = 0; i < N; i++) {
		for (j = 0; j < a[i];j++) {
			check=0;
			while (check == 0) {
				PA = gen_point(A_curve_a, A_curve_b, mod);
				if (check_point(PA, A_curve_a, A_curve_b, mod) == 0) cout<<"ellisoncurve-OK"<<endl;
				//cout << "inf?:" << PA.inf << ends;
				add_fp(mod, 1, mod, &add_temp);
				div_fp(add_temp, primes[i], mod, &bai);
				//cout << bai << endl;
				PA=l_to_r_binary_mont(PA, A_curve_a, A_curve_b, mod, bai);
				if (check_point(PA, A_curve_a, A_curve_b, mod) == 0) cout << "ellisoncurve-OK2" << endl;
				//cout << "inf(af scalar)?:" << PA.inf << endl;
				//cout << "PA:" << PA.x << ", " << PA.y << endl;
				if (PA.inf==false) {
					//cout << "check OK!" << endl;
					check = 1;
				}
			}
			//Pの位数はl[i]
			cout << " " << primes[i] << ends;
			//cout << "a, b;" << A_curve_a << ", " << A_curve_b << endl;
			Velu(PA, A_curve_a, A_curve_b, primes[i], mod);
			//cout << "AA:" << A_curve_a << ", " << A_curve_b << endl;
		}
		cout << ", " << ends;
	}
	cout << "\n\nAさんの公開情報:" << A_curve_a << ", " << A_curve_b << endl;
	/*** Aさんのstep1終了 ***/


	/*** Bさんのstep1 ***/
	mpz_class B_curve_a = curve_a, B_curve_b = curve_b;
	

	Point PB;

	for (i = 0; i < N; i++) {
		for (j = 0; j < b[i]; j++) {
			check = 0;
			while (check == 0) {
				PB = gen_point(B_curve_a, B_curve_b, mod);
				//if (check_point(P, mod) == 0) cout<<"OK"<<endl;
				//cout << "inf?:" << PB.inf << ends;
				add_fp(mod, 1, mod, &add_temp);
				div_fp(add_temp, primes[i], mod, &bai);
				//cout << bai << endl;
				PB = l_to_r_binary_mont(PB, B_curve_a, B_curve_b, mod, bai);
				//cout << "inf(af scalar)?:" << PB.inf << endl;
				//cout << "P:" << P.x << ", " << P.y << endl;
				if (PB.inf == 0) {
					//cout << "check OK!" << endl;
					check = 1;
				}
			}
			//Pの位数はl[i]
			cout << " " << primes[i] << ends;
			//cout << "P:" << P.x << ", " << P.y << endl;
			Velu(PB, B_curve_a, B_curve_b, primes[i], mod);
		}
		cout << ", " << ends;
	}
	cout << "\n\nBさんの公開情報:" << B_curve_a << ", " << B_curve_b << endl;
	/*** Bさんのstep1終了 ***/

	/*** Aさんのstep2 ***/
	//Bさんの公開情報を使う
	mpz_class BA_curve_a = B_curve_a, BA_curve_b = B_curve_b;
	Point PBA;
	for (i = 0; i < N; i++) {
		for (j = 0; j < a[i]; j++) {
			check = 0;
			while (check == 0) {
				PBA = gen_point(BA_curve_a, BA_curve_b, mod);
				add_fp(mod, 1, mod, &add_temp);
				div_fp(add_temp, primes[i], mod, &bai);
				PBA = l_to_r_binary_mont(PBA, BA_curve_a, BA_curve_b, mod, bai);
				//cout << "PBA:" << PBA.x << ", " << PBA.y << endl;
				if (PBA.inf == 0) {
					//cout << "check OK!" << endl;
					check = 1;
				}
			}
			//Pの位数はl[i]
			cout << " " << primes[i] << ends;
			Velu(PBA, BA_curve_a, BA_curve_b, primes[i], mod);
			//cout <<  BA_curve_a << ", " << BA_curve_b << endl;
		}
		cout << ", " << ends;
	}
	cout << "\n\nAさんの最終データ:" << BA_curve_a << ", " << BA_curve_b << endl;
	/*** Aさんのstep2終了 ***/

	/*** Bさんのstep2 ***/
	//Aさんの公開情報を使う
	mpz_class AB_curve_a = A_curve_a, AB_curve_b = A_curve_b;
	Point PAB;
	for (i = 0; i < N; i++) {
		for (j = 0; j < b[i]; j++) {
			check = 0;
			while (check == 0) {
				PAB = gen_point(AB_curve_a, AB_curve_b, mod);
				add_fp(mod, 1, mod, &add_temp);
				div_fp(add_temp, primes[i], mod, &bai);
				PAB = l_to_r_binary_mont(PAB, AB_curve_a, AB_curve_b, mod, bai);
				if (PAB.inf == 0) {
					//cout << "check OK!" << endl;
					check = 1;
				}
			}
			//Pの位数はl[i]
			cout << " " << primes[i] << ends;
			Velu(PAB, AB_curve_a, AB_curve_b, primes[i], mod);
		}
		cout << ", " << ends;
	}
	cout << "\n\nBさんの最終データ:" << AB_curve_a << ", " << AB_curve_b << endl;
	/*** Aさんのstep2終了 ***/

	/*** A,Bの最終データの比較 ***/
	if (BA_curve_a == AB_curve_a && BA_curve_b == AB_curve_b) {
		cout << "CSIDH OK!!!" << endl;
	}
	else {
		cout << "CSIDH error" << endl;
	}
}

//CSIDHに使う素数の生成
mpz_class gen_primes()
{
	int i;
	mpz_class p, qq;
	qq = 4;
	for (i = 0; i < N; i++) {
		qq *= primes[i];
	}
	p = qq - 1;
	cout << p << endl;
	string bit;
	bit = p.get_str(2);
	cout << bit.size() << endl;
	return p;
}


