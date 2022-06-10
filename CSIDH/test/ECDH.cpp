#include "fp.h"

//右向きバイナリ法	楕円曲線の点p, パラメータa, Modulus, スカラー倍n
Point l_to_r_binary(Point& p, const mpz_class& a, const mpz_class& mod, const mpz_class& n)
{
	string bit;
	size_t bit_size;
	Point result, temp_point;
	result.inf = true;
	temp_point = copy_point(p);

	bit = n.get_str(2);
	bit_size = bit.size();

	//cout << "n=" << n << " = 0b" << bit << endl;

	result = p;

	for (int i = bit_size - 2; i >= 0; i--) {
		result = ec_add(result, result, a, mod);	//double

		if (bit[bit_size - i - 1] == '1') {
			result = ec_add(result, temp_point, a, mod);	//add
		}
	}
	return result;
}

//右向きバイナリ法チェック
void binary_check()
{
	mpz_class a = 1, mod, n;
	Point p, result;

	mod = 4093;
	cout << "p=" << mod << "\n" << endl;

	p.x = 2116; p.y = 4041; p.inf = false;
	result.inf = true;

	/*ループで確認 5000回*/
	for (n = 1; n <= 1000; n++) {
		result = l_to_r_binary(p, a, mod, n);
		if (result.inf) {
			cout << "Inf\n" << endl;
		}
		else {
			//cout << result.x << ", " << result.y << endl; cout << endl;
			cout << n << "P=[" << result.x << ", " << result.y << "]\n" << endl;
		}
	}
}


/*
* PARI/GPで試したECDH
256 bit ECDH
num=59443744980303287241239278097475230676395746758186821689120627187548319115458
G=[59095769176739621529737602749134185345320410125049059797711125909246444074505, 12514794149476622141303143643372751044504499175864172408561274462277184358990]
step 1
a=3750262606785314222042460036037476285899366197507318144014303672648591270701
step 2
GA=[4823588636300028222970756607341551981319475616573177364558822483171465199759, 26123008269383518719219866915448003855623917031026125690655277568924419271292]
step 3
b=46777701863925597189893824885909005514091783047908907721185076849708549960107
step 4
GB=[50950403907435387503162792544571685202870181193070868225146587443963300026911, 26009154553004092222386587500868753242904729701293569298867843569969303448011]
step 5
KA=27158317184322598896534584519113203787354732534795022582990980112225341821022
step 6
KB=27158317184322598896534584519113203787354732534795022582990980112225341821022
*/
void ecdh()
{
	size_t n = 256;
	mpz_class mod;
	cout << n << "bit ECDH" << endl;

	mod.set_str("59443744980303287241239278097475230676879741160801108473218702094315841735873", 10);

	mpz_class a = 5, b = 0;
	Point p;
	p.x.set_str("59095769176739621529737602749134185345320410125049059797711125909246444074505", 10);
	p.y.set_str("12514794149476622141303143643372751044504499175864172408561274462277184358990", 10);

	//step1
	cout << "step1" << endl;
	mpz_class a_key;
	a_key.set_str("3750262606785314222042460036037476285899366197507318144014303672648591270701", 10);
	cout << "a_key: " << a_key << endl;

	//step2
	cout << "step2" << endl;
	Point GA;
	GA = l_to_r_binary(p, a, mod, a_key);
	cout << "GA:" << GA.x << ", " << GA.y << endl;

	//step3
	cout << "step3" << endl;
	Point GB;
	mpz_class b_key;
	b_key.set_str("46777701863925597189893824885909005514091783047908907721185076849708549960107", 10);
	cout << "b_key: " << b_key << endl;

	//step4
	cout << "step4" << endl;
	GB = l_to_r_binary(p, a, mod, b_key);
	cout << GB.x << ", " << GB.y << endl;

	//step5
	cout << "step5" << endl;
	Point KA;
	KA = l_to_r_binary(GB, a, mod, a_key);
	cout << KA.x << ", " << KA.y << endl;

	//step6
	cout << "step6" << endl;
	Point KB;
	KB = l_to_r_binary(GA, a, mod, b_key);
	cout << KB.x << ", " << KB.y << endl;

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
