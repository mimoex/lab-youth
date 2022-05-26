#include <iostream>
#include <mpirxx.h>
#include <random>

using namespace std;

const size_t primeNum = 74;

const int prime[primeNum] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587 };


int main(void)
{
    int i;
    mpz_class qq, p;

    /*** 乱数生成 ***/
    random_device rd;
    default_random_engine eng(rd());
    uniform_int_distribution<int> distr(0, 10);  //0から10までの乱数を生成 //論文では-5から5?

    cout << "[" << ends;
    for (i = 0; i < primeNum; i++) {
        cout << prime[i] <<" " << ends;
    }
    cout<<"]\n" << endl;
    
    qq = 4;
    for (i = 0; i < primeNum; i++) {
        qq *= prime[i];
    }
    p = qq - 1;
    cout << "p=" << p << endl;

    /*** 素数判定 ***/
    cout << "isPrime? :" <<ends;
    if (mpz_probab_prime_p(p.get_mpz_t(), 50) == 1) {
        cout << "True" << endl;
    }
    else {
        cout << "False" << endl;
    }

    /*** 鍵生成 ***/

    int a[primeNum], b[primeNum]; //Aさんの鍵,Bさんの鍵

    for (i = 0; i < primeNum; i++) {
        a[i] = distr(eng);
        b[i] = distr(eng);

        
    }
    cout << "Aさんの鍵:" << ends;
    for (i = 0; i < primeNum; i++) {
        cout << a[i] <<" "<< ends;
    }
    cout << endl;

    cout << "Bさんの鍵:" << ends;
    for (i = 0; i < primeNum; i++) {
        cout << b[i] << " " << ends;
    }
    cout << endl;
        

    return 0;
}