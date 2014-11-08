/* Main.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/01/2014
 * Last Modified    : Always...
 * 
 *      Since ProtoNTT is really just a library, there really is no need for a
 *  main method. But it exists to develop and test the library.
 * 
 *  When it's finished, I may put a small UI here that will run and verify
 *  test multiplications using random data.
 * 
 */

#include "source/Modulus.h"

#include "source/PrimeSets/Primes3.h"
#include "source/PrimeSets/Primes5.h"
#include "source/PrimeSets/Primes7.h"
#include "source/PrimeSets/Primes9.h"

using namespace ProtoNTT;

#include <iostream>
using std::cout;
using std::endl;


void print(const uint64_t* A,size_t L){
    if (L == 0){
        std::cout << "{}" << std::endl;
        return;
    }

    std::cout << "{";

    L--;
    for (size_t c = 0; c < L; c++){
        std::cout << (A[c]) << ",";
    }
    std::cout << A[L] << "}" << std::endl;
}

int main(){

    //  Test a multiply modulus.
    {
        ModularRing p(7097673012735901697);
        TwiddleFactor Wf = p.make_twiddle(4614278974170858164); //  Forward twiddle
        TwiddleFactor Wi = p.make_twiddle(5800385079451287434); //  Inverse twiddle

        uint64_t x = 5441918181069802676;
        cout << "x            = " << x << endl;
        x = p.mulmod(x,Wf);
        cout << "x * W        = " << x << endl;
        x = p.mulmod(x,Wi);
        cout << "x * W * W^-1 = " << x << endl;
    }

    //  Test modulus construction.
    {
        for (int c = 0; c < 3; c++){
            Modulus p(p3m1,c);
        }
        for (int c = 0; c < 5; c++){
            Modulus p(p5m1,c);
        }
        for (int c = 0; c < 7; c++){
            Modulus p(p7m1,c);
        }
        for (int c = 0; c < 9; c++){
            Modulus p(p9m1,c);
        }
    }

    //  Test a tiny 4-point square convolution.
    {
        Modulus p(p3m1,0);

        int k = 2;
        size_t L = (size_t)1 << k;
        p.make_tables(k);

        uint64_t T[] = {5, 3, 7, 8};
        cout << "T   = "; print(T,L);

        p.forward(T);
        p.inverse_fmul(T,T);
        for (size_t c = 0; c < L; c++){
            T[c] = p.scale_down(k,T[c]);
        }

        cout << "T^2 = "; print(T,L);
    }






    system("pause");
}
