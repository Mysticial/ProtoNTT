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

#include "source/Transforms.h"
#include "source/ModulusSet.h"

#include "source/PrimeSets/Primes3.h"
#include "source/PrimeSets/Primes5.h"
#include "source/PrimeSets/Primes7.h"
#include "source/PrimeSets/Primes9.h"

using namespace ProtoNTT;

#include <string.h>
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

    //  Test a 8-point square convolution.
    {
        Modulus p(p3m1,0);

        int k = 3;
        size_t L = (size_t)p.multiplier << k;
        p.make_tables(k);

        uint64_t T[] = {6, 9, 0, 3, 2, 4, 7, 7};
        cout << "T   = "; print(T,L);

        transform_forward(p,k,T);
        transform_inverse_fmul(p,k,T,T);
        for (size_t c = 0; c < L; c++){
            T[c] = p.scale_down(k,T[c]);
        }

        cout << "T^2 = "; print(T,L);
    }

    //  Test a 24-point square convolution.
    {
        Modulus p(p3m3,0);

        int k = 3;
        size_t L = (size_t)p.multiplier << k;
        p.make_tables(k);

        uint64_t T[] = {
            4, 2, 8, 6, 8, 4, 2, 5, 1, 6, 8, 5,
            1, 7, 8, 3, 8, 8, 5, 0, 1, 0, 3, 5
        };
        cout << "T   = "; print(T,L);

        transform_forward(p,k,T);
        transform_inverse_fmul(p,k,T,T);
        for (size_t c = 0; c < L; c++){
            T[c] = p.scale_down(k,T[c]);
        }

        cout << "T^2 = "; print(T,L);
    }

    //  Test a tiny 20-point square convolution.
    {
        Modulus p(p3m5,0);

        int k = 3;
        size_t L = (size_t)p.multiplier << k;
        p.make_tables(k);

        uint64_t T[] = {
            9, 4, 5, 0, 5, 7, 4, 7, 8, 5,
            8, 1, 3, 2, 4, 0, 6, 4, 3, 9,
            4, 7, 5, 0, 0, 2, 9, 0, 4, 0,
            3, 7, 7, 0, 7, 8, 0, 6, 0, 8
        };
        cout << "T   = "; print(T,L);

        transform_forward(p,k,T);
        transform_inverse_fmul(p,k,T,T);
        for (size_t c = 0; c < L; c++){
            T[c] = p.scale_down(k,T[c]);
        }

        cout << "T^2 = "; print(T,L);
    }

    //  Test a tiny 28-point square convolution.
    {
        Modulus p(p3m7,0);

        int k = 3;
        size_t L = (size_t)p.multiplier << k;
        p.make_tables(k);

        uint64_t T[] = {
            0, 8, 7, 4, 4, 8, 2, 3, 4, 5, 5, 6, 4, 3,
            9, 9, 6, 2, 4, 3, 9, 3, 8, 7, 3, 6, 7, 3,
            4, 8, 3, 2, 2, 9, 5, 0, 0, 1, 6, 9, 9, 5,
            5, 2, 2, 9, 2, 0, 4, 3, 6, 5, 5, 7, 9, 5
        };
        cout << "T   = "; print(T,L);

        transform_forward(p,k,T);
        transform_inverse_fmul(p,k,T,T);
        for (size_t c = 0; c < L; c++){
            T[c] = p.scale_down(k,T[c]);
        }

        cout << "T^2 = "; print(T,L);
    }

    //  Test a CRT.
    {
        ModulusSet set(p9m1);

        uint64_t R[] = {
            1840612861089888724ULL, 211220548242431273ULL,
            10742799559211267726ULL, 10056489111922224330ULL
        };
        uint64_t T[9];
        uint64_t carry[5];
        memset(carry,0,sizeof(carry));

        cout << "start = "; print(R,4);

        set.start_block(R,T,1);
        memset(R,0,sizeof(R));
        set.finish_block(0,carry,R,T,1);

        cout << "end   = "; print(R,4);
    }






    return system("pause");
}
