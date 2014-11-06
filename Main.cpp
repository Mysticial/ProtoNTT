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

#include "source/ModularIntrinsics.h"

#include "source/Primes3.h"
#include "source/Primes5.h"
#include "source/Primes7.h"
#include "source/Primes9.h"

using namespace ProtoNTT;

#include <iostream>
using std::cout;
using std::endl;

int main(){

    //  Test a multiply modulus.

    ModularRing p(7097673012735901697);
    TwiddleFactor Wf = p.make_twiddle(4614278974170858164); //  Forward twiddle
    TwiddleFactor Wi = p.make_twiddle(5800385079451287434); //  Inverse twiddle

    uint64_t x = 5441918181069802676;
    cout << "x            = " << x << endl;
    x = p.mulmod(x,Wf);
    cout << "x * W        = " << x << endl;
    x = p.mulmod(x,Wi);
    cout << "x * W * W^-1 = " << x << endl;




    system("pause");
}
