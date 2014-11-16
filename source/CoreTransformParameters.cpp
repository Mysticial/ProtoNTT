/* CoreTransformParameters.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/06/2014
 * Last Modified    : 11/15/2014
 * 
 */

#include <iostream>
#include "CoreTransformParameters.h"
#include "TestFramework.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
CoreTransformParameters::CoreTransformParameters(uint64_t cbitlen){
    //  Makes parameters suitable for performing a convolution at least
    //  as long as cbitlen.

    //  This uses a very basic heuristic to select the parameters. They are
    //  definitely not optimal.

    //  Always use 3 primes for now. Once the benchmark framework is done, use
    //  it to find the thresholds for switching between the # of primes.
    primes = 3;

    uint64_t scale = (primes - 1) * 32;

    //  Scale up to the largest power-of-two that is too small to handle the
    //  requested convolution length.
    int try_k = Modulus::BASE_K;
    while ((scale << (try_k + 1)) < cbitlen){
        try_k++;
    }

    uint64_t try_cbitlen;
    do{
        //  Use: 5 * 2^(k - 2)
        multiplier = 5;
        k = try_k - 2;
        try_cbitlen = scale * multiplier << k;
        if (k >= Modulus::BASE_K && try_cbitlen >= cbitlen)
            break;

        //  Use: 3 * 2^(k - 1)
        multiplier = 3;
        k = try_k - 1;
        try_cbitlen = scale * multiplier << k;
        if (k >= Modulus::BASE_K && try_cbitlen >= cbitlen)
            break;

        //  Use: 7 * 2^(k - 2)
        multiplier = 7;
        k = try_k - 2;
        try_cbitlen = scale * multiplier << k;
        if (k >= Modulus::BASE_K && try_cbitlen >= cbitlen)
            break;

        //  Use: 1 * 2^(k + 1)
        multiplier = 1;
        k = try_k + 1;
        try_cbitlen = scale * multiplier << k;
    }while (false);

    update_cbitlen();
}
void CoreTransformParameters::print() const{
    std::cout << "Transform Parameters:" << std::endl;
    std::cout << "    Number of Primes   = " << primes << std::endl;
    std::cout << "    Transform Length   = " << multiplier << " * 2^" << k << std::endl;
    std::cout << "    Convolution Length = "; print_bits(cbitlen);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
