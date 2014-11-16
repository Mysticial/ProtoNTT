/* CoreTransformParameters.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/06/2014
 * Last Modified    : 11/15/2014
 * 
 *      This class represents the mathematical set of parameters that describe
 *  the transform for this multiplication algorithm.
 * 
 *      This class provides no information about the run-time parameters that
 *  are to be used to perform such a transform.
 * 
 */

#pragma once
#ifndef _ntt_CoreTransformParameters_H
#define _ntt_CoreTransformParameters_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include <stdint.h>
#include "Internals/Modulus.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class CoreTransformParameters{
    static const uint64_t MAX_CBITLEN;

protected:
    int primes;
    int multiplier;
    int k;

    uint64_t cbitlen;

public:
    CoreTransformParameters(uint64_t cbitlen);
    CoreTransformParameters(int primes,int multiplier,int k)
        : primes(primes)
        , multiplier(multiplier)
        , k(k)
    {
        //  Make parameters using the specified parameters.
        update_cbitlen();
    }

    void print() const;

    uint64_t get_mbitlen() const{
        //  Returns the approximate memory bit-length for this transform size.
        //  This is mainly used as a heuristic to generate lookup tables.
        uint64_t length = (uint64_t)multiplier << k;
        return primes * 64 * length;
    }
    uint64_t get_cbitlen() const{
        return cbitlen;
    }

private:
    void update_cbitlen(){
        uint64_t length = (uint64_t)multiplier << k;
        cbitlen = (primes - 1) * 32 * length;
    }
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
