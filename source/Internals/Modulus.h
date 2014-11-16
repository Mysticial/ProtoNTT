/* Modulus.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/01/2014
 * Last Modified    : 11/09/2014
 * 
 *      This class represents an NTT modulus. It encapsulates all metadata
 *  about the modulus including the precomputed twiddle factors.
 */

#pragma once
#ifndef _ntt_Modulus_H
#define _ntt_Modulus_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include "../PrimeSets/PrimeSet.h"
#include "ModularIntrinsics.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class Modulus : public ModularRing{
public:
    static const int BASE_K = 2;

public:
    //  These fields are public for easy external access. But they should not be
    //  modified outside of the class.
    int multiplier;
    int max_k;

    //  Transform length = multiplier * 2^k
    //  Maximum length = multiplier * 2^max_k

    TwiddleFactor one;  //  1
    TwiddleFactor word; //  2^64 mod p

    //  Size of the current twiddle factor table.
    int table_k = 0;

    //  Precomputed Twiddle Factor Tables
    //  Let w(k) be the root of unity: w^(multiplier*2^k) mod p = 1
    //  forward_table[k][i] = w(k)^i
    //  inverse_table[k][i] = w(k)^-i
    std::vector<std::vector<TwiddleFactor>> forward_table;
    std::vector<std::vector<TwiddleFactor>> inverse_table;

    //  scaling_factors[i] = (m * 2^i)^-1 mod p
    std::vector<TwiddleFactor> scaling_factors;

    //  Table Stubs: First entry in each table.
    //  stub_f[i]^( multiplier*2^k) mod p = 1
    //  stub_i[i]^(-multiplier*2^k) mod p = 1
    std::vector<TwiddleFactor> stub_f;
    std::vector<TwiddleFactor> stub_i;

private:
    //  End-Point Transforms: size = multiplier * 2^BASE_K
    void (*fp_forward)(const Modulus& p,uint64_t* T);
    void (*fp_inverse_fmul)(const Modulus& p,uint64_t* T,const uint64_t* A);

public:
    Modulus(const PrimeSet& set,int index);
    void make_tables(int k);

public:
    uint64_t power(const TwiddleFactor& x,uint64_t pow) const;
    uint64_t generate(const uint64_t* T,size_t L) const;

public:
    FORCE_INLINE void forward(uint64_t* T) const{
        //  Perform forward transform of length: m * 2^BASE_K
        fp_forward(*this,T);
    }
    FORCE_INLINE void inverse_fmul(uint64_t* T,const uint64_t* A) const{
        //  Perform pointwise products with A.
        //  Perform inverse transform of length: m * 2^BASE_K
        fp_inverse_fmul(*this,T,A);
    }
    FORCE_INLINE uint64_t scale_down(int k,uint64_t x) const{
        //  Compute: x * (m * 2^k)^-1 mod p
        return mulmod(x,scaling_factors[k]);
    }

private:
    //  Internal Helpers
    void make_scaling_factors(uint64_t scaler,int max_k);
    void make_table_stubs(uint64_t primitive_root);
    void sanity_check() const;
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
