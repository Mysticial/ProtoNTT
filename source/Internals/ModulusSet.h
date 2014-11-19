/* ModulusSet.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/03/2014
 * Last Modified    : 11/09/2014
 * 
 *      This class represents a set of moduli that are used for the NTT algorithm.
 * 
 */

#pragma once
#ifndef _ntt_ModulusSet_H
#define _ntt_ModulusSet_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include "../PrimeSets/PrimeSet.h"
#include "Modulus.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class ModulusSet{
public:
    static const int MAX_PRIMES = 9;

    const int multiplier;
    const int factors;
    const int max_k;
    const int operand_words;
    const int carry_words;

private:
    std::vector<Modulus> p;
    const PrimeSet* set;

    int table_k = 0;

public:
    ModulusSet(const PrimeSet& set);
    void ensure_tables(int k);
    uint64_t table_bytes(int k) const;
    int get_table_k() const{
        return table_k;
    }

    void forward(int k,uint64_t* T) const;
    void inverse_fmul(int k,uint64_t* T,const uint64_t* A) const;

    //  Generate moduli from a block of raw data and write it to the transform.
    void start_block(const uint64_t* R,uint64_t* T,size_t stride) const;

    //  Scale down, construct CRT, perform carryout, write to destination.
    void finish_block(int k,uint64_t* carry,uint64_t* R,const uint64_t* T,size_t stride) const;

    //  Return the modulus object at the specified index.
    const Modulus& operator[](int index) const{
        return p[index];
    }

private:
    void CRT_build(int k,uint64_t P[MAX_PRIMES + 2],const uint64_t* W,size_t stride) const;
    void CRT_reduce1(uint64_t P[MAX_PRIMES + 2]) const;
    void CRT_reduce2(uint64_t P[MAX_PRIMES + 1]) const;
    void carryout(uint64_t* carry,uint64_t* T,uint64_t* P) const;

    void operator=(const ModulusSet&) = delete;
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
