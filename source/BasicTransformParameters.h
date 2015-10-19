/* BasicTransformParameters.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/06/2014
 * Last Modified    : 11/15/2014
 * 
 *      This class represents the run-time parameters for performing a specific
 *  transform for this algorithm.
 * 
 */

#pragma once
#ifndef _ntt_BasicTransformParameters_H
#define _ntt_BasicTransformParameters_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include "CoreTransformParameters.h"
#include "Internals/ModulusSet.h"
#include "TwiddleTable.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class BasicTransformParameters : public CoreTransformParameters{
public:
    //  These are public for easy access. But they should not be modified
    //  externally.

    uint64_t Tsize; //  Transform size in bytes.
    uint64_t Psize; //  Scratch memory in bytes.

private:
    ModulusSet* set;

public:
    BasicTransformParameters(uint64_t cwordlen, TwiddleTable& tables)
        : CoreTransformParameters(cwordlen * 64)
    {
        set_modset(tables);
        check_k();
        update_sizes();
    }
    BasicTransformParameters(int primes, int multiplier, int k, TwiddleTable& tables)
        : CoreTransformParameters(primes, multiplier, k)
    {
        set_modset(tables);
        check_k();
        update_sizes();
    }

    static void time_benchmark(int primes, int multiplier, int k, double seconds = 4.0, int table_reduction = 0);

    void bench_multiply(int table_reduction = 0);
    void bench_multiply(size_t AL, size_t BL, int table_reduction = 0);
    void test() const;

    void print() const;
    void ensure_tables(uint64_t cwordlen = 0){
        //  Ensures that the tables are large enough so that a convolution
        //  length of "cbitlen" does not require on-the-fly twiddle generation.
        uint64_t cbitlen = cwordlen * 64;

        if (cbitlen == 0)
            cbitlen = CoreTransformParameters::cbitlen;

        int k = TwiddleTable::get_table_k(cbitlen, *set);
        set->ensure_tables(k);
    }
    uint64_t table_bytes(uint64_t cwordlen = 0) const{
        //  Returns the # of bytes that calling "ensure_tables()" will require.
        uint64_t cbitlen = cwordlen * 64;

        if (cbitlen == 0)
            cbitlen = CoreTransformParameters::cbitlen;

        int k = TwiddleTable::get_table_k(cbitlen, *set);
        return set->table_bytes(k);
    }

    void sqr(uint64_t* C, const uint64_t* A, size_t AL) const;
    void mul(uint64_t* C, const uint64_t* A, size_t AL, const uint64_t* B, size_t BL) const;

private:
    void set_modset(TwiddleTable& tables);
    void update_sizes(){
        Tsize = (uint64_t)primes * sizeof(uint64_t) * multiplier << k;
        Psize = 0;
    }
    void check_k() const{
        if (k < Modulus::BASE_K)
            throw "Transform length is too small.";
        if (k > set->max_k)
            throw "Maximum transform length exceeded for this prime/multiplier combination.";
    }

    void raw_to_NTT(uint64_t* T, const uint64_t* R, size_t RL) const;
    void NTT_to_raw(const uint64_t* T, uint64_t* R, size_t RL) const;
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
