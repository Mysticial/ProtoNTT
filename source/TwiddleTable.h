/* TwiddleTable.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/06/2014
 * Last Modified    : 11/06/2014
 * 
 *      This is the master twiddle factor table used to store all the
 *  twiddle factors from all the transform modes supported by this algorithm.
 * 
 */

#pragma once
#ifndef _ntt_TwiddleTable_H
#define _ntt_TwiddleTable_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include "PrimeSets/Primes3.h"
#include "PrimeSets/Primes5.h"
#include "PrimeSets/Primes7.h"
#include "PrimeSets/Primes9.h"
#include "Internals/ModulusSet.h"
#include "CoreTransformParameters.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class TwiddleTable{
public:
    //  I should probably refactor this into a 2D array.
    ModulusSet p3m1, p3m3, p3m5, p3m7;
    ModulusSet p5m1, p5m3, p5m5, p5m7;
    ModulusSet p7m1, p7m3, p7m5, p7m7;
    ModulusSet p9m1, p9m3, p9m5, p9m7;

public:
    TwiddleTable()
        : p3m1(ProtoNTT::p3m1), p3m3(ProtoNTT::p3m3), p3m5(ProtoNTT::p3m5), p3m7(ProtoNTT::p3m7)
        , p5m1(ProtoNTT::p5m1), p5m3(ProtoNTT::p5m3), p5m5(ProtoNTT::p5m5), p5m7(ProtoNTT::p5m7)
        , p7m1(ProtoNTT::p7m1), p7m3(ProtoNTT::p7m3), p7m5(ProtoNTT::p7m5), p7m7(ProtoNTT::p7m7)
        , p9m1(ProtoNTT::p9m1), p9m3(ProtoNTT::p9m3), p9m5(ProtoNTT::p9m5), p9m7(ProtoNTT::p9m7)
    {
        populate_all_tables(Modulus::BASE_K);
    }
    void populate_all_tables(uint64_t cbitlen){
        //  Populates all the tables so that they can all handle a convolution
        //  length of "cbitlen" bits without on-the-fly twiddle generation.

        //  Warning: This may use a LOT of memory.

        p3m1.ensure_tables(get_table_k(cbitlen, p3m1));
        p3m3.ensure_tables(get_table_k(cbitlen, p3m3));
        p3m5.ensure_tables(get_table_k(cbitlen, p3m5));
        p3m7.ensure_tables(get_table_k(cbitlen, p3m7));

        p5m1.ensure_tables(get_table_k(cbitlen, p5m1));
        p5m3.ensure_tables(get_table_k(cbitlen, p5m3));
        p5m5.ensure_tables(get_table_k(cbitlen, p5m5));
        p5m7.ensure_tables(get_table_k(cbitlen, p5m7));

        p7m1.ensure_tables(get_table_k(cbitlen, p7m1));
        p7m3.ensure_tables(get_table_k(cbitlen, p7m3));
        p7m5.ensure_tables(get_table_k(cbitlen, p7m5));
        p7m7.ensure_tables(get_table_k(cbitlen, p7m7));

        p9m1.ensure_tables(get_table_k(cbitlen, p9m1));
        p9m3.ensure_tables(get_table_k(cbitlen, p9m3));
        p9m5.ensure_tables(get_table_k(cbitlen, p9m5));
        p9m7.ensure_tables(get_table_k(cbitlen, p9m7));
    }

    static int get_table_k(uint64_t cbitlen, const ModulusSet& set){
        //  Returns the table size k such that the table will be large enough to
        //  handle a convolution of cbitlen without resorting to expensive
        //  on-the-fly twiddle factor generation.
        int k = Modulus::BASE_K;
        while (CoreTransformParameters(set.factors, set.multiplier, k).get_cbitlen() < cbitlen){
            k++;
        }
        return k;
    }
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
