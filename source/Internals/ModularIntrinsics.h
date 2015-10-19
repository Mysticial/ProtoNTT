/* ModularIntrinsics.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/01/2014
 * Last Modified    : 11/05/2014
 * 
 *  These are all the performance-critical subroutines.
 * 
 */

#pragma once
#ifndef _ntt_ModularIntrinsics_H
#define _ntt_ModularIntrinsics_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Dependencies
#include "ArchIntrinsics.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct TwiddleFactor{
    uint64_t W;
    uint64_t Wp;    //  floor(W * 2^64 / p)
};
struct ModularRing{
protected:
    uint64_t prime;

    //  2^128 / p
    uint64_t iL;
    uint64_t iH;

public:
    ModularRing(uint64_t prime)
        : prime(prime)
    {
        invert64(iL, iH, prime);
    }

public:
    //  We need to keep the conditional trivial to have any chance that
    //  Visual Studio will be smart enough to generate a conditional move.
    FORCE_INLINE uint64_t reduce_p(uint64_t x) const{
        register uint64_t tmp = x - prime;
        return x >= prime ? tmp : x;
    }
    FORCE_INLINE uint64_t reduce_n(int64_t x) const{
        register uint64_t tmp = x + prime;
        return x < 0 ? tmp : x;
    }

public:
    FORCE_INLINE uint64_t mulmod(uint64_t x, const TwiddleFactor& W) const{
        //  Victor Shoup's multiply-modulus.

        //Conditions:
        //  -   p < 2^63
        //  -   W < p
        //  -   Wp = floor(W * 2^64 / p)

        return reduce_p(x * W.W - mulH(x, W.Wp) * prime);
    }
    FORCE_INLINE void butterfly2_forward(uint64_t& A, uint64_t& B) const{
        //  Non-twiddled butterfly with reductions at the end.
        uint64_t r0 = A + B;
        int64_t r1 = (int64_t)A - (int64_t)B;
        A = reduce_p(r0);
        B = reduce_n(r1);
    }
    FORCE_INLINE void butterfly2_inverse(uint64_t& A, uint64_t& B) const{
        //  Non-twiddled butterfly with reductions at the beginning.
        uint64_t r0 = reduce_p(A);
        uint64_t r1 = reduce_p(B);
        A = r0 + r1;
        B = r0 - r1 + prime;
    }
    FORCE_INLINE void butterfly2_forward(uint64_t& A, uint64_t& B, const TwiddleFactor& W) const{
        //  Forward butterfly with 2 reductions.
        uint64_t r0 = A + B;
        uint64_t r1 = mulmod(A - B + prime, W);
        A = reduce_p(r0);
        B = r1;
    }
    FORCE_INLINE void butterfly2_inverse(uint64_t& A, uint64_t& B, const TwiddleFactor& W) const{
        //  Inverse butterfly with 2 reductions.
        uint64_t r0 = reduce_p(A);
        uint64_t r1 = mulmod(B, W);
        A = r0 + r1;
        B = r0 - r1 + prime;
    }

    FORCE_INLINE TwiddleFactor make_twiddle(uint64_t W) const{
        TwiddleFactor out;
        out.W = W;
        out.Wp = get_Wp(W);
        return out;
    }

private:
    FORCE_INLINE uint64_t get_Wp(uint64_t W) const{
        register uint64_t q = mulH2x1(iL, iH, W);
        register uint64_t m = 0 - q*prime;
        register uint64_t tmp = q + 1;
        return m >= prime ? tmp : q;
    }
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
