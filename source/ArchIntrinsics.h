/* ArchIntrinsics.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/01/2014
 * Last Modified    : 11/05/2014
 * 
 *      This file contains all the bignum primitives.
 * 
 *  These are all extremely performance-critical and their implementations tend
 *  to be compiler-specific.
 */

#pragma once
#ifndef _ntt_ArchIntrinsics_H
#define _ntt_ArchIntrinsics_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include <stddef.h>
#include <stdint.h>

#ifdef _WIN32
//  Windows Includes
#include <Windows.h>
#include <intrin.h>
#else
//  Linux/GCC Includes
#endif

#if !(defined _WIN64) && !(defined __x86_64__)
#error "x64 is required."
#endif

namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
//  Visual Studio
#ifdef _WIN32
#define FORCE_INLINE inline __forceinline
FORCE_INLINE uint64_t mulH(uint64_t a,uint64_t b){
    uint64_t H;
    _umul128(a,b,&H);
    return H;
}
FORCE_INLINE void mulF(uint64_t& L,uint64_t& H,uint64_t a,uint64_t b){
    L = _umul128(a,b,&H);
}
FORCE_INLINE uint64_t mulH2x1(uint64_t aL,uint64_t aH,uint64_t b){
    //  (aL,aH) * b
    uint64_t H;
    _umul128(aL,b,&H);
    H += aH * b;
    return H;
}
inline void invert64(uint64_t& L,uint64_t& H,uint64_t p){
    //  (L,H) = 2^128 / p
    //  This function is not performance critical.
    double inv = 1. / p;
    double r = 79228162514264337593543950336. * inv;
    r *= 0.9999999999999964472863211994990706443786621093750;   //  1 - 1/2^48
    uint64_t qH = (uint64_t)r;

    uint64_t rL,rH;
    mulF(rL,rH,qH,p);

    L = 0 - (rL << 32);
    H = 0 - ((rL >> 32) | (rH << 32)) - 1;

    r = (double)L + (double)H  * 18446744073709551616.;
    r *= inv;

    //  This round-to-nearest will guarantee that we never under-estimate the
    //  quotient. But it may cause us to over-shoot the quotient by at most 1.
    uint64_t qL = (uint64_t)(r + 0.5);
    mulF(rL,rH,qL,p);

    char carry = _subborrow_u64(0,L,rL,&L);
    carry = _subborrow_u64(carry,H,rH,&H);

    //  Make correction
    if (carry)
        qL--;

    qH += qL >> 32;
    qL &= 0xffffffff;
    qL |= qH << 32;
    qH >>= 32;
    L = qL;
    H = qH;
}
////////////////////////////////////////////////////////////////////////////////
//  GCC
#else
#define FORCE_INLINE inline __attribute__ ((always_inline))
FORCE_INLINE uint64_t mulH(uint64_t a,uint64_t b){
    return (uint64_t)((unsigned __int128)a * b >> 64);
}
FORCE_INLINE void mulF(uint64_t& L,uint64_t& H,uint64_t a,uint64_t b){
    unsigned __int128 x = (unsigned __int128)a * b;
    L = (uint64_t)x;
    H = (uint64_t)(x >> 64);
}
FORCE_INLINE uint64_t mulH2x1(uint64_t aL,uint64_t aH,uint64_t b){
    //  (aL,aH) * b
    unsigned __int128 temp = ((unsigned __int128)aH << 64) | aL;
    temp *= b;
    return (uint64_t)(temp >> 64);
}
inline void invert64(uint64_t& L,uint64_t& H,uint64_t p){
    //  (L,H) = 2^128 / p
    //  This function is not performance critical.
    unsigned __int128 temp = 0;
    temp -= 1;
    temp /= p;
    L = (uint64_t)temp;
    H = (uint64_t)(temp >> 64);
}
#endif
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
