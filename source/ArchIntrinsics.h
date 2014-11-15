/* ArchIntrinsics.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/01/2014
 * Last Modified    : 11/09/2014
 * 
 *      This file contains all the bignum primitives.
 * 
 *  Most of these are extremely performance-critical and their implementations
 *  tend to be compiler-specific.
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
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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
FORCE_INLINE void add128(uint64_t& L,uint64_t& H,uint64_t x){
    char carry = _addcarry_u64(0,L,x,&L);
    _addcarry_u64(carry,H,0,&H);
}
FORCE_INLINE void muladd(uint64_t& L,uint64_t& H,uint64_t A,uint64_t B,uint64_t C){
    L = _umul128(A,B,&H);
    char carry = _addcarry_u64(0,L,C,&L);
    _addcarry_u64(carry,H,0,&H);
}
FORCE_INLINE void muladd2(uint64_t& L,uint64_t& H,uint64_t A,uint64_t B,uint64_t C,uint64_t D){
    L = _umul128(A,B,&H);
    char carry = _addcarry_u64(0,L,C,&L);
    _addcarry_u64(carry,H,0,&H);
    carry = _addcarry_u64(0,L,D,&L);
    _addcarry_u64(carry,H,0,&H);
}
FORCE_INLINE void add_and_carry(uint64_t* T,const uint64_t* A,size_t L){
    //Conditions:
    //  -   2 <= L

    char carry = _addcarry_u64(0,T[0],A[0],&T[0]);
    size_t c = 1;
    do{
        carry = _addcarry_u64(carry,T[c],A[c],&T[c]);
        c++;
    }while (c < L);

    if (carry == 0)
        return;

    while (++T[c++] == 0);
}
FORCE_INLINE void sub(uint64_t* T,const uint64_t* A,size_t L){
    //Conditions:
    //  -   2 <= L

    char carry = _subborrow_u64(0,T[0],A[0],&T[0]);
    size_t c = 1;
    do{
        carry = _subborrow_u64(carry,T[c],A[c],&T[c]);
        c++;
    }while (c < L);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
FORCE_INLINE void add128(uint64_t& L,uint64_t& H,uint64_t x){
    unsigned __int128 temp = ((unsigned __int128)H << 64) | L;
    temp += x;
    L = (uint64_t)temp;
    H = (uint64_t)(temp >> 64);
}
FORCE_INLINE void muladd(uint64_t& L,uint64_t& H,uint64_t A,uint64_t B,uint64_t C){
    unsigned __int128 temp = (unsigned __int128)A * B + C;
    L = (uint64_t)temp;
    H = (uint64_t)(temp >> 64);
}
FORCE_INLINE void muladd2(uint64_t& L,uint64_t& H,uint64_t A,uint64_t B,uint64_t C,uint64_t D){
    unsigned __int128 temp = (unsigned __int128)A * B + C + D;
    L = (uint64_t)temp;
    H = (uint64_t)(temp >> 64);
}
FORCE_INLINE void add_and_carry(uint64_t* T,const uint64_t* A,size_t L){
    //Conditions:
    //  -   2 <= L

    unsigned __int128 temp = (unsigned __int128)T[0] + A[0];
    T[0] = (uint64_t)temp;
    size_t c = 1;
    do{
        temp >>= 64;
        temp += T[c];
        temp += A[c];
        T[c] = (uint64_t)temp;
        c++;
    }while (c < L);

    temp >>= 64;
    if (temp == 0)
        return;

    while (++T[c++] == 0);
}
FORCE_INLINE void sub(uint64_t* T,const uint64_t* A,size_t L){
    //Conditions:
    //  -   2 <= L

    __int128 temp = (__int128)T[0] - A[0];
    T[0] = (uint64_t)temp;
    size_t c = 1;
    do{
        temp >>= 64;
        temp += T[c];
        temp -= A[c];
        T[c] = (uint64_t)temp;
        c++;
    }while (c < L);
}
#endif
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Shared
FORCE_INLINE uint64_t mul(uint64_t* T,const uint64_t* A,size_t L,uint64_t x){
    //  T = A * x

    //Conditions:
    //  -   2 <= L

    uint64_t carry;
    mulF(T[0],carry,A[0],x);
    size_t c = 1;
    do{
        muladd(T[c],carry,A[c],x,carry);
        c++;
    }while (c < L);
    return carry;
}
FORCE_INLINE uint64_t muladd(uint64_t* T,const uint64_t* A,size_t L,uint64_t x){
    //  T += A * x

    //Conditions:
    //  -   2 <= L

    uint64_t carry;
    muladd(T[0],carry,A[0],x,T[0]);
    size_t c = 1;
    do{
        muladd2(T[c],carry,A[c],x,carry,T[c]);
        c++;
    }while (c < L);
    return carry;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
