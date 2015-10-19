/* ModulusSet.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/03/2014
 * Last Modified    : 11/09/2014
 * 
 */

#include <string.h>
#include "Transforms.h"
#include "ModulusSet.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
ModulusSet::ModulusSet(const PrimeSet& set)
    : multiplier(set.m)
    , factors(set.factors)
    , max_k(set.max_k)
    , operand_words((factors - 1) / 2)
    , carry_words(factors - operand_words)
    , set(&set)
{
    if (factors > MAX_PRIMES)
        throw "Maximum # of primes exceeded.";

    p.reserve(factors);
    for (int c = 0; c < factors; c++){
        p.emplace_back(set, c);
    }
}
void ModulusSet::ensure_tables(int k){
    if (k <= table_k)
        return;

    for (auto& modulus : p){
        modulus.make_tables(k);
    }
    table_k = k;
}
uint64_t ModulusSet::table_bytes(int k) const{
    if (k < Modulus::BASE_K)
        k = Modulus::BASE_K;

    //  (k + 2) because...
    //  Each table is k - 1                     : 2^(k - 1)
    //  Each twiddle has 2 words. (W and Wp)    : 2^1
    //  Forward and Inverse tables.             : 2^1
    //  Geometric summation of all the tables.  : 2^1
    //  Total                                   : 2^(k + 2)

    return sizeof(uint64_t) * factors * multiplier << (k + 2);
}
void ModulusSet::forward(int k, uint64_t* T) const{
    size_t stride = (size_t)multiplier << k;
    for (int f = 0; f < factors; f++){
        transform_forward(p[f], k, T + f*stride);
    }
}
void ModulusSet::inverse_fmul(int k, uint64_t* T, const uint64_t* A) const{
    size_t stride = (size_t)multiplier << k;
    for (int f = 0; f < factors; f++){
        transform_inverse_fmul(p[f], k, T + f*stride, A + f*stride);
    }
}
void ModulusSet::start_block(const uint64_t* R, uint64_t* T, size_t stride) const{
    //Parameters:
    //  -   [R, operand_block)
    //  -   T is the transform array

    T[0] = p[0].generate(R, operand_words);
    int c = 1;
    do{
        T[c * stride] = p[c].generate(R, operand_words);
        c++;
    }while (c < factors);
}
void ModulusSet::finish_block(int k, uint64_t* carry, uint64_t* R, const uint64_t* T, size_t stride) const{
    //Parameters:
    //  -   [R, operand_block)
    //  -   [carry, carry_words)
    //  -   T is the transform array

    //  Construct CRT
    uint64_t P[MAX_PRIMES + 2];
    CRT_build(k, P, T, stride);
    CRT_reduce1(P);
    CRT_reduce2(P);

    //  Propagate carryout
    carryout(carry, R, P);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  CRT Construction
void ModulusSet::CRT_build(int k, uint64_t P[MAX_PRIMES + 2], const uint64_t* W, size_t stride) const{
    const uint64_t* CRT = set->E;

    //  First factor.
    {
        uint64_t word = p[0].scale_down(k, W[0]);
        P[factors + 0] = mul(P, CRT, factors, word);
        P[factors + 1] = 0;
    }

    //  Rest of the factors.
    int f = 1;
    do{
        CRT += factors;
        uint64_t word = p[f].scale_down(k, W[f * stride]);
        uint64_t carry = muladd(P, CRT, factors, word);
        add128(P[factors], P[factors + 1], carry);
        f++;
    }while (f < factors);
}
void ModulusSet::CRT_reduce1(uint64_t P[MAX_PRIMES + 2]) const{
    //  Yeah okay... We need to do a multi-precision modulus.
    //  The quotient can be slightly more than 64 bits large. So two
    //  reductions is unavoidable assuming word-size arithmetic.
    //  To make things easier, we'll use floating-point.

    //  Get a floating-point estimate of P.
    double f = P[factors + 1] * 18446744073709551616. + P[factors];

    //  Multiply by reciprocal.
    f *= set->R;

    //  Scale down by 2^32. Purturb downward to avoid overestimating.
    f *= (1. / 4294967296) * 0.9999999999999964472863211994990706443786621093750;   //  2^-32 * (1 - 2^-48)

    uint64_t quo = (int64_t)f;

    //  Multiply
    uint64_t S[MAX_PRIMES + 1];
    S[factors] = mul(S, set->P, factors, quo);

    //  Shift up by 32-bits.
    memmove((uint32_t*)S + 1, S, (factors * 2 + 1) * sizeof(uint32_t));
    memset(S, 0, sizeof(uint32_t));

    //  Subtract
    sub(P, S, factors + 1);
}
void ModulusSet::CRT_reduce2(uint64_t P[MAX_PRIMES + 1]) const{
    //  This is the second and final reduction.

    //  Get a floating-point estimate of P.
    double f = P[factors - 1] * (1. / 18446744073709551616.) + P[factors];

    //  Multiply by reciprocal.
    f *= set->R;

    //  No purturb is needed. In fact, this quotient will always be correct.
    //  The quotient can only be wrong if the fractional part is close to 0.5
    //  where it can be rounded the wrong way. But this is impossible because
    //  we have ensured that the coefficient will always have 1 bit to spare.
    //  Therefore, no final correction will be needed.
    uint64_t quo = (int64_t)(f + 0.5);

    //  Multiply
    uint64_t S[MAX_PRIMES];
    mul(S, set->P, factors, quo);

    //  Subtract
    sub(P, S, factors);
}
void ModulusSet::carryout(uint64_t* carry, uint64_t* T, uint64_t* P) const{
    //  Add carryout into accumulator P.
    add_and_carry(P, carry, factors - operand_words);

    //  Copy lower words into destination.
    memcpy(T, P, operand_words * sizeof(uint64_t));

    //  Copy upper words into the carry buffer.
    memcpy(carry, P + operand_words, carry_words * sizeof(uint64_t));
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
