/* Transform3.ipp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/02/2014
 * Last Modified    : 11/02/2014
 * 
 *      Transform Length = 3 * 2^k
 */

#include "Modulus.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <bool forward> FORCE_INLINE
void transform3(const Modulus& p,uint64_t T[3]){
    //  Inputs must be fully reduced.
    //  Outputs are fully reduced.

    uint64_t r0,r1,r2;
    r0 = T[0];
    r1 = T[1];
    r2 = T[2];

    uint64_t e1,e2;
    uint64_t rA,rB,rC;
    e1 = r1;
    e2 = r2;
    rA = r0;
    rB = r0;
    rC = r0;

    {
        rA = p.reduce_p(rA + r1);
        rA = p.reduce_p(rA + r2);
        T[0] = rA;
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][1]
            : p.inverse_table[0][1];
        uint64_t bw,cw;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        rB = p.reduce_p(rB + bw);
        rC = p.reduce_p(rC + cw);
    }

    T[1] = p.reduce_n((int64_t)rB - (int64_t)e2);
    T[2] = p.reduce_n((int64_t)rC - (int64_t)e1);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void transform3_forward(const Modulus& p,uint64_t T[3]){
    transform3<true>(p,T);
}
void transform3_inverse_fmul(const Modulus& p,uint64_t T[3],const uint64_t A[3]){
    T[0] = p.mulmod(T[0],p.make_twiddle(A[0]));
    T[1] = p.mulmod(T[1],p.make_twiddle(A[1]));
    T[2] = p.mulmod(T[2],p.make_twiddle(A[2]));
    transform3<false>(p,T);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void transform12_forward(const Modulus& p,uint64_t T[12]){
    const auto& tw2 = p.forward_table[2];
    const auto& tw1 = p.forward_table[1];
    {
        uint64_t r0,r1,r2,r3;
        r0 = T[0];
        r1 = T[3];
        r2 = T[6];
        r3 = T[9];
        p.butterfly2_forward(r0,r2);
        p.butterfly2_forward(r1,r3,tw2[3]);
        p.butterfly2_forward(r0,r1);
        p.butterfly2_forward(r2,r3);
        T[0] = r0;
        T[3] = r1;
        T[6] = r2;
        T[9] = r3;
    }
    for (int c = 1; c < 3; c++){
        uint64_t r0,r1,r2,r3;
        r0 = T[c + 0];
        r1 = T[c + 3];
        r2 = T[c + 6];
        r3 = T[c + 9];
        p.butterfly2_forward(r0,r2,tw2[c + 0]);
        p.butterfly2_forward(r1,r3,tw2[c + 3]);
        p.butterfly2_forward(r0,r1,tw1[c + 0]);
        p.butterfly2_forward(r2,r3,tw1[c + 0]);
        T[c + 0] = r0;
        T[c + 3] = r1;
        T[c + 6] = r2;
        T[c + 9] = r3;
    }

    transform3_forward(p,T + 0);
    transform3_forward(p,T + 3);
    transform3_forward(p,T + 6);
    transform3_forward(p,T + 9);
}
void transform12_inverse_fmul(const Modulus& p,uint64_t T[12],const uint64_t A[12]){
    transform3_inverse_fmul(p,T + 0,A + 0);
    transform3_inverse_fmul(p,T + 3,A + 3);
    transform3_inverse_fmul(p,T + 6,A + 6);
    transform3_inverse_fmul(p,T + 9,A + 9);

    const auto& tw2 = p.inverse_table[2];
    const auto& tw1 = p.inverse_table[1];
    {
        uint64_t r0,r1,r2,r3;
        r0 = T[0];
        r1 = T[3];
        r2 = T[6];
        r3 = T[9];
        p.butterfly2_inverse(r0,r1);
        p.butterfly2_inverse(r2,r3);
        p.butterfly2_inverse(r0,r2);
        p.butterfly2_inverse(r1,r3,tw2[3]);
        T[0] = r0;
        T[3] = r1;
        T[6] = r2;
        T[9] = r3;
    }
    for (int c = 1; c < 3; c++){
        uint64_t r0,r1,r2,r3;
        r0 = T[c + 0];
        r1 = T[c + 3];
        r2 = T[c + 6];
        r3 = T[c + 9];
        p.butterfly2_inverse(r0,r1,tw1[c + 0]);
        p.butterfly2_inverse(r2,r3,tw1[c + 0]);
        p.butterfly2_inverse(r0,r2,tw2[c + 0]);
        p.butterfly2_inverse(r1,r3,tw2[c + 3]);
        T[c + 0] = r0;
        T[c + 3] = r1;
        T[c + 6] = r2;
        T[c + 9] = r3;
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
