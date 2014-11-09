/* Transform5.ipp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/02/2014
 * Last Modified    : 11/02/2014
 * 
 *      Transform Length = 5 * 2^k
 */

#include "Modulus.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <bool forward> FORCE_INLINE
void transform5(const Modulus& p,uint64_t T[5]){
    //  Inputs must be fully reduced.
    //  Outputs are fully reduced.

    uint64_t r0,r1,r2,r3,r4;
    r0 = T[0];
    r1 = T[1];
    r2 = T[2];
    r3 = T[3];
    r4 = T[4];

    uint64_t e1,e2,e3,e4;
    uint64_t rA,rB,rC,rD,rE;
    e1 = r1;
    e2 = r2;
    e3 = r3;
    e4 = r4;
    rA = r0;
    rB = r0;
    rC = r0;
    rD = r0;
    rE = r0;

    {
        rA = p.reduce_p(rA + r1);
        rA = p.reduce_p(rA + r2);
        rA = p.reduce_p(rA + r3);
        rA = p.reduce_p(rA + r4);
        T[0] = rA;
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][1]
            : p.inverse_table[0][1];
        uint64_t bw,cw,dw,ew;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        dw = p.mulmod(r3,w1);
        ew = p.mulmod(r4,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        e3 = p.reduce_p(e3 + dw);
        e4 = p.reduce_p(e4 + ew);
        rB = p.reduce_p(rB + bw);
        rC = p.reduce_p(rC + dw);
        rD = p.reduce_p(rD + cw);
        rE = p.reduce_p(rE + ew);
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][2]
            : p.inverse_table[0][2];
        uint64_t bw,cw,dw,ew;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        dw = p.mulmod(r3,w1);
        ew = p.mulmod(r4,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        e3 = p.reduce_p(e3 + dw);
        e4 = p.reduce_p(e4 + ew);
        rB = p.reduce_p(rB + cw);
        rC = p.reduce_p(rC + bw);
        rD = p.reduce_p(rD + ew);
        rE = p.reduce_p(rE + dw);
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][3]
            : p.inverse_table[0][3];
        uint64_t bw,cw,dw,ew;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        dw = p.mulmod(r3,w1);
        ew = p.mulmod(r4,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        e3 = p.reduce_p(e3 + dw);
        e4 = p.reduce_p(e4 + ew);
        rB = p.reduce_p(rB + dw);
        rC = p.reduce_p(rC + ew);
        rD = p.reduce_p(rD + bw);
        rE = p.reduce_p(rE + cw);
    }

    T[1] = p.reduce_n((int64_t)rB - (int64_t)e4);
    T[2] = p.reduce_n((int64_t)rC - (int64_t)e2);
    T[3] = p.reduce_n((int64_t)rD - (int64_t)e3);
    T[4] = p.reduce_n((int64_t)rE - (int64_t)e1);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void transform5_forward(const Modulus& p,uint64_t T[5]){
    transform5<true>(p,T);
}
void transform5_inverse_fmul(const Modulus& p,uint64_t T[5],const uint64_t A[5]){
    T[0] = p.mulmod(T[0],p.make_twiddle(A[0]));
    T[1] = p.mulmod(T[1],p.make_twiddle(A[1]));
    T[2] = p.mulmod(T[2],p.make_twiddle(A[2]));
    T[3] = p.mulmod(T[3],p.make_twiddle(A[3]));
    T[4] = p.mulmod(T[4],p.make_twiddle(A[4]));
    transform5<false>(p,T);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void transform20_forward(const Modulus& p,uint64_t T[20]){
    const auto& tw2 = p.forward_table[2];
    const auto& tw1 = p.forward_table[1];
    {
        uint64_t r0,r1,r2,r3;
        r0 = T[ 0];
        r1 = T[ 5];
        r2 = T[10];
        r3 = T[15];
        p.butterfly2_forward(r0,r2);
        p.butterfly2_forward(r1,r3,tw2[5]);
        p.butterfly2_forward(r0,r1);
        p.butterfly2_forward(r2,r3);
        T[ 0] = r0;
        T[ 5] = r1;
        T[10] = r2;
        T[15] = r3;
    }
    for (int c = 1; c < 5; c++){
        uint64_t r0,r1,r2,r3;
        r0 = T[c +  0];
        r1 = T[c +  5];
        r2 = T[c + 10];
        r3 = T[c + 15];
        p.butterfly2_forward(r0,r2,tw2[c + 0]);
        p.butterfly2_forward(r1,r3,tw2[c + 5]);
        p.butterfly2_forward(r0,r1,tw1[c + 0]);
        p.butterfly2_forward(r2,r3,tw1[c + 0]);
        T[c +  0] = r0;
        T[c +  5] = r1;
        T[c + 10] = r2;
        T[c + 15] = r3;
    }

    transform5_forward(p,T +  0);
    transform5_forward(p,T +  5);
    transform5_forward(p,T + 10);
    transform5_forward(p,T + 15);
}
void transform20_inverse_fmul(const Modulus& p,uint64_t T[20],const uint64_t A[20]){
    transform5_inverse_fmul(p,T +  0,A +  0);
    transform5_inverse_fmul(p,T +  5,A +  5);
    transform5_inverse_fmul(p,T + 10,A + 10);
    transform5_inverse_fmul(p,T + 15,A + 15);

    const auto& tw2 = p.inverse_table[2];
    const auto& tw1 = p.inverse_table[1];
    {
        uint64_t r0,r1,r2,r3;
        r0 = T[ 0];
        r1 = T[ 5];
        r2 = T[10];
        r3 = T[15];
        p.butterfly2_inverse(r0,r1);
        p.butterfly2_inverse(r2,r3);
        p.butterfly2_inverse(r0,r2);
        p.butterfly2_inverse(r1,r3,tw2[5]);
        T[ 0] = r0;
        T[ 5] = r1;
        T[10] = r2;
        T[15] = r3;
    }
    for (int c = 1; c < 5; c++){
        uint64_t r0,r1,r2,r3;
        r0 = T[c +  0];
        r1 = T[c +  5];
        r2 = T[c + 10];
        r3 = T[c + 15];
        p.butterfly2_inverse(r0,r1,tw1[c + 0]);
        p.butterfly2_inverse(r2,r3,tw1[c + 0]);
        p.butterfly2_inverse(r0,r2,tw2[c + 0]);
        p.butterfly2_inverse(r1,r3,tw2[c + 5]);
        T[c +  0] = r0;
        T[c +  5] = r1;
        T[c + 10] = r2;
        T[c + 15] = r3;
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
