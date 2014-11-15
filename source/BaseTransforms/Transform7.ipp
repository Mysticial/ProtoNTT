/* Transform7.ipp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/02/2014
 * Last Modified    : 11/02/2014
 * 
 *      Transform Length = 7 * 2^k
 */

#include "../Modulus.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <bool forward> FORCE_INLINE
void transform7(const Modulus& p,uint64_t T[7]){
    //  Inputs must be fully reduced.
    //  Outputs are fully reduced.

    uint64_t r0,r1,r2,r3,r4,r5,r6;
    r0 = T[0];
    r1 = T[1];
    r2 = T[2];
    r3 = T[3];
    r4 = T[4];
    r5 = T[5];
    r6 = T[6];

    uint64_t e1,e2,e3,e4,e5,e6;
    uint64_t rA,rB,rC,rD,rE,rF,rG;
    e1 = r1;
    e2 = r2;
    e3 = r3;
    e4 = r4;
    e5 = r5;
    e6 = r6;
    rA = r0;
    rB = r0;
    rC = r0;
    rD = r0;
    rE = r0;
    rF = r0;
    rG = r0;

    {
        rA = p.reduce_p(rA + r1);
        rA = p.reduce_p(rA + r2);
        rA = p.reduce_p(rA + r3);
        rA = p.reduce_p(rA + r4);
        rA = p.reduce_p(rA + r5);
        rA = p.reduce_p(rA + r6);
        T[0] = rA;
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][1]
            : p.inverse_table[0][1];
        uint64_t bw,cw,dw,ew,fw,gw;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        dw = p.mulmod(r3,w1);
        ew = p.mulmod(r4,w1);
        fw = p.mulmod(r5,w1);
        gw = p.mulmod(r6,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        e3 = p.reduce_p(e3 + dw);
        e4 = p.reduce_p(e4 + ew);
        e5 = p.reduce_p(e5 + fw);
        e6 = p.reduce_p(e6 + gw);
        rB = p.reduce_p(rB + bw);
        rC = p.reduce_p(rC + ew);
        rD = p.reduce_p(rD + fw);
        rE = p.reduce_p(rE + cw);
        rF = p.reduce_p(rF + dw);
        rG = p.reduce_p(rG + gw);
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][2]
            : p.inverse_table[0][2];
        uint64_t bw,cw,dw,ew,fw,gw;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        dw = p.mulmod(r3,w1);
        ew = p.mulmod(r4,w1);
        fw = p.mulmod(r5,w1);
        gw = p.mulmod(r6,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        e3 = p.reduce_p(e3 + dw);
        e4 = p.reduce_p(e4 + ew);
        e5 = p.reduce_p(e5 + fw);
        e6 = p.reduce_p(e6 + gw);
        rB = p.reduce_p(rB + cw);
        rC = p.reduce_p(rC + bw);
        rD = p.reduce_p(rD + dw);
        rE = p.reduce_p(rE + ew);
        rF = p.reduce_p(rF + gw);
        rG = p.reduce_p(rG + fw);
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][3]
            : p.inverse_table[0][3];
        uint64_t bw,cw,dw,ew,fw,gw;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        dw = p.mulmod(r3,w1);
        ew = p.mulmod(r4,w1);
        fw = p.mulmod(r5,w1);
        gw = p.mulmod(r6,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        e3 = p.reduce_p(e3 + dw);
        e4 = p.reduce_p(e4 + ew);
        e5 = p.reduce_p(e5 + fw);
        e6 = p.reduce_p(e6 + gw);
        rB = p.reduce_p(rB + dw);
        rC = p.reduce_p(rC + fw);
        rD = p.reduce_p(rD + bw);
        rE = p.reduce_p(rE + gw);
        rF = p.reduce_p(rF + cw);
        rG = p.reduce_p(rG + ew);
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][4]
            : p.inverse_table[0][4];
        uint64_t bw,cw,dw,ew,fw,gw;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        dw = p.mulmod(r3,w1);
        ew = p.mulmod(r4,w1);
        fw = p.mulmod(r5,w1);
        gw = p.mulmod(r6,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        e3 = p.reduce_p(e3 + dw);
        e4 = p.reduce_p(e4 + ew);
        e5 = p.reduce_p(e5 + fw);
        e6 = p.reduce_p(e6 + gw);
        rB = p.reduce_p(rB + ew);
        rC = p.reduce_p(rC + cw);
        rD = p.reduce_p(rD + gw);
        rE = p.reduce_p(rE + bw);
        rF = p.reduce_p(rF + fw);
        rG = p.reduce_p(rG + dw);
    }
    {
        const TwiddleFactor& w1 = forward
            ? p.forward_table[0][5]
            : p.inverse_table[0][5];
        uint64_t bw,cw,dw,ew,fw,gw;
        bw = p.mulmod(r1,w1);
        cw = p.mulmod(r2,w1);
        dw = p.mulmod(r3,w1);
        ew = p.mulmod(r4,w1);
        fw = p.mulmod(r5,w1);
        gw = p.mulmod(r6,w1);
        e1 = p.reduce_p(e1 + bw);
        e2 = p.reduce_p(e2 + cw);
        e3 = p.reduce_p(e3 + dw);
        e4 = p.reduce_p(e4 + ew);
        e5 = p.reduce_p(e5 + fw);
        e6 = p.reduce_p(e6 + gw);
        rB = p.reduce_p(rB + fw);
        rC = p.reduce_p(rC + gw);
        rD = p.reduce_p(rD + ew);
        rE = p.reduce_p(rE + dw);
        rF = p.reduce_p(rF + bw);
        rG = p.reduce_p(rG + cw);
    }

    T[1] = p.reduce_n((int64_t)rB - (int64_t)e6);
    T[2] = p.reduce_n((int64_t)rC - (int64_t)e3);
    T[3] = p.reduce_n((int64_t)rD - (int64_t)e2);
    T[4] = p.reduce_n((int64_t)rE - (int64_t)e5);
    T[5] = p.reduce_n((int64_t)rF - (int64_t)e4);
    T[6] = p.reduce_n((int64_t)rG - (int64_t)e1);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void transform7_forward(const Modulus& p,uint64_t T[7]){
    transform7<true>(p,T);
}
void transform7_inverse_fmul(const Modulus& p,uint64_t T[7],const uint64_t A[7]){
    T[0] = p.mulmod(T[0],p.make_twiddle(A[0]));
    T[1] = p.mulmod(T[1],p.make_twiddle(A[1]));
    T[2] = p.mulmod(T[2],p.make_twiddle(A[2]));
    T[3] = p.mulmod(T[3],p.make_twiddle(A[3]));
    T[4] = p.mulmod(T[4],p.make_twiddle(A[4]));
    T[5] = p.mulmod(T[5],p.make_twiddle(A[5]));
    T[6] = p.mulmod(T[6],p.make_twiddle(A[6]));
    transform7<false>(p,T);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void transform28_forward(const Modulus& p,uint64_t T[28]){
    const auto& tw2 = p.forward_table[2];
    const auto& tw1 = p.forward_table[1];
    {
        uint64_t r0,r1,r2,r3;
        r0 = T[ 0];
        r1 = T[ 7];
        r2 = T[14];
        r3 = T[21];
        p.butterfly2_forward(r0,r2);
        p.butterfly2_forward(r1,r3,tw2[7]);
        p.butterfly2_forward(r0,r1);
        p.butterfly2_forward(r2,r3);
        T[ 0] = r0;
        T[ 7] = r1;
        T[14] = r2;
        T[21] = r3;
    }
    for (int c = 1; c < 7; c++){
        uint64_t r0,r1,r2,r3;
        r0 = T[c +  0];
        r1 = T[c +  7];
        r2 = T[c + 14];
        r3 = T[c + 21];
        p.butterfly2_forward(r0,r2,tw2[c + 0]);
        p.butterfly2_forward(r1,r3,tw2[c + 7]);
        p.butterfly2_forward(r0,r1,tw1[c + 0]);
        p.butterfly2_forward(r2,r3,tw1[c + 0]);
        T[c +  0] = r0;
        T[c +  7] = r1;
        T[c + 14] = r2;
        T[c + 21] = r3;
    }

    transform7_forward(p,T +  0);
    transform7_forward(p,T +  7);
    transform7_forward(p,T + 14);
    transform7_forward(p,T + 21);
}
void transform28_inverse_fmul(const Modulus& p,uint64_t T[28],const uint64_t A[28]){
    transform7_inverse_fmul(p,T +  0,A +  0);
    transform7_inverse_fmul(p,T +  7,A +  7);
    transform7_inverse_fmul(p,T + 14,A + 14);
    transform7_inverse_fmul(p,T + 21,A + 21);

    const auto& tw2 = p.inverse_table[2];
    const auto& tw1 = p.inverse_table[1];
    {
        uint64_t r0,r1,r2,r3;
        r0 = T[ 0];
        r1 = T[ 7];
        r2 = T[14];
        r3 = T[21];
        p.butterfly2_inverse(r0,r1);
        p.butterfly2_inverse(r2,r3);
        p.butterfly2_inverse(r0,r2);
        p.butterfly2_inverse(r1,r3,tw2[7]);
        T[ 0] = r0;
        T[ 7] = r1;
        T[14] = r2;
        T[21] = r3;
    }
    for (int c = 1; c < 7; c++){
        uint64_t r0,r1,r2,r3;
        r0 = T[c +  0];
        r1 = T[c +  7];
        r2 = T[c + 14];
        r3 = T[c + 21];
        p.butterfly2_inverse(r0,r1,tw1[c + 0]);
        p.butterfly2_inverse(r2,r3,tw1[c + 0]);
        p.butterfly2_inverse(r0,r2,tw2[c + 0]);
        p.butterfly2_inverse(r1,r3,tw2[c + 7]);
        T[c +  0] = r0;
        T[c +  7] = r1;
        T[c + 14] = r2;
        T[c + 21] = r3;
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
