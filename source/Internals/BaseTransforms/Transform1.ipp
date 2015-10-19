/* Transform1.ipp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/02/2014
 * Last Modified    : 11/02/2014
 * 
 *      Transform Length = 1 * 2^k
 */

#include "../Modulus.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void transform4_forward(const Modulus& p, uint64_t T[4]){
    uint64_t r0, r1, r2, r3;
    r0 = T[0];
    r1 = T[1];
    r2 = T[2];
    r3 = T[3];
    p.butterfly2_forward(r0, r2);
    p.butterfly2_forward(r1, r3, p.forward_table[2][1]);
    p.butterfly2_forward(r0, r1);
    p.butterfly2_forward(r2, r3);
    T[0] = r0;
    T[1] = r1;
    T[2] = r2;
    T[3] = r3;
}
void transform4_inverse_fmul(const Modulus& p, uint64_t T[4], const uint64_t A[4]){
    uint64_t r0, r1, r2, r3;
    r0 = T[0];
    r1 = T[1];
    r2 = T[2];
    r3 = T[3];
    r0 = p.mulmod(r0, p.make_twiddle(A[0]));
    r1 = p.mulmod(r1, p.make_twiddle(A[1]));
    r2 = p.mulmod(r2, p.make_twiddle(A[2]));
    r3 = p.mulmod(r3, p.make_twiddle(A[3]));
    p.butterfly2_inverse(r0, r1);
    p.butterfly2_inverse(r2, r3);
    p.butterfly2_inverse(r0, r2);
    p.butterfly2_inverse(r1, r3, p.inverse_table[2][1]);
    T[0] = r0;
    T[1] = r1;
    T[2] = r2;
    T[3] = r3;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
