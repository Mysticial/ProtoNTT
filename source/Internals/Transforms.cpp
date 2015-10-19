/* Transforms.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/02/2014
 * Last Modified    : 11/06/2014
 * 
 */

#include "Transforms.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void transform_forward(const Modulus& p, int k, uint64_t* T){
    if (k == 2){
        p.forward(T);
        return;
    }

    //  Radix-2 Reduction
    size_t block = (size_t)p.multiplier << (k - 1);

    {
        //  Pull out first iteration with no twiddle factor.
        p.butterfly2_forward(T[0], T[block]);
    }

    if (k <= p.table_k){
        //  All required twiddle factors are in the table.
        size_t c = 1;
        const auto& tw = p.forward_table[k];
        do{
            p.butterfly2_forward(T[c], T[c + block], tw[c]);
            c++;
        }while (c < block);
    }else{
        //  Need to generate twiddle factors on the fly.
        //  This will be very slow. :(
        TwiddleFactor current = p.stub_f[k];
        TwiddleFactor increment = current;
        size_t c = 1;
        do{
            p.butterfly2_forward(T[c], T[c + block], current);
            current = p.make_twiddle(p.mulmod(current.W, increment));
            c++;
        }while (c < block);
    }

    //  Sub-transforms
    transform_forward(p, k - 1, T);
    transform_forward(p, k - 1, T + block);
}
void transform_inverse_fmul(const Modulus& p, int k, uint64_t* T, const uint64_t* A){
    if (k == 2){
        p.inverse_fmul(T, A);
        return;
    }

    //  Radix-2 Reduction
    size_t block = (size_t)p.multiplier << (k - 1);

    //  Sub-transforms
    transform_inverse_fmul(p, k - 1, T, A);
    transform_inverse_fmul(p, k - 1, T + block, A + block);

    {
        //  Pull out first iteration with no twiddle factor.
        p.butterfly2_inverse(T[0], T[block]);
    }

    if (k <= p.table_k){
        //  All required twiddle factors are in the table.
        const auto& tw = p.inverse_table[k];
        size_t c = 1;
        do{
            p.butterfly2_inverse(T[c], T[c + block], tw[c]);
            c++;
        }while (c < block);
    }else{
        //  Need to generate twiddle factors on the fly.
        //  This will be very slow. :(
        TwiddleFactor current = p.stub_i[k];
        TwiddleFactor increment = current;
        size_t c = 1;
        do{
            p.butterfly2_inverse(T[c], T[c + block], current);
            current = p.make_twiddle(p.mulmod(current.W, increment));
            c++;
        }while (c < block);
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
