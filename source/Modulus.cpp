/* Modulus.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/01/2014
 * Last Modified    : 11/05/2014
 * 
 */

#include "Modulus.h"
#include "Transform1.ipp"
//#include "Transform3.ipp"
//#include "Transform5.ipp"
//#include "Transform7.ipp"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Constructors
Modulus::Modulus(const PrimeSet& set,int index)
    : ModularRing(set.factor[index])
    , multiplier(set.m)
    , max_k(set.max_k)
{
    //  Set the function pointers for the end-point transforms.
    switch (multiplier){
        case 1:
            fp_forward = transform4_forward;
            fp_inverse_fmul = transform4_inverse_fmul;
            break;
//        case 3:
//            fp_forward = transform12_forward;
//            fp_inverse_fmul = transform12_inverse_fmul;
//            break;
//        case 5:
//            fp_forward = transform20_forward;
//            fp_inverse_fmul = transform20_inverse_fmul;
//            break;
//        case 7:
//            fp_forward = transform28_forward;
//            fp_inverse_fmul = transform28_inverse_fmul;
//            break;
        default:
            throw "Unsupported Multiplier";
    }

    make_scaling_factors(set.scaler[index],set.max_k);
    make_table_stubs(set.root_f[index]);
    sanity_check();
}
void Modulus::make_scaling_factors(uint64_t scaler,int max_k){
    TwiddleFactor base = make_twiddle(scaler);

    TwiddleFactor inv_2 = make_twiddle(mulmod(multiplier,base));
    uint64_t inv_m = mulmod(2,base);

    scaling_factors.clear();
    scaling_factors.emplace_back(make_twiddle(inv_m));
    for (int k = 1; k <= max_k; k++){
        inv_m = mulmod(inv_m,inv_2);
        scaling_factors.emplace_back(make_twiddle(inv_m));
    }
}
void Modulus::make_table_stubs(uint64_t primitive_root){
    uint64_t max_pow = (uint64_t)multiplier << max_k;

    stub_f.resize(max_k + 1);
    stub_i.resize(max_k + 1);

    int c = max_k;
    stub_f[c] = make_twiddle(primitive_root);
    stub_i[c] = make_twiddle(power(make_twiddle(primitive_root),max_pow - 1));

    while (c-- > 0){
        stub_f[c] = make_twiddle(mulmod(stub_f[c + 1].W,stub_f[c + 1]));
        stub_i[c] = make_twiddle(mulmod(stub_i[c + 1].W,stub_i[c + 1]));
    }
}
void Modulus::sanity_check() const{
    //  Check the root-of-unity.
    {
        const TwiddleFactor& root_of_unity = stub_f.back();
        uint64_t test_power = (uint64_t)multiplier << (max_k - 1);
        uint64_t actual = power(root_of_unity,test_power);
        uint64_t expected = prime - 1;
        if (actual != expected)
            throw "Root of unity does not satisfy: root^(m*2^(k - 1)) mod p = -1";
    }

    //  Check the scaling factor.
    {
        uint64_t actual = mulmod(multiplier,scaling_factors[0]);
        uint64_t expected = 1;
        if (actual != expected)
            throw "Scaling factor does not satisfy: scale * 2m mod p = 1";
    }
}
void Modulus::make_tables(int k){
    forward_table.clear();
    inverse_table.clear();

    if (k < BASE_K)
        k = BASE_K;

    {
        //  Pull out first table and pad with empty sub-tables.
        forward_table.emplace_back(std::vector<TwiddleFactor>());
        inverse_table.emplace_back(std::vector<TwiddleFactor>());
    }
    for (int c = 1; c <= k; c++){
        //  Make empty sub-tables.
        forward_table.emplace_back(std::vector<TwiddleFactor>());
        inverse_table.emplace_back(std::vector<TwiddleFactor>());

        auto& forward = forward_table.back();
        auto& inverse = inverse_table.back();

        //  Half a full cycle.
        size_t stop = (size_t)multiplier << (c - 1);
        forward.reserve(stop);
        inverse.reserve(stop);

        //  First twiddle
        {
            TwiddleFactor ONE = make_twiddle(1);
            forward.emplace_back(ONE);
            inverse.emplace_back(ONE);
        }
        for (size_t i = 1; i < stop; i++){
            uint64_t a = mulmod(forward[i - 1].W,stub_f[c]);
            uint64_t b = mulmod(inverse[i - 1].W,stub_i[c]);
            forward.emplace_back(make_twiddle(a));
            inverse.emplace_back(make_twiddle(b));
        }
    }
    {
        //  Build multiplier tables.
        auto& forward = forward_table[1];
        auto& inverse = inverse_table[1];

        for (int c = 0; c < multiplier; c++){
            uint64_t a = mulmod(forward[c].W,forward[c]);
            uint64_t b = mulmod(inverse[c].W,inverse[c]);
            forward_table[0].emplace_back(make_twiddle(a));
            inverse_table[0].emplace_back(make_twiddle(b));
        }
    }

    table_k = k;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Scalar Operations
uint64_t Modulus::power(const TwiddleFactor& x,uint64_t pow) const{
    if (pow == 0)
        return 1;

    const uint64_t TOP_BIT = 0x8000000000000000;

    int c = 64;
    while ((pow & TOP_BIT) == 0){
        c--;
        pow <<= 1;
    }
    pow <<= 1;

    uint64_t out = x.W;
    while (--c > 0){
        out = mulmod(out,make_twiddle(out));

        if (pow & TOP_BIT){
            out = mulmod(out,x);
        }

        pow <<= 1;
    }

    return out;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
