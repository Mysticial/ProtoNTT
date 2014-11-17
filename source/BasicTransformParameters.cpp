/* BasicTransformParameters.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/06/2014
 * Last Modified    : 11/15/2014
 * 
 */

#include <iostream>
#include <memory>
#include "BasicTransformParameters.h"
#include "TestFramework.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Misc.
void BasicTransformParameters::print() const{
    CoreTransformParameters::print();
    std::cout << "    Transform Bytes    = "; print_bytes(Tsize);
    std::cout << "    Scratch Bytes      = "; print_bytes(Psize);
}
void BasicTransformParameters::set_modset(TwiddleTable& tables){
    switch (primes){
        case 3:
            switch (multiplier){
                case 1: set = &tables.p3m1; break;
                case 3: set = &tables.p3m3; break;
                case 5: set = &tables.p3m5; break;
                case 7: set = &tables.p3m7; break;
                default: throw "Invalid multiplier.";
            }
            break;
        case 5:
            switch (multiplier){
                case 1: set = &tables.p5m1; break;
                case 3: set = &tables.p5m3; break;
                case 5: set = &tables.p5m5; break;
                case 7: set = &tables.p5m7; break;
                default: throw "Invalid multiplier.";
            }
            break;
        case 7:
            switch (multiplier){
                case 1: set = &tables.p7m1; break;
                case 3: set = &tables.p7m3; break;
                case 5: set = &tables.p7m5; break;
                case 7: set = &tables.p7m7; break;
                default: throw "Invalid multiplier.";
            }
            break;
        case 9:
            switch (multiplier){
                case 1: set = &tables.p9m1; break;
                case 3: set = &tables.p9m3; break;
                case 5: set = &tables.p9m5; break;
                case 7: set = &tables.p9m7; break;
                default: throw "Invalid multiplier.";
            }
            break;
        default:
            throw "Invalid # of primes.";
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Raw Conversion
void BasicTransformParameters::raw_to_NTT(uint64_t* T,const uint64_t* R,size_t RL) const{
    //  Convert from raw operand into NTT array.

    int operand_words = set->operand_words;

    size_t index_R = 0;
    size_t index_T = 0;

    size_t stride = (size_t)set->multiplier << k;
    if (RL*64 > cbitlen)
        throw "Convolution length is too small.";

    if (RL >= set->operand_words){
        size_t stop = RL - operand_words + 1;
        while (index_R < stop){
            set->start_block(R + index_R,T + index_T,stride);
            index_R += operand_words;
            index_T++;
        }
    }

    //  Take care of partial top word.
    if (index_R < RL){
        uint64_t P[ModulusSet::MAX_PRIMES];
        size_t left = RL - index_R;
        memcpy(P,R + index_R,left * sizeof(uint64_t));
        memset(P + left,0,(operand_words - left) * sizeof(uint64_t));
        set->start_block(P,T + index_T,stride);
        index_T++;
    }

    //  Zero-pad the rest of the transform.
    for (int f = 0; f < primes; f++){
        memset(T + index_T + f*stride,0,(stride - index_T) * sizeof(uint64_t));
    }
}
void BasicTransformParameters::NTT_to_raw(const uint64_t* T,uint64_t* R,size_t RL) const{
    //  Convert from NTT array to raw operand applying carryout.

    int operand_words = set->operand_words;

    size_t stride = (size_t)set->multiplier << k;
    if (RL*64 > cbitlen)
        throw "Convolution length is too small.";

    size_t index_R = 0;
    size_t index_T = 0;

    uint64_t carry[ModulusSet::MAX_PRIMES];
    memset(carry,0,sizeof(carry));

    if (RL >= operand_words){
        size_t stop = RL - operand_words + 1;
        while (index_R < stop){
            set->finish_block(k,carry,R + index_R,T + index_T,stride);
            index_R += operand_words;
            index_T++;
        }
    }

    //  Take care of partial top word.
    if (index_R < RL){
        uint64_t P[ModulusSet::MAX_PRIMES];
        set->finish_block(k,carry,P,T + index_T,stride);
        memcpy(R + index_R,P,(RL - index_R) * sizeof(uint64_t));
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Ready-to-go Multiplication
void BasicTransformParameters::sqr(uint64_t* C,const uint64_t* A,size_t AL) const{
    size_t TL = (size_t)multiplier << k;

    std::unique_ptr<uint64_t> T_uptr(new uint64_t[TL * primes]);
    uint64_t* T = T_uptr.get();

    raw_to_NTT(T,A,AL);         //  Convert 1st operand.
    set->forward(k,T);          //  Forward transform 1st operand.
    set->inverse_fmul(k,T,T);   //  Pointwise multiply and inverse transform.
    NTT_to_raw(T,C,2*AL);       //  Convert back and carryout.
}
void BasicTransformParameters::mul(uint64_t* C,const uint64_t* A,size_t AL,const uint64_t* B,size_t BL) const{
    size_t TL = (size_t)multiplier << k;

    std::unique_ptr<uint64_t> T_uptr(new uint64_t[TL * primes]);
    std::unique_ptr<uint64_t> U_uptr(new uint64_t[TL * primes]);
    uint64_t* T = T_uptr.get();
    uint64_t* U = U_uptr.get();

    raw_to_NTT(T,A,AL);         //  Convert 1st operand.
    set->forward(k,T);          //  Forward transform 1st operand.
    raw_to_NTT(U,B,BL);         //  Convert 1st operand.
    set->forward(k,U);          //  Forward transform 1st operand.
    set->inverse_fmul(k,T,U);   //  Pointwise multiply and inverse transform.
    NTT_to_raw(T,C,AL + BL);    //  Convert back and carryout.
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Other
void BasicTransformParameters::test() const{
    std::cout << "p = " << primes << ",  m = " << multiplier << ",  k = " << k << "   :   ";

    size_t L = cbitlen / 128;
    set->ensure_tables(k);

    //  Allocate Operands
    std::unique_ptr<uint64_t> O_uptr(new uint64_t[2*L]);
    uint64_t* A = O_uptr.get();
    uint64_t* B = A + L;
    random(A,L);
    random(B,L,L);

    //  Hash Operands
    uint64_t hashA = hash_compute(A,L);
    uint64_t hashB = hash_compute(B,L);

    //  Multiply
    mul(A,A,L,B,L);

    //  Hash Result
    uint64_t hashC = hash_compute(A,2*L);
    uint64_t hashD = hash_mul(hashA,hashB);

    if (hashC == hashD){
        std::cout << "Pass" << std::endl;
    }else{
        std::cout << "Fail" << std::endl;
        pause();
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
