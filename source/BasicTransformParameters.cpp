/* BasicTransformParameters.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/06/2014
 * Last Modified    : 11/15/2014
 * 
 */

#include <string.h>
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
void BasicTransformParameters::raw_to_NTT(uint64_t* T, const uint64_t* R, size_t RL) const{
    //  Convert from raw operand into NTT array.

    size_t operand_words = set->operand_words;

    size_t index_R = 0;
    size_t index_T = 0;

    size_t stride = (size_t)set->multiplier << k;
    if (RL*64 > cbitlen)
        throw "Convolution length is too small.";

    if (RL >= operand_words){
        size_t stop = RL - operand_words + 1;
        while (index_R < stop){
            set->start_block(R + index_R, T + index_T, stride);
            index_R += operand_words;
            index_T++;
        }
    }

    //  Take care of partial top word.
    if (index_R < RL){
        uint64_t P[ModulusSet::MAX_PRIMES];
        size_t left = RL - index_R;
        memcpy(P, R + index_R, left * sizeof(uint64_t));
        memset(P + left, 0, (operand_words - left) * sizeof(uint64_t));
        set->start_block(P, T + index_T, stride);
        index_T++;
    }

    //  Zero-pad the rest of the transform.
    for (int f = 0; f < primes; f++){
        memset(T + index_T + f*stride, 0, (stride - index_T) * sizeof(uint64_t));
    }
}
void BasicTransformParameters::NTT_to_raw(const uint64_t* T, uint64_t* R, size_t RL) const{
    //  Convert from NTT array to raw operand applying carryout.

    size_t operand_words = set->operand_words;

    size_t stride = (size_t)set->multiplier << k;
    if (RL*64 > cbitlen)
        throw "Convolution length is too small.";

    size_t index_R = 0;
    size_t index_T = 0;

    uint64_t carry[ModulusSet::MAX_PRIMES];
    memset(carry, 0, sizeof(carry));

    if (RL >= operand_words){
        size_t stop = RL - operand_words + 1;
        while (index_R < stop){
            set->finish_block(k, carry, R + index_R, T + index_T, stride);
            index_R += operand_words;
            index_T++;
        }
    }

    //  Take care of partial top word.
    if (index_R < RL){
        uint64_t P[ModulusSet::MAX_PRIMES];
        set->finish_block(k, carry, P, T + index_T, stride);
        memcpy(R + index_R, P, (RL - index_R) * sizeof(uint64_t));
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Ready-to-go Multiplication
void BasicTransformParameters::sqr(uint64_t* C, const uint64_t* A, size_t AL) const{
    size_t TL = (size_t)multiplier << k;

    std::unique_ptr<uint64_t> T_uptr(new uint64_t[TL * primes]);
    uint64_t* T = T_uptr.get();

    raw_to_NTT(T, A, AL);         //  Convert 1st operand.
    set->forward(k, T);          //  Forward transform 1st operand.
    set->inverse_fmul(k, T, T);   //  Pointwise multiply and inverse transform.
    NTT_to_raw(T, C, 2*AL);       //  Convert back and carryout.
}
void BasicTransformParameters::mul(uint64_t* C, const uint64_t* A, size_t AL, const uint64_t* B, size_t BL) const{
    size_t TL = (size_t)multiplier << k;

    std::unique_ptr<uint64_t> T_uptr(new uint64_t[TL * primes]);
    std::unique_ptr<uint64_t> U_uptr(new uint64_t[TL * primes]);
    uint64_t* T = T_uptr.get();
    uint64_t* U = U_uptr.get();

    raw_to_NTT(T, A, AL);         //  Convert 1st operand.
    set->forward(k, T);          //  Forward transform 1st operand.
    raw_to_NTT(U, B, BL);         //  Convert 1st operand.
    set->forward(k, U);          //  Forward transform 1st operand.
    set->inverse_fmul(k, T, U);   //  Pointwise multiply and inverse transform.
    NTT_to_raw(T, C, AL + BL);    //  Convert back and carryout.
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Benchmark
void BasicTransformParameters::time_benchmark(int primes, int multiplier, int k, double seconds, int table_reduction){
    TwiddleTable table;
    BasicTransformParameters tp(primes, multiplier, k, table);

    size_t L = tp.cbitlen / 128;
    size_t CL = 2*L;

    std::cout << "p = " << primes << ",  m = " << multiplier << ",  k = " << k << ", cbitlen = ";
    print_commas(tp.cbitlen); std::cout << "\t";

    //  Make Table
    size_t table_size = CL >> table_reduction;
    tp.ensure_tables(table_size);

    //  Generate Operands
    std::unique_ptr<uint64_t> O_uptr(new uint64_t[CL]);
    uint64_t* A = O_uptr.get();
    uint64_t* B = A + L;
    random(A, L);
    random(B, L, L);

    //  Allocate Transforms
    std::unique_ptr<uint64_t> T_uptr(new uint64_t[tp.Tsize / sizeof(uint64_t)]);
    std::unique_ptr<uint64_t> U_uptr(new uint64_t[tp.Tsize / sizeof(uint64_t)]);
    uint64_t* T = T_uptr.get();
    uint64_t* U = U_uptr.get();
    memset(T, 0, tp.Tsize);
    memset(U, 0, tp.Tsize);

    //  Benchmark
    uint64_t iterations = 0;
    double start = wall_clock();
    double total_time;
    do{
        tp.raw_to_NTT(T, A, L);           //  Convert 1st operand.
        tp.set->forward(k, T);           //  Forward transform 1st operand.
        tp.raw_to_NTT(U, B, L);           //  Convert 1st operand.
        tp.set->forward(k, U);           //  Forward transform 1st operand.
        tp.set->inverse_fmul(k, T, U);    //  Pointwise multiply and inverse transform.
        tp.NTT_to_raw(T, A, CL);          //  Convert back and carryout.
        iterations++;
    }while ((total_time = wall_clock() - start) < seconds);

    double seconds_per_iteration = total_time / (double)iterations;
    std::cout << seconds_per_iteration << "\tsecs/iteration" << std::endl;
}
void BasicTransformParameters::bench_multiply(int table_reduction){
    size_t L = cbitlen / 128;
    bench_multiply(L, L, table_reduction);
}
void BasicTransformParameters::bench_multiply(size_t AL, size_t BL, int table_reduction){
    std::cout << "Benchmark Multiplication" << std::endl;
    std::cout << std::endl;

    size_t CL = AL + BL;

    std::cout << "Length A   = "; print_words(AL);
    std::cout << "Length B   = "; print_words(BL);
    std::cout << "Length A*B = "; print_words(CL);
    std::cout << std::endl;
    print();
    std::cout << std::endl;

    size_t table_size = CL >> table_reduction;

    uint64_t Osize = CL * sizeof(uint64_t);
    uint64_t Msize = Tsize*2 + Psize;
    uint64_t Wsize = table_bytes(table_size);

    std::cout << "Operand Size  = "; print_bytes(Osize);
    std::cout << "Multiply Size = "; print_bytes(Msize);
    std::cout << "Table Size    = "; print_bytes(Wsize);
    std::cout << std::endl;
    std::cout << "Total Memory Required:  "; print_bytes(Osize + Msize + Wsize);
    std::cout << std::endl << std::endl;
    pause();
    std::cout << std::endl << std::endl;

    std::cout << "Constructing Twiddle Tables...    ";
    double time0 = wall_clock();
    ensure_tables(table_size);
    std::cout << wall_clock() - time0 << std::endl;

    std::cout << "Generating Operands...            ";
    double time1 = wall_clock();
    std::unique_ptr<uint64_t> O_uptr(new uint64_t[CL]);
    uint64_t* A = O_uptr.get();
    uint64_t* B = A + AL;
    random(A, AL);
    random(B, BL, AL);
    std::cout << wall_clock() - time1 << std::endl;

    std::cout << "Allocating Transforms...          ";
    double time2 = wall_clock();
    std::unique_ptr<uint64_t> T_uptr(new uint64_t[Tsize / sizeof(uint64_t)]);
    std::unique_ptr<uint64_t> U_uptr(new uint64_t[Tsize / sizeof(uint64_t)]);
    uint64_t* T = T_uptr.get();
    uint64_t* U = U_uptr.get();
    memset(T, 0, Tsize);
    memset(U, 0, Tsize);
    std::cout << wall_clock() - time2 << std::endl;

    std::cout << "Hashing Operands...               ";
    double time3 = wall_clock();
    uint64_t hashA = hash_compute(A, AL);
    uint64_t hashB = hash_compute(B, BL);
    std::cout << wall_clock() - time3 << std::endl;

    std::cout << std::endl;

    std::cout << "Generating Residuals A...         ";
    double time4 = wall_clock();
    raw_to_NTT(T, A, AL);
    std::cout << wall_clock() - time4 << std::endl;

    std::cout << "Forward Transform A...            ";
    double time5 = wall_clock();
    set->forward(k, T);
    std::cout << wall_clock() - time5 << std::endl;

    std::cout << "Generating Residuals B...         ";
    double time6 = wall_clock();
    raw_to_NTT(U, B, BL);
    std::cout << wall_clock() - time6 << std::endl;

    std::cout << "Forward Transform B...            ";
    double time7 = wall_clock();
    set->forward(k, U);
    std::cout << wall_clock() - time7 << std::endl;

    std::cout << "Inverse Transform...              ";
    double time8 = wall_clock();
    set->inverse_fmul(k, T, U);
    std::cout << wall_clock() - time8 << std::endl;

    std::cout << "Constructing CRT...               ";
    double time9 = wall_clock();
    NTT_to_raw(T, A, CL);
    double time10 = wall_clock();
    std::cout << time10 - time9 << std::endl;

    std::cout << std::endl;
    std::cout << "Total Multiply Time...            ";
    std::cout << time10 - time4 << std::endl;
    std::cout << std::endl;

    std::cout << "Hashing Result...                 ";
    double time11 = wall_clock();
    uint64_t hashC = hash_compute(A, CL);
    uint64_t hashD = hash_mul(hashA, hashB);
    std::cout << wall_clock() - time11 << std::endl;
    std::cout << std::endl;

    std::cout << "Expected Hash = " << hashD << std::endl;
    std::cout << "Actual Hash   = " << hashC << std::endl;
    std::cout << std::endl;
}
void BasicTransformParameters::test() const{
    std::cout << "p = " << primes << ",  m = " << multiplier << ",  k = " << k << "   :   ";

    size_t L = cbitlen / 128;
    set->ensure_tables(k);

    //  Allocate Operands
    std::unique_ptr<uint64_t> O_uptr(new uint64_t[2*L]);
    uint64_t* A = O_uptr.get();
    uint64_t* B = A + L;
    random(A, L);
    random(B, L, L);

    //  Hash Operands
    uint64_t hashA = hash_compute(A, L);
    uint64_t hashB = hash_compute(B, L);

    //  Multiply
    mul(A, A, L, B, L);

    //  Hash Result
    uint64_t hashC = hash_compute(A, 2*L);
    uint64_t hashD = hash_mul(hashA, hashB);

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
