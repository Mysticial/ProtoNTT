/* Menus.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/07/2014
 * Last Modified    : 11/15/2014
 * 
 */

#include <limits.h>
#include <iostream>
#include "BasicTransformParameters.h"
#include "Menus.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int prompt_table_reduction(){
    std::cout << std::endl;
    std::cout << "    Reduce Table Size (0 for maximum size)" << std::endl;
    std::cout << "    Reduce by: 2^";
    int table_reduction;
    std::cin >> table_reduction;
    std::cin.ignore(INT_MAX,'\n');
    return table_reduction;
}
void bench_equal_length(){
    std::cout << "Benchmark Equal Length Multiply:" << std::endl << std::endl;

    std::cout << "    Length of operand (in 64-bit words): ";
    size_t L;
    std::cin >> L;
    std::cin.ignore(INT_MAX,'\n');

    int table_reduction = prompt_table_reduction();
    std::cout << std::endl << std::endl;

    TwiddleTable table;
    BasicTransformParameters(2*L,table).bench_multiply(L,L,table_reduction);
}
void bench_different_length(){
    std::cout << "Benchmark Variable Length Multiply:" << std::endl << std::endl;

    std::cout << "    Length of operand A (in 64-bit words): ";
    size_t AL;
    std::cin >> AL;
    std::cin.ignore(INT_MAX,'\n');

    std::cout << "    Length of operand B (in 64-bit words): ";
    size_t BL;
    std::cin >> BL;
    std::cin.ignore(INT_MAX,'\n');

    int table_reduction = prompt_table_reduction();
    std::cout << std::endl << std::endl;

    TwiddleTable table;
    BasicTransformParameters(AL + BL,table).bench_multiply(AL,BL,table_reduction);
}
void bench_parameters(){
    std::cout << "Benchmark NTT Parameters:" << std::endl << std::endl;

    int primes;
    do{
        std::cout << "    # of primes (3, 5, 7, or 9): ";
        std::cin >> primes;
        std::cin.ignore(INT_MAX,'\n');
    }while (primes < 3 || primes > 9 || primes % 2 != 1);

    std::cout << std::endl;
    std::cout << "    Select a transform length: m * 2^k   (m must be 1, 3, 5, or 7)" << std::endl;

    int multiplier;
    do{
        std::cout << "        m = ";
        std::cin >> multiplier;
        std::cin.ignore(INT_MAX,'\n');
    }while (multiplier < 1 || multiplier > 7 || multiplier % 2 != 1);

    int k;
    do{
        std::cout << "        k = ";
        std::cin >> k;
        std::cin.ignore(INT_MAX,'\n');
    }while (k < 2 || k > 55);

    int table_reduction = prompt_table_reduction();
    std::cout << std::endl << std::endl;

    TwiddleTable table;
    BasicTransformParameters tp(primes,multiplier,k,table);
    size_t L = tp.get_cbitlen() / 128;
    tp.bench_multiply(L,L,table_reduction);
}
void run_integration_tests(){
    std::cout << "Running start-to-finish tests for all modes..." << std::endl << std::endl;

    int k = 10;
    TwiddleTable table;
    BasicTransformParameters(3,1,k,table).test();
    BasicTransformParameters(3,3,k,table).test();
    BasicTransformParameters(3,5,k,table).test();
    BasicTransformParameters(3,7,k,table).test();
    BasicTransformParameters(5,1,k,table).test();
    BasicTransformParameters(5,3,k,table).test();
    BasicTransformParameters(5,5,k,table).test();
    BasicTransformParameters(5,7,k,table).test();
    BasicTransformParameters(7,1,k,table).test();
    BasicTransformParameters(7,3,k,table).test();
    BasicTransformParameters(7,5,k,table).test();
    BasicTransformParameters(7,7,k,table).test();
    BasicTransformParameters(9,1,k,table).test();
    BasicTransformParameters(9,3,k,table).test();
    BasicTransformParameters(9,5,k,table).test();
    BasicTransformParameters(9,7,k,table).test();
    std::cout << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
