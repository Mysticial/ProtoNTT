/* Main.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/01/2014
 * Last Modified    : Always...
 * 
 *      Since ProtoNTT is really just a library, there really is no need for a
 *  main method. But it exists to provide a basic UI to run benchmarks and tests.
 * 
 */

#include "source/BasicTransformParameters.h"
#include "source/TestFramework.h"
#include "source/Menus.h"

using namespace ProtoNTT;

#include <limits.h>
#include <string.h>
#include <iostream>
using std::cout;
using std::endl;


void sample(){
    //  This is a small example that shows how to use ProtoNTT to perform a
    //  multiplication of two large numbers represented as little-endian arrays
    //  of 64-bit integers.

    //  Inputs: A and B with lengths of AL and BL respectively.

    uint64_t A[] = {6, 8, 8, 7, 1, 2, 6, 5, 6, 8};
    uint64_t B[] = {2, 2, 9, 4, 1, 7, 4, 6, 6, 1, 0, 9};

    const size_t AL = sizeof(A) / sizeof(uint64_t);
    const size_t BL = sizeof(B) / sizeof(uint64_t);

    //  Output: C with length AL + BL;
    const size_t CL = AL + BL;
    uint64_t C[CL];

    //  Make twiddle factor table.
    TwiddleTable table;

    //  Get transform parameters.
    BasicTransformParameters tp(CL,table);

    //  (Optional) Populate the table. The size is up to preference.
    //  A larger table will be able to run larger multiplications at full
    //  speed. But it will also require more memory.
    tp.ensure_tables(CL);

    //  Perform multiply.
    tp.mul(C,A,AL,B,BL);

    //  Print the result.
    print(C,CL);


    //  Notes:
    //      It is not necessary to have the table be large enough. If it is
    //  too small, it will still work correctly, but slower. The difference
    //  between having no table, and having a large enough table is about 3x
    //  in performance.
    //
    //      Naturally the twiddle table should be stored globally and used by
    //  all multiplications. But care must be taken to not modify the table
    //  concurrently as it is thread-safe at all.
    //
    //      An alternative to calling "ensure_tables()" at run-time is to call
    //  "table.populate_all_tables()" at the start so that it is not necessary
    //  to resize it during run-time. But be aware that this will build the
    //  tables for all 16 modes which could use a lot of memory.
}

int main(){
    try{
        std::cout << "ProtoNTT Tester v1.0.0" << std::endl;
        std::cout << std::endl;

        std::cout << "Select an option:" << std::endl;
        std::cout << "    0     Benchmark: Equal Length Operands" << std::endl;
        std::cout << "    1     Benchmark: Different Length Operands" << std::endl;
        std::cout << "    2     Benchmark: Custom NTT Parameters" << std::endl;
        std::cout << "    3     Run Integration Tests" << std::endl;
        std::cout << std::endl;

        while (1){
            std::cout << "choice = ";
            int choice;
            std::cin >> choice;
            std::cin.ignore(INT_MAX,'\n');
            std::cout << std::endl;
            switch (choice){
                case 0:
                    bench_equal_length();
                    break;
                case 1:
                    bench_different_length();
                    break;
                case 2:
                    bench_parameters();
                    break;
                case 3:
                    run_integration_tests();
                    break;
                default:
                    continue;
            }
            break;
        };
    }catch (const char* err){
        std::cout << "Error: " << err << std::endl;
    }


#ifdef _WIN32
    pause();
#endif
}
