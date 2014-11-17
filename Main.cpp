/* Main.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/01/2014
 * Last Modified    : Always...
 * 
 *      Since ProtoNTT is really just a library, there really is no need for a
 *  main method. But it exists to develop and test the library.
 * 
 *  When it's finished, I may put a small UI here that will run and verify
 *  test multiplications using random data.
 * 
 */

#include "source/BasicTransformParameters.h"
#include "source/TestFramework.h"

using namespace ProtoNTT;

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

int main(){

    run_integration_tests();


#ifdef _WIN32
    pause();
#endif
}
