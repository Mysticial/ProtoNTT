ProtoNTT
========

ProtoNTT is a basic implementation of the [Small Primes Number-Theoretic Transform (NTT) algorithm](http://www.apfloat.org/ntt.html) for multiplication of large integers.

ProtoNTT is an early prototype for the Small Primes NTT that was added to [y-cruncher v0.6.8](http://www.numberworld.org/y-cruncher/). This prototype has most of the low-effort/high-payoff optimizations. But it lacks all the difficult memory optimizations that are needed to make it viable in a serious bignum application. Nevertheless, ProtoNTT is still "reasonably" efficient and is open-sourced here for educational purposes.

**Build Instructions/Requirements:**<br>
 - Hardware: x64 is required. The program will not compile for 32-bit.
 - Compiler:
     - Windows: Visual Studio 2013 update 2 or later. (update 2 is needed for the add-with-carry intrinsics)
     - Linux: GCC 4.8 or later.

A Visual Studio project is included. The Linux version can be built by simply running `compile-gcc.sh`.

**Comparison with y-cruncher v0.6.8's NTT:**<br>

|Feature(s)                         |ProtoNTT      |y-cruncher v0.6.8                         |
|-----------------------------------|--------------|------------------------------------------|
|Prime Size                         |63-bit        |63-bit                                    |
|# of Primes                        |3, 5, 7, or 9 |3, 4, 5, 6, 7, 8, or 9                    |
|Transform Lengths                  |2<sup>k</sup>, 3·2<sup>k</sup>, 5·2<sup>k</sup>, 7·2<sup>k</sup>|2<sup>k</sup>, 3·2<sup>k</sup>, 5·2<sup>k</sup>, 7·2<sup>k</sup>|
|Maximum Convolution Length         |~`10^18` bits |~`10^18` bits                             |
|Parallelization                    |No            |Yes                                       |
|Swap Mode (Out-of-core Computation)|No            |Yes                                       |
|Supported Architectures            |x64           |Any 32-bit or 64-bit                      |
|Processor-Specific Optimizations   |x64           |x86, x64, SSE4.2, XOP, AVX2, AVX512-(F/DQ)|
|Butterfly Radix                    |Radix 2       |Generalized Bailey's 4-step               |
|Twiddle Factor Table Size          |O(N)          |O(1)                                      |

The only thing that prevents ProtoNTT from being a viable algorithm is that:
 - It uses much memory when all twiddle factors are precomputed. (4x the transform size)
 - It is too slow when no twiddle factors are precomputed. (3x slower)

The space-time trade-off curve for precomputing a subset of twiddle factors is particularly unfavorable for ProtoNTT. The cuts needed to reduce the table down to a reasonable size will easily slap on a performance penalty of 50% or more to an already slow algorithm.
