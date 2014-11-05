ProtoNTT
========

ProtoNTT is a work-in-progress implementation of the Small Primes Number-Theoretic Transform (NTT) algorithm for multiplication of large integers.

As its name implies, ProtoNTT is a prototype and is not intended for use in serious applications.
It is open-sourced here for educational purposes.

**Motivation:**<br>
ProtoNTT is yet another by-product of the [y-cruncher project](http://www.numberworld.org/y-cruncher/). As of 2014, the world record of 13.3 trillion digits of Pi is within an order-of-magnitude of the limit of the program (~90 trillion). So there will eventually be a need to implement a new algorithm that be able to handle these extreme sizes.

Past experience has shown that these large-scale multiply algorithms tend to quickly devolve into incredible works of mess once they are optimized, parallelized, and equipped with swap mode support. Typically, the only way to clean it up is a complete rewrite. This time I'll try to get it right the first time by using a throw-away prototype to get some foresight on the complications that will be encountered.

**(Intended) ProtoNTT Features:**<br>
 - Uses Victor Shoup's butterfly with 63-bit primes.
 - Supports transform sizes of: `2^k`, `3*2^k`, `5*2^k`, and `7*2^k`.
 - Supports any odd number of primes.
 - Maximum convolution length is upwards of `10^18` bits.
 - Basic CPU optimizations for x64. (x64 is required)
 - Basic micro-optimizations.
 - Basic radix-2 recursive FFT.

**Missing Features:**
 - No cache or memory optimizations.
 - No attempt to be efficient when twiddle factors are not precomputed.
 - No parallelization.
 - No vectorization.
 - No swap mode (out-of-core) support.

