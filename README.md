ProtoNTT
========

ProtoNTT is a basic implementation of the [Small Primes Number-Theoretic Transform (NTT) algorithm](http://www.apfloat.org/ntt.html) for multiplication of large integers.

As its name implies, ProtoNTT is a prototype. While it has all the basic low-effort/high-payoff CPU optimizations, it lacks the difficult-to-do memory optimizations that are needed to make it viable in a serious bignum application. Nevertheless, ProtoNTT is still "reasonably" efficient and is open-sourced here for educational purposes.

**Motivation:**<br>
ProtoNTT is an experiment to test the viability of the Small Primes NTT algorithm against the [y-cruncher project](http://www.numberworld.org/y-cruncher/).

Historically, the Small Primes NTT has been a popular last-resort algorithm for extremely large multiplications where Floating-Point FFT is no longer viable. It is als onotorious for being slow. But this has changed ever since Victor Shoup found an efficient way to perform a multiply-modulus using only 3 multiplications and no divisions.

**Results:**<br>

Unfortunately, ProtoNTT still turned out to be quite a disappointment. At best it barely comes within a factor of 2x against y-cruncher on Sandy Bridge (and worse on other processors). While ProtoNTT still has room for optimizations, overcoming a difference of > 2x seems a bit of a long-shot.

Nevertheless, y-cruncher will probably get a Small Primes NTT anyway as it's necessary to overcome its current limitation of approximately 90 trillion digits. At that size, performance matters less as it will probably be disk-bound.

**ProtoNTT Features:**<br>
 - Uses Victor Shoup's butterfly with 63-bit primes.
 - Supports transform sizes of: `2^k`, `3*2^k`, `5*2^k`, and `7*2^k`.
 - Supports 3, 5, 7, and 9 primes.
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

**Problems with ProtoNTT:**
 - **The size of the twiddle factor table is horrendous:**
     - Precomputing all the necessary twiddle factors requires 4x more memory than the transform itself. Proportionally, it needs 4x more than the Floating-Point FFT. In the FFT, you save a factor of 2 because forward/inverse twiddles are complex conjugates. In the NTT, you can't do that and you lose a factor of 2 because you also need `Wp = W*2^64 / p` for each twiddle. Precomputing only *some* of the twiddles leads to a space-time trade-off.
 - **The algorithm will not vectorize on current hardware:**
     - There is no vectorized 64 x 64 -> 128-bit multiply on x86. AVX512-DQ has `vpmullq` which gets the lower-half, but we still need the top half.
     - There is no 64-bit integer compare prior to SSE4.2. And then, SSE4.2 brings only a greater-than compare. AMD's XOP finally gets it right and AVX512-F is perfect since it also embedded masking designed specifically for conditional operations.
     - The CRT construction step requires bignum arithmetic. Again we need a vectorized 64 x 64 -> 128-bit multiply to be efficient. As a consolation prize, we have MULX and ADCX/ADOX.
 

