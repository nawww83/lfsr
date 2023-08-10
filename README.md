# Linear Feedback Shift Registers
Linear Feedback Shift Registers, LFSR, and its cryptographic applications.

Prime number $p$ and register length $m$ are used.

It contains some functions to find vectors of coefficients $K$ that provide periods $T_0 = {p}^{m} - 1$ and $T_1 = {p}^{m-1} - 1$.

**Cryptographic** 16-bit hash is calculated using XOR combination of two LFSR pairs: $p = 251$ and $p = 241$ with the periods $T_0$ and $T_1$ for each $p$. Also, 32-bit and 64-bit **cryptographic** hashes are added.

## Build
g++ main.cpp lfsr_hash.cpp -std=c++20 -msse4.1 -O3 -o lfsr
## Run
./lfsr

## Hash performance
Approx. 210 MB/s @ Intel i7-8565U CPU 4.2GHz, GCC 11.4