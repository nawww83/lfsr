# Linear Feedback Shift Registers
Linear Feedback Shift Registers, LFSR, and its cryptographic applications.

Prime number $p$ and register length $m$ are used.

You can calculate 16, 32, 64 and 128-bit **cryptographic** hashes.
Additionally, src code contains some functions to find vectors of coefficients $K$ that provide periods $T_0 = {p}^{m} - 1$ and $T_1 = {p}^{m-1} - 1$.

## Build
g++ main.cpp lfsr_hash.cpp -std=c++20 -msse4.1 -O3 -o lfsr
## Run
./lfsr

## Hash performance
Approx. 210 MB/s @ Intel i7-8565U CPU 4.2GHz, GCC 11.4