# Linear Feedback Shift Registers
LFSR and its applications. Prime number p and register length m are used.
It contains function find_max_period_polynomial() to find LFSR coefficients K that provide maximal period T = p^m - 1.
## Build
g++ main.cpp -std=c++17 -O3 -o lfsr
## Run
./lfsr
## Possible output
LFSR with modulo p: 131, length m: 4
Wait for max period T = p^m - 1 polynomial look up...
Found coefficients K: (7, 1, 0, 0)
Period T: 294'499'920
