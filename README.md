# Linear Feedback Shift Registers (LFSR) in cryptography
LFSR and its cryptographic applications.

Prime number $p$ and register length $m$ are used.

You can calculate 32, 64 and 128-bit **cryptographic** hashes with **salt**. Also, you can add 256-bit hash etc using similar extention approach as you need.
Additionally, src code contains some functions to find vectors of coefficients $K$ that provide periods $T_0 = {p}^{m} - 1$ and $T_1 = {p}^{m-1} - 1$.

## Build
g++ main.cpp lfsr_hash.cpp -std=c++20 -msse4.1 -O3 -o lfsr
## Run
./lfsr

## Hash performance
Approx. 210 MB/s @ Intel i7-8565U CPU 4.2GHz, GCC 11.4

## LFSR hash principles
Base unit for LFSR hash is LFSR register which can be interpreted as Finite State Machine (Moore machine if be more precise).

Classicaly, LFSR operates with binary symbols, $0$ and $1$, but nobody forbid us to use any $M$-ary symbol - integer numbers from $0$ to $M-1$.

Due to special mathematical property of prime numbers, it makes sense to use $M$ as a prime number: $2, 3, 5, 7, \ldots{}$ That property is the maximal period property when we calculate all integer powers by modulo $p$: ${g}^{i} \mod p, \quad i \in \mathbb{Z}$.

For example, let's $p=11$, $g=2$, $i=0...10$, so we have the following powers $1, 2, 4, 8, 5, 10, 9, 7, 3, 6, 1$.

We can see that all numbers are different (except the last), from $1$ to $10$, and ${g}^{0} = {g}^{p-1} \mod p = 1$ or equivalently $g^p - g = 0$.
The last expression is true for all $a$, but not all $a$ give us unique powers. The maximal number of unique powers is $p-1$. Some $a$ give us maximal period $T = p - 1$, other ones give us smaller periods, for example, $T/2$ or smaller.

The previous example is a $p$-ary LFSR register wtih unit length, $m=1$, but nobody forbid us to use any length $m$. For that case maximal period is $T = {p}^{m} - 1$, and we should find $m$ integer numbers $\left( g_0, \ldots{}, g_{m-1} \right)$, which will provide maximal period. The $m$ numbers is the coefficients of the generator polynomial $g(x) = x^m - {g}_{m-1} {x}^{m-1} - \ldots{} - g_0$ with degree $m$.

In general, we use $v$ as the last element of LFSR state $\vec s$, multiply generator vector $\vec g$ by the scalar $v$, and then add the latter to the 1-delayed LFSR state (padded by zero). So we can express LFSR loop in vector form
$$\vec s = \left( v \vec g + D[\vec s, 1] \right) \mod p,$$
$$v = {s}_{m-1}.$$

Higher symbols of $\vec s$, i.e. $s_m$, ${s}_{m+1},...$ can be ignored.

Here the delay operator $D[\vec v, 1]$ means $(v_0, v_1, v_2, ...) \rightarrow (0, v_0, v_1, v_2, ...)$ transformation for any vector $\vec v$.

As a rule, LFSR is initialized by unit state $\vec s = (1, 0, 0, ... , 0)$. Having $m$ cycles, LFSR will be in the state $\vec s = \vec g$ exactly, which is the same as generator coefficients. So, we consider that LFSR is **saturated** at that moment. When we continue cycles, we will observe some different states, and at some moment the current state will be equal to the initial state $\vec a$. What cycles we done between two equal states will determine LFSR period $T$.

If a generator $\vec g$ provides maximal period $T = T_{max} = p^m - 1$, then LFSR evolutes all possible states $\vec s$ except zero-state, and it doesn't matter what the initial state has been set. For a period $T < {T}_{\max}$ the initial state has some influence, but most of initial states will provide the fix period which can be called as the "main period".

It was numerically shown there are generators $\vec g$ which provide period $T_1 = {p}^{m-1} - 1$, wherein the maximal period is marked as $T_0 \equiv {T}_{\max}$. In general, we can write $T_q = {p}^{m-q} - 1$.

LFSR generator can have the input port, $a$. For that more general case we can rewrite LFSR loop
$$\vec s = \left( v \vec g + D[\vec s, a, 1] \right) \mod p,$$
$$v = {s}_{m-1}.$$

Here the delay operator $D[\vec v, a, 1]$ means $(v_0, v_1, v_2, ...) \rightarrow (a, v_0, v_1, v_2, ...)$ transformation for any vector $\vec v$ and scalar $a$.

Note, all discussed LFSR periods have sense for zero input. That periods can be called as "free periods".

If we have two LFSR generators, with periods $T_0$ and $T_1$ and common input, we can combine (mix) the final LFSR states into one state that can be interpreted as LFSR hash of the input. Two periods have Greatest Common Divisor $p-1$, so for zero input we will have the total period
$${{T_0 T_1} \over {p-1}} \approx p^{2m-2} \approx {T_1}^{2}.$$

For states combination, the bit-wise XOR operator is chosen. We chose $m=4$ with $p=251$ and $p=241$.

Two generators ${\vec g}^{(0)}$ and ${\vec g}^{(1)}$ are packed into one SIMD 128-bit vector in C++ code for better SIMD registers usage. Such generator is called as **LFSR pair**. For better performance, the input is processed by 16-bit samples.

We use two LFSR pairs, with $p=251$ and $p=241$. Final states of two pairs are XOR-ed as usual. So, the total free period of given LFSR hash is about ${251}^{6} \cdot {241}^{6} \approx 95$ bits, i.e. ${2}^{95} \approx {10}^{28}$.

To improve the crypto resilience, the first LFSR pair is driven by the original input, but the second pair is driven by reversed input vector. Also, the salt is added before and after input loop. The minimal size of the input vector is $2$ elements.

Bit scaling of LFSR hash is done by salt adding.