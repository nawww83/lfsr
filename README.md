# Linear Feedback Shift Registers (LFSR) in cryptography
LFSR and its cryptographic applications.

Prime number $p$ and register length $m$ are used.

1. You can get 32, 64 and 128-bit **cryptographic** hashes with the **salt**. Also, you can add 256-bit hash etc using similar extention approach as you need.

It was estimated the proposed LFSR hash is comparable to **SHA-512**. The proposed implementation is both for Little and Big Endianess.

Additionally, the source code contains some functions to find vectors of coefficients $K$ that provide periods $T_0 = {p}^{m} - 1$ and $T_1 = {p}^{m-1} - 1$.

2. You can generate $64$-bit **cryptographic** random numbers based on two LFSR pairs (**total** length $4m$) with small primes $p$: $17$, $19$, and **Sawtooth modulation** of that pairs with small periods: $7$, $11$, see Random Generators V1 and V3.

## Build
g++ main.cpp lfsr_hash.cpp -std=c++20 -msse4.1 -O3 -o lfsr

## Run
./lfsr

## Hash performance
Approx. 210 MB/s @ Intel i7-8565U CPU 4.2GHz, GCC 11.4

## Random number generator performance
V1: ~130 MB/s @ Period: ~77 bit/sample
V3: ~50MB/s @ Period: ~155 bit/sample
@ Intel i7-8565U CPU 4.2GHz, GCC 12.3
1 sample = 16 bit.

## LFSR hash principles
The base unit for LFSR hash is a LFSR register, which can be interpreted as a Finite State Machine (Moore machine to be precise).

Usually, LFSR operates on binary symbols, $0$ and $1$, but no one forbid us to use any $M$-ary symbols - integer from $0$ to $M-1$.

Due to special mathematical property of prime numbers, it makes sense to use prime numbers as $M$: $2, 3, 5, 7, \ldots{}$. This property is the maximum period of a sequence of integer powers of some number $g$ by modulo $p$: ${g}^{i} \mod p, \quad i \in \mathbb{Z}$.

For example, consider $p=11$, $g=2$, $i=0...10$, then we have the following sequence of powers: $1, 2, 4, 8, 5, 10, 9, 7, 3, 6, 1$.

We can see that all numbers are different (except the last one), from $1$ to $10$, and ${g}^{0} = {g}^{p-1} \mod p = 1$ or equivalently $g^p - g = 0$. The last expression is true for all $g$, but not all $g$ provides a maximum period. The largest number of unique numbers is $p-1$. Some $g$ give us a maximum period $T = p - 1$, some don't.

In the previous example, a unit length $p$-ary LFSR register is considered, $m=1$, but no one forbids considering a register of arbitrary length $m$. In this case, the maximum period is $T = {p}^{m} - 1$, and we should find suitable coefficients $\left( g_0, \ldots{}, g_{m-1} \right)$ so that the period is just that. These coefficients are the coefficients of the generator polynomial $g(x) = x^m - {g}_{m-1} {x}^{m-1} - \ldots{} - g_0$ of degree $m$.

The working cycle of LFSR generator is as follows. We use the scalar $v$, wich holds the value of the last cell in LFSR state $\vec s$. We multiply generator vector $\vec g$ by scalar $v$, and then add the 1-delayed LFSR state (padded by zero) to the product. So, we can express LFSR cycle in vector form
$$\vec s = \left( v \vec g + D[\vec s, 1] \right) \mod p,$$
$$v = {s}_{m-1}.$$

The leading elements of the vector $\vec s$, i.e. $s_m$, ${s}_{m+1},...$ can be ignored.

Here the delay operator $D[\vec v, 1]$ means $(v_0, v_1, v_2, ...) \rightarrow (0, v_0, v_1, v_2, ...)$ transformation for some vector $\vec v$.

As a rule, LFSR is initialized by unit state $\vec s = (1, 0, 0, ... , 0)$. After $m$ cycles, LFSR will be in state $\vec s = \vec g$, which is equal to the generator coefficients. We consider that LFSR is **saturated** at that moment. When we continue cycles, we will observe different states, and at some point the current state will be equal to the initial one $\vec g$. The number of completed cycles will determine the period $T$ of the generator.

If some generator $\vec g$ provides the maximum period $T = T_{max} = p^m - 1$, then LFSR goes through all states $\vec s$ except zero, and it doesn't matter what initial state was set. For periods $T < {T}_{\max}$, the initial state affects the period, however most initial states will provide a fixed period, which can be called the **main period**.

It was numerically shown, there are generators $\vec g$ which provide period $T_1 = {p}^{m-1} - 1$, while the maximum period is marked as $T_0 \equiv {T}_{\max}$. In general, we can write $T_q = {p}^{m-q} - 1$.

In general, LFSR generator has an input port, $a$. Then its work cycle will look like this:
$$\vec s = \left( v \vec g + D[\vec s, a, 1] \right) \mod p,$$
$$v = {s}_{m-1}.$$

Here the delay operator $D[\vec v, a, 1]$ means $(v_0, v_1, v_2, ...) \rightarrow (a, v_0, v_1, v_2, ...)$ transformation for any vector $\vec v$ and scalar $a$.

Note that all the discussed LFSR periods make sense at zero input. These periods may be called **free periods**.

If we have two LFSR generators with periods $T_0$ and $T_1$ and a common input, then we can combine (mix) their final states into one, thereby forming the required LFSR hash of the input data. The Greatest Common Divisor of two periods is $p-1$. Thus, with a zero input, we get the total period
$${{T_0 T_1} \over {p-1}} \approx p^{2m-2} \approx {T_1}^{2}.$$

A bit-wise XOR is chosen to combine the states.

Two generators ${\vec g}^{(0)}$ and ${\vec g}^{(1)}$ are packed into one SIMD 128-bit vector in C++ code for better SIMD registers usage. Such generator is called as **LFSR pair**. For better performance, the input is processed by 16-bit samples.

We use two LFSR pairs of length $m=4$, and primes $p_1 = 251$ and $p_2 = 241$. Final states of two pairs are XOR-ed as usual. So, the total free period of given LFSR hash is equal to $$\text{LCM} \left( {{(p_1^4 - 1)(p_1^3 - 1)} \over {p_1 - 1}}; {{(p_2^4 - 1)(p_2^3 - 1)} \over {p_2 - 1}} \right) = 68'604'332'454'972'583'997'112'000,$$ i.e. $\approx 86$ bits. Here LCM - least common multiple.

To improve the crypto resilience, the first LFSR pair is driven by the original input, but the second one is driven by reversed input vector. Also, the **salt** is added before/after input loop.

Bit scaling of LFSR hash is done by adding **salt**.

## Random number generation principles
We have two LFSR pairs, one has $p=17$, the other has $p=19$. For example, lets consider the first pair. We can control the LFSR generator using a **sawtooth generator**:
$$i = (i + 1) \mod q.$$

Here $i$ - the output of the sawtooth generator with period $q < p$. We start generation from some initial value $i = i_0$. It is better to choose the sawtooth period $q$ so that the LFSR period $T = p^m - 1$ is not divided by that period. In this case we will visit all possible $i$ when the **fixed non-zero** LFSR state is repeated, while the total period of the system will be maximum and equal to $(p^m-1) \cdot q$.

Observations have shown that driving the LFSR generator by a proper sawtooth generator results in a sequence of $q$ LFSR **pseudo-random periods** $T_j$ such that their sum is equal to the total period of the system:
$${T}_{total} = \sum\_{j=1}^{q} T_j = (p^m-1) \cdot q.$$

This allows the sawtooth generator to be reset when $q-1$ periods are reached, skipping the last one. Thus, we will obtain a certain random period of the "LFSR + Sawtooth" system under consideration, slightly less than the maximum:
$${T}_{sub} = \sum\_{j=1}^{q-1} T_j < (p^m-1) \cdot q.$$

This will make it possible to combine several such subsystems into a common system with an increasing period. The GCD of these periods will most likely be small or even equal to one. If necessary, you can select the initial state of the sawtooth generator so as to achieve a unit GCD.

The $64$-bit generator output is generated after $4$ work cycles, accumulating $16$-bit states obtained after XORing the cells of all LFSR generators.

The total period of the proposed LFSR generator controlled by the sawtooth generator is $77...79$ bits: the period is a random value with relatively small variance.

The average byte-wise chi-square value is $256$ as expected.
