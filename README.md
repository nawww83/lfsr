# Linear Feedback Shift Registers (LFSR) in cryptography
LFSR and its cryptographic applications.

Prime number $p$ and register length $m$ are used.

1. You can get 32, 64 and 128-bit **cryptographic** hashes with **salt**. Also, you can add 256-bit hash etc using similar extention approach as you need. It was estimated the proposed LFSR hash is comparable to SHA-512.

Additionally, src code contains some functions to find vectors of coefficients $K$ that provide periods $T_0 = {p}^{m} - 1$ and $T_1 = {p}^{m-1} - 1$.

2. You can generate 64-bit **cryptographic** random numbers based on two LFSR pairs (**total** length $4m$) with small primes $p$: 17, 19, and sawtooth modulation of that pairs with small periods: 7, 11.

## Build
g++ main.cpp lfsr_hash.cpp -std=c++20 -msse4.1 -O3 -o lfsr
## Run
./lfsr

## Hash performance
Approx. 210 MB/s @ Intel i7-8565U CPU 4.2GHz, GCC 11.4

## Random number generation performance
Approx. 120 MB/s @ Intel i7-8565U CPU 4.2GHz, GCC 11.4

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

We use two LFSR pairs of length $m=4$, and primes $p=251$ and $p=241$. Final states of two pairs are XOR-ed as usual. So, the total free period of given LFSR hash is about ${251}^{6} \cdot {241}^{6} \approx 95$ bits, i.e. ${2}^{95} \approx {10}^{28}$.

To improve the crypto resilience, the first LFSR pair is driven by the original input, but the second one is driven by reversed input vector. Also, the **salt** is added before/after input loop.

Bit scaling of LFSR hash is done by adding **salt**.

## Random number generation principles
We have two LFSR pairs, one has $p=17$, another has $p=19$. For example, lets consider the first pair. We can control an LFSR generator by some sawtooth generator:
$$i = (i + 1) \mod q.$$

Here $i$ - output of sawtooth generator with the period $q < p$. We start with some $i = i0$ - the initial sawtooth state. It is better to choose $q$ such the LFSR period $T = p^m - 1$ is not divisible by the sawtooth generator period $q$. In this case we will visit all possible $i$ when LFSR has some **fixed** non-zero state, and the total period will be maximal and equal to $(p^m-1)*q$.

Observations have shown that when we control LFSR by sawtooth generator we will have some sequence of LFSR periods $T[j]$, $j = [0..q-1]$ such that the sum of that periods is equal to the total period. The interesting property is the inner periods $T[i]$ are like some random numbers, so we can reset sawtooth generator after the fixed LFSR state was acheived $q-1$ times exactly, not $q$. In this case each sawtooth-controlled LFSR will have the period which is equal to the sum of $T[j]$, $j=[0..q-2]$. This period is slightly smaller than $(p^m-1)*q$. Each LFSR will have their own period; frequently, that periods have unit (or small) Greatest Common Divisor and, therefore, the total period will be equal to the multiplication of all periods. It is possible to choose the sawtooth initial states such the GCD will be equal to $1$.

At each step, the states of all LFSR generators are XORed and form $16$-bit output. After $4$ steps a $64$-bit number is formed as actual Random Generator output.
