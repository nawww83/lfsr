# Linear Feedback Shift Registers (LFSR) in cryptography
LFSR and its cryptographic applications.

Prime number $p$ and register length $m$ are used.

You can calculate 16, 32, 64 and 128-bit **cryptographic** hashes.
Additionally, src code contains some functions to find vectors of coefficients $K$ that provide periods $T_0 = {p}^{m} - 1$ and $T_1 = {p}^{m-1} - 1$.

## Build
g++ main.cpp lfsr_hash.cpp -std=c++20 -msse4.1 -O3 -o lfsr
## Run
./lfsr

## Hash performance
Approx. 210 MB/s @ Intel i7-8565U CPU 4.2GHz, GCC 11.4

## LFSR hash principles
Base unit for LFSR hash is LFSR register which can be interpreted as Finite State Machine (Moore machine if be more precise).

Classicaly, LFSR operates with binary symbols, 0 and 1, but nobody forbid us to use any M-ary symbol - integer numbers from 0 to M-1.

Due to special mathematical property of prime numbers, it makes sense to use $M$ as a prime number: $2, 3, 5, 7, \ldots{}$ That property is the maximal period property when we calculate all integer powers by modulo $p$: ${a}^{i} \mod p, \quad i \in \mathbb{Z}$.

For example, let's $p=11$, $a=2$, $i=0...10$, so we have the following powers $1, 2, 4, 8, 5, 10, 9, 7, 3, 6, 1$.

We can see that all numbers are different (except the last), from $1$ to $10$, and ${a}^{0} = {a}^{p-1} \mod p = 1$ or equivalently $a^p - a = 0$.
The last expression is true for all $a$, but not all $a$ give us unique powers. The maximal number of unique powers is $p-1$. Some $a$ give us maximal period $T = p - 1$, other ones give us smaller periods, for example, $T/2$ or smaller.

The previous example is a $p$-ary LFSR register wtih unit length, $m=1$, but nobody forbid us to use any length $m$. For that case maximal period is $T = {p}^{m} - 1$, and we should find $m$ integer numbers $({a}_{0},..., {a}_{m-1})$, which provide maximal period. The $m$ numbers is the coefficients of the generator polynomial $g(x)$ with degree $m$.

In general, we use $v$ as the last element of LFSR state $\vec s$, then multiply generator vector $\vec a$ by the scalar $v$, then add the latter to the left shifted by 1 element LFSR state (with zero insertion). Left shift can be interpreted as multiplication by modulo $p$ (the direct analogy is the binary left shift and multiplication by $2$ equivalence), so theoretically we can express LFSR loop in vector form
$$\begin{matrix}
\vec s =& (v \vec a + p\vec s) \mod p \newline
v =& {s}_{m-1}
\end{matrix}$$

High symbols of $\vec s$, i.e. $s_m$, ${s}_{m+1}, ...$ can be ignored.

As a rule, LFSR is initialized by unit state $\vec s = (1, 0, 0, ..., 0)$. Having $m$ cycles, LFSR will be in the state $\vec s = \vec a$ exactly, which is the same as generator coefficients. So, we consider that LFSR is **saturated** at that moment. When we continue cycles, we will observe some different states, and at some moment the current state will be equal to the initial state $\vec a$. What cycles we done between two equal states will determine LFSR period $T$.

To be continued...