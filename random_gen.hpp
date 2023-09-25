#pragma once

#include "lfsr.hpp"

#include <utility>
#include <cmath>
#include <numeric>
#include <list>

// #include <iostream>

namespace lfsr_rng {

constexpr int p = 23;

using LFSR_pair = lfsr8::LFSR_paired_2x4<p>;
using STATE = lfsr8::u16x8;
using u16 = lfsr8::u16;
using u32 = lfsr8::u32;
using u64 = lfsr8::u64;
using u128 = std::pair<lfsr8::u64, lfsr8::u64>;

// static constexpr STATE K1 = {7, 1, 2, 3, 12, 4, 14, 1};    // p=41:   T0 = p^4 - 1
// static constexpr STATE K2 = {7, 0, 3, 1, 11, 2, 3, 2};    // p=41:   T0 = p^4 - 1

static constexpr STATE K1 = {2, 2, 3, 7, 3, 2, 3, 9};    // p=23:   T0 = p^4 - 1
static constexpr STATE K2 = {3, 1, 4, 0, 3, 6, 10, 7};    // p=23:   T0 = p^4 - 1

// static constexpr STATE K = {3, 4, 3, 3, 3, 4, 2, 1};    // p=251:   T0 = p^4 - 1

// static constexpr STATE K = {3, 2, 2, 0, 7, 3, 3, 2};    // p=251:   T0 = p^4 - 1

static const std::list<lfsr8::u16> primes {7, 11, 13, 17};

struct gens {
    LFSR_pair gp1;
    LFSR_pair gp2;
    lfsr8::u16 ii1; // Sawtooth initial state 1
    lfsr8::u16 ii2; // --//-- 2
    lfsr8::u16 qq1; // Sawtooth period 1
    lfsr8::u16 qq2; // Sawtooth period 2
public:
    constexpr gens(): gp1(K1), gp2(K2) {}
    bool seed(STATE st) {
        gp1.set_state(st);
        gp2.set_state(st);
        for (int i=0; i<12; ++i) {
            gp1.next();
            gp2.next();
        }
        const auto ref1 = gp1.get_state();
        const auto ref2 = gp2.get_state();
        long long T1;
        long long T2;
        long long T3;
        long long T4;

        auto test251 = [&T1, &T2, &T3, &T4, ref1, ref2, this](u16 i01, u16 q1, u16 i02, u16 q2) {
            long long T = 1;
            const long long T0 = std::pow(p, 4) - 1;
            gp1.set_state(ref1);
            gp2.set_state(ref2);
            u16 i1 = i01;
            u16 i2 = i02;
            bool f1 = false;
            bool f2 = false;
            bool f3 = false;
            bool f4 = false;
            while (true) {
                gp1.next(i1^i2); // Joint sawtooth modulation
                gp2.next(i1^i2); // cross scheme to get random periods T0 < T < p*T0
                i1++; i1 %= q1;
                i2++; i2 %= q2;
                const auto t1 = gp1.get_state();
                const auto t2 = gp2.get_state();
                if (! f1) {
                    if ((t1[0] == ref1[0]) && (t1[1] == ref1[1]) && (t1[2] == ref1[2]) && (t1[3] == ref1[3]) && (i1 = i01) && (i2 = i02)) {
                        T1 = T;
                        f1 = true;
                    }
                }
                if (! f2) {
                    if ((t1[4] == ref1[4]) && (t1[5] == ref1[5]) && (t1[6] == ref1[6]) && (t1[7] == ref1[7]) && (i2 == i02) && (i1 = i01)) {
                        T2 = T;
                        f2 = true;
                    }
                }
                if (! f3) {
                    if ((t2[0] == ref2[0]) && (t2[1] == ref2[1]) && (t2[2] == ref2[2]) && (t2[3] == ref2[3]) && (i1 = i01) && (i2 = i02)) {
                        T3 = T;
                        f3 = true;
                    }
                }
                if (! f4) {
                    if ((t2[4] == ref2[4]) && (t2[5] == ref2[5]) && (t2[6] == ref2[6]) && (t2[7] == ref2[7]) && (i2 == i02) && (i1 = i01)) {
                        T4 = T;
                        f4 = true;
                    }
                }
                if (f1 && f2 && f3 && f4) {
                    break;
                }
                T += 1;
                if (2*T <= p*T0) { // up to p*T/2
                    continue;
                }
                T = 0; // Abnormal period: unachievable state
                break;
            }
            return T;
        };
        const long long T0 = std::pow(p, 4) - 1;
        // std::cout << " T0 " << T0 << std::endl;
        bool is_finded = false;
        for (auto qT1 : primes) {
            for (auto qT2 : primes) {
                if (qT1 == qT2) {
                    continue;
                }
                ii1 = 0; // disabled value
                ii2 = 0;
                qq1 = 0;
                qq2 = 0;
                int idx1 = -1;
                int idx2 = -1;
                for (int i1=1; i1<qT1; ++i1) {
                    for (int i2=1; i2<qT2; ++i2) {
                        T1 = 0; T2 = 0; T3 = 0; T4 = 0;
                        auto T = test251(i1, qT1, i2, qT2);
                        // std::cout << " test: i1: " << i1 << ", i2: " << i2 << ", T1: " << T1 << ", T2: " << T2 << ", T: " << T << ", GCD: " << std::gcd(T1, T2) << std::endl;
                        // std::cout << " test: i1: " << i1 << ", i2: " << i2 << ", T3: " << T3 << ", T4: " << T4 << ", T: " << T << ", GCD: " << std::gcd(T3, T4) << std::endl;
                        const auto gcd1 = std::gcd(T1, T2);
                        const auto gcd2 = std::gcd(T3, T4);
                        const auto gcd = std::gcd(gcd1, gcd2);
                        is_finded = ((T > 0) && (gcd < 3) && (T1 > T0) && (T2 > T0) && (T3 > T0) && (T4 > T0));
                        if (is_finded) {
                            idx2 = i2;
                            ii2 = i2;
                            qq2 = qT2;
                            break;
                        }
                    }
                    if (is_finded) {
                        idx1 = i1;
                        ii1 = i1;
                        qq1 = qT1;
                        break;
                    }
                }
                if (is_finded) {
                    // const auto gcd1 = std::gcd(T1, T2);
                    // const auto gcd2 = std::gcd(T3, T4);
                    // const auto gcd = std::gcd(gcd1, gcd2);
                    // std::cout << "FOUND: " << "T1 = " << T1 << ", T2 = " << T2 << ", T3 = " << T3 << ", T4 = " << T4 << ", idx1 = " << idx1 << ", idx2 = " << idx2 << 
                        // ", qq1: " << qq1 << ", qq2: " << qq2 << std::endl;
                    // double Tb = std::log2(T1) + std::log2(T2) + std::log2(T3) + std::log2(T4);
                    // std::cout << "Total period: " << std::floor(Tb*100.)*0.01 << " bits. GCD: " << gcd << std::endl;
                    break;
                }
                if (is_finded) {
                    break;
                }
            }
            if (is_finded) {
                break;
            }
        }
        // if (! is_finded) {
            // std::cout << "Not found!" << std::endl;
        // }
        gp1.set_state(ref1); // must: restore initial states
        gp2.set_state(ref2);
        return is_finded;
    }
    void next() {
        gp1.next(ii1^ii2); // must: the same operator as in the seed()
        gp2.next(ii1^ii2);
        ii1++; ii1 %= qq1;
        ii2++; ii2 %= qq2;
    }
    auto get_u32() {
        auto st1 = gp1.get_state();
        auto st2 = gp2.get_state();
        lfsr8::u32 hash1; // form low 4 bits in byte: p=23: 5 bits => shift right by 1 bits to keep 4 bits
        hash1  = (((st1[0] ^ st1[4]) % p) >> 1);
        hash1 <<= 8;
        hash1 |= (((st1[1] ^ st1[5]) % p) >> 1);
        hash1 <<= 8;
        hash1 |= (((st1[2] ^ st1[6]) % p) >> 1);
        hash1 <<= 8;
        hash1 |= (((st1[3] ^ st1[7]) % p) >> 1);

        lfsr8::u32 hash2; // form high 4 bits in byte: p=23: 5 bits => shift left by 4 bits to keep 4 bits
        hash2  = (((st2[0] ^ st2[4]) % p) << 4);
        hash2 <<= 8;
        hash2 |= (((st2[1] ^ st2[5]) % p) << 4);
        hash2 <<= 8;
        hash2 |= (((st2[2] ^ st2[6]) % p) << 4);
        hash2 <<= 8;
        hash2 |= (((st2[3] ^ st2[7]) % p) << 4);
        
        return hash1 ^ hash2;
    }
};


}