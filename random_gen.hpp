#pragma once

#include "lfsr.hpp"

#include <utility>
#include <cmath>
#include <numeric>
#include <list>


namespace lfsr_rng {

constexpr int p = 23;

using LFSR_pair = lfsr8::LFSR_paired_2x4<p>;
using STATE = lfsr8::u16x8;
using u16 = lfsr8::u16;
using u32 = lfsr8::u32;
using u64 = lfsr8::u64;
using u128 = std::pair<lfsr8::u64, lfsr8::u64>;


static constexpr STATE K1 = {2, 2, 3, 7, 3, 2, 3, 9};    // p=23:   T0 = p^4 - 1
static constexpr STATE K2 = {3, 1, 4, 0, 3, 6, 10, 7};    // p=23:   T0 = p^4 - 1

class MyPrimes {
public:
    static const std::list<lfsr8::u16> primes() {return {7, 11, 13, 17, 19};}
};

struct gens {
    LFSR_pair gp1;
    LFSR_pair gp2;
    lfsr8::u16 ii11;  // Sawtooth initial state 1
    lfsr8::u16 ii12;  //   --//-- 2
    lfsr8::u16 ii21;  // Sawtooth initial state 2
    lfsr8::u16 ii22;  //   --//-- 2
    lfsr8::u16 qq;   // Sawtooth period
    lfsr8::u16 ii01; // Sawtooth initial state 1
    lfsr8::u16 ii02; //   --//-- 2
    lfsr8::u16 iie11; // Sawtooth end state 1
    lfsr8::u16 iie12; //   --//-- 2
    lfsr8::u16 iie21; // Sawtooth end state 2
    lfsr8::u16 iie22; //   --//-- 2
    STATE r1;
    STATE r2;
    bool is_finded;
    bool _dummy[3]{};
public:
    constexpr gens(): gp1(K1), gp2(K2),
        ii11(0), ii12(0), ii21(0), ii22(0), qq(0),
        ii01(0), ii02(0),
        iie11(0), iie12(0), iie21(0), iie22(0),
        r1({}), r2({}), is_finded(false)
        {}
    bool is_succes() const {
        return is_finded;
    }
    auto get_ii1() const {
        return ii01;
    }
    auto get_ii2() const {
        return ii02;
    }
    auto get_iie11() const {
        return iie11;
    }
    auto get_iie12() const {
        return iie12;
    }
    auto get_iie21() const {
        return iie21;
    }
    auto get_iie22() const {
        return iie22;
    }
    auto get_qq() const {
        return qq;
    }
    void seed(STATE st) {
        gp1.set_state(st);
        gp2.set_state(st);
        for (int i=0; i<12; ++i) { // a non-zero salt
            gp1.next(5);
            gp2.next(13);
        }
        const auto ref1 = gp1.get_state();
        const auto ref2 = gp2.get_state();
        auto test251 = [ref1, ref2, this](u16 i01, u16 i02, u16 q) {
            long long T1 = 0;
            long long T2 = 0;
            long long T3 = 0;
            long long T4 = 0;
            long long T = 1;
            const long long T0 = std::pow(p, 4) - 1;
            gp1.set_state(ref1);
            gp2.set_state(ref2);
            u16 i1 = i01;
            u16 i2 = i02;
            while (true) {
                gp1.next(i1); //  sawtooth modulation
                gp2.next(i2); //  to get random periods T0 < T < q*T0
                i1++; i1 %= q;
                i2++; i2 %= q;
                if (gp1.is_state_low(ref1)) {
                    T1 = T;
                    iie11 = i1;
//                        std::cout << " T1 detected: " << "i1 = " << i1 << ", T1 = " << T1 << std::endl;
                }
                if (gp1.is_state_high(ref1)) {
                    T2 = T;
                    iie12 = i1;
//                        std::cout << " T2 detected: " << "i1 = " << i1 << ", T2 = " << T2 << std::endl;
                }
                if (gp2.is_state_low(ref2)) {
                    T3 = T;
                    iie21 = i2;
//                        std::cout << " T3 detected: " << "i2 = " << i2 << ", T3 = " << T3 << std::endl;
                }
                if (gp2.is_state_high(ref2)) {
                    T4 = T;
                    iie22 = i2;
//                        std::cout << " T4 detected: " << "i2 = " << i2 << ", T4 = " << T4 << std::endl;
                }
                T += 1;
                if (T < q*T0) {
                    continue;
                }
                break;
            }
            const auto gcd1 = std::gcd(T1, T2);
            const auto gcd2 = std::gcd(T3, T4);
            const auto gcd = std::gcd(gcd1, gcd2);
            // std::cout << " GCD: " << gcd << std::endl;
            return ((gcd < 2) && (T1 > T0) && (T2 > T0) && (T3 > T0) && (T4 > T0));
        };
//        const long long T0 = std::pow(p, 4) - 1;
        // std::cout << " T0 " << T0 << std::endl;
        is_finded = false;
        for (auto qT : MyPrimes::primes()) {
            ii01 = 0; // disabled value
            ii02 = 0;
            qq = 0;
            std::cout << " test: qT: " << qT << std::endl;
            for (int i1=1; i1<qT; ++i1) {
                for (int i2=1; i2<qT; ++i2) {
                    is_finded = test251(i1, i2, qT);
                    if (is_finded) {
                        ii01 = i1;
                        ii02 = i2;
                        qq = qT;
                        break;
                    }
                }
                if (is_finded) {
                    break;
                }
            }
            if (is_finded) {
                break;
            }
        }
        ii11 = ii01;
        ii12 = ii01;
        ii21 = ii02;
        ii22 = ii02;
        r1 = ref1;
        r2 = ref2;
        gp1.set_state(ref1); // must: restore initial states
        gp2.set_state(ref2);
    }
    void next() {
        gp1.next(ii11, ii12); // must: the same operator as in the seed()
        gp2.next(ii21, ii22);
        ii11++; ii11 %= qq;
        ii12++; ii12 %= qq;
        ii21++; ii21 %= qq;
        ii22++; ii22 %= qq;
        bool f1 = gp1.is_state_low(r1);
        bool f2 = gp1.is_state_high(r1);
        bool f3 = gp2.is_state_low(r2);
        bool f4 = gp2.is_state_high(r2);
        ii11 = (! f1) ? ii11 : ((ii11 != iie11) ? ii11 : ii01);
        ii12 = (! f2) ? ii12 : ((ii12 != iie12) ? ii12 : ii01);
        ii21 = (! f3) ? ii21 : ((ii21 != iie21) ? ii21 : ii02);
        ii22 = (! f4) ? ii22 : ((ii22 != iie22) ? ii22 : ii02);
    }
    auto get_u32() const {
        auto st1 = gp1.get_state();
        auto st2 = gp2.get_state();
        lfsr8::u32 hash1; // form low 4 bits in byte: p=23: 5 bits => shift right by 1 bits to keep 4 bits
        hash1  = (((st1[0] ^ st2[4]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[1] ^ st2[5]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[2] ^ st2[6]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[3] ^ st2[7]) % p) >> 1) & 255;

        lfsr8::u32 hash2; // form high 4 bits in byte: p=23: 5 bits => shift left by 4 bits to keep 4 bits
        hash2  = (((st2[0] ^ st1[4]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[1] ^ st1[5]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[2] ^ st1[6]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[3] ^ st1[7]) % p) << 4) & 255;

        return hash1 ^ hash2;
    }
};


}
