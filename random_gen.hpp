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


static constexpr STATE K1 = {2, 2, 3, 7, 3, 2, 3, 9};     // p=23:   T0 = p^4 - 1
static constexpr STATE K2 = {3, 1, 4, 0, 3, 6, 10, 7};    // p=23:   T0 = p^4 - 1

class MyPrimes {
public:
    static const std::list<lfsr8::u16> primes() {return {7, 11, 13, 17, 19};} // less than p
};

struct gens {
    LFSR_pair gp1;
    LFSR_pair gp2;
    std::array<long long, 8> T;
    lfsr8::u16 ii11;  // Sawtooth initial states
    lfsr8::u16 ii12;  //
    lfsr8::u16 ii21;  //
    lfsr8::u16 ii22;  //
    lfsr8::u16 qq;    // Sawtooth period
    lfsr8::u16 ii01;  //
    lfsr8::u16 ii02;  //
    lfsr8::u16 _dummy1;
    int is_finded;
public:
    constexpr gens(): gp1(K1), gp2(K2), T({}),
        ii11(0), ii12(0), ii21(0), ii22(0), qq(0),
        ii01(0), ii02(0), _dummy1(0),
        is_finded(0)
        {}
    bool is_succes() const {
        return (is_finded != 0);
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
        auto test251 = [ref1, ref2, this](u16 i01, u16 i02, u16 q, u16 dummy=0) -> int {
            T = {1, 1, 1, 1, 0, 0, 0, 0};
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
                T[4] = (! gp1.is_state_low(ref1)) ? T[4] : T[0];
                T[5] = (! gp1.is_state_high(ref1)) ? T[5] : T[1];
                T[6] = (! gp2.is_state_low(ref2)) ? T[6] : T[2];
                T[7] = (! gp2.is_state_high(ref2)) ? T[7] : T[3];
                T[0]++; T[1]++; T[2]++; T[3]++;
                if (T[0] < q*T0) {
                    continue;
                }
                break;
            }
            auto gcd = std::gcd(std::gcd(T[4], T[5]), std::gcd(T[6], T[7]));
            return ((gcd < 2) && (T[4] > T0) && (T[5] > T0) && (T[6] > T0) && (T[7] > T0)) ? 1 : 0;
        };
        is_finded = 0;
        for (auto qT : MyPrimes::primes()) {
            ii01 = 0;
            ii02 = 0;
            qq = 0;
            for (int i1=1; i1<qT; ++i1) {
                for (int i2=1; i2<qT; ++i2) {
                    is_finded = test251(i1, i2, qT);
                    if (is_finded != 0) {
                        ii01 = i1;
                        ii02 = i2;
                        qq = qT;
                        break;
                    }
                }
                if (is_finded != 0) {
                    break;
                }
            }
            if (is_finded != 0) {
                break;
            }
        }
        ii11 = ii01;
        ii12 = ii01;
        ii21 = ii02;
        ii22 = ii02;
        T[0] = 1; // reset counters
        T[1] = 1;
        T[2] = 1;
        T[3] = 1;
        gp1.set_state(ref1); // must: restore initial states
        gp2.set_state(ref2);
    }
    void next() {
        gp1.next(ii11, ii12); // must: the same operator as in the seed()
        gp2.next(ii21, ii22);
        ii11++; ii11 %= qq;  // Sawtooth
        ii12++; ii12 %= qq;
        ii21++; ii21 %= qq;
        ii22++; ii22 %= qq;
        ii11 = (T[0] != T[4]) ? ii11 : ii01;
        ii12 = (T[1] != T[5]) ? ii12 : ii01;
        ii21 = (T[2] != T[6]) ? ii21 : ii02;
        ii22 = (T[3] != T[7]) ? ii22 : ii02;
        T[0] = (T[0] != T[4]) ? T[0] : 0; // reset counters
        T[1] = (T[1] != T[5]) ? T[1] : 0;
        T[2] = (T[2] != T[6]) ? T[2] : 0;
        T[3] = (T[3] != T[7]) ? T[3] : 0;
        T[0]++; T[1]++; T[2]++; T[3]++;
    }
    lfsr8::u32 get_u32() const {
        return get_u32_and_or();
    }
    lfsr8::u32 get_u32_xor() const {
        auto st1 = gp1.get_state();
        auto st2 = gp2.get_state();
        lfsr8::u32 hash1; // form low 4 bits in byte: p=23: 5 bits => shift right by 1 bits to keep 4 bits
        hash1  = (((st1[0] ^ st2[4]) % p) >> 1) & 255;   // XOR + modulo p: more secure
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
    lfsr8::u32 get_u32_and_or() const { // more secure than XOR
        auto st1 = gp1.get_state();
        auto st2 = gp2.get_state();
        lfsr8::u32 hash1;
        hash1  = (((st1[0] & st2[4]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[1] & st2[5]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[2] & st2[6]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[3] & st2[7]) % p) >> 1) & 255;

        lfsr8::u32 hash2;
        hash2  = (((st2[0] | st1[4]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[1] | st1[5]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[2] | st1[6]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[3] | st1[7]) % p) << 4) & 255;

        lfsr8::u32 hash3;
        hash3  = (((st1[0] | st1[4]) % p) >> 1) & 255;
        hash3 <<= 8;
        hash3 |= (((st1[1] | st1[5]) % p) >> 1) & 255;
        hash3 <<= 8;
        hash3 |= (((st1[2] | st1[6]) % p) >> 1) & 255;
        hash3 <<= 8;
        hash3 |= (((st1[3] | st1[7]) % p) >> 1) & 255;

        lfsr8::u32 hash4;
        hash4  = (((st2[0] & st2[4]) % p) << 4) & 255;
        hash4 <<= 8;
        hash4 |= (((st2[1] & st2[5]) % p) << 4) & 255;
        hash4 <<= 8;
        hash4 |= (((st2[2] & st2[6]) % p) << 4) & 255;
        hash4 <<= 8;
        hash4 |= (((st2[3] & st2[7]) % p) << 4) & 255;

        return hash1 ^ hash2 ^ hash3 ^ hash4;
    }
     lfsr8::u64 get_u64() const { //
        auto st1 = gp1.get_state();
        auto st2 = gp2.get_state();
        lfsr8::u32 hash1;
        hash1  = (((st1[0] & st2[4]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[1] & st2[5]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[2] & st2[6]) % p) >> 1) & 255;
        hash1 <<= 8;
        hash1 |= (((st1[3] & st2[7]) % p) >> 1) & 255;

        lfsr8::u32 hash2;
        hash2  = (((st2[0] | st1[4]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[1] | st1[5]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[2] | st1[6]) % p) << 4) & 255;
        hash2 <<= 8;
        hash2 |= (((st2[3] | st1[7]) % p) << 4) & 255;

        lfsr8::u32 hash3;
        hash3  = (((st1[0] | st1[4]) % p) >> 1) & 255;
        hash3 <<= 8;
        hash3 |= (((st1[1] | st1[5]) % p) >> 1) & 255;
        hash3 <<= 8;
        hash3 |= (((st1[2] | st1[6]) % p) >> 1) & 255;
        hash3 <<= 8;
        hash3 |= (((st1[3] | st1[7]) % p) >> 1) & 255;

        lfsr8::u32 hash4;
        hash4  = (((st2[0] & st2[4]) % p) << 4) & 255;
        hash4 <<= 8;
        hash4 |= (((st2[1] & st2[5]) % p) << 4) & 255;
        hash4 <<= 8;
        hash4 |= (((st2[2] & st2[6]) % p) << 4) & 255;
        hash4 <<= 8;
        hash4 |= (((st2[3] & st2[7]) % p) << 4) & 255;

        //

        lfsr8::u32 hash5;
        hash5  = (((st1[0] & st1[4]) % p) >> 1) & 255;
        hash5 <<= 8;
        hash5 |= (((st1[1] & st1[5]) % p) >> 1) & 255;
        hash5 <<= 8;
        hash5 |= (((st1[2] & st1[6]) % p) >> 1) & 255;
        hash5 <<= 8;
        hash5 |= (((st1[3] & st1[7]) % p) >> 1) & 255;

        lfsr8::u32 hash6;
        hash6  = (((st2[0] | st2[4]) % p) << 4) & 255;
        hash6 <<= 8;
        hash6 |= (((st2[1] | st2[5]) % p) << 4) & 255;
        hash6 <<= 8;
        hash6 |= (((st2[2] | st2[6]) % p) << 4) & 255;
        hash6 <<= 8;
        hash6 |= (((st2[3] | st2[7]) % p) << 4) & 255;

        lfsr8::u32 hash7;
        hash7  = (((st1[0] | st2[4]) % p) >> 1) & 255;
        hash7 <<= 8;
        hash7 |= (((st1[1] | st2[5]) % p) >> 1) & 255;
        hash7 <<= 8;
        hash7 |= (((st1[2] | st2[6]) % p) >> 1) & 255;
        hash7 <<= 8;
        hash7 |= (((st1[3] | st2[7]) % p) >> 1) & 255;

        lfsr8::u32 hash8;
        hash8  = (((st2[0] & st1[4]) % p) << 4) & 255;
        hash8 <<= 8;
        hash8 |= (((st2[1] & st1[5]) % p) << 4) & 255;
        hash8 <<= 8;
        hash8 |= (((st2[2] & st1[6]) % p) << 4) & 255;
        hash8 <<= 8;
        hash8 |= (((st2[3] & st1[7]) % p) << 4) & 255;

        hash1 = hash1 ^ hash2 ^ hash3 ^ hash4;
        hash5 = hash5 ^ hash6 ^ hash7 ^ hash8;
        return (lfsr8::u64(hash1) << 32) | lfsr8::u64(hash5);
    }
};


}

