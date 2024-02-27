#pragma once


#include "lfsr.hpp"

#include <utility>
#include <cmath>
#include <numeric>
#include <array>


namespace lfsr_rng {

using u16 = lfsr8::u16;
using u32 = lfsr8::u32;
using u64 = lfsr8::u64;
using u128 = std::pair<lfsr8::u64, lfsr8::u64>;

static constexpr u16 p1 = 19;
static constexpr u16 p2 = 17;

using LFSR_pair_1 = lfsr8::LFSR_paired_2x4<p1>;
using LFSR_pair_2 = lfsr8::LFSR_paired_2x4<p2>;
using STATE = lfsr8::u16x8;

static constexpr STATE K1 = {9, 5, 2, 0, 4, 2, 2, 6};    // p=19
static constexpr STATE K2 = {3, 4, 2, 1, 6, 1, 2, 1};    // p=17
static constexpr std::array<u16, 2> primes {7, 11}; // Sawtooth periods: such that T = p^4 - 1 is not divisible by the primes

// class MyPrimes {
// public:
//     static const std::list<lfsr8::u16> primes() {return {7, 11};}
// };

struct gens {
    LFSR_pair_1 gp1;
    LFSR_pair_2 gp2;
    std::array<u32, 8> T;
    u16 ii11;  // Sawtooth initial states
    u16 ii12;  //
    u16 ii21;  //
    u16 ii22;  //
    u16 ii01;  // initial states
    u16 ii02;  //
    int is_finded;
public:
    constexpr gens(): gp1(K1), gp2(K2), T({}),
        ii11(0), ii12(0), ii21(0), ii22(0),
        ii01(0), ii02(0),
        is_finded(0)
        {}
    bool is_succes() const {
        return (is_finded != 0);
    }
    void seed(STATE st) {
        STATE tmp1 = st;
        STATE tmp2 = st;
        for (auto& el : tmp1) {
            el %= 16;
        }
        for (auto& el : tmp2) {
            el >>= 4;
        }
        gp1.set_state(tmp1);
        gp2.set_state(tmp2);
        u16 h1 = 1;
        for (const auto& el : tmp1) {
            h1 ^= el;
        }
        u16 h2 = 2;
        for (const auto& el : tmp2) {
            h2 ^= el;
        }
        u16 i1 = h1;
        u16 i2 = h2;
        for (int i=0; i<int(primes[0])*int(primes[1]); ++i) { // saturate LFSRs
            gp1.next(i1);
            gp2.next(i2);
            i1++; i2++;
            i1 %= primes[0];
            i2 %= primes[1];
        }
        const auto ref1 = gp1.get_state();
        const auto ref2 = gp2.get_state();
        auto test251 = [ref1, ref2, this](u16 i01, u16 i02) -> int {
            T = {1, 1, 1, 1, 0, 0, 0, 0};
            const u32 T01 = std::pow((long)p1, 4) - 1;
            const u32 T02 = std::pow((long)p2, 4) - 1;
            gp1.set_state(ref1);
            gp2.set_state(ref2);
            u16 i1 = i01;
            u16 i2 = i02;
            while (true) {
                gp1.next(i1); //  Sawtooth modulation
                gp2.next(i2); //  to get random periods T[i] such that sum of T[i] is equal to q*T0, where the q - Sawtooth period
                i1++; i2++;   //  The remainder mod(p^4 - 1 , q) is not zero => we achieve all indexes i in [0, q) when LFSR is in the same reference state
                i1 %= primes[0];
                i2 %= primes[1]; // We visit almost all i ecxept one => we set the restriction T < q*T0 below
                T[4] = (! gp1.is_state_low(ref1)) ? T[4] : ((T[0] < primes[0]*T01) ? T[0] : T[4]);
                T[5] = (! gp1.is_state_high(ref1)) ? T[5] : ((T[1] < primes[0]*T01) ? T[1] : T[5]);
                T[6] = (! gp2.is_state_low(ref2)) ? T[6] : ((T[2] < primes[1]*T02) ? T[2] : T[6]);
                T[7] = (! gp2.is_state_high(ref2)) ? T[7] : ((T[3] < primes[1]*T02) ? T[3] : T[7]);
                T[0]++; T[1]++;
                T[2]++; T[3]++;
                // All counters are out of the range
                if ((T[0] >= primes[0]*T01) && (T[1] >= primes[0]*T01) && (T[2] >= primes[1]*T02) && (T[3] >= primes[1]*T02)) {
                    break;
                }
            }
            auto gcd = std::gcd(std::gcd(T[4], T[5]), std::gcd(T[6], T[7]));
            return ((gcd < 2) && (T[4] > T01) && (T[5] > T01) && (T[6] > T02) && (T[7] > T02)) ? 1 : 0;
        };
        is_finded = 0;
        ii01 = 0;
        ii02 = 0;
        for (u16 i1=1; i1<primes[0]; ++i1) {
            for (u16 i2=1; i2<primes[1]; ++i2) {
                is_finded = test251(i1, i2);
                if (is_finded != 0) {
                    ii01 = i1;
                    ii02 = i2;
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
    u64 next_u64() {
        u64 x = 0;
        STATE st1{}; STATE st2{};
        //
        gp1.next(ii11, ii12); // must: the same operator as in the seed()
        gp2.next(ii21, ii22);
        ii11++; ii12++; ii21++; ii22++;
        ii11 %= primes[0];  // Sawtooth
        ii12 %= primes[0];
        ii21 %= primes[1];
        ii22 %= primes[1];
        ii11 = (T[0] != T[4]) ? ii11 : ii01;
        ii12 = (T[1] != T[5]) ? ii12 : ii01;
        ii21 = (T[2] != T[6]) ? ii21 : ii02;
        ii22 = (T[3] != T[7]) ? ii22 : ii02;
        T[0] = (T[0] != T[4]) ? T[0] : 0; // reset counters
        T[1] = (T[1] != T[5]) ? T[1] : 0;
        T[2] = (T[2] != T[6]) ? T[2] : 0;
        T[3] = (T[3] != T[7]) ? T[3] : 0;
        T[0]++; T[1]++; T[2]++; T[3]++;
        //
        st1 = gp1.get_state();
        st2 = gp2.get_state();
        st1[0] ^= (st1[4] ^ st2[0] ^ st2[4]);
        st1[1] ^= (st1[5] ^ st2[1] ^ st2[5]);
        st1[2] ^= (st1[6] ^ st2[2] ^ st2[6]);
        st1[3] ^= (st1[7] ^ st2[3] ^ st2[7]);
        st1[0] &= 15;
        st1[1] &= 15;
        st1[2] &= 15;
        st1[3] &= 15;
        x <<= 4;
        x |= st1[0];
        x <<= 4;
        x |= st1[1];
        x <<= 4;
        x |= st1[2];
        x <<= 4;
        x |= st1[3];
        //
        gp1.next(ii11, ii12); // must: the same operator as in the seed()
        gp2.next(ii21, ii22);
        ii11++; ii12++; ii21++; ii22++;
        ii11 %= primes[0];  // Sawtooth
        ii12 %= primes[0];
        ii21 %= primes[1];
        ii22 %= primes[1];
        ii11 = (T[0] != T[4]) ? ii11 : ii01;
        ii12 = (T[1] != T[5]) ? ii12 : ii01;
        ii21 = (T[2] != T[6]) ? ii21 : ii02;
        ii22 = (T[3] != T[7]) ? ii22 : ii02;
        T[0] = (T[0] != T[4]) ? T[0] : 0; // reset counters
        T[1] = (T[1] != T[5]) ? T[1] : 0;
        T[2] = (T[2] != T[6]) ? T[2] : 0;
        T[3] = (T[3] != T[7]) ? T[3] : 0;
        T[0]++; T[1]++; T[2]++; T[3]++;
        //
        st1 = gp1.get_state();
        st2 = gp2.get_state();
        st1[0] ^= (st1[4] ^ st2[0] ^ st2[4]);
        st1[1] ^= (st1[5] ^ st2[1] ^ st2[5]);
        st1[2] ^= (st1[6] ^ st2[2] ^ st2[6]);
        st1[3] ^= (st1[7] ^ st2[3] ^ st2[7]);
        st1[0] &= 15;
        st1[1] &= 15;
        st1[2] &= 15;
        st1[3] &= 15;
        x <<= 4;
        x |= st1[0];
        x <<= 4;
        x |= st1[1];
        x <<= 4;
        x |= st1[2];
        x <<= 4;
        x |= st1[3];
        //
        gp1.next(ii11, ii12); // must: the same operator as in the seed()
        gp2.next(ii21, ii22);
        ii11++; ii12++; ii21++; ii22++;
        ii11 %= primes[0];  // Sawtooth
        ii12 %= primes[0];
        ii21 %= primes[1];
        ii22 %= primes[1];
        ii11 = (T[0] != T[4]) ? ii11 : ii01;
        ii12 = (T[1] != T[5]) ? ii12 : ii01;
        ii21 = (T[2] != T[6]) ? ii21 : ii02;
        ii22 = (T[3] != T[7]) ? ii22 : ii02;
        T[0] = (T[0] != T[4]) ? T[0] : 0; // reset counters
        T[1] = (T[1] != T[5]) ? T[1] : 0;
        T[2] = (T[2] != T[6]) ? T[2] : 0;
        T[3] = (T[3] != T[7]) ? T[3] : 0;
        T[0]++; T[1]++; T[2]++; T[3]++;
        //
        st1 = gp1.get_state();
        st2 = gp2.get_state();
        st1[0] ^= (st1[4] ^ st2[0] ^ st2[4]);
        st1[1] ^= (st1[5] ^ st2[1] ^ st2[5]);
        st1[2] ^= (st1[6] ^ st2[2] ^ st2[6]);
        st1[3] ^= (st1[7] ^ st2[3] ^ st2[7]);
        st1[0] &= 15;
        st1[1] &= 15;
        st1[2] &= 15;
        st1[3] &= 15;
        x <<= 4;
        x |= st1[0];
        x <<= 4;
        x |= st1[1];
        x <<= 4;
        x |= st1[2];
        x <<= 4;
        x |= st1[3];
        //
        gp1.next(ii11, ii12); // must: the same operator as in the seed()
        gp2.next(ii21, ii22);
        ii11++; ii12++; ii21++; ii22++;
        ii11 %= primes[0];  // Sawtooth
        ii12 %= primes[0];
        ii21 %= primes[1];
        ii22 %= primes[1];
        ii11 = (T[0] != T[4]) ? ii11 : ii01;
        ii12 = (T[1] != T[5]) ? ii12 : ii01;
        ii21 = (T[2] != T[6]) ? ii21 : ii02;
        ii22 = (T[3] != T[7]) ? ii22 : ii02;
        T[0] = (T[0] != T[4]) ? T[0] : 0; // reset counters
        T[1] = (T[1] != T[5]) ? T[1] : 0;
        T[2] = (T[2] != T[6]) ? T[2] : 0;
        T[3] = (T[3] != T[7]) ? T[3] : 0;
        T[0]++; T[1]++; T[2]++; T[3]++;
        //
        st1 = gp1.get_state();
        st2 = gp2.get_state();
        st1[0] ^= (st1[4] ^ st2[0] ^ st2[4]);
        st1[1] ^= (st1[5] ^ st2[1] ^ st2[5]);
        st1[2] ^= (st1[6] ^ st2[2] ^ st2[6]);
        st1[3] ^= (st1[7] ^ st2[3] ^ st2[7]);
        st1[0] &= 15;
        st1[1] &= 15;
        st1[2] &= 15;
        st1[3] &= 15;
        x <<= 4;
        x |= st1[0];
        x <<= 4;
        x |= st1[1];
        x <<= 4;
        x |= st1[2];
        x <<= 4;
        x |= st1[3];
        //
        return x;
    }
};


}

