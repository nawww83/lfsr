#include "lfsr.hpp"


namespace lfsr_hash {

using LFSR251x4 = lfsr8::LFSR_paired_2x4<251>;
using LFSR241x4 = lfsr8::LFSR_paired_2x4<241>;
using STATE = lfsr8::u16x8;

static constexpr STATE K1 = {7, 1, 6, 0, 4, 1, 3, 2};    // p=251
static constexpr STATE K2 = {13, 2, 5, 10, 7, 0, 10, 1}; // p=241

struct salt {
    lfsr8::u16 s0;
    lfsr8::u16 s1;
    lfsr8::u16 s2;
    lfsr8::u16 s3;
};

static constexpr salt S0 {2, 3, 4, 7};


struct gens {
    LFSR251x4 g_251x4;
    LFSR241x4 g_241x4;
public:
    constexpr gens(): g_251x4(K1),
        g_241x4(K2) {}
    void reset() {
        g_251x4.set_unit_state();
        g_241x4.set_unit_state();
    }
    void process_input(const uint8_t* input, int n) {
        for (int i=0; i<7; ++i) {
            g_251x4.next(S0.s0);
            g_241x4.next(S0.s1);
        }
        for (int i=0; i<n; ++i) {
            g_251x4.next(input[i]);
            g_241x4.next(input[n-i-1]);
        }
        for (int i=0; i<6; ++i) {
            g_251x4.next(S0.s2);
            g_241x4.next(S0.s3);
        }
    }
    auto form_hash() {
        lfsr8::u16 hash;
        auto st1 = g_251x4.get_state();
        auto st2 = g_241x4.get_state();
        hash = (st1[0] ^ st1[4] ^ st2[0] ^ st2[4]) ^ (st1[2] ^ st1[6] ^ st2[2] ^ st2[6]);
        hash <<= 8;
        hash |= (st1[1] ^ st1[5] ^ st2[1] ^ st2[5]) ^ (st1[3] ^ st1[7] ^ st2[3] ^ st2[7]);
        return hash;
    }
};

lfsr8::u16 hash(const uint8_t* input, int n);


}