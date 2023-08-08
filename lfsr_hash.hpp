#include "lfsr.hpp"


namespace lfsr_hash {

using LFSR251x4 = lfsr8::LFSR<251, 4>;
using LFSR241x4 = lfsr8::LFSR<241, 4>;
using STATE = lfsr8::u16x8;

static constexpr STATE K1 = {3, 0, 0, 5};
static constexpr STATE K2 = {3, 7, 1, 7};
static constexpr STATE K3 = {7, 10, 1, 4};
static constexpr STATE K4 = {6, 1, 4, 2};

struct salt {
    lfsr8::u16 s0;
    lfsr8::u16 s1;
    lfsr8::u16 s2;
    lfsr8::u16 s3;
};

static constexpr salt S0 {2, 3, 4, 6};
static constexpr salt S1 {2, 1, 3, 7};


struct gens {
    LFSR251x4 g1_251x4;
    LFSR251x4 g2_251x4;
    LFSR241x4 g1_241x4;
    LFSR241x4 g2_241x4;
public:
    constexpr gens(): g1_251x4(K1),
        g2_251x4(K2),
        g1_241x4(K3),
        g2_241x4(K4) {}
    void reset() {
        g1_251x4.set_unit_state();
        g2_251x4.set_unit_state();
        g1_241x4.set_unit_state();
        g2_241x4.set_unit_state();
    }
    void process_input(const uint8_t* input, int n) {
        for (int i=0; i<5; ++i) {
            g1_251x4.next(S0.s0);
            g2_251x4.next(S0.s1);
            g1_241x4.next(S0.s2);
            g2_241x4.next(S0.s3);
        }
        for (int i=0; i<n; ++i) {
            g1_251x4.next(input[i]);
            g2_251x4.next(input[i]);
            g1_241x4.next(input[n-i-1]);
            g2_241x4.next(input[n-i-1]);
        }
        for (int i=0; i<6; ++i) {
            g1_251x4.next(S1.s0);
            g2_251x4.next(S1.s1);
            g1_241x4.next(S1.s2);
            g2_241x4.next(S1.s3);
        }
    }
    auto form_hash() {
        lfsr8::u16 hash;
        auto st11 = g1_251x4.get_state();
        auto st21 = g2_251x4.get_state();
        auto st12 = g1_241x4.get_state();
        auto st22 = g2_241x4.get_state();
        hash = st11[0] ^ st12[0] ^ st11[2] ^ st12[2];
        hash <<= 8;
        hash |= (st11[1] ^ st12[1] ^ st11[3] ^ st12[3]);
        return hash;
    }
};

lfsr8::u16 hash(const uint8_t* input, int n);


}