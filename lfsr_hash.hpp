#include "lfsr.hpp"


namespace lfsr_hash {

using LFSR251x4 = lfsr8::LFSR<251, 4>;
using LFSR241x4 = lfsr8::LFSR<241, 4>;
using STATE = lfsr8::u32x8;

static constexpr STATE K1 = {1, 3, 1, 1};
static constexpr STATE K2 = {2, 1, 1, 1};
static constexpr STATE K3 = {1, 1, 4, 1};
static constexpr STATE K4 = {7, 1, 5, 1};

struct salt {
    lfsr8::u32 s0;
    lfsr8::u32 s1;
    lfsr8::u32 s2;
    lfsr8::u32 s3;
};


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
    void salt(salt s, int q) {
        for (int i=0; i<q; ++i) {
            g1_251x4.next(s.s0);
            g2_251x4.next(s.s1);
            g1_241x4.next(s.s2);
            g2_241x4.next(s.s3);
        }
    }
    void do_main_loop(const uint8_t* input, int n) {
        for (int i=0; i<n; ++i) {
            g1_251x4.next(input[i]);
            g2_251x4.next(input[n-i-1]);
            g1_241x4.next(input[i]);
            g2_241x4.next(input[n-i-1]);
        }
    }
    auto get_hash() {
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