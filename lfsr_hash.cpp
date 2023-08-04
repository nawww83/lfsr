#include "lfsr_hash.hpp"


lfsr_hash::gens g;

lfsr8::u16 lfsr_hash::hash(const uint8_t* input, int n) {
    g.reset();
    constexpr salt s0 {2, 3, 4, 6};
    constexpr salt s1 {2, 1, 3, 7};
    g.salt(s0, 5);
    g.do_main_loop(input, n);
    g.salt(s1, 6);
    return g.get_hash();
}