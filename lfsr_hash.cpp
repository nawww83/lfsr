#include "lfsr_hash.hpp"


lfsr_hash::gens g;

lfsr_hash::u16 lfsr_hash::hash16(const uint8_t* input, int n) {
    g.reset();
    g.add_salt(S0);
    g.process_input(input, n);
    g.add_salt(S1);
    return g.form_hash16();
}

lfsr_hash::u32 lfsr_hash::hash32(const uint8_t* input, int n) {
    g.reset();
    g.add_salt(S1);
    g.process_input(input, n);
    g.add_salt(S0);
    u32 h1 = g.form_hash16();
    g.add_salt(S2);
    u32 h2 = g.form_hash16();
    return h1 | (h2 << 16);
}