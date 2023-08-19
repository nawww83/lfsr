#include "lfsr_hash.hpp"


lfsr_hash::gens g;

// lfsr_hash::u16 lfsr_hash::hash16(const uint8_t* input, int n) {
//     g.reset();
//     g.add_salt(S0);
//     g.process_input(input, n);
//     g.add_salt(S1);
//     return g.form_hash16();
// }

lfsr_hash::u32 lfsr_hash::hash32(const uint8_t* input, int n) {
    g.reset();
    g.add_salt(S1);
    g.process_input(input, n);
    g.add_salt(S0);
    u32 h = g.form_hash32();
    return h;
}

lfsr_hash::u64 lfsr_hash::hash64(const uint8_t* input, int n) {
    g.reset();
    g.add_salt(S0);
    g.process_input(input, n);
    g.add_salt(S1);
    u64 h1 = g.form_hash32();
    g.add_salt(S2);
    u64 h2 = g.form_hash32();
    return (h1 << 32) | h2;
}

lfsr_hash::u128 lfsr_hash::hash128(const uint8_t* input, int n) {
    g.reset();
    g.add_salt(S1);
    g.add_salt(S0);
    g.process_input(input, n);
    g.add_salt(S0);
    g.add_salt(S1);

    u64 h1 = g.form_hash32();
    g.add_salt(S3);
    u64 h2 = g.form_hash32();
    g.add_salt(S2);
    u64 h3 = g.form_hash32();
    g.add_salt(S4);
    u64 h4 = g.form_hash32();
    
    return {
        h1 | (h2 << 32),
        h3 | (h4 << 16)
    };
}