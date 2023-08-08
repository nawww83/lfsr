#include "lfsr_hash.hpp"


lfsr_hash::gens g;

lfsr8::u16 lfsr_hash::hash(const uint8_t* input, int n) {
    g.reset();
    g.process_input(input, n);
    return g.form_hash();
}