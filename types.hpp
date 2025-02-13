#pragma once

#include <cstdint>
#include <type_traits>
#include "lfsr.hpp"

using u64 = uint64_t;
using i64 = int64_t;
using u32 = uint32_t;
using u16 = uint16_t;

template <int p, int m>
using LFSR = lfsr8::LFSR<p, m>;

template <int p>
using LFSR_paired = lfsr8::LFSR_paired_2x4<p>;

template <int m>
using STATE = lfsr8::MType<m>::STATE;