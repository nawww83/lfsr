#pragma once

#include <iostream>
#include <cmath>
#include <cstring>
#include <utility>
#include <limits>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <bit>
#include <bitset>
#include <concepts>
#include <charconv>

#include "timer.hpp"
#include "io_utils.hpp"

#include "random_utils.hpp"
#include "types.hpp"

template<int x>
inline constexpr auto is_prime() {
	int d = 2;
	bool ok = true;
	while (d < x) {
		ok &= ((x % d) != 0);
		d += 1;
	}
	return ok;
}

inline constexpr u64 sqrt_int(u64 x) {
    return static_cast<i64>(std::sqrt(x));
}

inline auto factor(u64 x) {
    std::map<u64, int> res{};
    auto inner_loop = [&res](u64 d, u64& x) -> bool {
        bool x_changed = false;
        while ((x % d == 0) && (x > 1)) {
            x /= d;
            res[d]++;
            x_changed = true;
        }
        return x_changed;
    };
    constexpr auto max_d = [](u64 x) -> u64 {
        return sqrt_int(x);
    };
    inner_loop(2, x);
    u64 d = 3;
    u64 md = max_d(x);
    while (d <= md) {
        md = inner_loop(d, x) ? max_d(x) : md;
        d += 2;
        d += (d % 3 == 0) ? 2 : 0;
    }
    if (x > 1 || res.empty()) {
        res[x]++;
    }

    return res;
}


inline auto get_random_u16x8(int64_t offset) {
    const int64_t since_epoch_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    std::seed_seq g_rnd_sequence{since_epoch_ms & 255, (since_epoch_ms >> 8) & 255,
                                 (since_epoch_ms >> 16) & 255, (since_epoch_ms >> 24) & 255, offset};
    lfsr8::u16x8 st;
    g_rnd_sequence.generate(st.begin(), st.end());
    return st;
}

inline auto get_random_u32x4(int64_t offset) {
    const int64_t since_epoch_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    std::seed_seq g_rnd_sequence{since_epoch_ms & 255, (since_epoch_ms >> 8) & 255,
                                 (since_epoch_ms >> 16) & 255, (since_epoch_ms >> 24) & 255, offset};
    lfsr8::u32x4 st;
    g_rnd_sequence.generate(st.begin(), st.end());
    return st;
}

template <char sep=' ', int width=3>
inline void show_with_thousands_separator(std::integral auto x) {
    constexpr int max_digits = int(std::log2(2.) * std::popcount(-1ull) / std::log2(10.)) + 1;
    std::array<char, max_digits> str;
    std::array<char, max_digits + max_digits/width> str_th;
    if (auto [ptr, ec] = std::to_chars(str.data(), str.data() + str.size(), x); ec == std::errc()) {
        const int sign = (str[0] == '-') ? 1 : 0;
        const int digits = ptr - str.data();
        const int num_of_seps = (digits-1-sign) / width;
        const int new_len = num_of_seps + digits;
        char* c_ptr = str_th.data() + new_len - width;
        for (int i=0; i<num_of_seps; ++i) {
            std::memcpy(c_ptr, str.data() + digits - (i+1)*width, width);
            c_ptr -= sizeof(sep);
            *c_ptr = sep;
            c_ptr -= width;
        }
        const int reseque = (digits-sign) % width;
        const bool is_full_group = (reseque == 0);
        const int res_len = (!is_full_group) ? reseque + sign : width + sign;
        std::memcpy(str_th.data(), str.data(), res_len);
        std::cout << std::string_view(str_th.data(), str_th.data() + new_len) << '\n';
    }
}

template <typename T, typename Generator>
inline T GetQuasiGaussSample(Generator& generator) { // Sum of 4 uniform samples.
    uint64_t g_sample_u64 = generator.next_u64();
    T g_sample_f = 0;
    g_sample_f += T(g_sample_u64 & 0xffff) / T(65535);
    g_sample_u64 >>= 16;
    g_sample_f += T(g_sample_u64 & 0xffff) / T(65535);
    g_sample_u64 >>= 16;
    g_sample_f += T(g_sample_u64 & 0xffff) / T(65535);
    g_sample_u64 >>= 16;
    g_sample_f += T(g_sample_u64 & 0xffff) / T(65535);
    g_sample_f = T(0.5) * (g_sample_f - 2);
    return g_sample_f; // In range (-1, 1).
}

template <int p, int m>
inline void increment_state(STATE<m>& x) {
	int carry = 0;
	x[0]++;
	if (x[0] == p) {
		carry = 1;
		x[0] = 0;
	}
	for (int i=1; i<m; ++i) {
		x[i] += carry;
		if (x[i] == p) {
			carry = 1;
			x[i] = 0;
		} else {
			carry = 0;
		}
	};
}

template <int p, int m>
inline std::pair<bool, u64> calculate_period(LFSR<p, m>& g) {
	const u64 T0 = std::pow(p, m) - 1;
	const auto multipliers = factor(T0);
    std::set<u64> divisors{1}; // must be sorted!
    //
    for (const auto [t, c] : multipliers) {
        u64 T_ = t;
        std::set<u64> div_tmp{};
        for (int i=0; i<c; ++i) {
            div_tmp.insert(T_);
            T_ *= t;
        }
        std::set<u64> d_copy = divisors;
        for (auto d1 : d_copy) {
            for (auto d2 : div_tmp) {
                divisors.insert(d1*d2);
            }
        }
    }
    //
	const auto ref = g.get_state();
    u64 T = 0;
	for (const auto t : divisors) {
		while (T != t) {
			g.next();
			T++;
		}
		if (g.is_state(ref)) {
			return {true, T};
		}
	}
	return {false, 0};
}

template <int p, int m>
inline bool is_maximal_period(LFSR<p, m>& g, u64& T) {
	const u64 T0 = std::pow(p, m) - 1;
	const auto multipliers = factor(T0);
    std::set<u64> divisors{1}; // must be sorted!
    //
    for (const auto [t, c] : multipliers) {
        u64 T_ = t;
        std::set<u64> div_tmp{};
        for (int i=0; i<c; ++i) {
            div_tmp.insert(T_);
            T_ *= t;
        }
        std::set<u64> d_copy = divisors;
        for (auto d1 : d_copy) {
            for (auto d2 : div_tmp) {
                divisors.insert(d1*d2);
            }
        }
    }
    //
	T = 0;
	constexpr int num_of_bits = 64;
	g.set_unit_state();
	const auto ref = g.get_state();
	LFSR<p, m> tmp_reg = g;
	for (const u64 t : divisors) {
		const auto& bits_str = std::bitset<num_of_bits>(t).to_string();
		const int bit_width = bits_str.length();
		tmp_reg.set_unit_state();
		for (int i=0; i<bit_width; ++i) {
			if (bits_str[i] == '0') continue;
			g.set_unit_state();
			g.next();
			for (int c = 0; c < (num_of_bits-1-i); ++c) {
				g.square();
			}
			tmp_reg.mult_by(g.get_state());
		}
		if (tmp_reg.is_state(ref) && (t != T0)) {
			T = t;
			return false;
		}
		if (tmp_reg.is_state(ref) && (t == T0)) {
			T = t;
			return true;
		}
	}
	return false;
}

template <int p, int m>
inline auto calculate_sawtooth_period(LFSR<p, m>& g, int q, int& i0) {
	assert( q > 0 );
	u64 T = 1;
	const auto ref = g.get_state();
	i0 %= q;
	while (true) {
		g.next(i0);
		i0++;
		i0 %= q;
		if (! g.is_state(ref)) {
			T++;
			continue;
		}
		break;
	}
	return T;
}

template <int p, int m>
inline auto research_periods(LFSR<p, m>& g, u64 max_T, int iters) {
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	std::set<u64> Ts;
	const u64 T0 = std::pow(p, m) - 1;
	int iter = 0;
	while (true) {
		auto st = get_random_state<p, m>(r);
		g.set_state(st);
		const auto [is_divisor, T] = calculate_period<p, m>(g);
		if (T < 1) {
			std::cout << "Abnormal period detected! Skip it.\n";
			continue;
		}
		if (T > max_T) {
			break;
		}
		Ts.insert( T );
		if (T == T0) {
			break;
		}
		++iter;
		if (iter >= iters) {
			break;
		}
	}
	return Ts;
}

template <int p, int m>
auto find_T1_polynomial(u64& T) { // Period T1 = p^(m-1) - 1.
	STATE<m> K = {1};
	LFSR<p, m> g(K);
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	T = 1;
	const u64 T_ref = std::pow(p, m-1) - 1;
	while (true) {
		K = get_random_coeffs<p, m>(r);
		g.set_K(K);
		auto Ts = research_periods(g, T_ref, 30);
		if (Ts.contains(T_ref)) {
			T = T_ref;
			break;
		}
	}
	return K;
}

template <int p, int m>
auto find_T0_polynomial(u64& T) { // maximal period Tmax = T0 = p^m - 1.
	STATE<m> K = {1};
	LFSR<p, m> g(K);
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	T = 0;
	for (;;) {
		K = get_random_coeffs<p, m>(r);
		g.set_K(K);
		const auto is_divisible = is_maximal_period<p, m>(g, T);
        if (!is_divisible) {
            std::cout << " ... skipped coefficients K: (";
            for (int i=0; i<m-1; ++i) {
                std::cout << K[i] << ", ";
            }
            std::cout << K[m-1] << ")\n";
        } else {
			break;
		}
	}
	return K;
}

template <int p>
void test_next_back_inner_1() {
	STATE<4> K = {1, 0, 0, 0};
	LFSR<p, 4> g(K);
	using SAMPLE = lfsr8::MType<4>::SAMPLE;
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	const int iters = 256;
	auto sub_test = [&g, &K, &r](int saturation){
		int iter = 0;
		const int check_pos = 0;
		std::vector<SAMPLE> v{};
		assert(check_pos >= 0);
		assert(check_pos < 4);
		while (iter < iters) {
			K = rnd_n::get_random_coeffs<p, 4>(r);
			g.set_K(K);
			g.set_unit_state();
			for (int i=0; i<saturation; ++i) {
				g.next();
			}
			const int max_i = 32 + ((16383*iter)/iters);
			const auto ref = g.get_state();
			for (int i=0; i<max_i; ++i) {
				g.next(i % p); // mod p: keep linearity
				v.push_back( g.get_state()[check_pos] );
			}
			bool result = true;
			for (int i=0; i<max_i; ++i) {
				result &= (v.back() == g.get_state()[check_pos]);
				g.back((max_i - 1 - i) % p); // mod p: keep linearity
				v.pop_back();
			}
			assert(result);
			assert(g.is_state(ref));
			iter++;
		}
	};
	sub_test(8);
	sub_test(9);
	sub_test(12);
	sub_test(14);
}

template <int p>
void test_next_back_inner_2() {
	STATE<8> K = {1, 0, 0, 0, 0, 0, 0, 0};
	LFSR<p, 8> g(K);
	using SAMPLE = lfsr8::MType<8>::SAMPLE;
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	const int iters = 256;
	auto sub_test = [&g, &K, &r](int saturation){
		int iter = 0;
		const int check_pos = 0;
		std::vector<SAMPLE> v{};
		assert(check_pos >= 0);
		assert(check_pos < 8);
		while (iter < iters) {
			K = rnd_n::get_random_coeffs<p, 8>(r);
			g.set_K(K);
			g.set_unit_state();
			for (int i=0; i<saturation; ++i) {
				g.next();
			}
			const int max_i = 32 + ((16383*iter)/iters);
			const auto ref = g.get_state();
			for (int i=0; i<max_i; ++i) {
				g.next(i % p); // mod p: keep linearity
				v.push_back( g.get_state()[check_pos] );
			}
			bool result = true;
			for (int i=0; i<max_i; ++i) {
				result &= (v.back() == g.get_state()[check_pos]);
				g.back((max_i - 1 - i) % p); // mod p: keep linearity
				v.pop_back();
			}
			assert(result);
			assert(g.is_state(ref));
			iter++;
		}
	};
	sub_test(8);
	sub_test(9);
	sub_test(12);
	sub_test(14);
}

template <int p>
void test_next_back_inner_pair() {
	STATE<8> K = {1, 0, 0, 0, 1, 0, 0, 0};
	LFSR_paired<p> g(K);
	using SAMPLE = lfsr8::MType<8>::SAMPLE;
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	const int iters = 256;
	auto sub_test = [&g, &K, &r](int saturation){
		int iter = 0;
		const int check_pos = 0;
		std::vector<SAMPLE> v{};
		assert(check_pos >= 0);
		assert(check_pos < 8);
		while (iter < iters) {
			K = rnd_n::get_random_paired_coeffs<p>(r);
			g.set_K(K);
			g.set_unit_state();
			for (int i=0; i<saturation; ++i) {
				g.next();
			}
			const int max_i = 32 + ((16383*iter)/iters);
			const auto ref = g.get_state();
			for (int i=0; i<max_i; ++i) {
				g.next(i % p, i % p); // mod p: keep linearity
				v.push_back( g.get_state()[check_pos] );
			}
			bool result = true;
			for (int i=0; i<max_i; ++i) {
				result &= (v.back() == g.get_state()[check_pos]);
				g.back((max_i - 1 - i) % p, (max_i - 1 - i) % p); // mod p: keep linearity
				v.pop_back();
			}
			assert(result);
			assert(g.is_state_low(ref));
			assert(g.is_state_high(ref));
			iter++;
		}
	};
	sub_test(8);
	sub_test(9);
	sub_test(12);
	sub_test(14);
}

void test_next_back();

void test_state_increment();

void test_special();

void test_total_period();

void test_bias();

void test_random_generators();

void test_random_generator_next_back();

void find_lfsr_coefficients_T0_period();

template <int p, int m>
void test_some_poly(STATE<m> K) {
    const u64 T0 = std::pow(p, m) - 1;
	LFSR<p, m> g(K);
    g.set_K(K);
	u64 T;
    const auto is_divisor = is_maximal_period(g, T);
    assert(is_divisor);
	assert(T == T0);
}

void find_lfsr_coefficients_T1_period();

void test_lfsr_hash_coverage_1();

void test_lfsr_hash_coverage_2();

void test_lfsr_hash_coverage_3();

void test_lfsr_hash_coverage_4();

void test_square_of_lfsr();

void test_mult_by_of_lfsr();