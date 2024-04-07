
#include <iostream>
#include <cmath>
#include <cstring>
#include <set>
#include <map>
#include <limits>
#include <vector>
#include <algorithm>
#include <bit>
#include <concepts>
#include <charconv>

#include "lfsr_hash.hpp"
#include "random_gen.hpp"
#include "random_gen_2.hpp"
#include "random_gen_3.hpp"
#include "timer.hpp"
#include "io_utils.hpp"

#include "random_utils.hpp"
#include "types.hpp"


template <int p, int m>
void increment_state(STATE<m>& x) {
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

template<int x>
constexpr auto is_prime() {
	int d = 2;
	bool ok = true;
	while (d < x) {
		ok &= ((x % d) != 0);
		d += 1;
	}
	return ok;
}

template <char sep=' ', int width=3>
void show_with_thousands_separator(std::integral auto x) {
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

static timer_n::Timer timer;

template <int p, int m>
auto calculate_period(LFSR<p, m>& g) {
	long long T = 1;
	const long long T0 = std::pow(p, m) - 1;
	const auto ref = g.get_state();
	while (true) {
		g.next();
		if (g.is_state(ref)) {
			break;
		}
		T += 1;
		if (T <= T0) {
			continue;
		}
		T = 0; // Abnormal period: unachievable state
		break;
	}
	return T;
}

template <int p, int m>
auto calculate_sawtooth_period(LFSR<p, m>& g, int q, int& i0) {
	assert( q > 0 );
	long long T = 1;
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
auto research_periods(LFSR<p, m>& g, long long max_T, int iters) {
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	std::set<long long> Ts;
	const long long T0 = std::pow(p, m) - 1;
	int iter = 0;
	while (true) {
		auto st = get_random_state<p, m>(r);
		g.set_state(st);
		auto T = calculate_period(g);
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
auto find_T1_polynomial(long long& T) { // Period T1 = p^(m-1) - 1.
	STATE<m> K = {1};
	LFSR<p, m> g(K);
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	T = 1;
	const long long T_ref = std::pow(p, m-1) - 1;
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
auto find_T0_polynomial(long long& T) { // maximal period Tmax = T0 = p^m - 1.
	STATE<m> K = {1};
	LFSR<p, m> g(K);
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	T = 1;
	const long long T_ref = std::pow(p, m) - 1;
	while (T != T_ref) {
		K = get_random_coeffs<p, m>(r);
		g.set_K(K);
		g.set_state(K);
		T = calculate_period(g);
	}
	return K;
}

template <int p, int m>
void test_next_back_inner() {
	STATE<m> K = {1};
	LFSR<p, m> g(K);
	using SAMPLE = lfsr8::MType<m>::SAMPLE;
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	const int iters = 256;
	auto sub_test = [&g, &K, &r](int saturation){
		int iter = 0;
		const int check_pos = 0;
		std::vector<SAMPLE> v{};
		assert(check_pos >= 0);
		assert(check_pos < m);
		while (iter < iters) {
			K = rnd_n::get_random_coeffs<p, m>(r);
			g.set_K(K);
			g.set_unit_state();
			g.saturate(saturation);
			const int max_i = 32 + ((16383*iter)/iters);
			const auto ref = g.get_state();
			for (int i=0; i<max_i; ++i) {
				g.next(i % p); // mod p: keep linearity
				v.push_back( g.get_cell(check_pos) );
			}
			bool result = true;
			for (int i=0; i<max_i; ++i) {
				result &= (v.back() == g.get_cell(check_pos));
				g.back((max_i - 1 - i) % p); // mod p: keep linearity
				v.pop_back();
			}
			assert(result);
			assert(g.is_state(ref));
			iter++;
		}
	};
	sub_test(4);
	sub_test(5);
	sub_test(8);
	sub_test(10);
}

void test_next_back() {
    std::cout << "Wait for Next-Back test...\n";
    test_next_back_inner<19, 4>();
    std::cout << " All Ok! Completed.\n";
}

void test_state_increment() {
    std::cout << "Wait for state increment test...\n";
    STATE<4> x {};
    const STATE<4> zero_state {};
    long long T = 0;
    const long long T0 = std::pow(19, 4);
    while (1) {
        increment_state<19, 4>(x);
        T++;
        if (x == zero_state) {
            break;
        }
    }
    // show_with_thousands_separator(T0);
    assert(T == T0);
    std::cout << " All Ok! Completed.\n";
}

void test_special() {
    const STATE<4> K = {4, 2, 2, 6};
    const int T_sawtooth = 7; // a prime
    LFSR<19, 4> g(K);
    STATE<4> initial_state {8, 7, 15, 13}; // special non-zero state
    g.set_state(initial_state);
    int i0 = 0;
    while (1) {
        auto T = calculate_sawtooth_period(g, T_sawtooth, i0);
        (void)T;
        const auto st = g.get_state();
        assert(initial_state == st);
        if (i0 == 0) {
            break;
        }
    }
}

void test_total_period() {
    const int p = 11;
    const int m = 4;
    std::cout << "Wait for total period test...\n";	
    const STATE<m> zero_state {};
    const long long T0 = std::pow(p, m) - 1;
    const STATE<m> K = {5, 1, 4, 0}; // p = 11
    const int init_i = 1;
    const int T_sawtooth = 7; // a prime: (p^m - 1)  mod T_sawtooth is not zero
    LFSR<p, m> g(K);
    STATE<m> initial_state {}; // any non-zero state
    increment_state<p, m>(initial_state);
    long long sum_T = 0;
    timer.reset();
    while (1) {
        g.set_state(initial_state);
        int i0 = init_i;
        long long TT = 0;
        int c = 0;
        while (1) {
            auto T = calculate_sawtooth_period(g, T_sawtooth, i0);
            TT += T;
            c++;
            if (i0 == init_i) {
                break;
            }
        }
        sum_T += TT;
        increment_state<p, m>(initial_state);
        if (initial_state != zero_state) {
            continue;
        }
        break;
    }
    auto dt = timer.elapsed_ns();
    std::cout << "Sum of T: ";
    show_with_thousands_separator(sum_T);
    assert(sum_T == (T0*T0*(long long)(T_sawtooth)) - T0*(long long)(T_sawtooth) + T_sawtooth);
    std::cout << " All Ok! Completed. Elapsed: " << dt*1e-9 << " s.\n";
}

void find_lfsr_coefficients_T0_period() {
    constexpr int p = 13;
    constexpr int m = 4;
    std::cout << "LFSR with modulo p: " << p << ", length m: " << m << std::endl;
    {
        std::cout << "Wait for maximal period T0 = p^m - 1 polynomial look up...\n";
        long long T;
        timer.reset();
        STATE<m> K = find_T0_polynomial<p, m>(T);
        auto dt = 1.e-9*timer.elapsed_ns();
        std::cout << " Found coefficients K: (";
        for (int i=0; i<m-1; ++i) {
            std::cout << K[i] << ", ";
        }
        std::cout << K[m-1] << ")" << ", dt: " << dt << " sec.\n";
        std::cout << " Period T0: ";
        show_with_thousands_separator(T);
    }
}

void find_lfsr_coefficients_T1_period() {
    constexpr int p = 13;
    constexpr int m = 4;
    long long T;
    std::cout << "Wait for period T1 = p^(m-1) - 1 polynomial look up...\n";
    timer.reset();
    STATE<m> K = find_T1_polynomial<p, m>(T);
    auto dt = 1.e-9*timer.elapsed_ns();
    std::cout << " Found coefficients K: (";
    for (int i=0; i<m-1; ++i) {
        std::cout << K[i] << ", ";
    }
    std::cout << K[m-1] << ")" << ", dt: " << dt << " sec.\n";
    std::cout << " Period T1: ";
    show_with_thousands_separator(T);
}

void test_lfsr_hash_bench() {
    std::cout << "Wait for LFSR hashes benchmark...\n";
    const size_t N = 65536*16*16;
    auto v = std::vector<uint8_t>(N);
    std::cout << "Input array of " << N << " bytes is allocated.\n";
    
    timer.reset();
    auto hash32 = lfsr_hash::hash32(v.data(), N);
    auto dt2 = timer.elapsed_ns();
    double perf2 = (1.e3*N)/dt2;

    timer.reset();
    auto hash64 = lfsr_hash::hash64(v.data(), N);
    auto dt3 = timer.elapsed_ns();
    double perf3 = (1.e3*N)/dt3;

    timer.reset();
    auto hash128 = lfsr_hash::hash128(v.data(), N);
    auto dt4 = timer.elapsed_ns();
    double perf4 = (1.e3*N)/dt4;
    //
    std::cout << "LFSR hashes:\n";
    std::cout << " 32-bit hash:  " << std::hex << hash32 << std::dec << "\t\t\t\t\t\tperf: " << perf2 << " MB/s\n";
    std::cout << " 64-bit hash:  " << std::hex << hash64 << std::dec << "\t\t\t\t\tperf: " << perf3 << " MB/s\n";
    std::cout << " 128-bit hash: " << std::hex << hash128.first << ":" << hash128.second << std::dec << "\t\tperf: " << perf4 << " MB/s\n";
    std::cout << " All Ok! Completed.\n";
}

void test_lfsr_hash_coverage_1() {
    io_u::io_utils io;
    std::cout << "Wait for LFSR 32-bit hashes coverage test 1.1...\n";	
    std::set<lfsr8::u32> hashes;
    for (int i=0; i<256; i++) {
        const uint8_t x = i;
        hashes.insert(lfsr_hash::hash32(&x, 1));
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == 256);
    for (int i=0; i<256*256; i++) {
        const lfsr8::u16 x = i;
        hashes.insert(lfsr_hash::hash32((uint8_t*)&x, 2));
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == (256u + 65536u));
    for (size_t i=0; i<256*256*2; i++) {
        const lfsr8::u32 x = i;
        uint8_t b[4];
        io.copy_to_mem_32(x, b, 4); // LE guaranteed
        const auto hash = lfsr_hash::hash32(b, 3);
        hashes.insert( hash );
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == (256u + 65536u + 65536u*2));
    for (size_t i=0; i<256*256; i++) {
        const lfsr8::u32 x = i;
        uint8_t b[4];
        io.copy_to_mem_32(x, b, 4); // LE guaranteed
        const auto hash = lfsr_hash::hash32(b, 4);
        hashes.insert( hash );
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == (256u + 65536u + 65536u*2 + 65536u));
    std::cout << " All Ok! Completed.\n";
}

void test_lfsr_hash_coverage_2() {
    io_u::io_utils io;
    std::cout << "Wait for LFSR 64-bit hashes coverage test 1.2...\n";	
    std::set<lfsr8::u64> hashes;
    for (int i=0; i<256; i++) {
        const uint8_t x = i;
        hashes.insert( lfsr_hash::hash64((uint8_t*)&x, 1) );
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == 256);
    for (size_t i=0; i<256*256; i++) {
        const lfsr8::u16 x = i;
        hashes.insert( lfsr_hash::hash64((uint8_t*)&x, 2) );
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == (256u + 65536u));
    for (size_t i=0; i<256ull*256ull*256ull; i++) {
        const lfsr8::u32 x = i;
        uint8_t b[4];
        io.copy_to_mem_32(x, b, 4); // LE guaranteed
        hashes.insert( lfsr_hash::hash64(b, 3) );
        hashes.insert( lfsr_hash::hash64(b, 4) );
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == (256ull + 65536ull + 2ull*256ull*256ull*256ull));
    std::cout << " All Ok! Completed.\n";
}

void test_lfsr_hash_coverage_3() {
    std::cout << "Wait for LFSR 128-bit hashes coverage test 1.3...\n";	
    std::set<std::pair<lfsr8::u64, lfsr8::u64>> hashes;
    for (int i=0; i<256; i++) {
        const uint8_t x = i;
        hashes.insert( lfsr_hash::hash128((uint8_t*)&x, 1) );
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == 256);
    for (size_t i=0; i<256*256; i++) {
        const lfsr8::u16 x = i;
        hashes.insert( lfsr_hash::hash128((uint8_t*)&x, 2) );
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == (256u + 65536u));
    for (size_t i=0; i<256ull*256ull; i++) {
        const auto x = std::array<lfsr8::u64, 32> {i};
        for (int i=3; i<=256; ++i) {
            hashes.insert( lfsr_hash::hash128((uint8_t*)&x, i) );
        }
    }
    std::cout << hashes.size() << std::endl;
    assert(hashes.size() == (256ull + 65536ull + (256ull - 3ull + 1ull)*65536ull));
    std::cout << " All Ok! Completed.\n";
}
	
void test_lfsr_hash_coverage_4() {
    std::cout << "Wait for LFSR hashes coverage test 2...\n";
    const size_t N = 8192;
    std::vector<uint8_t> v(N);
    std::set<lfsr8::u32> hashes;
    size_t s = 0;
    for (int val=0; val<256; val++) {
        v.assign(N, (uint8_t)val);
        assert(v.size() == N);
        for (size_t i=1; i<=N; i++) {
            hashes.insert( lfsr_hash::hash32(v.data(), i) );
        }
    }
    s += hashes.size();
    std::cout << " s = " << s << std::endl;
    // assert(N*256ull == s);
    assert(s >= 2096661);
    // when N = 8192 we can see the LFSR hash is comparable with SHA-512 by the size of hashes set 's = hashes.size()'
    // sha-512 (but it's low 32-bit) has s = 2'096'661
    // Current LFSR has s = 2'096'662 (~comparable)
    std::cout << " All Ok! Completed.\n";
}

void test_random_generators() {
	// Random generator infinite test
	#define gen_version 3
	#if gen_version == 1
		lfsr_rng::gens g;
	#endif
	#if gen_version == 2
		lfsr_rng_2::gens g;
	#endif
	#if gen_version == 3
		lfsr_rng_3::Generators g;
	#endif
	rnd_n::GeometricDistribution<int> r(0.3);
	r.seed();
	#if gen_version == 1 || gen_version == 3
		auto state_conversion = [](lfsr8::u32x4 st) {
			lfsr8::u16x8 st1;
			st1[0] = st[0];
			st1[1] = st[0] >> 16;
			st1[2] = st[1];
			st1[3] = st[1] >> 16;

			st1[4] = st[2];
			st1[5] = st[2] >> 16;
			st1[6] = st[3];
			st1[7] = st[3] >> 16;

			return st1;
		};
	#endif
	//
	long long c = 0;
	#if gen_version == 1 || gen_version == 3
		long long skeep = 0;
	#endif
	double ave_dt = 0;
	double ave_perf = 0;
	double max_perf = 0;
	double ave_var_dt = 0;
	double min_dt = std::numeric_limits<double>::infinity();
	double max_dt = 0;
	while (true) {
		std::cout << std::endl;
		auto st = rnd_n::get_random_u32x4<4>(r);
		timer.reset();
		#if gen_version == 1 || gen_version == 3
			const auto st_c = state_conversion(st);
			g.seed(st_c);
		#endif
		#if gen_version == 2
			g.seed(st);
		#endif
		double dt = timer.elapsed_ns();
		//
		#if gen_version == 1 || gen_version == 3
			if (! g.is_succes()) {
				skeep++; // it is better to achieve zero skeeps
				std::cout << "Skipped!\n";
				continue;
			}
		#endif
		//
		c++;
		ave_dt += (dt - ave_dt) / (1.*c);
		ave_var_dt += (dt*dt - ave_var_dt) / (1.*c);
		min_dt = std::min(min_dt, dt);
		max_dt = std::max(max_dt, dt);
		#if gen_version == 1
			const double T01 = std::pow(lfsr_rng::p1, 4) - 1;
			const double T02 = std::pow(lfsr_rng::p2, 4) - 1;
			const double T_bits = std::log2(g.T[4]) + std::log2(g.T[5]) + std::log2(g.T[6]) + std::log2(g.T[7]);
			cout << "Counter: " << c << ", skeep: " << skeep << ", ave dt: " << ave_dt*1e-9 << " s, rms dt: " << std::sqrt(ave_var_dt - ave_dt*ave_dt)*1e-9 <<
				", max dt: " << max_dt*1.e-09 << ", min dt: " << min_dt*1e-9 << "; " << g.ii01 << " : " << g.ii02 << " : " << 
				g.T[4]/T01 << " : " << g.T[5]/T01 << " : " << g.T[6]/T02 << " : " << g.T[7]/T02 << ", T bits: " << T_bits << endl;
		#endif
		#if gen_version == 3
			const double T_bits = g.period();
			std::cout << "Counter: " << c << ", skeep: " << skeep << ", ave dt: " << ave_dt*1e-9 << " s, rms dt: " << std::sqrt(ave_var_dt - ave_dt*ave_dt)*1e-9 <<
				", max dt: " << max_dt*1.e-09 << ", min dt: " << min_dt*1e-9 << "; T bits: " << T_bits << std::endl;
		#endif
		// Byte-wise chi-square test
		{
			static std::vector<double> frequencies(256);
			auto chisq = [](const std::vector<double>& f, double E) {
				double chi2 = 0;
				for (const auto& el : f) {
					chi2 += (el - E)*(el - E);
				}
				return chi2 / E;
			};
			frequencies.assign(256, 0.);
			static double ave_chi2 = 0;
			static double max_chi2 = 0;
			const int N = 16384;
			for (int i=0; i<N; ++i) {
		 		lfsr8::u64 x = g.next_u64();
				frequencies[(x >> 0) & 255]++;
				frequencies[(x >> 8) & 255]++;
				frequencies[(x >> 16) & 255]++;
				frequencies[(x >> 24) & 255]++;
				frequencies[(x >> 32) & 255]++;
				frequencies[(x >> 40) & 255]++;
				frequencies[(x >> 48) & 255]++;
				frequencies[(x >> 56) & 255]++;
			}
			double chi2 =  chisq(frequencies, 8*N/256.);
			ave_chi2 += (chi2 - ave_chi2) / (1.*c);
			max_chi2 = std::max(chi2, max_chi2);
			std::cout << "Byte-wise chi-square: " << chi2 << ", N: " << 8*N << " bytes. Ave chi2: " << ave_chi2 << ", max chi2: " << max_chi2 << std::endl;
		}
		//
		auto measure_time = [&g, ave_perf, max_perf](size_t n) {
			timer.reset();
			lfsr8::u64 h = 0;
			for (size_t i=0; i<n; ++i) {
				h ^= g.next_u64();
			}
			double dt = timer.elapsed_ns();
			std::cout << "Hash u64: " << std::hex << h << std::dec;
			double perf = 8. * 1000. * double(n) / dt; // 16 bit = 2 bytes; 64 bit = 8 bytes
			std::cout << ", Random Generator performance: " << perf << ", on average: " << ave_perf << ", max: " << max_perf << " MB/s" << std::endl;
			return perf;
		};
		auto perf = measure_time(10000000);
		ave_perf += (perf - ave_perf) / (1.*c);
		max_perf = std::max(perf, max_perf);
	}
}
