#include <iostream>
#include <random>
#include <cmath>
#include <cstring>
#include <set>
#include <map>
#include <limits>
#include <vector>

#include "lfsr_hash.hpp"
#include "random_gen.hpp"

#include "timer.hpp"


constexpr int p = 241;
constexpr int m = 4;

using LFSR = lfsr8::LFSR<p, m>;
using STATE = lfsr8::MType<m>::STATE;


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

template <typename T>
std::string format_with_commas(T x) {
	auto s = std::to_string(x);
	int n = s.length() - 3;
	const int stop = (x >= 0) ? 0 : 1;
	while (n > stop) {
		s.insert(n, "'");
		n -= 3;
	}
	return s;
}

static std::random_device rd;


template <typename T>
class GeometricDistribution {
private:
    std::mt19937 gen;
    std::geometric_distribution<T> dist;
public:
    GeometricDistribution(double p): gen(rd()), dist(p) {}
	GeometricDistribution(GeometricDistribution&& other) noexcept: dist(std::move(other.dist)), gen(std::move(other.gen)) {}
    T operator()() { return dist(gen); }
	void seed() {
		auto t = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    	std::seed_seq seq {t & 255, (t >> 8) & 255, (t >> 16) & 255, (t >> 24) & 255};
    	gen.seed( seq );
	}
};

template<class Generator>
auto get_random_state(Generator& g) {
	STATE st {};
	while (1) {
		int sum = 0;
		for (int i=0; i<m; ++i) {
			st[i] = g() % p;
			sum += st[i];
		}
		if (sum != 0) {
			break;
		}
	}
	return st;
}

template<class Generator>
auto get_random_coeffs(Generator& g) {
	STATE st {};
	while (1) {
		for (int i=0; i<m; ++i) {
			st[i] = g() % p;
		}
		if (st[0] != 0) {
			break;
		}
	}
	return st;
}

template<class Generator>
auto get_random_u32x4(Generator& g) {
	lfsr8::u32x4 st {};
	for (int i=0; i<m; ++i) {
		st[i] = g() & 255;
		st[i] <<= 8;
		st[i] |= g() & 255;
		st[i] <<= 8;
		st[i] |= g() & 255;
		st[i] <<= 8;
		st[i] |= g() & 255;
	}
	return st;
}

auto calculate_period(LFSR& g) {
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

// auto calculate_sawtooth_period(LFSR& g, int q, int i0=0) {
// 	assert( q > 0 );
// 	const long long T0 = std::pow(p, m) - 1;
// 	long long T = 1;
// 	long long prev_T = 0;
// 	const auto ref = g.get_state();
// 	int i = i0;
// 	while (true) {
// 		g.next(i);
// 		i++;
// 		i %= q;
// 		if (! g.is_state(ref)) {
// 			T++;
// 			continue;
// 		}
// 		long long curr_T = T - prev_T;
// 		std::cout << " sub T: " << format_with_commas(curr_T) << ", i: " << i << std::endl;
// 		prev_T = T;
// 		if (i == i0) {
// 			break;
// 		}
// 		T++;
// 		if (T > (long long)(q)*T0) {
// 			T = 0; // Abnormal period: unachievable state
// 			break;	
// 		}
// 	}
// 	return T;
// }

auto research_periods(LFSR& g, long long max_T, int iters) {
	GeometricDistribution<int> r(0.3);
	r.seed();
	std::set<long long> Ts;
	const long long T0 = std::pow(p, m) - 1;
	int iter = 0;
	while (true) {
		auto st = get_random_state(r);
		g.set_state(st);
		auto T = calculate_period(g);
		if (T < 1) {
			std::cout << "Abnormal period detected! Skip it." << std::endl;
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

auto find_T1_polynomial(long long& T) {
	STATE K = {1};
	LFSR g(K);
	GeometricDistribution<int> r(0.3);
	r.seed();
	T = 1;
	const long long T_ref = std::pow(p, m-1) - 1;
	while (true) {
		K = get_random_coeffs(r);
		g.set_K(K);
		auto Ts = research_periods(g, T_ref, 30);
		if (Ts.contains(T_ref)) {
			T = T_ref;
			break;
		}
	}
	return K;
}

auto find_T0_polynomial(long long& T) {
	STATE K = {1};
	LFSR g(K);
	GeometricDistribution<int> r(0.3);
	r.seed();
	T = 1;
	const long long T_ref = std::pow(p, m) - 1;
	while (T != T_ref) {
		K = get_random_coeffs(r);
		g.set_K(K);
		g.set_state(K);
		T = calculate_period(g);
	}
	return K;
}

void test_next_back() {
	STATE K = {1};
	LFSR g(K);
	GeometricDistribution<int> r(0.3);
	r.seed();
	const int iters = 256;
	auto sub_test = [&g, &K, &r](int saturation){
		int iter = 0;
		const int check_pos = 0;
		std::vector<lfsr8::MType<m>::SAMPLE> v{};
		assert(check_pos >= 0);
		assert(check_pos < m);
		while (iter < iters) {
			K = get_random_coeffs(r);
			g.set_K(K);
			g.set_unit_state();
			g.saturate(saturation);
			const int max_i = 32 + ((16383*iter)/iters);
			const STATE ref = g.get_state();
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


static timer_n::Timer timer;

int main() {
	using namespace std;
	
	static_assert(is_prime<p>());
	/*
	cout << "LFSR with modulo p: " << p << ", length m: " << m << endl;
	{
		cout << "Wait for maximal period T0 = p^m - 1 polynomial look up..." << endl;
		long long T;
		timer.reset();
		STATE K = find_T0_polynomial(T);
		auto dt = 1.e-9*timer.elapsed_ns();
		cout << " Found coefficients K: (";
		for (int i=0; i<m-1; ++i) {
			cout << K[i] << ", ";
		}
		cout << K[m-1] << ")" << ", dt: " << dt << " sec." << endl;
		auto T_str = format_with_commas(T);
		cout << " Period T0: " << T_str << endl;
	}
	*/
	//
	/*
	{
		long long T;
		cout << "Wait for period T1 = p^(m-1) - 1 polynomial look up..." << endl;
		timer.reset();
		STATE K = find_T1_polynomial(T);
		auto dt = 1.e-9*timer.elapsed_ns();
		cout << " Found coefficients K: (";
		for (int i=0; i<m-1; ++i) {
			cout << K[i] << ", ";
		}
		cout << K[m-1] << ")" << ", dt: " << dt << " sec." << endl;
		auto T_str = format_with_commas(T);
		cout << " Period T1: " << T_str << endl;
	}
	*/
	//
	{
		cout << "Wait for Next-Back test..." << endl;
		test_next_back();
		cout << " All Ok! Completed." << endl;
	}
	// {
	// 	const STATE K = {13, 2, 5, 10};  // T = 241^4 - 1
	// 	const STATE initial_state = {1, 1, 0, 0};
	// 	LFSR g(K);
	// 	g.set_state(initial_state);
	// 	cout << "Initial state: ";
	// 	for (int i=0; i<m; ++i) {
	// 		cout << g.get_cell(i) << ", ";
	// 	}
	// 	cout << endl;
	// 	cout << "Wait for period calculation..." << endl;
	// 	auto T = calculate_sawtooth_period(g, 7);
	// 	cout << " T: " << format_with_commas(T) << endl;
	// }
	// LFSR hashes test
	{
		const int N = 65536*16*16;
		auto v = std::vector<uint8_t>(N);
		cout << "Input array of " << N << " bytes is allocated." << endl;
		
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
		cout << "LFSR hashes:" << endl;
		cout << " 32-bit hash:  " << std::hex << hash32 << std::dec << "\t\t\t\t\t\tperf: " << perf2 << " MB/s" << endl;
		cout << " 64-bit hash:  " << std::hex << hash64 << std::dec << "\t\t\t\t\tperf: " << perf3 << " MB/s" << endl;
		cout << " 128-bit hash: " << std::hex << hash128.first << ":" << hash128.second << std::dec << "\t\tperf: " << perf4 << " MB/s" << endl;
	}
	{
		std::map<lfsr8::u32, int> hashes;
		for (int i=0; i<256; i++) {
			const uint8_t x = i;
			hashes[lfsr_hash::hash32(&x, 1)] = i;
		}
		cout << hashes.size() << endl;
		assert(hashes.size() == 256);
		for (int i=0; i<256*256; i++) {
			const lfsr8::u16 x = i;
			hashes[lfsr_hash::hash32((uint8_t*)&x, 2)] = i;
		}
		cout << hashes.size() << endl;
		assert(hashes.size() == (256u + 65536u));
	}
	//
	{
		const long N = 1024*8; // the more N is passed, the better the LFSR hash
		std::vector<uint8_t> v(N);
		std::map<lfsr8::u32, int> hashes;
		cout << "N: " << N << ", wait test..." << endl;
		long s = 0;
		for (int val=0; val<256; val++) {
			hashes.clear();
			v.assign(N, (uint8_t)val);
			assert(v.size() == N);
			for (int i=1; i<=N; i++) {
				hashes[lfsr_hash::hash32(v.data(), i)] = i;
			}
			// cout << hashes.size() << endl;
			// assert(hashes.size() > N - 2);
			s += hashes.size();
		}
		assert((N*256 - s) == 0);
		cout << " => Passed." << endl;
	}
	// Random generator test
	lfsr_rng::gens g;
	GeometricDistribution<int> r(0.3);
	r.seed();
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
	//
	long long c = 0;
	long long skeep = 0;
	double ave_dt = 0;
	double ave_perf = 0;
	double ave_var_dt = 0;
	double min_dt = std::numeric_limits<double>::infinity();
	double max_dt = 0;
	while (true) {
		cout << endl;
		STATE st = get_random_u32x4(r);
		const auto st_c = state_conversion(st);
		timer.reset();
		g.seed(st_c);
		double dt = timer.elapsed_ns();
		//
		if (! g.is_succes()) {
			skeep++; // it is better to get zero skeeps
			std::cout << "Skipped! " << std::endl;
			continue;
		}
		//
		c++;
		ave_dt += (dt - ave_dt) / (1.*c);
		ave_var_dt += (dt*dt - ave_var_dt) / (1.*c);
		min_dt = std::min(min_dt, dt);
		max_dt = std::max(max_dt, dt);
		const double T01 = std::pow(lfsr_rng::p1, 4) - 1;
		const double T02 = std::pow(lfsr_rng::p2, 4) - 1;
		const double T_bits = std::log2(g.T[4]) + std::log2(g.T[5]) + std::log2(g.T[6]) + std::log2(g.T[7]);
		cout << "Counter: " << c << ", skeep: " << skeep << ", ave dt: " << ave_dt*1e-9 << " s, rms dt: " << std::sqrt(ave_var_dt - ave_dt*ave_dt)*1e-9 <<
			 ", max dt: " << max_dt*1.e-09 << ", min dt: " << min_dt*1e-9 << "; " << g.ii01 << " : " << g.ii02 << " : " << 
			 g.T[4]/T01 << " : " << g.T[5]/T01 << " : " << g.T[6]/T02 << " : " << g.T[7]/T02 << ", T bits: " << T_bits << endl;
		
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
			const int N = 16384;
			for (int i=0; i<N; ++i) {
		 		lfsr8::u64 x = g.next_u64();
				// lfsr8::u16 x = g.get_u16();
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
			cout << "Byte-wise chi-square: " << chi2 << ", N: " << 8*N << " bytes. Ave chi2: " << ave_chi2 << endl;
		}
		//
		// auto pretty_print = [&g](int n) {
		// 	for (int i=0; i<n; ++i) {
		// 		g.next();
		// 		cout << std::hex << g.get_u16() << std::dec << ((i < (n-1)) ? ", " : "");
		// 	}
		// 	cout << endl;
		// };
		// cout << "First 32 16-bit random numbers: ";
		// pretty_print(32);
		auto measure_time = [&g, ave_perf](int n) {
			timer.reset();
			lfsr8::u64 h = 0;
			for (int i=0; i<n; ++i) {
				h ^= g.next_u64();
				// h ^= g.get_u16();
			}
			double dt = timer.elapsed_ns();
			cout << "Hash u64: " << std::hex << h << std::dec;
			double perf = 8. * 1000. * double(n) / dt; // 16 bit = 2 bytes; 64 bit = 8 bytes
			cout << ", Random Generator performance: " << perf << ", on average: " << ave_perf << " MB/s" << endl;
			return perf;
		};
		auto perf = measure_time(10000000);
		// measure_time(10000);
		// measure_time(100000);
		//
		ave_perf += (perf - ave_perf) / (1.*c);
	}
	
    return 0;
}
