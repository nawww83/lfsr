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


constexpr int p = 251;
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
	while (1) {
		lfsr8::u64 sum = 0;
		for (int i=0; i<m; ++i) {
			st[i] = g() & 255;
			st[i] <<= 8;
 			st[i] |= g() & 255;
			st[i] <<= 8;
 			st[i] |= g() & 255;
			st[i] <<= 8;
 			st[i] |= g() & 255;
			sum += st[i];
		}
		if (sum != 0) {
			break;
		}
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

auto research_periods(LFSR& g, long long max_T, int iters) {
	GeometricDistribution<int> r(0.3);
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
	const int iters = 256;
	auto sub_test = [&g, &K, &r](int saturation){
		int iter = 0;
		while (iter < iters) {
			//
			K = get_random_coeffs(r);
			g.set_K(K);

			g.set_unit_state();
			g.saturate(saturation);
			const STATE ref = g.get_state();
			for (int i=0; i<32 + (31*iter/iters); ++i) {
				g.next();
			}
			for (int i=0; i<32 + (31*iter/iters); ++i) {
				g.back();
			}
			bool test = g.is_state(ref);
			assert(test);
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
	//
	/*
	{
		cout << "Wait for Next-Back test..." << endl;
		test_next_back();
		cout << " Completed." << endl;
	}
	*/
	
	//
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
	{
		cout << "Wait test..." << endl;
		const long N = 1024*8;
		std::vector<uint8_t> v(N);
		std::map<lfsr8::u32, int> hashes;
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
	/*
	lfsr_rng::gens g;
	GeometricDistribution<int> r(0.3);

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

	int c = 0;
	int skeep = 0;
	double ave_dt = 0.;
	double ave_var_dt = 0.;
	double min_dt = std::numeric_limits<double>::infinity();
	double max_dt = 0.;
	while (1) {
		c++;
		STATE st = get_random_u32x4(r);
		auto st_c = state_conversion(st);
		timer.reset();
		g.seed(st_c);
		double dt = timer.elapsed_ns();
		ave_dt += (dt - ave_dt) / (1.*c);
		ave_var_dt += (dt*dt - ave_var_dt) / (1.*c);
		min_dt = std::min(min_dt, dt);
		max_dt = std::max(max_dt, dt);
		//
		if (! g.is_succes()) { // if not good, we must change initial state (seed), for ex., increment it etc
			skeep++; // it is better to get zero skeeps
			continue;
		}

		cout << "Counter: " << c << ", skeep: " << skeep << ", ave dt: " << ave_dt*1e-9 << " s, rms dt: " << std::sqrt(ave_var_dt - ave_dt*ave_dt)*1e-9 <<
			 ", max dt: " << max_dt*1e-9 << ", min dt: " << min_dt*1e-9 << endl;
		cout << "First 16 random numbers: ";
		auto pretty_print = [&g](int n) {
			for (int i=0; i<n; ++i) {
				g.next();
				cout << std::hex << g.get_u32() << std::dec << ((i < (n-1)) ? ", " : "");
			}
			cout << endl;
		};
		pretty_print(16);
		//
		auto measure_time = [&g](int n) {
			lfsr8::u32 h = 0;
			timer.reset();
			for (int i=0; i<n; ++i) {
				g.next();
				h ^= g.get_u32(); // 4 bytes = u32 value
			}
			double dt = timer.elapsed_ns();
			cout << "Total hash: " << std::hex << h << std::dec;
			cout << ", random generator performance: " << 4*1e+3*double(n)/dt << " MB/s" << endl;
		};
		measure_time(150000);
		cout << endl;
	}
	*/
	
    return 0;
}
