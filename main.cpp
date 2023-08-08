#include <iostream>
#include <random>
#include <cmath>
#include <cstring>
#include <set>

#include "lfsr_hash.hpp"

#include "timer.hpp"


constexpr int p = 251;
constexpr int m = 4;

using LFSR = lfsr8::LFSR<p, m>;
using STATE = lfsr8::u16x8;


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

template <typename T>
class GeometricDistribution {
private:
    std::random_device rd;
    std::mt19937 gen;
    std::geometric_distribution<T> dist;
public:
    GeometricDistribution(double p): gen(rd()), dist(p) {}
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

auto calculate_period(LFSR& g) {
	long long T = 1;
	const long long T0 = std::pow(p, m) - 1;
	const auto ref = g.get_state();
	while (true) {
		g.next();
		if (g.is_state(ref))
			break;
		T += 1;
		if (T > T0) {
			T = 0; // Abnormal period: unachievable state
			break;
		}
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



int main() {
	using namespace std;

	timer_n::Timer timer;
	
	static_assert(is_prime<p>());
	
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
	//
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
	//
	//
	{
		cout << "Wait for Next-Back test..." << endl;
		test_next_back();
		cout << " Completed." << endl;
	}
	//
	{
		const int N = 65536*16;
		auto v = new uint8_t[N];
		assert(v != nullptr);
		cout << "Input array of " << N << " bytes is allocated." << endl;
		timer.reset();
		auto hash = lfsr_hash::hash(v, N);
		auto dt = timer.elapsed_ns();
		double perf = (1.e3*N)/dt;
		delete [] v;
		cout << "LFSR 16 bit hash: " << std::hex << hash << std::dec << ", perf: " << perf << " MB/s." << endl;
	}

    return 0;
}
