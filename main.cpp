#include <iostream>
#include <random>
#include <cmath>
#include <cstring>
#include "lfsr.hpp"

constexpr int p = 131;
constexpr int m = 4;

using LFSR = lfsr8::LFSR<p, m>;
using STATE = lfsr8::u32x8;

template<int x>
constexpr auto is_prime() {
	int d = 2;
	bool ok = true;
	while (d <= int(std::sqrt(x))) {
		ok &= ((x % d) != 0);
		d += 1;
	}
	return ok;
}

auto add_separators(std::string& value, char thousandSep = '\'') {
    int len = value.length();
    int dlen = 3;
    while (len > dlen) {
        value.insert(len - dlen, 1, thousandSep);
        dlen += 4;
        len += 1;
    }
    return value;
}

auto calculate_period(LFSR& g) {
	long long T = 1;
	// const auto old_state = g.get_state();
	g.set_unit_state();
	g.saturate();
	const auto ref = g.get_state();
	while (true) {
		g.next();
		if (g.is_state(ref))
			break;
		T += 1;
	}
	// g.set_state(old_state);
	return T;
}

auto find_max_period_polynomial(long long& T) {
	STATE K = {1};
	LFSR g(K);
	
	std::random_device rd;
    std::mt19937 gen(rd());
    std::geometric_distribution<> d(0.3);
	T = 0;
	const long long T_ref = std::pow(p, m) - 1;
	while (T != T_ref) {
		for (int i=0; i<m; ++i) {
			K[i] = d(gen) % p;
		}
		g.set_K(K);
		T = calculate_period(g);
	}
	return K;
}


int main() {
	using namespace std;
	
	static_assert(is_prime<p>());
	cout << "LFSR with modulo p: " << p << ", length m: " << m << endl;
	cout << "Wait for max period T = p^m - 1 polynomial look up..." << endl;
	long long T;
	STATE K = find_max_period_polynomial(T);
	cout << "Found coefficients K: (";
	for (int i=0; i<m-1; ++i) {
		cout << K[i] << ", ";
	}
	cout << K[m-1] << ")" << endl;
	auto T_str = std::to_string(T);
	cout << "Period T: " << add_separators(T_str) << endl;

    return 0;
}
