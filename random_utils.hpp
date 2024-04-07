#pragma once

#include <random>

#include "types.hpp"

namespace rnd_n {

inline auto G_SEED = std::random_device{}();

template <typename T>
class GeometricDistribution {
private:
    std::mt19937 gen;
    std::geometric_distribution<T> dist;
public:
    GeometricDistribution(double p): gen(G_SEED), dist(p) {}
	GeometricDistribution(GeometricDistribution&& other) noexcept: dist(std::move(other.dist)), gen(std::move(other.gen)) {}
    T operator()() { return dist(gen); }
	void seed() {
		auto t = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    	std::seed_seq seq {t & 255, (t >> 8) & 255, (t >> 16) & 255, (t >> 24) & 255};
    	gen.seed( seq );
	}
};

template<int p, int m, class Generator>
inline auto get_random_state(Generator& g) {
	STATE<m> st {};
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

template<int p, int m, class Generator>
inline auto get_random_coeffs(Generator& g) {
	STATE<m> st {};
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

template<int m, class Generator>
inline auto get_random_u32x4(Generator& g) {
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

}