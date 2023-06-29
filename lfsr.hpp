#include <cstdint>
#include <cassert>
#include <array>

namespace lfsr8 {
using u32 = uint32_t;
using u32x8 = std::array<u32, 8>;

static auto is_zero(u32x8 state) {
	bool zf = true;
	for (int i=0; i<8; ++i) {
		zf &= (state[i] == 0);
	}
	return zf;
}

template <int p, int m>
class LFSR {
	static_assert(m <= 8);
	static_assert(m > 0);
	static_assert(p < 256);
	static_assert(p > 1);
public:
    	constexpr LFSR(u32x8 K): m_K(K) {};
	void set_state(u32x8 state) {
		m_state = state;
		v = m_state[m-1];
	}
   	void set_unit_state() {
		m_state = {1};
		v = m_state[m-1];
	}
	void set_K(u32x8 K) {
		m_K = K;
	}
	void next(u32 input=0) {
	    for (int i=m-1; i>0; i--) {
	        m_state[i] = (m_state[i-1] + v*m_K[i]) % p;
		}
		m_state[0] = (input + v*m_K[0]) % p;
		v = m_state[m-1];
	}
	void saturate(int q=m) {
		assert(! is_zero(m_state));
		assert(q >= m);
		for (int i=0; i<q; ++i) {
			next();
		}
	}
	bool is_state(u32x8 state) const {
		return (state == m_state);
	}
	auto get_state() const {
		return m_state;
	}
private:
	u32x8 m_state {};
	u32x8 m_K {};
	u32 v {};
};

}
