#include <cstdint>
#include <cassert>
#include <array>

namespace lfsr8 {
using u32 = uint32_t;
using u16 = uint16_t;
using u16x8 = std::array<u16, 8>;
using u32x8 = std::array<u32, 8>;


template <typename T>
static auto is_zero(T state) {
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
	constexpr LFSR(u16x8 K): m_K(K) {
		m_calc_inv_K0();
	};
	void set_state(u16x8 state) {
		m_state = state;
		m_v = m_state[m-1];
	}
   	void set_unit_state() {
		m_state = {1};
		m_v = m_state[m-1];
	}
	void set_K(u16x8 K) {
		m_K = K;
		m_calc_inv_K0();
	}
	void next(u16 input=0) {
	    for (int i=m-1; i>0; i--) {
	        m_state[i] = (m_state[i-1] + m_v*m_K[i]) % p;
		}
		m_state[0] = (input + m_v*m_K[0]) % p;
		m_v = m_state[m-1];
	}
	void back(u16 input=0) {
	    m_v = (m_inv_K0*(m_state[0] - input + p)) % p;
		for (int i=0; i<m-1; i++) {
	        m_state[i] = (m_state[i+1] - m_v*m_K[i+1] + p*p) % p;
		}
		m_state[m-1] = m_v;
	}
	void saturate(int q=m) {
		assert(! is_zero(m_state));
		assert(q >= m);
		for (int i=0; i<q; ++i) {
			next();
		}
	}
	bool is_state(u16x8 state) const {
		return (state == m_state);
	}
	auto get_state() const {
		return m_state;
	}
private:
	u16x8 m_state {};
	u16x8 m_K {};
	u16 m_v {};
	u16 m_inv_K0 {};
	void m_calc_inv_K0() {
		const u16 x = m_K[0];
		assert(x != 0);
		m_inv_K0 = 1;
		while (true) {
			if (((x*m_inv_K0) % p) == 1) {
				break;
			}
			m_inv_K0++;
		}
	}
};

}
