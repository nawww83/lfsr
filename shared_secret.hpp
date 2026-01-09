#pragma once

#include "types.hpp"
#include "random_utils.hpp"
#include <cmath>

namespace shared_secret_n
{

static rnd_n::GeometricDistribution<int> g_rnd(0.3);

/**
 * @brief Генератор общего секрета.
 */
template <int prime_modulo, int register_length>
class SharedSecret
{
    using _LFSR = LFSR<prime_modulo, register_length>;
private:
    /**
     * @brief Генератор g(x) стороны 1.
     */
    _LFSR mGp1;

    /**
     * @brief Генератор g(x) стороны 2.
     */
    _LFSR mGp2;

    /**
     * @brief Состояние v(x) стороны 1, которое формирует общий секрет.
     */
    STATE<register_length> mSt1;

    /**
     * @brief Состояние v(x) стороны 2, которое формирует общий секрет.
     */
    STATE<register_length> mSt2;
public:
    /**
     * @brief Конструктор с параметром.
     * @param g Коэффициенты генераторного полинома g(x).
     */
    constexpr explicit SharedSecret<prime_modulo, register_length>(STATE<register_length> g)
        : mGp1{g}
        , mGp2{g}
    {
        g_rnd.seed();
        mSt1 = rnd_n::get_random_state<prime_modulo, register_length>(g_rnd);
        mSt2 = rnd_n::get_random_state<prime_modulo, register_length>(g_rnd);
    }

    /**
     * @brief Сгенерировать общий секрет стороной 1.
     */
    STATE<register_length> GenerateSecretBySide1()
    {
        const u32 T_max = std::pow(prime_modulo, register_length) - 1; // Максимальный период генератора.
        constexpr long pm = prime_modulo * register_length;
        const long M = pm + rnd_n::get_random_long_positive(0) % (T_max - 2*pm + 1);
        const long N = pm + rnd_n::get_random_long_positive(0) % (T_max - 2*pm + 1);
        
        STATE<register_length> tmp_state = mSt2; // Временная переменная для хранения состояния.
        mGp2.set_unit_state(); // x^0.
        mGp2.next(); // x^1.
        mGp2.power_by(N); // Пара next() и power_by(q) соответствуют next() в цикле q раз. Так вычисляется x^q.
        mGp2.mult_by(tmp_state); // x^q * v2(x)

        mGp1.set_state(mGp2.get_state());
        
        mGp1.mult_by(mSt1); // *= v1(x)
        tmp_state = mGp1.get_state();
        mGp1.set_unit_state();
        mGp1.next();
        mGp1.power_by(M);
        mGp1.mult_by(tmp_state);

        mGp2.set_state(mGp1.get_state());
    
        mGp2.mult_by(mSt2);
        tmp_state = mGp2.get_state();
        mGp2.set_unit_state();
        mGp2.next();
        mGp2.power_by(T_max - N);
        mGp2.mult_by(tmp_state);
        
        mGp1.set_state(mGp2.get_state());

        mGp1.mult_by(mSt1);
        tmp_state = mGp1.get_state();
        mGp1.set_unit_state();
        mGp1.next();
        mGp1.power_by(T_max - M);
        mGp1.mult_by(tmp_state);
    
        return mGp1.get_state();
    }

    /**
     * @brief Сгенерировать общий секрет стороной 2.
     */
    STATE<register_length> GenerateSecretBySide2()
    {
        const u32 T_max = std::pow(prime_modulo, register_length) - 1; // Максимальный период генератора.
        constexpr long pm = prime_modulo * register_length;
        const long P = pm + rnd_n::get_random_long_positive(0) % (T_max - 2*pm + 1);
        const long Q = pm + rnd_n::get_random_long_positive(0) % (T_max - 2*pm + 1);

        STATE<register_length> tmp_state = mSt1;
        mGp1.set_unit_state();
        mGp1.next();
        mGp1.power_by(P);
        mGp1.mult_by(tmp_state);

        mGp2.set_state(mGp1.get_state());
        
        mGp2.mult_by(mSt2);
        tmp_state = mGp2.get_state();
        mGp2.set_unit_state();
        mGp2.next();
        mGp2.power_by(Q);
        mGp2.mult_by(tmp_state);

        mGp1.set_state(mGp2.get_state());
    
        mGp1.mult_by(mSt1);
        tmp_state = mGp1.get_state();
        mGp1.set_unit_state();
        mGp1.next();
        mGp1.power_by(T_max - P);
        mGp1.mult_by(tmp_state);
        
        mGp2.set_state(mGp1.get_state());

        mGp2.mult_by(mSt2);
        tmp_state = mGp2.get_state();
        mGp2.set_unit_state();
        mGp2.next();
        mGp2.power_by(T_max - Q);
        mGp2.mult_by(tmp_state);
    
        return mGp2.get_state();
    }
};
    
} // shared_secret_n