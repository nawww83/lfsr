#pragma once

#include "types.hpp"
#include "random_utils.hpp"
#include <cmath>

namespace shared_secret_n
{

/**
 * @brief Генератор общего секрета.
 */
template <int prime_modulo, int register_length>
class SharedSecret
{
    using _LFSR = LFSR<prime_modulo, register_length>;
    static constexpr long PM = prime_modulo * register_length;
private:
    /**
     * @brief Максимальный период генератора.
     */
    static const u32 T_MAX = std::pow(prime_modulo, register_length) - 1;

    /**
     * @brief Случайная степень x. Задаёт x^M mod g(x) - некоторое случайное состояние.
     */
    long M;

    /**
     * @brief LFSR генератор g(x). Должен иметь максимальный период.
     */
    _LFSR mGenerator;

    /**
     * @brief Случайное состояние v(x) стороны, которое формирует общий секрет.
     */
    const STATE<register_length> mRandomState;
public:
    /**
     * @brief Конструктор с параметром.
     * @param g Коэффициенты генераторного полинома g(x).
     */
    constexpr explicit SharedSecret<prime_modulo, register_length>(STATE<register_length> g)
        : mGenerator{g}
        , mRandomState{rnd_n::get_random_state<prime_modulo, register_length>(std::move(rnd_n::GeometricDistribution<int>{0.3}))}
    {}

    /**
     * @brief Инициализировать генератор общего секрета.
     */
    void Init()
    {
       M = PM + rnd_n::get_random_long_positive(0) % (T_MAX - 2 * PM + 1);
       mGenerator.set_unit_state();
    }

    /**
     * @brief Обновить состояние.
     * @param state Входящее состояние.
     * @details Делается путём умножения текущего состояния генератора на входящее состояние.
     */
    void SetState(STATE<register_length> state)
    {
        mGenerator.set_state(state);
    }

    /**
     * @brief Получить текущее состояние генератора.
     */
    STATE<register_length> GetState() const
    {
        return mGenerator.get_state();
    }

    /**
     * @brief Сделать шаг вперёд.
     */
    void Forward()
    {
        STATE<register_length> mTempState;
        mTempState = mGenerator.get_state();
        mGenerator.set_unit_state();
        mGenerator.power_by(M);
        mGenerator.mult_by(mTempState);
        mGenerator.mult_by(mRandomState);
    }

    /**
     * @brief Сделать шаг назад.
     */
    void Backward()
    {
        STATE<register_length> mTempState;
        mTempState = mGenerator.get_state();
        mGenerator.set_unit_state();
        mGenerator.power_by(T_MAX - M);
        mGenerator.mult_by(mTempState);
        mGenerator.mult_by(mRandomState);
    }
};
    
} // shared_secret_n