#pragma once

#include "types.hpp"
#include "random_utils.hpp"
#include <cmath>

namespace shared_secret_n
{

/**
 * @brief Вспомогательный класс для двухфазного трехрауднового алгоритма генерации общего секрета.
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
     * @brief Случайное ненулевое состояние-маска для одной фазы алгоритма генерации общего секрета.
     */
    STATE<register_length> mMaskState;

    /**
     * @brief Обратное состояние-маска.
     */
    STATE<register_length> mInverseMaskState;

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
        , mRandomState{rnd_n::get_random_state<prime_modulo, register_length>( std::move(rnd_n::GeometricDistribution<int>{0.3}) )}
    {}

    /**
     * @brief Инициализировать генератор общего секрета перед началом новой фазы.
     */
    void Init()
    {
        auto g = rnd_n::GeometricDistribution<int>{0.3};
        for (;;) 
        {
            mMaskState = rnd_n::get_random_state<prime_modulo, register_length>( g );
            if (mMaskState != mRandomState) break;
        }
        mInverseMaskState = _LFSR::inverse_of(mMaskState, mGenerator.get_generator_coeffs());
        mGenerator.set_unit_state();
    }

    /**
     * @brief Установить состояние.
     * @param state Входящее состояние.
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
        mGenerator.mult_by(mMaskState);
        mGenerator.mult_by(mRandomState);
    }

    /**
     * @brief Сделать шаг назад.
     */
    void Backward()
    {
        mGenerator.mult_by(mInverseMaskState);
        mGenerator.mult_by(mRandomState);
    }
};
    
} // shared_secret_n