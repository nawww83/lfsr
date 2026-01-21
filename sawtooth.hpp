#pragma once

#include <array>
#include <cstddef>

namespace saw_n
{
    template <typename T, size_t N>
    inline void sawtooth_forward(std::array<T, N> &v, const std::array<int, N> &p)
    {
        for (size_t i = 0; auto &el : v)
        {
            el++;
            el %= static_cast<T>(p.at(i++));
        }
    }

    template <typename T, size_t N>
    inline void sawtooth_backward(std::array<T, N> &v, const std::array<int, N> &p)
    {
        for (size_t i = 0; auto &el : v)
        {
            el += static_cast<T>(p.at(i++));
        }
        for (size_t i = 0; auto &el : v)
        {
            el--;
            el %= static_cast<T>(p.at(i++));
        }
    }
}