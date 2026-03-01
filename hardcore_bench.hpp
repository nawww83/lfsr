#include "lfsr.hpp"

#include <string>
#include <iostream>  // Для вывода результатов (std::cout)
#include <chrono>    // Для замера времени в наносекундах (std::chrono)
#include <vector>    // Для хранения семплов и сортировки (медиана)
#include <iomanip>   // Для красивого форматирования таблицы (std::setw)
#include <algorithm> // Для std::sort и std::max/min

namespace hard_bench {

int run_all();

// --- Hardcore Benchmark Tool ---
inline uint64_t rdtsc() { return __rdtsc(); }

template <typename Func>
void run_bench(const std::string &name, uint64_t iterations, Func &&func)
{
    for (int i = 0; i < 1000; ++i)
        func(); // Warmup
    std::vector<uint64_t> samples;
    auto start = std::chrono::high_resolution_clock::now();
    for (int s = 0; s < 100; ++s)
    {
        uint64_t t1 = rdtsc();
        for (uint64_t i = 0; i < iterations; ++i)
            func();
        samples.push_back((rdtsc() - t1) / iterations);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::sort(samples.begin(), samples.end());
    double avg_ns = (std::chrono::duration<double>(end - start).count() * 1e9) / (iterations * 100);
    std::cout << std::left << std::setw(12) << name << " | Median: " << std::setw(3)
              << samples[50] << " ticks | Avg: " << std::fixed << std::setprecision(2)
              << avg_ns << " ns | Speed: " << (1000.0 / avg_ns) << " Mop/s\n";
}

// Вспомогательная функция для скалярного Next (без SSE)
template <int p, int m>
void scalar_next(typename lfsr8::MType<m>::STATE &state, const typename lfsr8::MType<m>::STATE &K)
{
    using SAMPLE = typename lfsr8::MType<m>::SAMPLE;
    const SAMPLE m_v = state[m - 1];
    for (int i = m - 1; i > 0; i--)
    {
        state[i] = (state[i - 1] + m_v * K[i]) % (SAMPLE)p;
    }
    state[0] = (m_v * K[0]) % (SAMPLE)p; // упростим для теста без input
}
}