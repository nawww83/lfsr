#pragma once

#include "lfsr.hpp"

#include <string>
#include <iostream>  // Для вывода результатов (std::cout)
#include <chrono>    // Для замера времени в наносекундах (std::chrono)
#include <vector>    // Для хранения семплов и сортировки (медиана)
#include <iomanip>   // Для красивого форматирования таблицы (std::setw)
#include <algorithm> // Для std::sort и std::max/min

// --- Hardcore Benchmark Tool ---
#include <stdint.h>
#include <thread>

#ifdef _WIN32
#include <windows.h>
#include <intrin.h>
#else
#include <x86intrin.h>
#include <sched.h>
#endif

#include <cpuid.h> // Для GCC/Clang. В MSVC используется __cpuid

namespace hard_bench
{

inline bool has_invariant_tsc()
{
    unsigned int eax, ebx, ecx, edx;

    // Сначала проверяем, поддерживает ли CPU расширенные функции (80000007h)
    if (__get_cpuid(0x80000000, &eax, &ebx, &ecx, &edx) && eax >= 0x80000007)
    {
        // Вызываем функцию Advanced Power Management
        __get_cpuid(0x80000007, &eax, &ebx, &ecx, &edx);
        // Проверяем 8-й бит EDX (Invariant TSC)
        return (edx & (1 << 8)) != 0;
    }
    return false;
}

class CycleTimer
{
    double ghz;

public:
    CycleTimer()
    {
        fix_to_core();
        ghz = calibrate();
    }

    // Привязка к первому доступному ядру
    void fix_to_core()
    {
#ifdef _WIN32
        SetThreadAffinityMask(GetCurrentThread(), 1);
#else
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(0, &cpuset);
        sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
#endif
    }

    // Автоматическое определение частоты Invariant TSC
    double calibrate()
    {
        auto t_start = std::chrono::steady_clock::now();
        uint64_t c_start = __rdtsc();

        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        uint64_t c_end = __rdtsc();
        auto t_end = std::chrono::steady_clock::now();

        std::chrono::duration<double> elapsed = t_end - t_start;
        return (c_end - c_start) / elapsed.count() / 1e9;
    }

    inline uint64_t start() const
    {
        _mm_lfence();
        uint64_t tsc = __rdtsc();
        _mm_lfence();
        return tsc;
    }

    inline uint64_t stop(unsigned int &cpu_id) const
    {
        // rdtscp гарантирует выполнение всего кода ДО неё
        uint64_t tsc = __rdtscp(&cpu_id);
        _mm_lfence();
        return tsc;
    }

    double to_ms(uint64_t cycles) const { return (double)cycles / (ghz * 1e6); }
    double to_us(uint64_t cycles) const { return (double)cycles / (ghz * 1e3); }
    double to_ns(uint64_t cycles) const { return (double)cycles / (ghz * 1); }
    double get_ghz() const { return ghz; }
};

template <typename Func>
inline void run_bench(const std::string &name, uint64_t iterations, Func &&func)
{
    CycleTimer timer;
    unsigned int cpu_start, cpu_end;

    // std::cout << "Detected TSC Frequency: " << timer.get_ghz() << " GHz\n";

    // Замер
    __rdtscp(&cpu_start); // Берем ID ядра в начале

    double ave_delay_ns = 0.;
    std::vector<uint64_t> samples;
    samples.reserve(100);
    for (int s = 0; s < 100; ++s)
    {
        uint64_t t1 = timer.start();
        for (uint64_t i = 0; i < iterations; ++i)
        {
            func();
#ifdef _MSC_VER
            _ReadWriteBarrier();
#else
            asm volatile("" : : : "memory");
#endif
        }
        uint64_t t2 = timer.stop(cpu_end);
        samples.push_back(uint64_t(double(t2 - t1) / iterations + 0.5));
        ave_delay_ns += timer.to_ns(t2 - t1);
    }

    ave_delay_ns /= (100. * iterations);

    if (cpu_start != cpu_end)
        std::cout << "Warning: Context switch detected!\n";

    std::sort(samples.begin(), samples.end());
    std::cout << std::left << std::setw(12) << name << " | Median: " << std::setw(3)
              << samples[50] << " ticks | Avg: " << std::fixed << std::setprecision(2)
              << ave_delay_ns << " ns | Speed: " << (1000.0 / ave_delay_ns) << " Mop/s\n";
}

void run_all();

}
