#include "hardcore_bench.hpp"

namespace hard_bench
{
    void run_all()
    {
        std::cout << "Run the hardcore benchmark\n...";
        // if (has_invariant_tsc())
        // {
        //     std::cout << "TSC стабилен (Invariant). Замеры будут точными!\n";
        // }
        // else
        // {
        //     std::cout << "TSC плавает. Использовать RDTSC для замера времени опасно.\n";
        // }
        // {
        //     run_bench("Pure XOR", 1000000, [&]()
        //               {
        //         static volatile int x = 1;
        //         x ^= 0xDEADBEEF; });
        //
        // }
        {
            int p = 65521;
            lfsr8::u32x4 K4 = {65521 - 7, 65521 - 1, 0, 0};
            using STATE = typename lfsr8::MType<4>::STATE;
            STATE state{1};

            run_bench("Local Var(m=4)", 1000000, [&]()
                      { 
            using SAMPLE = typename lfsr8::MType<4>::SAMPLE;
            const SAMPLE m_v = state[4 - 1];
            for (int i = 4 - 1; i > 0; i--)
            {
                state[i] = (state[i - 1] + m_v * K4[i]) % (SAMPLE)p;
            }
            state[0] = (0 + m_v * K4[0]) % (SAMPLE)p;
            asm volatile("" : : "g"(&state) : ); });
        }
        {
            lfsr8::u32x4 K4 = {65521 - 7, 65521 - 1, 0, 0};
            lfsr8::LFSR<65521, 4> lfsr4(K4);
            lfsr4.set_unit_state();

            std::cout << "--- COMPARISON m=4 ---\n";

            // 1. Тест Next (m=4)
            run_bench("Next (m=4)", 1000000, [&]()
                      {
            lfsr4.next();
            asm volatile("" : : "g"(lfsr4.get_cell(0)) : ); });

            // 2. Тест Back (m=4)
            run_bench("Back (m=4)", 1000000, [&]()
                      {
            lfsr4.back();
            asm volatile("" : : "g"(lfsr4.get_cell(0)) : ); });
        }
        {
            lfsr8::u16x8 K8 = {1, 3, 0, 2, 0, 1, 1, 4};
            lfsr8::LFSR<251, 8> lfsr8(K8);
            lfsr8.set_unit_state();

            std::cout << "--- COMPARISON m=8 ---\n";

            // 1. Тест Next (m=8)
            run_bench("Next (m=8)", 1000000, [&]()
                      {
            lfsr8.next();
            asm volatile("" : : "g"(lfsr8.get_cell(0)) : ); });

            // 2. Тест Back (m=8)
            run_bench("Back (m=8)", 1000000, [&]()
                      {
            lfsr8.back();
            asm volatile("" : : "g"(lfsr8.get_cell(0)) : ); });
        }
        {
            const int p = 251;
            lfsr8::u16x8 K = {1, 3, 0, 2, 1, 3, 0, 2}; // Симметричные ключи для теста

            // Парный LFSR (2x4)
            lfsr8::LFSR_paired_2x4<p> paired(K);
            paired.set_unit_state();

            std::cout << "--- PAIR BENCHMARK ---\n";

            run_bench("Paired LFSR Next(i, j)", 1000000, [&]()
                      {
        paired.next(1, 2); // Два разных входа
        asm volatile("" : : "g"(paired.get_state()[0]) : ); });

            run_bench("Paired LFSR Back(i, j)", 1000000, [&]()
                      {
        paired.back(1, 2);
        asm volatile("" : : "g"(paired.get_state()[0]) : ); });
        }
        {
            lfsr8::u32x4 K4 = {65521 - 7, 65521 - 1, 0, 0};
            lfsr8::LFSR<65521, 4> lfsr4(K4);
            lfsr4.set_unit_state();

            std::cout << "--- LFSR (m=4) Multiply BENCHMARK ---\n";
            run_bench("Multiply", 100000, [&]()
                      { lfsr4.mult_by(K4); asm volatile("" : : "g"(lfsr4.get_cell(0)) : ); });
            run_bench("Square", 100000, [&]()
                      { lfsr4.square(); asm volatile("" : : "g"(lfsr4.get_cell(0)) : ); });
        }
        {
            lfsr8::u16x8 K8 = {1, 3, 0, 2, 0, 1, 1, 4};
            lfsr8::LFSR<251, 8> lfsr8(K8);
            lfsr8.set_unit_state();

            std::cout << "--- LFSR (m=8) Multiply BENCHMARK ---\n";
            run_bench("Multiply", 100000, [&]()
                      { lfsr8.mult_by(K8); asm volatile("" : : "g"(lfsr8.get_cell(0)) : ); });
            run_bench("Square", 100000, [&]()
                      { lfsr8.square(); asm volatile("" : : "g"(lfsr8.get_cell(0)) : ); });
        }
    }
}