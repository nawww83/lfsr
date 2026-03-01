#include "hardcore_bench.hpp"

namespace hard_bench {


int run_all()
{
    {
        lfsr8::u32x4 K4 = {1, 3, 0, 2};
        lfsr8::LFSR<251, 4> lfsr4(K4);
        lfsr4.set_unit_state();

        lfsr8::u32x4 state_scalar = {1, 0, 0, 0};

        std::cout << "--- COMPARISON m=4 (Whiskey Lake i7-8565U) ---\n";

        // 1. Тест SSE версии (m=4)
        run_bench("SSE Next (m=4)", 1000000, [&]()
                  {
            lfsr4.next();
            asm volatile("" : : "g"(lfsr4.get_cell(0)) : "memory"); });

        // 2. Тест Скалярной версии (m=4)
        run_bench("Scalar Next (m=4)", 1000000, [&]()
                  {
            scalar_next<251, 4>(state_scalar, K4);
            asm volatile("" : : "g"(state_scalar[0]) : "memory"); });
    }
    {

        lfsr8::u16x8 K = {1, 3, 0, 2, 0, 1, 1, 4};
        lfsr8::LFSR<251, 8> lfsr(K);
        lfsr.set_unit_state();

        std::cout << "--- LFSR HARDCORE BENCHMARK (Whiskey Lake i7-8565U) ---\n";
        run_bench("NextStep", 1000000, [&]()
                  { lfsr.next(); asm volatile("" : : "g"(lfsr.get_cell(0)) : "memory"); });
        run_bench("BackStep", 1000000, [&]()
                  { lfsr.back(); asm volatile("" : : "g"(lfsr.get_cell(0)) : "memory"); });
        run_bench("Multiply", 100000, [&]()
                  { lfsr.mult_by(K); asm volatile("" : : "g"(lfsr.get_cell(0)) : "memory"); });
        run_bench("Square", 100000, [&]()
                  { lfsr.square(); asm volatile("" : : "g"(lfsr.get_cell(0)) : "memory"); });
    }

    {
        const int p = 251;
        lfsr8::u16x8 K = {1, 3, 0, 2, 1, 3, 0, 2}; // Симметричные ключи для теста

        // 1. Одиночный LFSR (m=8)
        lfsr8::LFSR<p, 8> single(K);
        single.set_unit_state();

        // 2. Парный LFSR (2x4)
        lfsr8::LFSR_paired_2x4<p> paired(K);
        paired.set_unit_state();

        std::cout << "--- PAIR VS SINGLE BENCHMARK (i7-8565U) ---\n";

        run_bench("Single LFSR Next(m=8)", 1000000, [&]()
                  {
        single.next(1);
        asm volatile("" : : "g"(single.get_cell(0)) : "memory"); });

        run_bench("Paired LFSR Next(1 i)", 1000000, [&]()
                  {
        paired.next(1); // Один вход (дублируется или маскируется)
        asm volatile("" : : "g"(paired.get_state()[0]) : "memory"); });

        run_bench("Paired LFSR Next(2 i)", 1000000, [&]()
                  {
        paired.next(1, 2); // Два разных входа
        asm volatile("" : : "g"(paired.get_state()[0]) : "memory"); });

        run_bench("Paired LFSR Back(2 i)", 1000000, [&]()
                  {
        paired.back(1, 2);
        asm volatile("" : : "g"(paired.get_state()[0]) : "memory"); });
    }

    return 0;
}
}