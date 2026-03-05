#include <iostream>
#include <cassert>

#include "io_utils.hpp"
#include "tests.hpp"
#include "hardcore_bench.hpp"
#include "types.hpp"

int main() {
	using std::cout;
	using std::endl;

	io_u::io_utils io;

	cout << "Welcome!\n";
    cout << "Endianess: " << (io.is_little_endian() ? "LE\n" : "Not LE\n");
    cout << "Endianess: " << (io.is_big_endian() ? "BE\n" : "Not BE\n");
    assert( io.is_little_endian() ^ io.is_big_endian() );

	using namespace tests;
	
	// {
	// 	std::cout << "Wait for direct LFSR periods calculation...\n";
	// 	STATE<4> K11 = {3, 2, 2, 0}; 
	// 	STATE<4> K21 = {7, 0, 4, 3};
	// 	LFSR<251, 4> lfsr11{K11};
	// 	LFSR<241, 4> lfsr21{K21};

	// 	STATE<4> K12 = {9, 1, 0, 4};
	// 	STATE<4> K22 = {7, 3, 2, 0};
	// 	LFSR<251, 4> lfsr12{K12};
	// 	LFSR<241, 4> lfsr22{K22};
	// 	{
	// 		const u64 T_max = lfsr8::safe_ipow<u64>(251, 4) - 1;
	// 		u64 T = 1;
	// 		lfsr11.set_unit_state();
	// 		lfsr12.set_unit_state();
	// 		for (;;)
	// 		{
	// 			lfsr11.next();
	// 			lfsr12.next();
	// 			if (lfsr11.is_state({1}) || lfsr12.is_state({1})) break;
	// 			if (T > T_max) break;
	// 			T++;
	// 		}
	// 		std::cout << "T = " << T << '\n';
	// 		assert(T == T_max);

	// 	}
	// 	{
	// 		const u64 T_max = lfsr8::safe_ipow<u64>(241, 4) - 1;
	// 		u64 T = 1;
	// 		lfsr21.set_unit_state();
	// 		lfsr22.set_unit_state();
	// 		for (;;)
	// 		{
	// 			lfsr21.next();
	// 			lfsr22.next();
	// 			if (lfsr21.is_state({1}) || lfsr22.is_state({1})) break;
	// 			if (T > T_max) break;
	// 			T++;
	// 		}
	// 		std::cout << "T = " << T << '\n';
	// 		assert(T == T_max);

	// 	}
	// 	std::cout << "Ok\n";
	// }
	test_lfsr_hash_coverage_1();
	// test_lfsr_hash_coverage_2(); // You should have >= 16 GB memory to run this.
	test_lfsr_hash_coverage_3();
	test_lfsr_hash_coverage_4();

	// test_bias8bit();
	// test_bias16bit();
	test_next_back();
	// test_state_increment();
	// test_special();
	test_total_period();
	
	// test_random_generators();

	test_random_generator_next_back();
	find_lfsr_coefficients_T0_period();

	test_square_of_lfsr();
	test_mult_by_of_lfsr();

	test_power_of_lfsr();

	// test_lfsr_hash_bench();

	// test_debug();

	hard_bench::run_all();
    
	return 0;
}
