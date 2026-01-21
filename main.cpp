#include <iostream>
#include <cassert>

#include "io_utils.hpp"
#include "tests.hpp"

int main() {
	using std::cout;
	using std::endl;

	io_u::io_utils io;

	cout << "Welcome!\n";
    cout << "Endianess: " << (io.is_little_endian() ? "LE\n" : "Not LE\n");
    cout << "Endianess: " << (io.is_big_endian() ? "BE\n" : "Not BE\n");
    assert( io.is_little_endian() ^ io.is_big_endian() );

	using namespace tests;

	test_shared_secret_generation();
	
	// test_lfsr_hash_coverage_1();
	// test_lfsr_hash_coverage_2(); // You should have >= 16 GB memory to run this.
	// test_lfsr_hash_coverage_3();
	// test_lfsr_hash_coverage_4();

	// test_bias8bit();
	// test_bias16bit();
	// test_next_back();
	// test_state_increment();
	// test_special();
	// test_total_period();
	
	// test_random_generators();

	// test_random_generator_next_back();
	// find_lfsr_coefficients_T0_period();

	// test_square_of_lfsr();
	// test_mult_by_of_lfsr();

	// test_power_of_lfsr();
    
	return 0;
}
