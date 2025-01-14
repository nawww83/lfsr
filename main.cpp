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
	
	// test_lfsr_hash_coverage_1();

	// test_lfsr_hash_coverage_2();

	// test_lfsr_hash_coverage_3();

	// test_lfsr_hash_coverage_4();

	test_next_back();
	test_state_increment();
	test_special();
	// test_total_period();
	// test_random_generators();
	test_random_generator_next_back();
	// find_lfsr_coefficients_T0_period();
	test_some_poly_1();
	test_some_poly_2();
	test_some_poly_3();
	test_some_poly_4();
	test_some_poly_5();
	test_some_poly_6();
	test_some_poly_7();
	test_some_poly_8();
    
	return 0;
}
