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

	test_bias16bit();
	// test_next_back();
	// test_state_increment();
	// test_special();
	// test_total_period();
	// test_random_generators();
	// test_random_generator_next_back();
	// find_lfsr_coefficients_T0_period();
	
	// test_some_poly<19, 4>(STATE<4>{9, 5, 2, 0});
	// test_some_poly<19, 4>(STATE<4>{4, 2, 2, 6});
	// test_some_poly<17, 4>(STATE<4>{3, 4, 2, 1});
	// test_some_poly<17, 4>(STATE<4>{6, 1, 2, 1});
	// test_some_poly<17, 4>(STATE<4>{3, 2, 3, 4});
	// test_some_poly<17, 4>(STATE<4>{6, 2, 0, 7});
	// test_some_poly<13, 4>(STATE<4>{2, 3, 1, 1});
	// test_some_poly<13, 4>(STATE<4>{2, 0, 1, 7});

	// {
	// 	u64 T;
	// 	LFSR<2, 8> g ({1, 0, 0, 0, 1, 1, 0, 1});
	// 	assert(!is_maximal_period(g, T));
	// 	assert(T == 51);
	// }
	// {
	// 	u64 T;
	// 	LFSR<2, 8> g ({1, 0, 1, 1, 0, 1, 0, 0});
	// 	assert(is_maximal_period(g, T));
	// 	assert(T == 255);
	// }

	// test_square_of_lfsr();
	// test_mult_by_of_lfsr();
    
	return 0;
}
