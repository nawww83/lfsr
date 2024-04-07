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
	
	test_next_back();
	test_state_increment();
	test_special();
	test_total_period();
	// test_random_generators();
    
	return 0;
}
