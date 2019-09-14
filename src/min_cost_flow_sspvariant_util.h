#pragma once

#include <iostream>


#ifdef _MSC_VER
#define NORETURN __declspec(noreturn)
#else
#define NORETURN [[noreturn]]
#endif


NORETURN void unsolvable() {
	std::cout << "The instance is unsolvable." << std::endl;
	exit(0);
}
