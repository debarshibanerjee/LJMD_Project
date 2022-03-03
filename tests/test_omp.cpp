#include "omp.h"
#include "utilities.h"
#include "gtest/gtest.h"

TEST(omp_test, test1) {
#pragma omp parallel
	printf("Number of OpenMP Threads: %d\n", omp_get_num_threads());
}
