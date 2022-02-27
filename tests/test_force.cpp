#include "prototypes.h"
#include "structures.h"
#include "gtest/gtest.h"

class ForceTest : public ::testing::TestWithParam<double> {
   protected:
	mdsys_t* sys;
	int eps = GetParam();

	void SetUp() {
		sys = new mdsys_t;
		sys->natoms = 2;
		sys->epsilon = eps;
		sys->sigma = 1.0;
		sys->box = 10.0;

		sys->mpicomm = MPI_COMM_WORLD;

		sys->rx = new double[2]();
		sys->ry = new double[2]();
		sys->rz = new double[2]();

		sys->fx = new double[2];
		sys->fy = new double[2];
		sys->fz = new double[2];

		sys->cx = new double[2];
		sys->cy = new double[2];
		sys->cz = new double[2];

		sys->rx[0] = -1.0;
		sys->rx[1] = 1.0;
	}

	void TearDown() {
		delete[] sys->rx;
		delete[] sys->ry;
		delete[] sys->rz;

		delete[] sys->fx;
		delete[] sys->fy;
		delete[] sys->fz;

		delete[] sys->cx;
		delete[] sys->cy;
		delete[] sys->cz;

		delete sys;
	}
};

TEST_P(ForceTest, test1) {
	ASSERT_NE(sys, nullptr);
	ASSERT_DOUBLE_EQ(sys->natoms, 2);
	ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
	ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
	ASSERT_DOUBLE_EQ(sys->ry[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->ry[1], 0.0);
	ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);

	sys->rcut = 0.5;

	force(sys);

	EXPECT_DOUBLE_EQ(sys->fx[0], 0.0);
	EXPECT_DOUBLE_EQ(sys->fx[1], 0.0);
	EXPECT_DOUBLE_EQ(sys->fy[0], 0.0);
	EXPECT_DOUBLE_EQ(sys->fy[1], 0.0);
	EXPECT_DOUBLE_EQ(sys->fz[0], 0.0);
	EXPECT_DOUBLE_EQ(sys->fz[1], 0.0);
	EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}

INSTANTIATE_TEST_SUITE_P(ForceTest_parametric, ForceTest, ::testing::Values(0.0, 0.5, 1.0));

int main(int argc, char* argv[]) {
	int result = 0, size, rank;

	::testing::InitGoogleTest(&argc, argv);
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
	if (rank != 0) {
		delete listeners.Release(listeners.default_result_printer());
	}
	result = RUN_ALL_TESTS();
	MPI_Finalize();

	return result;
}
