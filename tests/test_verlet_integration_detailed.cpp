#include "prototypes.h"
#include "structures.h"
#include "gtest/gtest.h"

class Verlet_Integration_Test_Detailed
	: public ::testing::TestWithParam<std::tuple<double, double>> {
   protected:
	mdsys_t* sys;

	void SetUp() {
		sys = new mdsys_t;
		sys->natoms = 2;
		sys->mass = 1.0;

		sys->rx = new double[2]();
		sys->ry = new double[2]();
		sys->rz = new double[2]();

		sys->vx = new double[2]();
		sys->vy = new double[2]();
		sys->vz = new double[2]();

		sys->fx = new double[2]();
		sys->fy = new double[2]();
		sys->fz = new double[2]();
	}

	void TearDown() {
		delete[] sys->rx;
		delete[] sys->ry;
		delete[] sys->rz;

		delete[] sys->vx;
		delete[] sys->vy;
		delete[] sys->vz;

		delete[] sys->fx;
		delete[] sys->fy;
		delete[] sys->fz;

		delete sys;
	}
};

TEST_P(Verlet_Integration_Test_Detailed, test1) {
	ASSERT_NE(sys, nullptr);

	sys->dt = std::get<0>(GetParam());

	sys->ry[0] = 1.0;
	sys->ry[1] = -1.0;

	sys->fx[0] = std::get<1>(GetParam());
	sys->fx[1] = std::get<1>(GetParam());

	verlet_vel_propagation(sys);

	double v_inc = 0.5 * std::get<0>(GetParam()) * std::get<1>(GetParam()) / mvsq2e;
	double r_inc = v_inc * std::get<0>(GetParam());

	ASSERT_DOUBLE_EQ(sys->rx[0], r_inc);
	ASSERT_DOUBLE_EQ(sys->ry[0], 1.0);
	ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->vx[0], v_inc);
	ASSERT_DOUBLE_EQ(sys->vy[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->vz[0], 0.0);

	ASSERT_DOUBLE_EQ(sys->rx[1], r_inc);
	ASSERT_DOUBLE_EQ(sys->ry[1], -1.0);
	ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);
	ASSERT_DOUBLE_EQ(sys->vx[1], v_inc);
	ASSERT_DOUBLE_EQ(sys->vy[1], 0.0);
	ASSERT_DOUBLE_EQ(sys->vz[1], 0.0);
}

TEST_P(Verlet_Integration_Test_Detailed, test2) {
	ASSERT_NE(sys, nullptr);

	sys->dt = std::get<0>(GetParam());

	sys->ry[0] = 1.0;
	sys->ry[1] = -1.0;

	sys->fx[0] = std::get<1>(GetParam());
	sys->fx[1] = std::get<1>(GetParam());

	verlet_vel_update(sys);

	double v_inc = 0.5 * std::get<0>(GetParam()) * std::get<1>(GetParam()) / mvsq2e;

	ASSERT_DOUBLE_EQ(sys->rx[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->ry[0], 1.0);
	ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->vx[0], v_inc);
	ASSERT_DOUBLE_EQ(sys->vy[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->vz[0], 0.0);

	ASSERT_DOUBLE_EQ(sys->rx[1], 0.0);
	ASSERT_DOUBLE_EQ(sys->ry[1], -1.0);
	ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);
	ASSERT_DOUBLE_EQ(sys->vx[1], v_inc);
	ASSERT_DOUBLE_EQ(sys->vy[1], 0.0);
	ASSERT_DOUBLE_EQ(sys->vz[1], 0.0);
}

INSTANTIATE_TEST_SUITE_P(Verlet_Integration_Test_Detailed_parametric,
						 Verlet_Integration_Test_Detailed,
						 ::testing::Combine(::testing::Values(0.0, 0.5, 1.0),
											::testing::Values(0.0, 1.0, 2.0)));

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
