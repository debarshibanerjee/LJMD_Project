#include "prototypes.h"
#include "structures.h"
#include "gtest/gtest.h"

class MPITestEnv : public ::testing::Environment {
   public:
	int nprocs_arg;
	int nprocs;

	explicit MPITestEnv(const std::string& cli_arg) { nprocs_arg = std::stoi(cli_arg); }

	virtual ~MPITestEnv(){};

   protected:
	virtual void SetUp() {
		char** argv;
		int argc = 0;
		int mpistatus = MPI_Init(&argc, &argv);
		ASSERT_FALSE(mpistatus);

		nprocs = get_mpi_procs();
		ASSERT_EQ(nprocs_arg, nprocs);

		int proc_id;
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	}

	virtual void TearDown() {
		int mpistatus = MPI_Finalize();
		ASSERT_FALSE(mpistatus);
	}

   public:
	static int get_mpi_procs() {
		int procs;
		MPI_Comm_size(MPI_COMM_WORLD, &procs);
		return procs;
	}

	static int get_mpi_rank() {
		int id;
		MPI_Comm_rank(MPI_COMM_WORLD, &id);
		return id;
	}
};

class MPITest : public ::testing::Test {
   protected:
	int nprocs;
	void SetUp() { nprocs = MPITestEnv::get_mpi_procs(); }
};

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
		sys->nthreads = 1;
		sys->mpicomm = MPI_COMM_WORLD;
		int mpisize, mpirank;
		sys->mpirank = MPITestEnv::get_mpi_rank();
		sys->nsize = MPITestEnv::get_mpi_procs();

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

TEST_P(ForceTest, shorter_range) {
	ASSERT_NE(sys, nullptr);
	ASSERT_DOUBLE_EQ(sys->natoms, 2);
	ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
	ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
	ASSERT_DOUBLE_EQ(sys->ry[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->ry[1], 0.0);
	ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);

	sys->rcut = 1.0;

	force_omp_simple(sys);

	EXPECT_DOUBLE_EQ(sys->fx[0], 0.0);
	EXPECT_DOUBLE_EQ(sys->fx[1], 0.0);
	EXPECT_DOUBLE_EQ(sys->fy[0], 0.0);
	EXPECT_DOUBLE_EQ(sys->fy[1], 0.0);
	EXPECT_DOUBLE_EQ(sys->fz[0], 0.0);
	EXPECT_DOUBLE_EQ(sys->fz[1], 0.0);
	EXPECT_DOUBLE_EQ(sys->epot, 0.0);
}

TEST_P(ForceTest, longer_range) {
	ASSERT_NE(sys, nullptr);
	ASSERT_DOUBLE_EQ(sys->natoms, 2);
	ASSERT_DOUBLE_EQ(sys->rx[0], -1.0);
	ASSERT_DOUBLE_EQ(sys->rx[1], 1.0);
	ASSERT_DOUBLE_EQ(sys->ry[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->ry[1], 0.0);
	ASSERT_DOUBLE_EQ(sys->rz[0], 0.0);
	ASSERT_DOUBLE_EQ(sys->rz[1], 0.0);

	sys->rcut = 4.0;

	force_omp_simple(sys);

	double exp_epot = -eps * 63.0 / 1024.0;
	double exp_ff = -eps * 93.0 / 512.0;

	EXPECT_DOUBLE_EQ(sys->fx[0], -exp_ff);
	EXPECT_DOUBLE_EQ(sys->fx[1], exp_ff);
	EXPECT_DOUBLE_EQ(sys->fy[0], 0.0);
	EXPECT_DOUBLE_EQ(sys->fy[1], 0.0);
	EXPECT_DOUBLE_EQ(sys->fz[0], 0.0);
	EXPECT_DOUBLE_EQ(sys->fz[1], 0.0);
	EXPECT_DOUBLE_EQ(sys->epot, exp_epot);
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

	return 0;
}
