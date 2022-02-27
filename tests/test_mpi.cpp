#include "gtest/gtest.h"
#include <mpi.h>
#include <string>

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

TEST_F(MPITest, size_rank) {
	int me = MPITestEnv::get_mpi_rank();
	if (me == 0) {
		std::cout << "Root " << me << std::endl;
		std::cout << "Total processors running = " << nprocs << std::endl;
	}
}

TEST_F(MPITest, broadcast_check) {
	int me = MPITestEnv::get_mpi_rank();
	int root = 0;
	int buffer;
	if (me == 0) {
		buffer = 3;
	}
	MPI_Bcast(&buffer, 1, MPI_INT, root, MPI_COMM_WORLD);
	ASSERT_EQ(3, buffer);
	if (me == 0) {
		std::cout << "Original Buffer on root (rank=0) is " << 3 << std::endl;
	} else {
		std::cout << "I am rank " << me << " and the buffer I have been broadcast is = " << buffer
				  << std::endl;
	}
}

TEST_F(MPITest, reduce_check) {
	int me = MPITestEnv::get_mpi_rank();
	int sendbuf = 1;
	int recvbuf = 0;
	MPI_Reduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (me == 0) {
		ASSERT_EQ(nprocs, recvbuf);
		std::cout << "Result after Reduction = " << recvbuf << std::endl;
	} else {
		ASSERT_EQ(0, recvbuf);
	}
}

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
