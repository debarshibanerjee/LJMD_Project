// unit test example with test fixture
#include "gtest/gtest.h"
#include "structures.h"
#include "prototypes.h"

class VerletTest: public ::testing::Test {

protected:

    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->mass = 1.0;
        sys->rx = new double[2];
        sys->vx = new double[2];
        sys->fx = new double[2];
        sys->rx[0] = -1.0;
        sys->rx[1] = 1.0;
        sys->vx[0] = 0.0;
        sys->vx[1] = 0.0;
        sys->fx[0] = 1.0;
        sys->fx[1] = 0.2;
    }

    void TearDown()
        {
            delete[] sys->rx;
            delete[] sys->vx;
            delete[] sys->fx;

            delete sys;
        }
};

TEST_F(VerletTest, step1)
{
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0],-1.0);
    ASSERT_DOUBLE_EQ(sys->vx[0],0.0);
    verlet_1(sys);
    ASSERT_DOUBLE_EQ(sys->rx[0],-0.5);
    ASSERT_DOUBLE_EQ(sys->vx[0], 0.5);
}

TEST_F(VerletTest, step2)
{
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0],-1.0);
    ASSERT_DOUBLE_EQ(sys->vx[0],0.0);
    verlet_2(sys);
    ASSERT_DOUBLE_EQ(sys->rx[0],-1.0);
    ASSERT_DOUBLE_EQ(sys->vx[0], 0.5);
}


int main(int argc, char* argv[]) {
    int result = 0,size,rank;

    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
    if (rank != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }
    result = RUN_ALL_TESTS();
    MPI_Finalize();

    return result;
}
