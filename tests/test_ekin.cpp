#include "gtest/gtest.h"
#include "structures.h"
#include "prototypes.h"

class EkinTest : public ::testing::Test {

      protected:
        mdsys_t *sys;

        void SetUp() {
                sys = new mdsys_t;
                sys->natoms = 2;
                sys->mass = 1.0;
                sys->vx = new double[2];
                sys->vy = new double[2];
                sys->vz = new double[2];

                sys->vx[0] = 0.0;
                sys->vx[1] = 0.0;
                sys->vy[0] = 1.0;
                sys->vy[1] = -1.0;
                sys->vz[0] = 0.0;
                sys->vz[1] = 0.0;
        }

        void TearDown() {
                delete[] sys->vx;
                delete[] sys->vy;
                delete[] sys->vz;

                delete sys;
        }
};

TEST(EkinTestEmpty, empty) {
        mdsys_t *sys = new mdsys_t;
        sys->natoms = 0;

        ekin(sys);

        ASSERT_DOUBLE_EQ(sys->ekin, 0.0);
        ASSERT_DOUBLE_EQ(sys->temp, 0.0);

        delete sys;
}

TEST_F(EkinTest, test1) {
        ASSERT_NE(sys, nullptr);
        ASSERT_DOUBLE_EQ(sys->vy[0], 1.0);
        ASSERT_DOUBLE_EQ(sys->vy[1], -1.0);

        ekin(sys);
        double exp_ekin = mvsq2e;
        double exp_temp = 2 * exp_ekin / (3 * kboltz);

        ASSERT_DOUBLE_EQ(sys->ekin, exp_ekin);
        ASSERT_DOUBLE_EQ(sys->temp, exp_temp);
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
