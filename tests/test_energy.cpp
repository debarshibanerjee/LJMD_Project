// unit test example with test fixture
#include "gtest/gtest.h"
#include "prototypes.h"
#include "structures.h"
class EnergyTest: public ::testing::Test {

protected:

    mdsys_t *sys;

    void SetUp()
    {    
	char** argv;
	int argc=0;
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->mass = 1.0;
        sys->vx = new double[2];
        sys->vy = new double[2];
        sys->vz = new double[2];
	sys->vx[0] = 1.0;
        sys->vx[1] = 1.0;
        sys->vy[0] = 0.0;
        sys->vy[1] = 2.0;
        sys->vz[0] = 1.0;
        sys->vz[1] = 0.2;
    }

    void TearDown()
        {
            delete[] sys->vx;
            delete[] sys->vy;
            delete[] sys->vz;
            delete sys;
        }

};


TEST_F(EnergyTest, step1)
{
    ASSERT_NE(sys,nullptr);
    ekin(sys);
    ASSERT_DOUBLE_EQ(sys->ekin, 8413.0019125978852);
    ASSERT_DOUBLE_EQ(sys->temp, 2822387.797772584);
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
