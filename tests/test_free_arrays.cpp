#include "structures.h"
#include "utilities.h"
#include "gtest/gtest.h"

TEST(ArrayFreeingTest, test) {
	mdsys_t sys;
	sys.natoms = 2;
	allocate_sys_arrays(&sys);
	ASSERT_NE(sys.cx, nullptr);
	ASSERT_NE(sys.cy, nullptr);
	ASSERT_NE(sys.cz, nullptr);
	ASSERT_NE(sys.rx, nullptr);
	ASSERT_NE(sys.ry, nullptr);
	ASSERT_NE(sys.rz, nullptr);
	ASSERT_NE(sys.vx, nullptr);
	ASSERT_NE(sys.vy, nullptr);
	ASSERT_NE(sys.vz, nullptr);
	ASSERT_NE(sys.fx, nullptr);
	ASSERT_NE(sys.fy, nullptr);
	ASSERT_NE(sys.fz, nullptr);
	free_sys_arrays(&sys);
	ASSERT_EQ(sys.cx, nullptr);
	ASSERT_EQ(sys.cy, nullptr);
	ASSERT_EQ(sys.cz, nullptr);
	ASSERT_EQ(sys.rx, nullptr);
	ASSERT_EQ(sys.ry, nullptr);
	ASSERT_EQ(sys.rz, nullptr);
	ASSERT_EQ(sys.vx, nullptr);
	ASSERT_EQ(sys.vy, nullptr);
	ASSERT_EQ(sys.vz, nullptr);
	ASSERT_EQ(sys.fx, nullptr);
	ASSERT_EQ(sys.fy, nullptr);
	ASSERT_EQ(sys.fz, nullptr);
}
