#include "prototypes.h"
#include "gtest/gtest.h"

TEST(test_output, output) {
	mdsys_t sys;
	FILE *erg, *traj;

	sys.natoms = 23;
	sys.nfi = 137;
	sys.ekin = 153.371;
	sys.epot = 371.153;
	sys.temp = 273.15;

	sys.rx = (double*)malloc(sys.natoms * sizeof(double));
	sys.ry = (double*)malloc(sys.natoms * sizeof(double));
	sys.rz = (double*)malloc(sys.natoms * sizeof(double));

	for (int i = 0; i < sys.natoms; ++i) {
		sys.rx[i] = (double)i / 100. + 1.0;
		sys.ry[i] = (double)i / 100. + 2.0;
		sys.rz[i] = (double)i / 100. + 3.0;
	}

	erg = fopen("ergfile.txt", "w");
	traj = fopen("trajfile.txt", "w");

	output(&sys, erg, traj);

	fclose(erg);
	fclose(traj);

	free(sys.rx);
	free(sys.ry);
	free(sys.rz);

	// Verify Energy File
	double etot;
	mdsys_t sys_1;
	erg = fopen("ergfile.txt", "r");
	fscanf(erg, "%d %lf %lf %lf %lf", &(sys_1.nfi), &(sys_1.temp), &(sys_1.ekin),
		   &(sys_1.epot), &etot);
	fclose(erg);

	ASSERT_EQ(sys_1.nfi, 137);
	ASSERT_DOUBLE_EQ(sys_1.temp, 273.15);
	ASSERT_DOUBLE_EQ(sys_1.ekin, 153.371);
	ASSERT_DOUBLE_EQ(sys_1.epot, 371.153);
	ASSERT_DOUBLE_EQ(etot, sys_1.ekin + sys_1.epot);
	ASSERT_DOUBLE_EQ(etot, sys.ekin + sys.epot);

	// Verify Trajectory File
	mdsys_t sys_2;
	traj = fopen("trajfile.txt", "r");
	fscanf(traj, "%d\n nfi=%d etot=%lf\n", &(sys_2.natoms), &(sys_2.nfi), &etot);

	ASSERT_EQ(sys_2.natoms, 23);
	ASSERT_EQ(sys_2.nfi, 137);
	ASSERT_DOUBLE_EQ(etot, sys.ekin + sys.epot);

	sys_2.rx = (double*)malloc(sys_2.natoms * sizeof(double));
	sys_2.ry = (double*)malloc(sys_2.natoms * sizeof(double));
	sys_2.rz = (double*)malloc(sys_2.natoms * sizeof(double));

	for (int i = 0; i < sys.natoms; ++i) {
		fscanf(traj, "Ar %lf %lf %lf\n", &(sys_2.rx[i]), &(sys_2.ry[i]),
			   &(sys_2.rz[i]));

		ASSERT_DOUBLE_EQ(sys_2.rx[i], (double)i / 100. + 1.0);
		ASSERT_DOUBLE_EQ(sys_2.ry[i], (double)i / 100. + 2.0);
		ASSERT_DOUBLE_EQ(sys_2.rz[i], (double)i / 100. + 3.0);
	}

	fclose(traj);

	free(sys_2.rx);
	free(sys_2.ry);
	free(sys_2.rz);
}

