#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <mpi.h>
/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
	MPI_Comm mpicomm;
	int nsize, mpirank, nthreads;
	int natoms, nfi, nsteps;
	double dt, mass, epsilon, sigma, box, rcut;
	double ekin, epot, temp;
	double *rx, *ry, *rz;
	double *vx, *vy, *vz;
	double *fx, *fy, *fz;
	double *cx, *cy, *cz;
};
typedef struct _mdsys mdsys_t;

#endif
