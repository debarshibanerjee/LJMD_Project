#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <mpi.h>

/* structure to organize atoms in discrete cells */
struct _cell {
	int natoms;
	int x_pos, y_pos, z_pos;
	int* idxlist;
};
typedef struct _cell cell_t;

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
	MPI_Comm mpicomm;
	int nsize, mpirank, nthreads;
	int *plist, *cell_counter;
	cell_t* clist;
	int nfi, nsteps, natoms, ncells, ncel_d, npairs;
	double dt, mass, epsilon, sigma, box, rcut;
	double ekin, epot, temp;
	double *rx, *ry, *rz;
	double *vx, *vy, *vz;
	double *fx, *fy, *fz;
	double *cx, *cy, *cz;
	double De, a_m, re;
};
typedef struct _mdsys mdsys_t;

#endif
