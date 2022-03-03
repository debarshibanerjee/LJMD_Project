#include "prototypes.h"

const double kboltz = 0.0019872067;             /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; 	/* m*v^2 in kcal/mol */
//const double inv_mvsq2e= 0.00041839999;
/* velocity verlet propagation step*/
void verlet_vel_propagation(mdsys_t* sys) {
	int i;
	double inv_mass = 1/sys->mass;
	double inv_mvsq2e = 1/mvsq2e;
	/* first part: propagate velocities by half and positions by full step */
#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(i)
#endif
	for (i = 0; i < sys->natoms; ++i) {
		sys->vx[i] += 0.5 * sys->dt * inv_mvsq2e * sys->fx[i] * inv_mass; //some extra changes
		sys->vy[i] += 0.5 * sys->dt * inv_mvsq2e * sys->fy[i] * inv_mass;
		sys->vz[i] += 0.5 * sys->dt * inv_mvsq2e * sys->fz[i] * inv_mass;
		sys->rx[i] += sys->dt * sys->vx[i];
		sys->ry[i] += sys->dt * sys->vy[i];
		sys->rz[i] += sys->dt * sys->vz[i];
	}
}

/* velocity verlet update step*/
void verlet_vel_update(mdsys_t* sys) {
	int i;
	double inv_mass = 1/sys->mass;
	double inv_mvsq2e = 1/mvsq2e;
	/* second part: propagate velocities by another half step */
#ifdef _OPENMP
	#pragma omp parallel for default(shared) private(i)
#endif
	for (i = 0; i < sys->natoms; ++i) {
		sys->vx[i] += 0.5 * sys->dt * inv_mvsq2e * sys->fx[i] * inv_mass;
		sys->vy[i] += 0.5 * sys->dt * inv_mvsq2e * sys->fy[i] * inv_mass;
		sys->vz[i] += 0.5 * sys->dt * inv_mvsq2e * sys->fz[i] * inv_mass;
	}
}

void velverlet(mdsys_t* sys) {
	/* initial velocity propagation */
	verlet_vel_propagation(sys);

	/* compute forces and potential energy */
	force(sys);

	/* update the velocities */
	verlet_vel_update(sys);
}

void verlet_test_1(mdsys_t* sys) {
	int i;

	for (i = 0; i < sys->natoms; ++i) {
		sys->vx[i] += 0.5 * sys->fx[i] / sys->mass;
		sys->rx[i] += sys->vx[i];
	}
}

void verlet_test_2(mdsys_t* sys) {
	int i;
	for (i = 0; i < sys->natoms; ++i) {
		sys->vx[i] += 0.5 * sys->fx[i] / sys->mass;
	}
}

