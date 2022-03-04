#include "prototypes.h"

const double kboltz = 0.0019872067;             /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; 	/* m*v^2 in kcal/mol */

/* velocity verlet propagation step*/
void verlet_vel_propagation(mdsys_t* sys) {
	/* first part: propagate velocities by half and positions by full step */
        const double factor = 0.5 * sys->dt / (mvsq2e * sys->mass);
        int i;
        //const double delta_t = sys->dt;
	for (i = 0; i < sys->natoms; ++i) {
		sys->vx[i] += factor * sys->fx[i];
		sys->vy[i] += factor * sys->fy[i];
		sys->vz[i] += factor * sys->fz[i];
		sys->rx[i] += sys->vx[i] * sys->dt;
		sys->ry[i] += sys->dt * sys->vy[i];
		sys->rz[i] += sys->vz[i] * sys->dt;
	}
}

/* velocity verlet update step*/
void verlet_vel_update(mdsys_t* sys) {
	/* second part: propagate velocities by another half step */
        const double factor = 0.5 * sys->dt / (mvsq2e * sys->mass);
        int i;
	for (i = 0; i < sys->natoms; ++i) {
		sys->vx[i] += sys->fx[i] * factor;
		sys->vy[i] += factor * sys->fy[i];
		sys->vz[i] += sys->fz[i] * factor;
	}
}

void velverlet(mdsys_t* sys) {

	/* initial velocity propagation */
        
        //verlet_vel_propagation(&factor, sys);
	verlet_vel_propagation(sys);

	/* compute forces and potential energy */

        force_optimized_with3LawN(sys);

        //force_optimized(sys);
     
	//force(sys);

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

