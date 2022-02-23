#include "prototypes.h"
#include "variables.h"

const double kboltz = 0.0019872067;             /* boltzman constant in kcal/mol/K */
const double mvsq2e = 2390.05736153349; /* m*v^2 in kcal/mol */

/* velocity verlet propagation step*/
void verlet_vel_propagation(mdsys_t* sys) {

	int i;
	/* first part: propagate velocities by half and positions by full step */
#ifdef _OPENMP
#	pragma omp parallel for default(shared)
#endif
	for (i = 0; i < sys->natoms; ++i) {
		sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
		sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
		sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
		sys->rx[i] += sys->dt * sys->vx[i];
		sys->ry[i] += sys->dt * sys->vy[i];
		sys->rz[i] += sys->dt * sys->vz[i];
	}
}


/* velocity verlet update step*/
void verlet_vel_update(mdsys_t* sys) {
	int i;
	/* second part: propagate velocities by another half step */
#ifdef _OPENMP
#	pragma omp parallel for default(shared)
#endif
	for (i = 0; i < sys->natoms; ++i) {
		sys->vx[i] += 0.5 * sys->dt / mvsq2e * sys->fx[i] / sys->mass;
		sys->vy[i] += 0.5 * sys->dt / mvsq2e * sys->fy[i] / sys->mass;
		sys->vz[i] += 0.5 * sys->dt / mvsq2e * sys->fz[i] / sys->mass;
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

void verlet_1(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->fx[i] / sys->mass;
        sys->rx[i] += sys->vx[i];
    }
}


void verlet_2(mdsys_t *sys)
{
    int i;

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5 * sys->fx[i] / sys->mass;
    }
}




