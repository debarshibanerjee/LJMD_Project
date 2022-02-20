#include "prototypes.h"

/* compute forces */
void force(mdsys_t* sys) {
	double rsq, ffac;
	double rx, ry, rz;
	int i, j;

	/* zero energy and forces */
	sys->epot = 0.0;
	azzero(sys->fx, sys->natoms);
	azzero(sys->fy, sys->natoms);
	azzero(sys->fz, sys->natoms);

	double c12 = 4.0 * sys->epsilon * pow(sys->sigma, 12.0);
	double c6 = 4.0 * sys->epsilon * pow(sys->sigma, 6.0);
	double rcsq = sys->rcut * sys->rcut;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, j, rx, ry, rz, ffac, r1x, r1y, r1z, f1x, f1y, f1z, rsq, rinv, r6) reduction(+:epot)
#endif
	for (i = 0; i < (sys->natoms); ++i) {
		for (j = 0; j < (sys->natoms); ++j) {
			/* particles have no interactions with themselves */
			if (i == j)
				continue;

			/* get distance between particle i and j */
			rx = pbc(sys->rx[i] - sys->rx[j], 0.5 * sys->box);
			ry = pbc(sys->ry[i] - sys->ry[j], 0.5 * sys->box);
			rz = pbc(sys->rz[i] - sys->rz[j], 0.5 * sys->box);
			rsq = (rx * rx + ry * ry + rz * rz);

			/* compute force and energy if within cutoff */
			if (r < sys->rcut) {
				ffac = -4.0 * sys->epsilon *
					   (-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

				sys->epot += 0.5 * 4.0 * sys->epsilon *
							 (pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

				sys->fx[i] += rx / r * ffac;
				sys->fy[i] += ry / r * ffac;
				sys->fz[i] += rz / r * ffac;
			}
		}
	}
}
