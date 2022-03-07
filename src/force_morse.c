#include "prototypes.h"
#include "utilities.h"

void force_morse(mdsys_t* sys) {
	double epot = 0.0;
	int tid;

	double rcsq;
	double coeff = 2 * sys->a_m * sys->De;

	rcsq = sys->rcut * sys->rcut;

	MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
	MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
	MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);

	//	zero energy and forces
	/* azzero(sys->fx, sys->natoms); */
	/* azzero(sys->fy, sys->natoms); */
	/* azzero(sys->fz, sys->natoms); */
	/* azzero(sys->cx, sys->natoms); */
	/* azzero(sys->cy, sys->natoms); */
	/* azzero(sys->cz, sys->natoms); */
	azzero(sys->cx, sys->nthreads * sys->natoms);
	azzero(sys->cy, sys->nthreads * sys->natoms);
	azzero(sys->cz, sys->nthreads * sys->natoms);

#ifdef _OPENMP
	#pragma omp parallel private(tid)
#endif
	{
		double rx, ry, rz;
		double *cx, *cy, *cz;
		double r, rsq, ffac;
		int grids = sys->nsize * sys->nthreads;
		int i, ii, j;

#ifdef _OPENMP
		int tid = omp_get_thread_num();
#else
		int tid = 0;
#endif

		cx = sys->cx + (tid * sys->natoms);
		cy = sys->cy + (tid * sys->natoms);
		cz = sys->cz + (tid * sys->natoms);

		for (i = 0; i < sys->natoms - 1; i += grids) {
			ii = i + sys->mpirank * sys->nthreads + tid;
			if (ii > (sys->natoms) - 1)
				break;
			for (int j = ii + 1; j < (sys->natoms); ++j) {
				/* get distance between particle i and j */
				rx = pbc(sys->rx[ii] - sys->rx[j], 0.5 * sys->box);
				ry = pbc(sys->ry[ii] - sys->ry[j], 0.5 * sys->box);
				rz = pbc(sys->rz[ii] - sys->rz[j], 0.5 * sys->box);
				rsq = rx * rx + ry * ry + rz * rz;

				/* compute force and energy if within cutoff */
				if (rsq < rcsq) {
					r = sqrt(rsq);
					double rdiff = r - sys->re;
					double a_rdiff = sys->a_m * rdiff;
					double morse_pot_tmp = exp(-a_rdiff);
					double morse_pot =
						sys->De * (morse_pot_tmp * morse_pot_tmp - 2 * morse_pot_tmp);
					double ffac =
						coeff * (exp(-2 * a_rdiff) - exp(-a_rdiff)) / rsq;
#pragma omp atomic
					epot += (morse_pot);

					cx[ii] += rx * ffac;
					cy[ii] += ry * ffac;
					cz[ii] += rz * ffac;

					cx[j] -= rx * ffac;
					cy[j] -= ry * ffac;
					cz[j] -= rz * ffac;
				}
			}
		}
#ifdef _OPENMP
	#pragma omp barrier	 // sync everything
#endif
		i = 1 + (sys->natoms / sys->nthreads);
		int fromidx = tid * i;
		int toidx = fromidx + i;
		if (toidx > sys->natoms)
			toidx = sys->natoms;
		for (i = 1; i < sys->nthreads; ++i) {
			int offs = i * sys->natoms;
			for (int j = fromidx; j < toidx; ++j) {
				sys->cx[j] += sys->cx[offs + j];
				sys->cy[j] += sys->cy[offs + j];
				sys->cz[j] += sys->cz[offs + j];
			}
		}
	}

	MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
}

