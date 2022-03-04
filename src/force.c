#include "prototypes.h"
#include "utilities.h"

void force(mdsys_t* sys) {
	double epot = 0.0;
	int tid;

	double sigma, sigma6;
	double c6, c12, rcsq;

	sigma = sys->sigma;
	sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;

	c6 = 4.0 * sys->epsilon * sigma6;
	c12 = 4.0 * sys->epsilon * sigma6 * sigma6;
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
		double rsq, ffac;
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
					double r6, rsqinv;
					rsqinv = 1.0 / rsq;
					r6 = rsqinv * rsqinv * rsqinv;
					ffac = 12.0 * c12 * r6 * r6 * rsqinv - 6.0 * c6 * r6 * rsqinv;
#pragma omp atomic
					epot += (c12 * r6 * r6 - c6 * r6);
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

void force_omp_simple(mdsys_t* sys) {
	double epot = 0.0;
	double sigma = sys->sigma;
	double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
	double c6 = 4.0 * sys->epsilon * sigma6;
	double c12 = 4.0 * sys->epsilon * sigma6 * sigma6;
	double rcsq = sys->rcut * sys->rcut;

	/* zero energy and forces */
	sys->epot = 0.0;
	azzero(sys->fx, sys->natoms);
	azzero(sys->fy, sys->natoms);
	azzero(sys->fz, sys->natoms);
	azzero(sys->cx, sys->natoms);
	azzero(sys->cy, sys->natoms);
	azzero(sys->cz, sys->natoms);

	MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
	MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
	MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);

	double* rx = sys->rx;
	double* ry = sys->ry;
	double* rz = sys->rz;

	double* cx = sys->cx;
	double* cy = sys->cy;
	double* cz = sys->cz;

#ifdef _OPENMP
	#pragma omp parallel for shared(sys, rx, ry, rz, rcsq, c6, c12) reduction(+ : epot, cx[:sys->natoms], cy[:sys->natoms], cz[:sys->natoms])
#endif
	for (int i = 0; i < (sys->natoms - 1); i += sys->nsize) {
		int ii = i + sys->mpirank;
		if (ii < (sys->natoms - 1)) {
			for (int j = ii + 1; j < (sys->natoms); ++j) {
				/* get distance between particle i and j */
				double prx = pbc(rx[ii] - rx[j], 0.5 * sys->box);
				double pry = pbc(ry[ii] - ry[j], 0.5 * sys->box);
				double prz = pbc(rz[ii] - rz[j], 0.5 * sys->box);
				double rsq = (prx * prx + pry * pry + prz * prz);

				/* compute force and energy if within cutoff */
				if (rsq < rcsq) {
					double rsqinv = 1.0 / rsq;
					double r6 = rsqinv * rsqinv * rsqinv;
					double ffac = 12.0 * c12 * r6 * r6 * rsqinv - 6.0 * c6 * r6 * rsqinv;
					epot += c12 * r6 * r6 - c6 * r6;

					/* #ifdef _OPENMP */
					/* 	#pragma omp atomic */
					/* 				sys->cx[ii] += prx * ffac; */
					/* 	#pragma omp atomic */
					/* 				sys->cy[ii] += pry * ffac; */
					/* 	#pragma omp atomic */
					/* 				sys->cz[ii] += prz * ffac; */
					/* 	#pragma omp atomic */
					/* 				sys->cx[j] -= prx * ffac; */
					/* 	#pragma omp atomic */
					/* 				sys->cy[j] -= pry * ffac; */
					/* 	#pragma omp atomic */
					/* 				sys->cz[j] -= prz * ffac; */
					/* #else */
					cx[ii] += prx * ffac;
					cy[ii] += pry * ffac;
					cz[ii] += prz * ffac;
					cx[j] -= prx * ffac;
					cy[j] -= pry * ffac;
					cz[j] -= prz * ffac;
					/* #endif */
				}
			}
		}
	}

	MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
	MPI_Reduce(&epot, &(sys->epot), 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
}

