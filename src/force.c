#include "utilities.h"
#include "prototypes.h"
/* compute forces */
void force(mdsys_t* sys) {
	double r, ffac;
	double rx, ry, rz;
	int i, j;

	/* zero energy and forces */
	sys->epot = 0.0;
	azzero(sys->fx, sys->natoms);
	azzero(sys->fy, sys->natoms);
	azzero(sys->fz, sys->natoms);
	azzero(sys->cx, sys->natoms);
       	azzero(sys->cy, sys->natoms); 
	azzero(sys->cz, sys->natoms);

//	MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
//	MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
//	MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);

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
			r = sqrt(rx * rx + ry * ry + rz * rz);

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


/* compute forces with third Newton Law applied */
void force_optimized_with3LawN(mdsys_t *sys) {
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms)-1; ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {

                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r+ 6*pow(sys->sigma/r,6.0)/r);

                sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0) -pow(sys->sigma/r,6.0));

                sys->fx[i] += rx/r*ffac;
                sys->fx[j] -= rx/r*ffac;
                sys->fy[i] += ry/r*ffac;
                sys->fy[j] -= ry/r*ffac;
                sys->fz[i] += rz/r*ffac;
                sys->fz[j] -= rz/r*ffac;
            }
        }
    }
}



/* compute forces with third Newton Law applied, more fucntion modifications*/
void force_optimized_with3LawN_more_opt(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    double c12, c6, rcsq, s6;
    s6 = sys->sigma*sys->sigma*sys->sigma*sys->sigma*sys->sigma*sys->sigma;
    c12=4.0*sys->epsilon*s6*s6;
    c6 =4.0*sys->epsilon*s6;
    rcsq = sys->rcut * sys->rcut;
    for(i=0; i < (sys->natoms)-1; ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                double r6,rinv; 
                rinv=1.0/(r*r); 
                r6=rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6;
                sys->epot  += r6*(c12*r6 - c6);
                sys->fx[i] += rx/r*ffac;
                sys->fx[j] -= rx/r*ffac;
                sys->fy[i] += ry/r*ffac;
                sys->fy[j] -= ry/r*ffac;
                sys->fz[i] += rz/r*ffac;
                sys->fz[j] -= rz/r*ffac;
            }
        }
    }
}



