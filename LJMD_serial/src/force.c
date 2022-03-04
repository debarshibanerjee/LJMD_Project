#include "prototypes.h"
#include "utilities.h"




/* compute forces */
void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
        for(j=0; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);

            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r + 6*pow(sys->sigma/r,6.0)/r);

                sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0) -pow(sys->sigma/r,6.0));

                sys->fx[i] += rx/r*ffac;
                sys->fy[i] += ry/r*ffac;
                sys->fz[i] += rz/r*ffac;
            }
        }
    }
}




/* compute forces with third Newton Law applied, more fucntion modifications*/
void force_optimized(mdsys_t* sys) {
	double r, ffac;
	double rx, ry, rz;
	int i, j;
	double factor;
    	double A_x, A_y, A_z;
    	double r6, rinv_sq;
        double c12, c6, s6;
	double half_box;
	/* zero energy and forces */
	sys->epot = 0.0;
	azzero(sys->fx, sys->natoms);
	azzero(sys->fy, sys->natoms);
	azzero(sys->fz, sys->natoms);
	half_box = 0.5 * sys->box;
	double sigma = sys->sigma;
	s6 = sigma * sigma * sigma * sigma * sigma * sigma;
	c6 = 4.0 * sys->epsilon * s6;
	c12 = c6 * s6;
	for (i = 0; i < (sys->natoms); ++i) {
		for (j = 0; j < (sys->natoms); ++j) {
			/* particles have no interactions with themselves */
			if (i == j)
				continue;
			/* get distance between particle i and j */
			rx = pbc(sys->rx[i] - sys->rx[j], half_box);
			ry = pbc(sys->ry[i] - sys->ry[j], half_box);
			rz = pbc(sys->rz[i] - sys->rz[j], half_box);
			r = sqrt(rx * rx + ry * ry + rz * rz);

			/* compute force and energy if within cutoff */
			if (r < sys->rcut) {
				
				rinv_sq= 1.0 / (r * r);
				r6 = rinv_sq * rinv_sq * rinv_sq;
				sys->epot +=  (c12 * r6  - c6 ) * r6 * 0.5;
				ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6;
				factor = ffac * rinv_sq;
            			A_x = factor * rx;
            			A_y = factor * ry;
            			A_z = factor * rz;
            			sys->fx[i] += A_x;
				sys->fy[i] += A_y;
				sys->fz[i] += A_z;

				
			}
		}
	}
}



/* compute forces with third Newton Law applied */
void force_with3LawN(mdsys_t* sys) {
	double r, ffac;
	double rx, ry, rz;
	int i, j;

	/* zero energy and forces */
	sys->epot = 0.0;
	azzero(sys->fx, sys->natoms);
	azzero(sys->fy, sys->natoms);
	azzero(sys->fz, sys->natoms);

	for (i = 0; i < (sys->natoms) - 1; ++i) {
		for (j = i + 1; j < (sys->natoms); ++j) {
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
				ffac = -4.0 * sys->epsilon *(-12.0 * pow(sys->sigma / r, 12.0) / r + 6 * pow(sys->sigma / r, 6.0) / r);

				sys->epot += 0.5 * 4.0 * sys->epsilon *(pow(sys->sigma / r, 12.0) - pow(sys->sigma / r, 6.0));

				sys->fx[i] += rx / r * ffac;
				sys->fx[j] -= rx / r * ffac;
				sys->fy[i] += ry / r * ffac;
				sys->fy[j] -= ry / r * ffac;
				sys->fz[i] += rz / r * ffac;
				sys->fz[j] -= rz / r * ffac;
			}
		}
	}
}

/* compute forces with third Newton Law applied, more fucntion modifications*/
void force_optimized_with3LawN(mdsys_t* sys) {
	double rx, ry, rz;
	int i, j;
        
        /* factors useful for optimization */
        
        const double rcsq = sys->rcut * sys->rcut;
        const double half_box= 0.5 * sys->box;
        const double epsilon = sys->epsilon;
	const double sigma = sys->sigma;
        const int N_atoms = sys->natoms;
        const double s6 = sigma * sigma * sigma * sigma * sigma * sigma;
        const double c6 = 4.0 * epsilon * s6;
	const double c12 = c6 * s6;
        //double factor;
    	double r6, rsq, rinv_sq, ffac;
        double A_x, A_y, A_z;

	/* zero energy and forces */
	sys->epot = 0.0;
        
	azzero(sys->fx, N_atoms);
	azzero(sys->fy, N_atoms);
	azzero(sys->fz, N_atoms);
	
	for (i = 0; i < (N_atoms - 1); ++i) {
		for (j = i + 1; j < N_atoms; ++j) {
			/* particles have no interactions with themselves */
			if (i == j)
				continue;
			/* get distance between particle i and j */
			rx = pbc(sys->rx[i] - sys->rx[j], half_box);
			ry = pbc(sys->ry[i] - sys->ry[j], half_box);
			rz = pbc(sys->rz[i] - sys->rz[j], half_box);
			rsq = (rx * rx + ry * ry + rz * rz);
                        rinv_sq = 1.0 / (rsq);

			/* compute force and energy if within cutoff */
			if (rsq < rcsq) {
				r6 = rinv_sq * rinv_sq * rinv_sq;
                                //var1 = c12 * r6 * r6;
                                //var2 = c6 * r6;
				sys->epot += r6 * (r6* c12 - c6);
				ffac = 6.0* r6 * ( r6 * c12 * 2.0 - c6 ) * rinv_sq;
				//factor = ffac * rinv_sq;                                
            		        A_x = ffac * rx ;
            			sys->fx[i] += A_x;
				sys->fx[j] -= A_x;
                                A_y = ffac * ry ;
				sys->fy[i] += A_y;
				sys->fy[j] -= A_y;
                                A_z = ffac * rz ;
				sys->fz[i] += A_z;
				sys->fz[j] -= A_z;
				
			}
		}
	}
}

