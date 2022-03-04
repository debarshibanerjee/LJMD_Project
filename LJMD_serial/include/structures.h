#ifndef STRUCTURES_H
#define STRUCTURES_H

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
	int nsize;
	int nfi, nsteps, natoms;
	double dt, mass, epsilon, box, sigma, rcut;
	double ekin, temp, epot;
        double *fx, *fy, *fz;
	double *vx, *vy, *vz;
        double *rx, *ry, *rz;
        
};
typedef struct _mdsys mdsys_t;

#endif
