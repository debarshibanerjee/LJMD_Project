#include "utilities.h"
#include <sys/time.h>
#include <time.h>

double wallclock()
{
    struct timeval t; //struct defined in time.h
    gettimeofday(&t,0);
    return ((double) t.tv_sec) + 1.0e-6*((double) t.tv_usec);
}

/* helper function: sleep for some time */
void doublesleep(double t)
{
    struct timespec ts;
    ts.tv_sec = (time_t) t;
    ts.tv_nsec = (long)((t-(double)ts.tv_sec)*1000000000.0);
    nanosleep(&ts, NULL);
}


/* helper function: zero out an array */
void azzero(double* d, const int n) {
	int i;
	for (i = 0; i < n; ++i) {
		d[i] = 0.0;
	}
}

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2) {
	while (x > boxby2)
		x -= 2.0 * boxby2;
	while (x < -boxby2)
		x += 2.0 * boxby2;
	return x;
}


void allocate_sys_arrays ( mdsys_t * const sys ) {

  // allocate coordinates arrays
  sys->rx = (double *) malloc( sys->natoms * sizeof(double) );
  sys->ry = (double *) malloc( sys->natoms * sizeof(double) );
  sys->rz = (double *) malloc( sys->natoms * sizeof(double) );

  // allocate velocity arrays
  sys->vx = (double *) malloc( sys->natoms * sizeof(double) );
  sys->vy = (double *) malloc( sys->natoms * sizeof(double) );
  sys->vz = (double *) malloc( sys->natoms * sizeof(double) );

  // allocate forces arrays
  sys->fx = (double *) malloc( sys->natoms * sizeof(double) );
  sys->fy = (double *) malloc( sys->natoms * sizeof(double) );
  sys->fz = (double *) malloc( sys->natoms * sizeof(double) );

#ifdef USE_MPI
  // allocate support array fo centers
  sys->cx = (double *) malloc( sys->natoms * sizeof(double) );
  sys->cy = (double *) malloc( sys->natoms * sizeof(double) );
  sys->cz = (double *) malloc( sys->natoms * sizeof(double) );
#endif //USE_MPI

  // allocate forces arrays
  sys->fx = (double *) malloc( sys->natoms * sizeof(double) );
  sys->fy = (double *) malloc( sys->natoms * sizeof(double) );
  sys->fz = (double *) malloc( sys->natoms * sizeof(double) );
  

}

void free_sys_arrays ( mdsys_t * const sys ) {

  // free coordinates
  free( sys->rx );
  free( sys->ry );
  free( sys->rz );
  
  // free velocities
  free( sys->vx );
  free( sys->vy );
  free( sys->vz );
  
  
  // free forces
  free( sys->fx );
  free( sys->fy );
  free( sys->fz );

#ifdef USE_MPI
  // free support
  free( sys->cx );
  free( sys->cy );
  free( sys->cz );
#endif //USE_MPI
  

}
