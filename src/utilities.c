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
