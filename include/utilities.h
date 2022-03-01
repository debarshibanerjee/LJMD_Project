#ifndef UTILITIES_H
#define UTILITIES_H

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "structures.h"


#ifdef __cplusplus
extern "C"
{
#endif

//utilities
extern double wallclock();
extern void azzero(double *d, const int n);
extern double pbc(double x, const double boxby2);
extern void doublesleep(double t);
extern void allocate_sys_arrays ( mdsys_t * const sys );
extern void free_sys_arrays ( mdsys_t * const sys );

#ifdef __cplusplus
}
#endif

#endif

