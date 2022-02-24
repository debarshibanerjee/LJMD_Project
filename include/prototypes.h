#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include "structures.h"
#include "variables.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

//I/O function
extern int get_a_line(FILE *fp, char *buf);
extern void output(mdsys_t *sys, FILE *erg, FILE *traj);


//physical functions

extern void force(mdsys_t *sys);
extern void force_optimized_with3LawN(mdsys_t *sys)
extern void force_optimized_with3LawN_more_opt(mdsys_t *sys);
extern void ekin(mdsys_t *sys);
extern void verlet_vel_update(mdsys_t *sys);
extern void verlet_vel_propagation(mdsys_t *sys);
extern void velverlet(mdsys_t *sys);
extern void verlet_1(mdsys_t *sys);
extern void verlet_2(mdsys_t *sys);


#ifdef __cplusplus
}
#endif

#endif

