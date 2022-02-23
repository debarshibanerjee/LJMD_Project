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



#ifdef __cplusplus
extern "C"
{
#endif

//I/O function
int get_a_line(FILE *fp, char *buf);
void output(mdsys_t *sys, FILE *erg, FILE *traj);



//physical functions

void force(mdsys_t *sys);
void ekin(mdsys_t *sys);
void verlet_vel_update(mdsys_t *sys);
void verlet_vel_propagation(mdsys_t *sys);
static void velverlet(mdsys_t *sys);
extern void verlet_1(mdsys_t *sys);
extern void verlet_2(mdsys_t *sys);


#ifdef __cplusplus
}
#endif

#endif

