#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "structures.h"
#include "variables.h"

//I/O function
int get_a_line(FILE *fp, char *buf);
void output(mdsys_t *sys, FILE *erg, FILE *traj);

//utilities
double wallclock();
void azzero(double *d, const int n);
double pbc(double x, const double boxby2);

//physical functions

void force(mdsys_t *sys);
void ekin(mdsys_t *sys);
void verlet_vel_update(mdsys_t *sys);
void verlet_vel_propagation(mdsys_t *sys);
static void velverlet(mdsys_t *sys);




#endif

