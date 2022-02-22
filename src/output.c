<<<<<<< HEAD
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "structures.h"
=======
>>>>>>> 41feee22558ca29c77c5d3768789fc2e45ad37bd
#include "prototypes.h"

/* append data to output. */
static void output(mdsys_t* sys, FILE* erg, FILE* traj) {
	int i;

<<<<<<< HEAD
    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
    for (i=0; i<sys->natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
    }
=======
	printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot,
		   sys->ekin + sys->epot);
	fprintf(erg, "% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin,
			sys->epot, sys->ekin + sys->epot);
	fprintf(traj, "%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin + sys->epot);
	for (i = 0; i < sys->natoms; ++i) {
		fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
	}
>>>>>>> 41feee22558ca29c77c5d3768789fc2e45ad37bd
}
