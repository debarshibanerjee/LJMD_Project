#include "prototypes.h"
#include "variables.h"

/* compute kinetic energy */


void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    double prefactor = 0.5*mvsq2e*sys->mass;
    int N_atoms = sys->natoms;
    double factor_temp= (2.0/((3.0*N_atoms-3.0)*kboltz));
    for (i=0; i<N_atoms; ++i) {
        sys->ekin += prefactor*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = factor_temp*sys->ekin;
}
