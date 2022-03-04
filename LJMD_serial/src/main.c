/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */
//#include "variables.h"
#include "prototypes.h"
#include "utilities.h"


#define FLAG_OPTIMIZATION yes
#define PRINT_OPTION

/* main */
int main(int argc, char **argv)
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double t_start;

    printf("LJMD version %3.1f\n", LJMD_VERSION);

    t_start = wallclock();

    /* read input file */
    populate_data(stdin, &line, &restfile, &trajfile, &ergfile, &sys, &nprint);
    

    /* allocate memory */

    allocate_sys_arrays(&sys);
    

    /* read restart */
    fp=fopen(restfile,"r");
    const int N_atoms = sys.natoms;
    if(fp) {
        
        for ( i=0; i< N_atoms; ++i) {
            fscanf(fp,"%lf%lf%lf", sys.rx + i, sys.ry + i, sys.rz + i);
        }
        for ( i=0; i< N_atoms; ++i) {
            fscanf(fp,"%lf%lf%lf", sys.vx + i, sys.vy + i, sys.vz + i);
        }
        fclose(fp);
        
        //azzero(sys.fx, N_atoms);
        //azzero(sys.fy, N_atoms);
        //azzero(sys.fz, N_atoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* initialize forces and energies.*/
    sys.nfi=0;
    
  
    force_optimized_with3LawN(&sys);

    //force_optimized(&sys);
   
    //force(&sys);
   

    ekin(&sys);

    erg = fopen(ergfile,"w");
    traj= fopen(trajfile,"w");

    printf("Startup time: %10.3fs\n", wallclock()-t_start);
    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /* reset timer */
    t_start = wallclock();

    /**************************************************/
    /* main MD loop */

    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        #ifdef PRINT_OPTION
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);
        #endif

        /* propagate system and recompute energies */
        velverlet(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);
    fclose(erg);
    fclose(traj);

    free_sys_arrays(&sys);


    return 0;
}
