/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */
//#include "variables.h"
#include "prototypes.h"
#include "utilities.h"

/* main */
int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int nprint, i;
	char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
	FILE *fp, *traj, *erg;
	mdsys_t sys;
	double t_start, t_pre_output;
	double t_IO = 0;

	t_start = wallclock();

	/*set MPI parameters*/
	sys.mpirank = rank;
	sys.nsize = size;
	sys.mpicomm = MPI_COMM_WORLD;

#ifdef _OPENMP
	#pragma omp parallel
	sys.nthreads = omp_get_num_threads();
#else
	sys.nthreads = 1;
#endif

	if (sys.mpirank == 0) {
		printf("LJMD version %3.1f\n", LJMD_VERSION);
		printf("Number of OpenMP Threads: %d\n", sys.nthreads);

		populate_data(stdin, &line, &restfile, &trajfile, &ergfile, &sys, &nprint);
		/* read input file */

		// Morse parameters for Argon
#ifdef MORSE
		sys.De = sys.epsilon;		   // well depth
		sys.re = 1.12246 * sys.sigma;  // 2^(1/6) * sigma
		sys.a_m = 1.05;				   // guesstimation
#endif
	}
	/*Sending data to all MPI processors*/
	MPI_Bcast(&(sys.natoms), 1, MPI_INT, 0, sys.mpicomm);
	MPI_Bcast(&(sys.nsteps), 1, MPI_INT, 0, sys.mpicomm);
	MPI_Bcast(&(sys.mass), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.epsilon), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.sigma), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.box), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.rcut), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.dt), 1, MPI_DOUBLE, 0, sys.mpicomm);
#ifdef MORSE
	MPI_Bcast(&(sys.De), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.re), 1, MPI_DOUBLE, 0, sys.mpicomm);
	MPI_Bcast(&(sys.a_m), 1, MPI_DOUBLE, 0, sys.mpicomm);
#endif

	if (sys.mpirank == 0) {
		printf("Communication time: %10.3fs\n", wallclock() - t_start);
	}

	int ncell_perdim;
	int cells = sys.box / sys.rcut;
	ncell_perdim = (cells % 2 == 0) ? cells + 1 : cells;
	sys.ncel_d = ncell_perdim;
	sys.ncells = ncell_perdim * ncell_perdim * ncell_perdim;

	/* allocate memory */
	allocate_sys_arrays(&sys);

	/* cell and pair list creation*/
	cell_localization(&sys);
	pairlist_creation(&sys);

	/* read restart */
	if (sys.mpirank == 0) {
		fp = fopen(restfile, "r");
		if (fp) {
			for (i = 0; i < sys.natoms; ++i) {
				fscanf(fp, "%lf%lf%lf", sys.rx + i, sys.ry + i, sys.rz + i);
			}
			for (i = 0; i < sys.natoms; ++i) {
				fscanf(fp, "%lf%lf%lf", sys.vx + i, sys.vy + i, sys.vz + i);
			}
			fclose(fp);
			azzero(sys.fx, sys.natoms);
			azzero(sys.fy, sys.natoms);
			azzero(sys.fz, sys.natoms);
		} else {
			perror("cannot read restart file");
			return 3;
		}
	}

	/* initialize forces and energies.*/
	MPI_Barrier(sys.mpicomm);
	sys.nfi = 0;

	/* Indexing atoms by cell location, only working for serial version */
	if(sys.mpirank==0)
	       	ordering_atoms(&sys);

	/*calling the Force function*/
#ifdef SIMPLE
	if (sys.mpirank == 0)
		printf("Using simpler force function.\n");
	force_omp_simple(&sys);	 // simpler OMP+MPI implementation
#elif MORSE
	if (sys.mpirank == 0)
		printf("Using Morse Potential Function.\n");
	force_morse(&sys);
#else
	if (sys.mpirank == 0)
		printf("Using default force function.\n");
	force(&sys);  // better OMP+MPI implementation
#endif

	ekin(&sys);

	if (sys.mpirank == 0) {
		erg = fopen(ergfile, "w");
		traj = fopen(trajfile, "w");

		printf("Startup time: %10.3fs\n", wallclock() - t_start);
		printf("Starting simulation with %d atoms for %d steps.\n", sys.natoms, sys.nsteps);
		printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
		output(&sys, erg, traj);
	}

	/* reset timer */
	t_start = wallclock();

	/**************************************************/
	/* main MD loop */
	for (sys.nfi = 1; sys.nfi <= sys.nsteps; ++sys.nfi) {
		/* write output, if requested */
		if (sys.mpirank == 0) {
			if ((sys.nfi % nprint) == 0) {
				t_pre_output = wallclock();
				output(&sys, erg, traj);
				t_IO += (wallclock() - t_pre_output);
			}
		}

		/* propagate system and recompute energies */
		velverlet(&sys);
		ekin(&sys);
	}
	/**************************************************/

	/* clean up: close files, free memory */
	if (sys.mpirank == 0) {
		printf("Simulation Done. Run time: %10.3fs\n", wallclock() - t_start);
		printf("I/O time : %10.3fs\n", t_IO);
		fclose(erg);
		fclose(traj);
	}
	free_sys_arrays(&sys);

	MPI_Finalize();
	return 0;
}
