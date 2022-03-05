#include "utilities.h"
#include <sys/time.h>
#include <time.h>

double wallclock() {
	struct timeval t;  // struct defined in time.h
	gettimeofday(&t, 0);
	return ((double)t.tv_sec) + 1.0e-6 * ((double)t.tv_usec);
}

/* helper function: sleep for some time */
void doublesleep(double t) {
	struct timespec ts;
	ts.tv_sec = (time_t)t;
	ts.tv_nsec = (long)((t - (double)ts.tv_sec) * 1000000000.0);
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

void ordering_atoms(mdsys_t* const sys){
	int counter, cell_index;
	int cell_x,cell_y,cell_z;
	double local_posx,local_posy,local_posz;
	double cell_size= sys->box/sys->ncel_d;
	int half_box = sys->box*0.5;
	int n=sys->ncel_d;

	for(int i=0;i<sys->ncells;i++) sys->cell_counter[i]=0;

	for(int at=0;at<sys->natoms;at++){
	        local_posx= pbc(sys->rx[at],half_box);	
	        local_posy= pbc(sys->ry[at],half_box);	
	        local_posz= pbc(sys->rz[at],half_box);	

		cell_x=(local_posx+half_box)/cell_size;
		cell_y=(local_posy+half_box)/cell_size;
		cell_z=(local_posz+half_box)/cell_size;

		cell_index= cell_z+n*cell_y+n*n*cell_x;

		counter=sys->cell_counter[cell_index];
		sys->clist[cell_index].idxlist[counter]=at;
		sys->cell_counter[cell_index]+=1;
	}
	for(int i=0;i<sys->ncells;i++) sys->clist[i].natoms=sys->cell_counter[i]-1;
}

void cell_localization(mdsys_t* const sys){
	int counter=0;
	int n=sys->ncel_d;
	for(int i=0;i<n;i++){
		for(int j=0; j<n;j++){
			for(int k=0;k<n;k++){
				sys->clist[counter].x_pos=i;
				sys->clist[counter].y_pos=j;
				sys->clist[counter].z_pos=k;
				counter++;
			}
		}
	}
}

void pairlist_creation(mdsys_t* const sys){
	int counter=0,cell_index;
	int x_pos,y_pos,z_pos;
	int near_xpos,near_ypos,near_zpos;
	int n=sys->ncel_d;
	int N=sys->ncells;
	for(int c=0;c<N;c++){
		x_pos=sys->clist[c].x_pos;
		y_pos=sys->clist[c].y_pos;
		z_pos=sys->clist[c].z_pos;
		for(int i=-1;i<2;i++){
			for(int j=-1; j<2;j++){
				for(int k=-1;k<2;k++){
					near_xpos=x_pos+i;
					near_ypos=y_pos+j;
					near_zpos=z_pos+k;
					cell_index=near_zpos+n*near_ypos+n*n*near_xpos;
					if(near_xpos>n-1||near_xpos<0) continue;
					if(near_ypos>n-1||near_ypos<0) continue;
					if(near_zpos>n-1||near_zpos<0) continue;
					sys->plist[2*counter]=c;
					sys->plist[2*counter+1]=cell_index;
					counter++;
				}
			}
		}
	}
	sys->npairs = counter;
}

void allocate_sys_arrays(mdsys_t* const sys) {
	// allocate coordinates arrays
	sys->rx = (double*)malloc(sys->natoms * sizeof(double));
	sys->ry = (double*)malloc(sys->natoms * sizeof(double));
	sys->rz = (double*)malloc(sys->natoms * sizeof(double));

	// allocate velocity arrays
	sys->vx = (double*)malloc(sys->natoms * sizeof(double));
	sys->vy = (double*)malloc(sys->natoms * sizeof(double));
	sys->vz = (double*)malloc(sys->natoms * sizeof(double));

	// allocate forces arrays
	sys->fx = (double*)malloc(sys->natoms * sizeof(double));
	sys->fy = (double*)malloc(sys->natoms * sizeof(double));
	sys->fz = (double*)malloc(sys->natoms * sizeof(double));

	// allocate support array for forces
	sys->cx = (double*)malloc(sys->nthreads * sys->natoms * sizeof(double));
	sys->cy = (double*)malloc(sys->nthreads * sys->natoms * sizeof(double));
	sys->cz = (double*)malloc(sys->nthreads * sys->natoms * sizeof(double));

	// allocate cells
	sys->clist = (cell_t*)malloc(sys->ncells * sizeof(cell_t));

	// allocate cell_counter
	sys->cell_counter = (int*)malloc(sys->ncells * sizeof(int));

	// allocate pair_list
	sys->plist= (int*)malloc(2*26*(sys->ncells) * sizeof(int));

	// allocate atom indices
	for(int i=0;i<sys->ncells;i++)
		sys->clist[i].idxlist=(int*)malloc(sys->natoms * sizeof(int));
}

void free_sys_arrays(mdsys_t* const sys) {
	// free coordinates
	free(sys->rx);
	free(sys->ry);
	free(sys->rz);
	sys->rx = NULL;
	sys->ry = NULL;
	sys->rz = NULL;

	// free velocities
	free(sys->vx);
	free(sys->vy);
	free(sys->vz);
	sys->vx = NULL;
	sys->vy = NULL;
	sys->vz = NULL;

	// free forces
	free(sys->fx);
	free(sys->fy);
	free(sys->fz);
	sys->fx = NULL;
	sys->fy = NULL;
	sys->fz = NULL;

	// free support arrays
	free(sys->cx);
	free(sys->cy);
	free(sys->cz);
	sys->cx = NULL;
	sys->cy = NULL;
	sys->cz = NULL;

	// free cells allocations
	free(sys->plist);
	free(sys->clist);
	free(sys->cell_counter);
}
